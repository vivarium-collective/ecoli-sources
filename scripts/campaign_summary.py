#!/usr/bin/env python3
"""Generate a campaign summary HTML report from a synced sms-api pilot dir.

Reads:
  - nextflow/workflow_config.json (or _local.json) for variant -> dataset mapping
  - nextflow/nextflow_workdirs/*/* for task status, durations, error reasons
  - history/experiment_id=*/variant=*/seed=*/.../generation=*/agent_id=*/*.pq for sim outputs

Writes:
  - analyses/campaign_summary/index.html  — master report
  - analyses/campaign_summary/status.tsv  — per-task table
"""

import json
import os
import re
import sys
from pathlib import Path

import duckdb
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots

PILOT = Path(__file__).resolve().parent
CONFIG = json.loads((PILOT / "nextflow" / "workflow_config_local.json").read_text())
EXP_ID = CONFIG["experiment_id"]
PARCA_VARIANTS = CONFIG["parca_variants"]
WORKDIRS = PILOT / "nextflow" / "nextflow_workdirs"
OUTDIR = PILOT / "analyses" / "campaign_summary"
OUTDIR.mkdir(parents=True, exist_ok=True)

# How many sims/seeds would run *per variant* if its ParCa succeeded.
SEEDS_PER_VARIANT = int(CONFIG.get("n_init_sims", 1))
GENS_PER_VARIANT = int(CONFIG.get("generations", 1))


def _load_perturbation_map() -> dict[str, str]:
    """Build dataset_id -> human-readable perturbation label.

    Reads the campaign sidecar (referenced by sensitivity_overview's
    `campaign_sidecar`) to map opaque hashed dataset_ids back to the
    perturbation parameters that produced them. Generated IDs are emitted
    in the same order as the cartesian product of `param_grid`, so we can
    zip them.
    """
    sidecar_path = None
    sens = (CONFIG.get("analysis_options", {})
                  .get("multivariant", {})
                  .get("sensitivity_overview", {}))
    if "campaign_sidecar" in sens:
        sidecar_path = Path(sens["campaign_sidecar"])
    if sidecar_path is None or not sidecar_path.exists():
        return {}

    sidecar = json.loads(sidecar_path.read_text())
    spec = sidecar.get("spec", {})
    operator = spec.get("operator", "?")
    param_grid = spec.get("param_grid", {})
    generated_ids = sidecar.get("generated_dataset_ids", [])

    # Cartesian product of param_grid in the same key order Python iterates.
    import itertools
    keys = list(param_grid.keys())
    combos = list(itertools.product(*[param_grid[k] for k in keys]))

    label_map: dict[str, str] = {}
    for did, combo in zip(generated_ids, combos):
        kv = ", ".join(f"{k}={v}" for k, v in zip(keys, combo))
        label_map[did] = f"{operator}({kv})"
    # Source dataset → "baseline"
    src = spec.get("source_dataset_id")
    if src:
        label_map[src] = "baseline (no perturbation)"
    return label_map


def _last_meaningful(text: str, n: int = 1) -> str:
    """Return the last n non-blank lines, joined by ' / '."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    return " / ".join(lines[-n:]) if lines else ""


_REAL_ERROR_PATTERNS = (
    "ValueError:", "RuntimeError:", "KeyError:", "TypeError:", "AttributeError:",
    "AssertionError:", "FileNotFoundError:", "ImportError:",
)


def _extract_python_error(text: str) -> str:
    """Pull the last Python exception line, ignoring the post-failure aws-cli noise."""
    last = ""
    for ln in text.splitlines():
        s = ln.strip()
        if any(p in s for p in _REAL_ERROR_PATTERNS):
            last = s
    return last or _last_meaningful(text, n=1)


def _scan_workdirs() -> list[dict]:
    """One row per Nextflow task workdir."""
    rows = []
    for hh_dir in sorted(WORKDIRS.iterdir()):
        if not hh_dir.is_dir():
            continue
        for task_dir in sorted(hh_dir.iterdir()):
            if not task_dir.is_dir():
                continue
            sh = task_dir / ".command.sh"
            if not sh.exists():
                continue
            sh_text = sh.read_text()
            # Identify the process by first python/script reference
            if "runscripts/parca.py" in sh_text:
                kind = "runParca"
                m = re.search(r"parca_config_(\d+)\.json", sh_text)
                parca_id = int(m.group(1)) if m else -1
                tag = f"parca_{parca_id}"
                meta = {"parca_id": parca_id}
            elif "runscripts/create_variants.py" in sh_text:
                kind = "createVariants"
                m = re.search(r"--parca-id\s+(\d+)", sh_text)
                parca_id = int(m.group(1)) if m else -1
                tag = f"createVariants_{parca_id}"
                meta = {"parca_id": parca_id}
            elif "ecoli_master_sim.py" in sh_text:
                kind = "sim"
                m_v = re.search(r"--variant\s+(\d+)", sh_text)
                m_s = re.search(r"--seed\s+(\d+)", sh_text)
                m_g = re.search(r"--lineage_seed\s+(\d+)", sh_text)
                v = int(m_v.group(1)) if m_v else -1
                s = int(m_s.group(1)) if m_s else -1
                lineage = int(m_g.group(1)) if m_g else -1
                tag = f"sim_v{v}_s{s}"
                meta = {"variant": v, "seed": s, "lineage_seed": lineage}
            elif "metadata.json" in sh_text and "fsspec" in sh_text:
                kind = "mergeVariantMetadata"
                tag = "mergeVariantMetadata"
                meta = {}
            else:
                kind = "other"
                tag = task_dir.name[:8]
                meta = {}

            # Status & timing
            exitcode_file = task_dir / ".exitcode"
            exit_code = exitcode_file.read_text().strip() if exitcode_file.exists() else "?"
            begin_file = task_dir / ".command.begin"
            log_file = task_dir / ".command.log"
            try:
                t_begin = begin_file.stat().st_mtime if begin_file.exists() else None
                t_end = exitcode_file.stat().st_mtime if exitcode_file.exists() else None
                duration_s = (t_end - t_begin) if (t_begin and t_end) else None
            except OSError:
                duration_s = None

            # Failure reason — look for a real Python exception, not the
            # post-failure aws-cli message that always trails a failed task.
            err_text = ""
            if log_file.exists():
                err_text = _extract_python_error(log_file.read_text(errors="replace"))
            elif (task_dir / ".command.err").exists():
                err_text = _extract_python_error((task_dir / ".command.err").read_text(errors="replace"))

            rows.append({
                "kind": kind,
                "tag": tag,
                "exit_code": exit_code,
                "duration_s": round(duration_s, 1) if duration_s is not None else None,
                "error": err_text[:200] if exit_code != "0" else "",
                "workdir": f"{hh_dir.name}/{task_dir.name[:8]}",
                **meta,
            })
    return rows


def _build_status(rows: list[dict]) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Per-parca and per-sim status frames."""
    df = pl.DataFrame(rows)

    # Per-parca: dataset_id + outcome
    parca_rows = df.filter(pl.col("kind") == "runParca")
    parca_status = []
    for i, pv in enumerate(PARCA_VARIANTS):
        match = parca_rows.filter(pl.col("parca_id") == i)
        if len(match) == 0:
            parca_status.append({
                "parca_id": i,
                "dataset_id": pv["rnaseq_basal_dataset_id"],
                "exit_code": "?", "duration_s": None, "error": "(no workdir)",
            })
            continue
        r = match.row(0, named=True)
        parca_status.append({
            "parca_id": i,
            "dataset_id": pv["rnaseq_basal_dataset_id"],
            "exit_code": r["exit_code"],
            "duration_s": r["duration_s"],
            "error": r["error"],
        })
    parca_df = pl.DataFrame(parca_status)

    sim_df = df.filter(pl.col("kind") == "sim").select(
        "variant", "seed", "lineage_seed", "exit_code", "duration_s", "error", "workdir"
    ).sort("variant", "seed")
    return parca_df, sim_df


def _aa_synthesis_fold_change() -> tuple[pl.DataFrame, set[str]]:
    """
    Build a comprehensive enzyme fold-change table including:
      - Each AA's forward synthesis enzymes (summed)            label: "<AA> (synth)"
      - Each AA's reverse / degradation enzymes (summed)        label: "<AA> (deg)"
      - Each individual explicitly-adjusted RNA from
        rna_expression_adjustments.tsv                          label: "<gene> ★"

    Returns (fc_df, highlighted_labels).
    fc_df columns: row_label, group, variant, log2_fc, count.
    """
    import csv, json as _json, pickle
    import numpy as np

    repo = Path("/home/youdonotexist/code/omics-vEcoli")
    pathways_tsv = repo / "reconstruction/ecoli/flat/amino_acid_pathways.tsv"
    adjust_tsv = repo / "reconstruction/ecoli/flat/adjustments/rna_expression_adjustments.tsv"

    # ------- AA pathways: forward + reverse enzymes per AA -------
    aa_forward: dict[str, list[str]] = {}
    aa_reverse: dict[str, list[str]] = {}
    with pathways_tsv.open() as f:
        for row in csv.DictReader(
            (ln for ln in f if not ln.startswith("#")), delimiter="\t"
        ):
            aa = row.get("Amino acid")
            if not aa:
                continue
            for col, target in (("Enzymes", aa_forward), ("Reverse enzymes", aa_reverse)):
                try:
                    target[aa] = _json.loads(row.get(col, ""))
                except Exception:
                    target[aa] = []

    # ------- Explicitly adjusted RNAs (RNA cistron id -> gene comment) -------
    adjusted_rna: dict[str, str] = {}  # cistron_id -> short label
    if adjust_tsv.exists():
        with adjust_tsv.open() as f:
            for row in csv.DictReader(
                (ln for ln in f if not ln.startswith("#")), delimiter="\t"
            ):
                cistron = (row.get("name") or "").strip().strip('"')
                comment = (row.get("_comments") or "").strip().strip('"')
                if not cistron:
                    continue
                # Comment looks like "metA, homoserine O-succinyltransferase; ..."
                gene = comment.split(",", 1)[0].strip() if comment else cistron
                adjusted_rna[cistron] = gene

    # ------- sim_data lookups -------
    with open(PILOT / "parca_0" / "kb" / "simData.cPickle", "rb") as f:
        sd = pickle.load(f)
    arr = sd.process.translation.monomer_data.struct_array
    monomer_ids = [str(x) for x in arr["id"]]
    monomer_idx = {m: i for i, m in enumerate(monomer_ids)}
    cistron_for_monomer = {str(m): str(c) for m, c in zip(arr["id"], arr["cistron_id"])}
    monomers_for_cistron: dict[str, list[str]] = {}
    for m, c in cistron_for_monomer.items():
        monomers_for_cistron.setdefault(c, []).append(m)
    complexation = sd.process.complexation

    def _expand(complex_or_monomer: str) -> list[int]:
        """Return monomer indices for a complex (recursing once into sub-complexes)
        or a single monomer."""
        token = complex_or_monomer if "[" in complex_or_monomer else f"{complex_or_monomer}[c]"
        if token in monomer_idx:
            return [monomer_idx[token]]
        try:
            sub = complexation.get_monomers(token)
        except Exception:
            return []
        out: list[int] = []
        for m in sub.get("subunitIds", []):
            m = str(m)
            if m in monomer_idx:
                out.append(monomer_idx[m])
            else:
                try:
                    sub2 = complexation.get_monomers(m)
                    for m2 in sub2.get("subunitIds", []):
                        m2 = str(m2)
                        if m2 in monomer_idx:
                            out.append(monomer_idx[m2])
                except Exception:
                    pass
        return out

    # Build the row spec — list of (label, group, monomer_indices, is_highlight)
    row_specs: list[tuple[str, str, list[int], bool]] = []

    for aa in aa_forward.keys() | aa_reverse.keys():
        fwd_idxs: set[int] = set()
        for cplx in aa_forward.get(aa, []):
            fwd_idxs.update(_expand(cplx))
        if fwd_idxs:
            # ★-highlight if any of the underlying cistrons are in the adjusted list
            cistrons = {cistron_for_monomer.get(monomer_ids[i]) for i in fwd_idxs}
            highlight = bool(cistrons & adjusted_rna.keys())
            row_specs.append((f"{aa} (synth)", "AA synthesis", sorted(fwd_idxs), highlight))

        rev_idxs: set[int] = set()
        for cplx in aa_reverse.get(aa, []):
            rev_idxs.update(_expand(cplx))
        if rev_idxs:
            cistrons = {cistron_for_monomer.get(monomer_ids[i]) for i in rev_idxs}
            highlight = bool(cistrons & adjusted_rna.keys())
            row_specs.append((f"{aa} (deg)", "AA degradation", sorted(rev_idxs), highlight))

    # Individual rows for each explicitly-adjusted RNA — always highlighted.
    for cistron, gene in sorted(adjusted_rna.items(), key=lambda kv: kv[1].lower()):
        idxs: list[int] = []
        for m in monomers_for_cistron.get(cistron, []):
            if m in monomer_idx:
                idxs.append(monomer_idx[m])
        if idxs:
            row_specs.append((f"{gene} ({cistron})", "Adjusted RNA", sorted(set(idxs)), True))

    # ------- Pull mean monomer counts per (variant) from parquet -------
    glob = (
        f"{PILOT}/history/experiment_id={EXP_ID}/variant=*/lineage_seed=*/"
        f"generation=*/agent_id=*/*.pq"
    )
    con = duckdb.connect()
    sql = f"""
    SELECT variant, idx, avg(monomer) AS avg_count
    FROM (
      SELECT variant, generate_subscripts(listeners__monomer_counts, 1) AS idx,
             unnest(listeners__monomer_counts) AS monomer
      FROM read_parquet('{glob}', hive_partitioning=true)
    )
    GROUP BY variant, idx
    """
    long = con.execute(sql).pl()
    variants = sorted(set(long["variant"].to_list()))
    n_mono = max(long["idx"].max() or 0, len(monomer_ids))
    mono_avg: dict[int, np.ndarray] = {}
    for v in variants:
        arr_v = np.zeros(n_mono, dtype=np.float64)
        for r in long.filter(pl.col("variant") == v).iter_rows(named=True):
            arr_v[int(r["idx"]) - 1] = r["avg_count"]
        mono_avg[v] = arr_v

    if 0 not in mono_avg:
        return pl.DataFrame({"row_label": [], "group": [], "variant": [], "log2_fc": []}), set()

    out_rows = []
    highlighted_labels: set[str] = set()
    for label, group, idxs, highlight in row_specs:
        if highlight:
            highlighted_labels.add(label)
        baseline = float(mono_avg[0][idxs].sum())
        for v in variants:
            count = float(mono_avg[v][idxs].sum())
            log2_fc = float(np.log2((count + 1.0) / (baseline + 1.0)))
            out_rows.append({
                "row_label": label, "group": group,
                "variant": v, "log2_fc": log2_fc, "count": count,
            })
    return pl.DataFrame(out_rows), highlighted_labels


def _query_mass_data() -> pl.DataFrame:
    """Pull cell_mass and instantaneous_growth_rate per (variant, seed, time)."""
    glob = (
        f"{PILOT}/history/experiment_id={EXP_ID}/variant=*/lineage_seed=*/"
        f"generation=*/agent_id=*/*.pq"
    )
    con = duckdb.connect()
    sql = f"""
    SELECT
      variant,
      lineage_seed AS seed,
      generation,
      time,
      listeners__mass__cell_mass        AS cell_mass_fg,
      listeners__mass__dry_mass         AS dry_mass_fg,
      listeners__mass__protein_mass     AS protein_mass_fg,
      listeners__mass__rna_mass         AS rna_mass_fg,
      listeners__mass__instantaneous_growth_rate AS growth_rate
    FROM read_parquet('{glob}', hive_partitioning=true)
    ORDER BY variant, lineage_seed, time
    """
    return con.execute(sql).pl()


def _make_html(
    parca_df: pl.DataFrame,
    sim_df: pl.DataFrame,
    mass: pl.DataFrame,
    aa_fc: pl.DataFrame,
    aa_highlight: set[str],
) -> str:
    """Render the report HTML."""
    # Variant -> rich label that includes the perturbation if any. The raw
    # dataset_id is opaque (e.g. "...__add_log_normal_noise__e8e44a96") so we
    # fold in the sigma/seed from the campaign sidecar where available.
    perturbation_map = _load_perturbation_map()

    def _short_perturbation(label: str) -> str:
        # "add_log_normal_noise(sigma=0.2, seed=0)" -> "noise σ=0.2,s=0"
        if label.startswith("add_log_normal_noise("):
            inside = label[len("add_log_normal_noise("):].rstrip(")")
            inside = inside.replace("sigma", "σ").replace(", ", ",").replace("seed=", "s=")
            return f"noise {inside}"
        if label == "baseline (no perturbation)":
            return "baseline"
        return label

    def _rich_label(parca_id: int) -> str:
        ds = PARCA_VARIANTS[parca_id]["rnaseq_basal_dataset_id"]
        pert = perturbation_map.get(ds)
        if pert and pert != "(extra dataset)":
            return f"v{parca_id}: {_short_perturbation(pert)}"
        return f"v{parca_id}: {ds}"

    label_map = {i: _rich_label(i) for i in range(len(PARCA_VARIANTS))}
    mass = mass.with_columns(
        pl.col("variant").map_elements(lambda v: f"v{v}: {label_map.get(v, '?')}", return_dtype=pl.String).alias("label"),
        (pl.col("variant").cast(pl.String) + "/s" + pl.col("seed").cast(pl.String)).alias("trace"),
    )

    # 1. Cell-mass-over-time, color by variant, dashed by seed
    fig_mass = px.line(
        mass.to_pandas(),
        x="time", y="cell_mass_fg",
        color="label", line_dash="seed",
        title="Cell mass (fg) vs time — surviving variants × seeds",
        labels={"time": "time (s)", "cell_mass_fg": "cell mass (fg)", "label": "variant"},
    )
    # Tall enough that all 15 (variants × seeds) legend entries show
    # without scrolling. Keep the legend on the right (Plotly default)
    # so its grouping by variant remains visually obvious.
    fig_mass.update_layout(
        height=750,
        hovermode="x unified",
    )

    # 2. Mass-fraction stack (median across seeds, per variant)
    components = ["protein_mass_fg", "rna_mass_fg", "dry_mass_fg"]
    stack_data = mass.group_by(["label", "time"]).agg(
        [pl.col(c).median().alias(c) for c in components]
    ).sort("label", "time").to_pandas()
    fig_components = make_subplots(
        rows=1, cols=3, subplot_titles=["protein (fg)", "RNA (fg)", "dry mass (fg)"],
    )
    # Stable color per variant — Plotly's auto-coloring assigns sequential
    # colors per add_trace call, so without an explicit map each subplot
    # would get a fresh palette and the legend would only line up with
    # the last subplot.
    palette = px.colors.qualitative.Plotly
    sorted_labels = sorted(stack_data["label"].unique())
    color_for = {lbl: palette[i % len(palette)] for i, lbl in enumerate(sorted_labels)}
    for col_i, comp in enumerate(components):
        for label in sorted_labels:
            sub = stack_data[stack_data["label"] == label]
            fig_components.add_trace(
                go.Scatter(
                    x=sub["time"], y=sub[comp], name=label,
                    legendgroup=label, showlegend=(col_i == 0),
                    line={"color": color_for[label]},
                ),
                row=1, col=col_i + 1,
            )
    fig_components.update_layout(height=400, title="Mass components over time (median across seeds)")

    # 3. Growth rate distribution per variant (boxplot)
    fig_gr = px.box(
        mass.filter(pl.col("growth_rate").is_finite()).to_pandas(),
        x="label", y="growth_rate",
        title="Instantaneous growth rate distribution per variant",
        labels={"label": "variant", "growth_rate": "growth rate"},
    )
    fig_gr.update_layout(height=400, xaxis_tickangle=-20)

    # 3b. Enzyme fold-change heatmap (vs baseline = variant 0). Three
    # row groups: AA-synthesis (sum of forward enzymes per AA), AA-degradation
    # (sum of reverse), and individual explicitly-adjusted RNAs. Rows whose
    # underlying cistron is in rna_expression_adjustments.tsv are ★-marked.
    fig_aa = None
    if aa_fc.height:
        # Order rows: AA synthesis first, then AA degradation, then adjusted RNAs.
        # Within each group, sort alphabetically by row_label.
        group_order = {"AA synthesis": 0, "AA degradation": 1, "Adjusted RNA": 2}

        # Pivot to wide. Polars 1.x: pivot uses positional/column args.
        wide = aa_fc.with_columns(
            pl.col("variant").map_elements(
                lambda v: label_map.get(int(v), f"v{v}"), return_dtype=pl.String
            ).alias("variant_label"),
        ).pivot(
            values="log2_fc", index=["row_label", "group"], on="variant_label"
        ).with_columns(
            pl.col("group").map_elements(
                lambda g: group_order.get(g, 99), return_dtype=pl.Int64
            ).alias("_grp_ord"),
        ).sort("_grp_ord", "row_label")

        row_labels = wide["row_label"].to_list()
        groups = wide["group"].to_list()
        var_cols = [c for c in wide.columns if c not in ("row_label", "group", "_grp_ord")]
        order_lookup = {label_map[i]: i for i in range(len(PARCA_VARIANTS))}
        var_cols_sorted = sorted(var_cols, key=lambda c: order_lookup.get(c, 999))
        z = wide.select(var_cols_sorted).to_numpy()

        # Bold/star the highlighted rows.
        ytext = [
            f"<b>{lbl} ★</b>" if lbl in aa_highlight else lbl
            for lbl in row_labels
        ]

        # Use group-name suffix in the y label so the three sections are visible.
        ytext_grouped = []
        last_grp = None
        for lbl, grp in zip(ytext, groups):
            prefix = ""
            if grp != last_grp:
                prefix = f"<i>[{grp}]</i> "
                last_grp = grp
            ytext_grouped.append(prefix + lbl)

        fig_aa = go.Figure(go.Heatmap(
            z=z,
            x=var_cols_sorted,
            y=ytext_grouped,
            colorscale="RdBu_r", zmid=0,
            colorbar={"title": "log2(FC vs baseline)"},
            hovertemplate="row=%{y}<br>variant=%{x}<br>log2FC=%{z:.3f}<extra></extra>",
        ))
        fig_aa.update_layout(
            title=("Enzyme abundance — log2 fold change vs baseline (variant 0). "
                   f"★ = explicitly adjusted in rna_expression_adjustments.tsv "
                   f"({len(aa_highlight)} highlighted, {len(row_labels)} rows total)"),
            height=max(700, 18 * len(row_labels) + 200),
            xaxis_tickangle=-30,
            yaxis={"title": ""},
            margin={"l": 220},
        )

    # 4. Combined parca + sim status — left-outer join on dataset_id.
    # One row per (dataset_id, sim_seed). For datasets whose ParCa failed,
    # there are no sim rows, so we emit a single row with sim_* fields blank.
    def _perturbation(dataset_id: str) -> str:
        return perturbation_map.get(dataset_id, "(extra dataset)")

    parca_for_join = parca_df.rename({
        "exit_code": "parca_exit",
        "duration_s": "parca_dur_s",
        "error": "parca_error",
    }).with_columns(
        pl.col("dataset_id").map_elements(
            _perturbation, return_dtype=pl.String
        ).alias("perturbation"),
    )

    # Per-sim row keyed by dataset_id (resolve via variant -> dataset_id).
    # Use the *raw* dataset_id from the campaign config — NOT the rich
    # display label — so the join key matches the parca side.
    raw_dataset_id_map = {
        i: pv["rnaseq_basal_dataset_id"] for i, pv in enumerate(PARCA_VARIANTS)
    }
    sim_for_join = sim_df.with_columns(
        pl.col("variant").map_elements(
            lambda v: raw_dataset_id_map.get(int(v), "?"), return_dtype=pl.String
        ).alias("dataset_id"),
    ).rename({
        "exit_code": "sim_exit",
        "duration_s": "sim_dur_s",
        "error": "sim_error",
    }).select("dataset_id", "seed", "sim_exit", "sim_dur_s", "sim_error")

    combined = parca_for_join.join(sim_for_join, on="dataset_id", how="left").sort(
        "parca_id", "seed"
    ).select(
        "parca_id", "perturbation", "dataset_id",
        "parca_exit", "parca_dur_s", "parca_error",
        "seed", "sim_exit", "sim_dur_s", "sim_error",
    )

    # Build HTML rows. Color-code by:
    #  - red if parca failed (no sim ever ran)
    #  - red if sim failed
    #  - green only if both parca + sim succeeded
    #  - gray if parca succeeded but no sim row (shouldn't happen here)
    # Group seeds visually by rowspanning the parca-level columns across
    # all sim-seed rows belonging to the same parca run.
    parca_cols = ("parca_id", "perturbation", "dataset_id",
                  "parca_exit", "parca_dur_s", "parca_error")
    sim_cols = ("seed", "sim_exit", "sim_dur_s", "sim_error")
    header = "<tr>" + "".join(
        f"<th>{c}</th>" for c in (*parca_cols, *sim_cols)
    ) + "</tr>"

    # Group rows by parca_id so we can emit rowspan on the parca columns.
    rows_dicts = combined.to_dicts()
    by_parca: dict[int, list[dict]] = {}
    for r in rows_dicts:
        by_parca.setdefault(r["parca_id"], []).append(r)

    body_rows: list[str] = []
    for parca_id in sorted(by_parca.keys()):
        group = by_parca[parca_id]
        n = len(group)
        for i, r in enumerate(group):
            parca_ok = str(r["parca_exit"]) == "0"
            sim_ok = str(r.get("sim_exit", "")) == "0"
            if not parca_ok:
                css = "fail"
            elif r["sim_exit"] is None:
                css = "warn"
            elif sim_ok:
                css = "ok"
            else:
                css = "fail"
            cells = []
            if i == 0:
                # Emit parca-level cells with rowspan covering the whole group.
                for c in parca_cols:
                    v = r[c]
                    cells.append(
                        f'<td rowspan="{n}" class="grouped">'
                        f'{"" if v is None else v}</td>'
                    )
            for c in sim_cols:
                v = r[c]
                cells.append(f"<td>{'' if v is None else v}</td>")
            body_rows.append(f'<tr class="{css}">{"".join(cells)}</tr>')
    combined_html = f'<table class="tbl">{header}{"".join(body_rows)}</table>'

    n_parca_ok = parca_df.filter(pl.col("exit_code") == "0").height
    n_sim_ok = sim_df.filter(pl.col("exit_code") == "0").height
    n_parca = parca_df.height
    n_sim = sim_df.height
    # If every ParCa had succeeded, this many sim tasks would have run.
    n_sim_max = n_parca * SEEDS_PER_VARIANT * GENS_PER_VARIANT

    body = f"""
<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Campaign summary — {EXP_ID}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 1400px; margin: 2em auto; padding: 0 1em; color: #222; }}
  h1, h2 {{ border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }}
  table.tbl {{ border-collapse: collapse; font-size: 13px; margin: 0.5em 0; width: 100%; }}
  table.tbl th, table.tbl td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: left; vertical-align: top; }}
  table.tbl th {{ background: #f5f5f5; }}
  table.tbl tr.ok td   {{ background: #e8ffe8; }}
  table.tbl tr.fail td {{ background: #ffe8e8; }}
  table.tbl tr.warn td {{ background: #f5f5f5; }}
  table.tbl td.grouped {{ background: #fafafa; border-right: 2px solid #888; }}
  .stat {{ display: inline-block; padding: 0.2em 0.7em; margin-right: 0.5em; background: #eef; border-radius: 4px; }}
</style>
</head><body>

<h1>Campaign summary — {EXP_ID}</h1>
<p>
  <span class="stat">ParCa: {n_parca_ok}/{n_parca} succeeded</span>
  <span class="stat">Sim: {n_sim_ok}/{n_sim} of attempted succeeded
    &nbsp;<small>({n_sim} were attempted out of {n_sim_max} planned —
    {n_parca - n_parca_ok} variants skipped because their ParCa failed;
    plan was {n_parca} variants × {SEEDS_PER_VARIANT} seeds × {GENS_PER_VARIANT} gens)</small>
  </span>
</p>

<h2>Task status (ParCa + Sim)</h2>
<p>One row per Nextflow task. ParCa rows are model-fitting per RNA-seq input;
sim rows are cell-cycle simulations on a surviving variant. Green = exit 0,
red = failed. Duration is wall-clock seconds.</p>
{combined_html}

<h2>Cell mass over time</h2>
<div id="mass"></div>

<h2>Mass components over time (median across seeds)</h2>
<div id="components"></div>

<h2>Growth rate distribution</h2>
<div id="gr"></div>

<h2>Amino-acid synthesis enzyme abundance — log2 fold change vs baseline</h2>
<p>Per amino acid (rows): the model-summed protein count of all forward-synthesis
enzymes for that AA, averaged across all (seed × time) for the variant, then
log2-folded against variant 0 (the unperturbed baseline). Cells must synthesize
all 20 AAs in M9 minus AAs, so any cross-variant signal here is the perturbation
biting. ★-marked rows are AAs whose synthesis enzymes have explicit
RNA-expression overrides in <code>rna_expression_adjustments.tsv</code>.</p>
<div id="aa"></div>

<script>
  Plotly.newPlot('mass', {fig_mass.to_json()});
  Plotly.newPlot('components', {fig_components.to_json()});
  Plotly.newPlot('gr', {fig_gr.to_json()});
  {f"Plotly.newPlot('aa', {fig_aa.to_json()});" if fig_aa else "/* no AA fold-change data */"}
</script>
</body></html>
"""
    return body


def main() -> int:
    if not WORKDIRS.exists():
        print(f"ERROR: workdirs not found at {WORKDIRS}", file=sys.stderr)
        return 1
    rows = _scan_workdirs()
    parca_df, sim_df = _build_status(rows)
    mass = _query_mass_data()
    aa_fc, aa_highlight = _aa_synthesis_fold_change()

    # Save TSVs
    parca_df.write_csv(OUTDIR / "parca_status.tsv", separator="\t")
    sim_df.write_csv(OUTDIR / "sim_status.tsv", separator="\t")
    if aa_fc.height:
        aa_fc.write_csv(OUTDIR / "aa_synthesis_fold_change.tsv", separator="\t")
    print(f"wrote {OUTDIR / 'parca_status.tsv'}")
    print(f"wrote {OUTDIR / 'sim_status.tsv'}")
    print(f"wrote {OUTDIR / 'aa_synthesis_fold_change.tsv'} ({aa_fc.height} rows)")

    html = _make_html(parca_df, sim_df, mass, aa_fc, aa_highlight)
    out_html = OUTDIR / "index.html"
    out_html.write_text(html)
    print(f"wrote {out_html}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
