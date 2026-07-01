"""
Overlay-union loader for the validation-data subsystem.

Mirrors the ParCa overlay pattern (vEcoli ``wholecell/io/ingestion.py``'s
``ingest_rnaseq_manifest_with_overlays``, unique on ``dataset_id``) for the
*validation* flow: union a primary validation bundle with any overlay bundles,
unique on ``canonical_key``. The reusable loader lives here; overlay repos ship
only their own ``validation_data/`` and are unioned in at load time.

Two claim-table conventions are resolved transparently to one in-memory band:

* **per-observable file** (public ``ScalarClaimSchema``): one file per
  ``(condition, observable)``, e.g. ``data/basal/o2_uptake.tsv``; the measured
  band is the ``kind == measured`` rows.
* **tidy assay table** (private ``ScalarObservationSchema``): one table holding
  many observables across cultivations, e.g. ``fermentation.tsv``; the bundle row
  carries a ``(cultivation_id, observable)`` filter selecting the slot's rows.

Overlay bundle paths come from ``overlay_paths=`` or the
``$ECOLI_SOURCES_VALIDATION_OVERLAYS`` environment variable (``;`` preferred,
``:`` legacy; cloud ``scheme://`` URIs preserved) — the validation analogue of
``$ECOLI_SOURCES_OVERLAYS``.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from ecoli_sources import VALIDATION_BUNDLE_PATH

# Scalar claim-table schemas this loader resolves to a measured band. Vector
# schemas (ReactionFlux, ProteinAbundance, …) are handled by their own readers.
_SCALAR_SCHEMAS = {"ScalarClaimSchema", "ScalarObservationSchema"}

VALIDATION_OVERLAYS_ENV = "ECOLI_SOURCES_VALIDATION_OVERLAYS"


# ---------------------------------------------------------------------------
# Overlay path parsing (mirrors vEcoli ingestion._split_overlay_separator)
# ---------------------------------------------------------------------------

def _split_overlay_separator(raw: str) -> list[str]:
    """Split a separator-delimited overlay list. ``;`` is preferred (URI-safe);
    ``:`` is the legacy separator, with ``scheme://`` fragments reattached so a
    colon-split does not tear ``s3://bucket/...`` apart."""
    raw = raw.strip()
    if not raw:
        return []
    if ";" in raw:
        return [p.strip() for p in raw.split(";") if p.strip()]
    parts = raw.split(":")
    out: list[str] = []
    for part in parts:
        if out and out[-1].endswith("/") and part.startswith("/"):
            # reattach a <scheme>://… that the colon-split tore apart
            out[-1] = out[-1] + ":" + part
        elif out and out[-1] in ("s3", "gs", "gcs", "http", "https", "file"):
            out[-1] = out[-1] + ":" + part
        else:
            out.append(part)
    return [p.strip() for p in out if p.strip()]


def validation_overlay_paths(env: str = VALIDATION_OVERLAYS_ENV) -> list[Path]:
    """Overlay validation-bundle paths from ``$ECOLI_SOURCES_VALIDATION_OVERLAYS``."""
    return [Path(p) for p in _split_overlay_separator(os.environ.get(env, ""))]


# ---------------------------------------------------------------------------
# Band extraction from one bundle row (both conventions)
# ---------------------------------------------------------------------------

def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def _measured_band(df: pd.DataFrame, *, source_col: str) -> dict:
    """From a claim/observation frame (already filtered to one slot), pull the
    measured band + optional theoretical-max ceiling + provenance."""
    meas = df[df["kind"] == "measured"]
    tmax = df[df["kind"] == "theoretical_max"] if "kind" in df else df.iloc[0:0]
    tmax_row = tmax.loc[tmax["value"].idxmin()] if len(tmax) else None
    unc = (meas["uncertainty"] if "uncertainty" in meas
           else pd.Series([None] * len(meas)))
    return {
        "measured": [float(v) for v in meas["value"]],
        "measured_unc": [(float(u) if pd.notna(u) else None) for u in unc],
        "theoretical_max": (float(tmax_row["value"]) if tmax_row is not None else None),
        "theoretical_source": (str(tmax_row[source_col]) if tmax_row is not None
                               and source_col in tmax_row else None),
        "sources": [str(s) for s in meas.get(source_col, [])],
        "units": (str(meas["units"].iloc[0]) if len(meas) and "units" in meas else None),
    }


def _row_band(row: dict, root: Path) -> dict:
    """Resolve one validation-bundle row to a measured band, handling both the
    per-observable-file and tidy-assay-table conventions."""
    table = _read_tsv(root / str(row["source_path"]))
    cultivation_id = row.get("cultivation_id")
    observable = row.get("observable")
    is_tidy = (pd.notna(cultivation_id) if cultivation_id is not None else False) and (
        "observable" in table.columns)
    if is_tidy:
        sub = table[table["cultivation_id"].astype(str) == str(cultivation_id)]
        if pd.notna(observable):
            sub = sub[sub["observable"].astype(str) == str(observable)]
        band = _measured_band(sub, source_col="replicate")
        band["cultivation_id"] = str(cultivation_id)
        band["observable"] = (str(observable) if pd.notna(observable) else None)
        # derived-tier temporal provenance (constant across replicates of a slot)
        for col in ("phase", "window"):
            vals = ([v for v in sub[col].dropna().unique()] if col in sub.columns else [])
            band[col] = (str(vals[0]) if vals else None)
    else:
        band = _measured_band(table, source_col="source_id")
        band["cultivation_id"] = None
        # observable inferred from the canonical_key tail (``<condition>__<obs>``)
        key = str(row["canonical_key"])
        band["observable"] = key.split("__", 1)[1] if "__" in key else key
    band["canonical_key"] = str(row["canonical_key"])
    band["schema_name"] = str(row.get("schema_name") or "")
    band.setdefault("phase", None)    # per-observable convention carries no temporal tier
    band.setdefault("window", None)
    return band


# ---------------------------------------------------------------------------
# Union loader
# ---------------------------------------------------------------------------

def load_scalar_observations(
    primary_bundle: Path | str | None = None,
    overlay_paths: list[Path | str] | None = None,
    *,
    include_primary: bool = True,
    prefix: str | None = None,
) -> dict[str, dict]:
    """Union scalar validation slots across the primary + overlay bundles.

    Returns ``{canonical_key: band}`` where ``band`` carries ``measured`` (the
    band), ``measured_unc``, ``theoretical_max``/``theoretical_source``,
    ``sources``, ``units``, ``observable``, ``cultivation_id``, ``schema_name``.
    Only scalar slots (``ScalarClaimSchema`` / ``ScalarObservationSchema``) are
    returned; vector slots are skipped. ``canonical_key`` must be unique across
    the union — a collision raises ``ValueError`` (mirrors the ParCa overlay
    ``dataset_id`` contract). ``prefix`` optionally filters to keys starting with
    it (e.g. a ``"<cultivation_id>__"`` prefix)."""
    if primary_bundle is None:
        primary_bundle = VALIDATION_BUNDLE_PATH
    if overlay_paths is None:
        overlay_paths = validation_overlay_paths()

    bundles: list[Path] = []
    if include_primary:
        bundles.append(Path(primary_bundle))
    bundles.extend(Path(p) for p in overlay_paths)

    out: dict[str, dict] = {}
    for bundle_path in bundles:
        manifest = _read_tsv(bundle_path)
        root = bundle_path.parent
        for _, r in manifest.iterrows():
            row = r.to_dict()
            schema = str(row.get("schema_name") or "")
            if schema not in _SCALAR_SCHEMAS:
                continue
            key = str(row["canonical_key"])
            if prefix is not None and not key.startswith(prefix):
                continue
            if key in out:
                raise ValueError(
                    f"duplicate validation canonical_key '{key}' across the "
                    f"primary+overlay union (second seen in {bundle_path})")
            out[key] = _row_band(row, root)
    return out


def observations_for_cultivation(
    cultivation_id: str,
    primary_bundle: Path | str | None = None,
    overlay_paths: list[Path | str] | None = None,
    *,
    include_primary: bool = True,
) -> dict[str, dict]:
    """Scalar observations for one cultivation, keyed by ``observable``.

    Convenience view for a ``vs_experiment`` card: returns
    ``{observable: band}`` for the given ``cultivation_id`` (the shape the basal
    ``vs_literature`` render consumes, sourced from a cultivation instead of
    literature)."""
    allk = load_scalar_observations(
        primary_bundle, overlay_paths, include_primary=include_primary)
    return {b["observable"]: b for b in allk.values()
            if b.get("cultivation_id") == cultivation_id and b.get("observable")}


def load_cultivations(
    primary_registry: Path | str | None = None,
    overlay_registries: list[Path | str] | None = None,
) -> dict[str, dict]:
    """Union the cultivation registry across primary + overlays, keyed by
    ``cultivation_id`` (unique across the union or ``ValueError``). Registry
    paths default to ``cultivations.tsv`` beside each validation bundle."""
    registries: list[Path] = []
    if primary_registry is not None:
        registries.append(Path(primary_registry))
    elif (VALIDATION_BUNDLE_PATH.parent / "cultivations.tsv").exists():
        registries.append(VALIDATION_BUNDLE_PATH.parent / "cultivations.tsv")
    if overlay_registries is None:
        overlay_registries = [p.parent / "cultivations.tsv"
                              for p in validation_overlay_paths()]
    registries.extend(Path(p) for p in overlay_registries if Path(p).exists())

    out: dict[str, dict] = {}
    for reg in registries:
        for _, r in _read_tsv(reg).iterrows():
            cid = str(r["cultivation_id"])
            if cid in out:
                raise ValueError(f"duplicate cultivation_id '{cid}' across registries ({reg})")
            out[cid] = {k: (None if pd.isna(v) else v) for k, v in r.to_dict().items()}
    return out
