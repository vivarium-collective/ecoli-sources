"""
Perturbation operators for RNA-seq TPM datasets.

Each operator is a pure function that maps a TPM DataFrame
(``gene_id`` / ``tpm_mean`` / optional ``tpm_std``) to a perturbed TPM
DataFrame of the same shape. Operators are deterministic given their RNG
seed and parameters.

The ``make_perturbation_variant`` driver applies an operator to a source
dataset, writes the resulting TSV, and appends a provenance row to the
samples manifest. The ``generate_campaign`` helper iterates the driver
over a parameter grid for bulk dataset generation.

This is the Group 1 (expression profile) operator library from
``.claude/plans/dataset-sensitivity-exploration.md`` Part 3. Operators for
Groups 2-5 (physiology, kinetics, turnover, regulation) will land in
separate modules as those groups come online.
"""

from __future__ import annotations

import hashlib
import json
import os
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from schemas import RnaseqSamplesManifestSchema, RnaseqTpmTableSchema


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _drop_std(df: pd.DataFrame) -> pd.DataFrame:
    """Drop tpm_std if present — it's not meaningful after a mean perturbation."""
    if "tpm_std" in df.columns:
        return df.drop(columns="tpm_std")
    return df


def _renormalize_to_million(df: pd.DataFrame) -> pd.DataFrame:
    """Rescale tpm_mean so it sums to 1e6 (the TPM convention)."""
    total = df["tpm_mean"].sum()
    if total <= 0:
        return df
    df = df.copy()
    df["tpm_mean"] = df["tpm_mean"] * 1e6 / total
    return df


def _hash_params(params: dict) -> str:
    """8-char hash of a params dict for filenames and ids."""
    blob = json.dumps(params, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:8]


# ---------------------------------------------------------------------------
# Operators (Group 1: expression profile)
# ---------------------------------------------------------------------------


def add_log_normal_noise(
    tpm: pd.DataFrame,
    sigma: float,
    seed: int,
    renormalize: bool = False,
) -> pd.DataFrame:
    """
    Multiply each ``tpm_mean`` by ``exp(N(0, sigma))``.

    Log-normal noise is the natural scalar-invariant perturbation for
    TPM: it preserves positivity, respects order-of-magnitude structure,
    and is symmetric in log space. ``sigma == 0`` is identity.

    Args:
        sigma: log-space std dev (``sigma=0.2`` ≈ ±20% multiplicative).
        seed: RNG seed.
        renormalize: rescale so tpm_mean sums to 1e6 after perturbation.
    """
    rng = np.random.default_rng(seed)
    out = _drop_std(tpm).copy()
    factors = np.exp(rng.normal(0.0, sigma, size=len(out)))
    out["tpm_mean"] = out["tpm_mean"].astype(float) * factors
    return _renormalize_to_million(out) if renormalize else out


def scale_gene_set(
    tpm: pd.DataFrame,
    gene_ids: Iterable[str],
    factor: float,
    renormalize: bool = False,
) -> pd.DataFrame:
    """
    Multiply ``tpm_mean`` for a named set of genes by ``factor``.

    Inverse of the manual `rna_expression_adjustments` mechanism: tests
    whether parca can absorb explicit subsystem-level rescaling.
    """
    gene_set = set(gene_ids)
    out = _drop_std(tpm).copy()
    mask = out["gene_id"].isin(gene_set)
    out.loc[mask, "tpm_mean"] = out.loc[mask, "tpm_mean"].astype(float) * factor
    return _renormalize_to_million(out) if renormalize else out


def zero_genes(
    tpm: pd.DataFrame,
    gene_ids: Iterable[str],
) -> pd.DataFrame:
    """
    Set ``tpm_mean`` to 0 for the given genes (full knockout).

    Useful for single-gene sensitivity: which genes, when silenced, break
    the model? Notably overlaps with the adjustment-flagged gene set.
    Unlike :func:`make_gene_exclusion_variant`, rows are retained (value
    set to 0) so downstream joins don't drop them.
    """
    gene_set = set(gene_ids)
    out = _drop_std(tpm).copy()
    out.loc[out["gene_id"].isin(gene_set), "tpm_mean"] = 0.0
    return out


def interpolate_datasets(
    tpm_a: pd.DataFrame,
    tpm_b: pd.DataFrame,
    alpha: float,
) -> pd.DataFrame:
    """
    Linear interpolation: ``(1 - alpha) * tpm_a + alpha * tpm_b``.

    Restricted to genes present in both inputs. ``alpha=0`` returns
    tpm_a (intersection-of-genes); ``alpha=1`` returns tpm_b. Use to
    walk a continuous path between a working and a failing dataset and
    locate the break point.
    """
    if not 0.0 <= alpha <= 1.0:
        raise ValueError(f"alpha must be in [0, 1]; got {alpha}")
    merged = (
        _drop_std(tpm_a)
        .rename(columns={"tpm_mean": "tpm_a"})
        .merge(
            _drop_std(tpm_b).rename(columns={"tpm_mean": "tpm_b"}),
            on="gene_id",
            how="inner",
        )
    )
    merged["tpm_mean"] = (1 - alpha) * merged["tpm_a"] + alpha * merged["tpm_b"]
    return merged[["gene_id", "tpm_mean"]]


def quantile_match(
    tpm_source: pd.DataFrame,
    tpm_target: pd.DataFrame,
) -> pd.DataFrame:
    """
    Rank-preserving rescale of ``tpm_source`` to the distribution of
    ``tpm_target``.

    For each gene in the source, its rank is preserved but its value is
    replaced by the target value at the same rank. Separates
    "distribution shape matters" from "gene identity matters": if
    quantile-matching a failing dataset onto a working one rescues it,
    the dynamic range / outlier structure was to blame (not which genes
    are high).

    Restricted to genes present in both inputs, truncated to the min of
    the two lengths.
    """
    common = set(tpm_source["gene_id"]) & set(tpm_target["gene_id"])
    src = _drop_std(tpm_source).loc[tpm_source["gene_id"].isin(common)].copy()
    tgt = _drop_std(tpm_target).loc[tpm_target["gene_id"].isin(common)].copy()

    src = src.sort_values("tpm_mean").reset_index(drop=True)
    tgt_sorted_values = np.sort(tgt["tpm_mean"].astype(float).values)

    n = min(len(src), len(tgt_sorted_values))
    src = src.iloc[:n].copy()
    src["tpm_mean"] = tgt_sorted_values[:n]
    return src[["gene_id", "tpm_mean"]].reset_index(drop=True)


def exclude_genes(
    tpm: pd.DataFrame,
    gene_ids: Iterable[str],
) -> pd.DataFrame:
    """
    Remove rows for the given genes entirely (Annabelle's gene-exclusion
    approach from ``post_processing.make_gene_exclusion_variant``).

    Unlike :func:`zero_genes` which keeps rows at ``tpm_mean=0``, this
    **drops** them from the table. The difference matters downstream: parca's
    ``rnaseq_fill_missing_genes_from_ref`` fills absent genes from the
    reference, so excluded genes get reference-level expression rather than
    zero. Use this to test whether specific genes are load-bearing for the
    parca fit.
    """
    gene_set = set(gene_ids)
    out = _drop_std(tpm).copy()
    return out[~out["gene_id"].isin(gene_set)].reset_index(drop=True)


def drop_and_fill(
    tpm: pd.DataFrame,
    fraction: float,
    seed: int,
    fill_value: float = 0.0,
) -> pd.DataFrame:
    """
    Randomly replace the ``tpm_mean`` of ``fraction`` genes with
    ``fill_value``.

    Stress test for the ``rnaseq_fill_missing_genes_from_ref`` code path.
    With ``fill_value=0``, the parca's fill-from-ref logic (if enabled)
    will recover values from the reference dataset; we measure how much
    coverage loss parca can tolerate before the P-solve or downstream
    fits fail.
    """
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"fraction must be in [0, 1]; got {fraction}")
    rng = np.random.default_rng(seed)
    out = _drop_std(tpm).copy()
    n_drop = int(round(fraction * len(out)))
    if n_drop == 0:
        return out
    drop_idx = rng.choice(len(out), size=n_drop, replace=False)
    out.iloc[drop_idx, out.columns.get_loc("tpm_mean")] = fill_value
    return out


# ---------------------------------------------------------------------------
# Dispatcher + driver
# ---------------------------------------------------------------------------


#: Operator registry. Unary operators (one source dataset) only; binary
#: operators are called through :func:`make_binary_perturbation_variant`.
UNARY_OPERATORS: dict[str, Callable[..., pd.DataFrame]] = {
    "add_log_normal_noise": add_log_normal_noise,
    "scale_gene_set": scale_gene_set,
    "zero_genes": zero_genes,
    "exclude_genes": exclude_genes,
    "drop_and_fill": drop_and_fill,
}

BINARY_OPERATORS: dict[str, Callable[..., pd.DataFrame]] = {
    "interpolate_datasets": interpolate_datasets,
    "quantile_match": quantile_match,
}


def _load_source_tpm(manifest: pd.DataFrame, dataset_id: str, data_dir: str) -> pd.DataFrame:
    row = manifest.loc[manifest["dataset_id"] == dataset_id]
    if row.empty:
        raise KeyError(f"dataset_id {dataset_id!r} not found in manifest")
    path = os.path.join(data_dir, row.iloc[0]["file_path"])
    return pd.read_csv(path, sep="\t")


def _append_manifest_row(
    manifest: pd.DataFrame,
    source_row: pd.Series,
    new_dataset_id: str,
    new_file_path: str,
    new_description: str,
    operator: str,
    operator_params: dict,
    seed: int | None,
) -> pd.DataFrame:
    new_row = source_row.to_dict()
    new_row["dataset_id"] = new_dataset_id
    new_row["file_path"] = new_file_path
    new_row["dataset_description"] = new_description
    new_row["parent_dataset_id"] = source_row["dataset_id"]
    new_row["operator"] = operator
    # Drop seed from params JSON — it gets its own column.
    params_for_json = {k: v for k, v in operator_params.items() if k != "seed"}
    new_row["operator_params_json"] = json.dumps(params_for_json, sort_keys=True, default=str)
    new_row["seed"] = seed
    return pd.concat([manifest, pd.DataFrame([new_row])], ignore_index=True)


def make_perturbation_variant(
    manifest: pd.DataFrame,
    source_dataset_id: str,
    operator: str,
    operator_params: dict,
    data_dir: str,
    new_dataset_id: str | None = None,
    subdir: str = "perturbations",
    manifest_filename: str = "manifest.tsv",
    description_prefix: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """
    Apply a unary operator to a source dataset, write the perturbed TSV,
    and append a provenance-tagged row to the manifest.

    Returns the updated manifest and the new ``dataset_id``.

    If ``new_dataset_id`` is ``None``, it is derived from the source id,
    operator name, and a short hash of the operator parameters.
    """
    if operator not in UNARY_OPERATORS:
        raise KeyError(
            f"Unknown unary operator {operator!r}; "
            f"known: {sorted(UNARY_OPERATORS)}"
        )
    fn = UNARY_OPERATORS[operator]

    if source_dataset_id not in set(manifest["dataset_id"]):
        raise ValueError(f"source_dataset_id {source_dataset_id!r} not in manifest")

    source_row = manifest.loc[manifest["dataset_id"] == source_dataset_id].iloc[0]
    source_tpm = _load_source_tpm(manifest, source_dataset_id, data_dir)

    perturbed = fn(source_tpm, **operator_params)

    # Validate the perturbed table before writing.
    perturbed = RnaseqTpmTableSchema.validate(perturbed)

    seed = operator_params.get("seed")

    if new_dataset_id is None:
        params_hash = _hash_params(operator_params)
        new_dataset_id = f"{source_dataset_id}__{operator}__{params_hash}"
    if new_dataset_id in set(manifest["dataset_id"]):
        raise ValueError(f"new_dataset_id {new_dataset_id!r} already in manifest")

    # Write the TSV
    rel_dir = os.path.join(subdir, operator, source_dataset_id)
    abs_dir = os.path.join(data_dir, rel_dir)
    os.makedirs(abs_dir, exist_ok=True)
    rel_path = os.path.join(rel_dir, f"{new_dataset_id}.tsv")
    abs_path = os.path.join(data_dir, rel_path)
    perturbed.to_csv(abs_path, sep="\t", index=False)

    description = source_row["dataset_description"]
    if description_prefix is None:
        description_prefix = f"[{operator}] "
    new_description = f"{description_prefix}{description}"

    updated_manifest = _append_manifest_row(
        manifest, source_row,
        new_dataset_id=new_dataset_id,
        new_file_path=rel_path,
        new_description=new_description,
        operator=operator,
        operator_params=operator_params,
        seed=seed,
    )

    RnaseqSamplesManifestSchema.validate(updated_manifest)
    updated_manifest.to_csv(
        os.path.join(data_dir, manifest_filename), sep="\t", index=False
    )

    return updated_manifest, new_dataset_id


def make_binary_perturbation_variant(
    manifest: pd.DataFrame,
    source_dataset_ids: tuple[str, str],
    operator: str,
    operator_params: dict,
    data_dir: str,
    new_dataset_id: str | None = None,
    subdir: str = "perturbations",
    manifest_filename: str = "manifest.tsv",
    description_prefix: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """
    Variant of :func:`make_perturbation_variant` for binary operators
    (``interpolate_datasets``, ``quantile_match``).

    ``source_dataset_ids`` is an ordered pair (``source_a``, ``source_b``).
    The provenance row records ``parent_dataset_id = source_a`` and
    encodes ``source_b`` in ``operator_params_json``.
    """
    if operator not in BINARY_OPERATORS:
        raise KeyError(
            f"Unknown binary operator {operator!r}; "
            f"known: {sorted(BINARY_OPERATORS)}"
        )
    fn = BINARY_OPERATORS[operator]
    a_id, b_id = source_dataset_ids

    for did in source_dataset_ids:
        if did not in set(manifest["dataset_id"]):
            raise ValueError(f"source dataset {did!r} not in manifest")

    tpm_a = _load_source_tpm(manifest, a_id, data_dir)
    tpm_b = _load_source_tpm(manifest, b_id, data_dir)
    perturbed = RnaseqTpmTableSchema.validate(fn(tpm_a, tpm_b, **operator_params))

    params_for_record = {**operator_params, "source_b": b_id}
    if new_dataset_id is None:
        params_hash = _hash_params(params_for_record)
        new_dataset_id = f"{a_id}__{operator}__{b_id}__{params_hash}"
    if new_dataset_id in set(manifest["dataset_id"]):
        raise ValueError(f"new_dataset_id {new_dataset_id!r} already in manifest")

    rel_dir = os.path.join(subdir, operator, a_id)
    abs_dir = os.path.join(data_dir, rel_dir)
    os.makedirs(abs_dir, exist_ok=True)
    rel_path = os.path.join(rel_dir, f"{new_dataset_id}.tsv")
    perturbed.to_csv(os.path.join(data_dir, rel_path), sep="\t", index=False)

    source_row = manifest.loc[manifest["dataset_id"] == a_id].iloc[0]
    description = source_row["dataset_description"]
    if description_prefix is None:
        description_prefix = f"[{operator}({b_id})] "
    new_description = f"{description_prefix}{description}"

    updated_manifest = _append_manifest_row(
        manifest, source_row,
        new_dataset_id=new_dataset_id,
        new_file_path=rel_path,
        new_description=new_description,
        operator=operator,
        operator_params=params_for_record,
        seed=operator_params.get("seed"),
    )

    RnaseqSamplesManifestSchema.validate(updated_manifest)
    updated_manifest.to_csv(
        os.path.join(data_dir, manifest_filename), sep="\t", index=False
    )

    return updated_manifest, new_dataset_id


def generate_campaign(
    manifest: pd.DataFrame,
    source_dataset_id: str,
    operator: str,
    param_grid: list[dict],
    data_dir: str,
    **kwargs,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Run :func:`make_perturbation_variant` over a list of parameter dicts
    and return the fully-updated manifest plus the generated dataset ids.
    """
    ids = []
    current = manifest
    for params in param_grid:
        current, did = make_perturbation_variant(
            current, source_dataset_id, operator, params, data_dir, **kwargs
        )
        ids.append(did)
    return current, ids
