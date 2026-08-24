"""
Pandera schemas for experimental RNA-seq data.

Canonical format: one file per condition, two required columns (gene_id, tpm_mean)
plus optional tpm_std. Sample metadata is stored in a separate manifest so that
condition semantics (strain, media, is_basal, etc.) are not encoded in column
headers.
"""

import pandera.pandas as pa


# ---------------------------------------------------------------------------
# TPM table: one file per condition (gene_id, tpm_mean [, tpm_std])
# ---------------------------------------------------------------------------

RnaseqTpmTableSchema = pa.DataFrameSchema(
    name="rnaseq_tpm_table",
    columns={
        "gene_id": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description="Gene identifier; must match reference gene set (e.g. EcoCyc id like EG10001).",
        ),
        "tpm_mean": pa.Column(
            float,
            nullable=False,
            checks=[
                pa.Check.greater_than_or_equal_to(0),
                # TODO: add check to ensure tpm_mean sums to 1 million (perhaps just give a warning if not)
            ],
            description="Mean TPM (transcripts per million) for this gene in this condition.",
        ),
        "tpm_std": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=[
                pa.Check.greater_than_or_equal_to(0),
            ],
            description="Optional: standard deviation of TPM across replicates.",
        ),
    },
    strict="filter",  # allow extra columns but validate required ones
    coerce=True,
    description=(
        "RNA-seq TPM table. One file per sample/condition. "
        "Required columns: gene_id, tpm_mean. Optional: tpm_std."
    ),
)


# ---------------------------------------------------------------------------
# Sample manifest: maps sample_id to file path and metadata
# ---------------------------------------------------------------------------

RnaseqSamplesManifestSchema = pa.DataFrameSchema(
    name="rnaseq_samples_manifest",
    columns={
        "dataset_id": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description="Unique identifier for this dataset, corresponding to the file name of the TPM table.",
        ),
        "dataset_description": pa.Column(
            dtype=str,
            nullable=False,
            description="Description of the dataset (e.g. 'exp96546: MG1655 in M9 glucose, average of 3h and 4h timepoints').",
        ),
        "file_path": pa.Column(
            dtype=str,
            nullable=False,
            description="Path to the TPM table file (relative to manifest or absolute).",
        ),
        "data_source": pa.Column(
            dtype=str,
            nullable=False,
            description="Source of the data (e.g. 'Ginkgo', 'PNNL').",
        ),
        "data_source_experiment_id": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Experiment identifier from the data source (e.g. 'exp96546').",
        ),
        "data_source_date": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Date of the experiment from the data source (e.g. '2026-01-01').",
        ),
        "strain": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Optional: strain descriptor (e.g. 'MG1655 rph+').",
        ),
        "condition": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Optional: cultivation condition descriptor (e.g. 'Modified_M9_N_Fe').",
        ),
        # --- Provenance columns for synthetic / perturbed datasets ---
        "parent_dataset_id": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional: for synthetic datasets, the dataset_id this one "
                "was derived from. Null for primary (non-derived) datasets."
            ),
        ),
        "operator": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional: name of the perturbation operator that produced "
                "this dataset (e.g. 'add_log_normal_noise')."
            ),
        ),
        "operator_params_json": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional: JSON-serialized dict of operator parameters for "
                "reproducibility (excluding the RNG seed, which has its own column)."
            ),
        ),
        "seed": pa.Column(
            dtype="Int64",  # nullable integer
            nullable=True,
            required=False,
            description="Optional: RNG seed used by the operator (if stochastic).",
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Manifest of RNA-seq datasets. Primary columns: dataset_id, "
        "dataset_description, file_path, data_source (+ optional data_source_* "
        "/ strain / condition metadata). Provenance columns "
        "(parent_dataset_id, operator, operator_params_json, seed) are used "
        "for synthetic / perturbed datasets; null for primary sources."
    ),
)


# ---------------------------------------------------------------------------
# Differential-expression result table (one file per contrast)
# ---------------------------------------------------------------------------
#
# Promoted alongside the TPM tables an overlay derives FROM them, so a consumer
# can ask which genes moved and by how much, rather than only seeing the
# post-overlay expression. The table is deliberately LOSSLESS: every tested
# transcript is kept, with no significance filter applied, because the
# thresholds are an analysis choice belonging to the consumer. Use
# ``select_significant_genes`` to apply them — see the warning there about
# operation order.

DeseqResultTableSchema = pa.DataFrameSchema(
    name="deseq_result_table",
    columns={
        "transcript_id": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description=(
                "Transcript identifier as tested, one row per transcript. This is the "
                "table's grain: several transcripts may map to one gene."
            ),
        ),
        "gene_id": pa.Column(
            dtype=str,
            nullable=True,
            description=(
                "Gene identifier the transcript maps to (e.g. an EcoCyc id like EG10001), "
                "matching the reference gene set used by the TPM tables. NULLABLE and "
                "NOT unique: a heterologous transcript has no reference gene, and several "
                "transcripts can share one gene."
            ),
        ),
        "log2_fold_change": pa.Column(
            float,
            nullable=True,
            description=(
                "log2 fold change of the contrast, in the direction recorded in the "
                "dataset's provenance. Nullable: a transcript can be tested and yield no "
                "estimate."
            ),
        ),
        "padj": pa.Column(
            float,
            nullable=True,
            checks=[
                pa.Check.in_range(0.0, 1.0),
            ],
            description=(
                "Multiple-testing-adjusted p-value. Nullable: independent filtering "
                "leaves padj undefined for some transcripts, and that is meaningfully "
                "different from a non-significant value."
            ),
        ),
        "base_mean": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=[pa.Check.greater_than_or_equal_to(0)],
            description="Optional: mean normalized count across samples in the contrast.",
        ),
        "lfc_se": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=[pa.Check.greater_than_or_equal_to(0)],
            description="Optional: standard error of the log2 fold change.",
        ),
    },
    strict="filter",  # allow extra annotation columns; validate the contract
    coerce=True,
    description=(
        "Differential-expression results for ONE contrast, one row per tested "
        "transcript, unfiltered by significance. Required: transcript_id, gene_id, "
        "log2_fold_change, padj. The contrast's groups, sample composition and "
        "caller settings are NOT columns here — they belong to the dataset's "
        "provenance record, so that one file means one comparison."
    ),
)


def select_significant_genes(
    results,
    *,
    padj_thresh: float = 0.1,
    lfc_thresh: float = 1.0,
):
    """Canonical reading of a ``DeseqResultTableSchema`` table: the significant genes.

    Returns a gene-keyed frame, indexed by ``gene_id``, holding one row per gene.

    ⚠ **The order of operations is part of the contract, and reversing it silently
    changes the answer.** The steps are:

    1. drop rows with no ``padj``, no ``log2_fold_change``, or no ``gene_id`` —
       a transcript that was not testable, or that maps to no reference gene,
       is not a result about a gene;
    2. **then** apply the significance thresholds;
    3. **then** collapse transcripts to one row per ``gene_id``, keeping the
       transcript with the largest ``|log2_fold_change|``.

    Collapsing before filtering would pick a gene's most-extreme transcript and
    then test *that* for significance, which is a different question and yields a
    different gene set. Step 1 is likewise not a tidy-up: dropping rows without a
    gene id changes the count.

    This function exists so that consumers do not each re-derive those steps.
    Callers wanting a different definition should say so explicitly rather than
    reimplementing this one slightly differently.

    Args:
        results: a DataFrame conforming to ``DeseqResultTableSchema``.
        padj_thresh: keep rows with ``padj`` strictly below this.
        lfc_thresh: keep rows with ``|log2_fold_change|`` strictly above this.

    Returns:
        DataFrame indexed by ``gene_id``, one row per significant gene, carrying
        the columns of the selected transcript.
    """
    frame = results.dropna(subset=["padj", "log2_fold_change", "gene_id"])
    significant = frame[
        (frame["padj"] < padj_thresh)
        & (frame["log2_fold_change"].abs() > lfc_thresh)
    ].copy()
    significant["_abs_lfc"] = significant["log2_fold_change"].abs()
    collapsed = (
        significant.sort_values("_abs_lfc", ascending=False)
        .drop_duplicates("gene_id")
        .drop(columns="_abs_lfc")
        .set_index("gene_id")
    )
    return collapsed
