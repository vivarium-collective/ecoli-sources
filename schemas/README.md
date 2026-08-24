# Experimental data schemas

Pandera schemas for the ParCa data integration layer. These define the **canonical formats** for experimental inputs so that ingestion is validated and datasets can be substituted like-for-like.

## RNA-seq TPM tables

**Canonical format:** One file per condition. Two required columns:

| Column     | Type  | Required | Description |
|------------|--------|----------|-------------|
| `gene_id`  | string | yes      | Gene identifier; must match reference gene set (e.g. EcoCyc id like `EG10001`). Must be unique (one row per gene). |
| `tpm_mean` | float  | yes      | Mean TPM (transcripts per million) for this gene in this condition. Must be ≥ 0. |
| `tpm_std`  | float  | no       | Optional: standard deviation of TPM across replicates. Must be ≥ 0 if present. |

- File naming suggestion: `rnaseq_<dataset>_<sample_id>_tpm.tsv`
- Extra columns (e.g. gene symbol) are allowed; the schema uses `strict="filter"` and validates only the columns above.

### Validation

```python
import pandas as pd
from reconstruction.ecoli.experimental_data.schemas import RnaseqTpmTableSchema

df = pd.read_csv("path/to/rnaseq_exp96546_MG1655_M9_tpm.tsv", sep="\t")
validated = RnaseqTpmTableSchema.validate(df)
```

## RNA-seq sample manifest

Maps datasets to TPM table paths and metadata. One row per dataset. Used by ParCa config and for QC.

| Column                     | Type  | Required | Description |
|----------------------------|--------|----------|-------------|
| `dataset_id`               | string | yes      | Unique identifier for this dataset; typically matches the TPM table file name. |
| `dataset_description`      | string | yes      | Description of the dataset (e.g. "exp96546: MG1655 in M9 glucose, average of 3h and 4h timepoints"). |
| `file_path`                | string | yes      | Path to the TPM table file (relative to manifest or absolute). |
| `data_source`              | string | yes      | Source of the data (e.g. "Ginkgo", "PNNL"). |
| `data_source_experiment_id`| string | no       | Experiment identifier from the data source (e.g. "exp96546"). |
| `data_source_date`         | string | no       | Date of the experiment from the data source (e.g. "2026-01-01"). |
| `strain`                   | string | no       | Strain descriptor (e.g. "MG1655 rph+"). |
| `condition`                | string | no       | Cultivation condition descriptor (e.g. "Modified_M9_N_Fe"). |

Extra columns are allowed (`strict="filter"`).

### Validation

```python
import pandas as pd
from reconstruction.ecoli.experimental_data.schemas import RnaseqSamplesManifestSchema

manifest = pd.read_csv("path/to/rnaseq_samples.tsv", sep="\t")
validated = RnaseqSamplesManifestSchema.validate(manifest)
```

## Other vEcoli input tables

In addition to RNA-seq, this package validates the parca-time parameter and
adjustment tables that will migrate out of `omics-vEcoli/reconstruction/ecoli/flat/`.
Grouped by file:

| Schema | Source TSV(s) |
|---|---|
| `AdjustmentValueSchema` | `flat/adjustments/{rna_expression,rna_deg_rates,translation_efficiencies,protein_deg_rates}_adjustments.tsv` |
| `AminoAcidPathwayAdjustmentSchema` | `flat/adjustments/amino_acid_pathways.tsv` |
| `BalancedTranslationEfficiencyGroupSchema` | `flat/adjustments/balanced_translation_efficiencies.tsv` |
| `RelativeMetaboliteConcentrationChangeSchema` | `flat/adjustments/relative_metabolite_concentrations_changes.tsv` |
| `GrowthRateDependentParametersSchema` | `flat/growth_rate_dependent_parameters.tsv` |
| `DryMassCompositionSchema` | `flat/dry_mass_composition.tsv` |
| `RnaHalfLivesSchema` | `flat/rna_half_lives.tsv` |
| `ProteinHalfLivesMeasuredSchema` | `flat/protein_half_lives_measured.tsv` |
| `ProteinHalfLivesPulsedSilacSchema` | `flat/protein_half_lives_pulsed_silac.tsv` |
| `ProteinHalfLivesNEndRuleSchema` | `flat/protein_half_lives_n_end_rule.tsv` |
| `TranslationEfficiencySchema` | `flat/translation_efficiency.tsv` |
| `TranscriptionFactorsSchema` | `flat/transcription_factors.tsv` |
| `FoldChangesSchema` | `flat/fold_changes.tsv` |
| `PpgppRegulationSchema` | `flat/ppgpp_regulation.tsv` |

All schemas use `strict="filter"` (extra columns allowed) and `coerce=True`.
Source TSVs with `#`-prefixed comment lines should be read with
`pd.read_csv(..., sep="\t", comment="#")`.

### CLI

```bash
uv run python -m schemas.validate <schema_name> <path/to/file.tsv>
```

## Not yet schematized

Core biology tables (`genes.tsv`, `rnas.tsv`, `proteins.tsv`, `metabolites.tsv`,
reactions, transcription units, curated `_added/_modified/_removed` edits,
and the `condition/`, `mass_fractions/`, `rna_seq_data/` subdirectories)
are planned for a second phase.

## Differential-expression result tables

**Canonical format:** One file per **contrast**, one row per **tested transcript**, unfiltered by significance.

| Column             | Type   | Required | Description |
|--------------------|--------|----------|-------------|
| `transcript_id`    | string | yes      | The table's grain. Unique. |
| `gene_id`          | string | yes      | Reference gene id the transcript maps to. **Nullable and NOT unique** — a heterologous transcript maps to no reference gene, and several transcripts can share one gene. |
| `log2_fold_change` | float  | yes      | Nullable: a transcript can be tested and yield no estimate. |
| `padj`             | float  | yes      | Adjusted p-value, 0–1. Nullable: independent filtering leaves it undefined for some transcripts, which is **not** the same as "not significant". |
| `base_mean`        | float  | no       | Mean normalized count. |
| `lfc_se`           | float  | no       | Standard error of the log2 fold change. |

The stored table is deliberately **lossless** — no significance filter is applied, because padj/lfc cutoffs are an analysis choice belonging to the consumer.

The contrast's identity — groups compared, sample composition, caller settings — is **not** in the table. It belongs to the dataset's provenance record, so that one file means exactly one comparison.

### ⚠ Mapping a producer's output onto this schema is not a pure rename

A differential-expression tool's output commonly carries **several gene-identifier columns**, and one of them may itself be named `gene_id` while holding a symbol or locus tag rather than the reference id. Mapping the wrong one **validates cleanly** and silently changes which transcripts collapse together — so the gene count changes with no error at any point. `strict="filter"` then drops the unmapped columns, removing the evidence.

Map whichever column matches the TPM tables' `gene_id`, whatever the producer called it, and verify against the reference gene set.

### Reading one

```python
from schemas import DeseqResultTableSchema, select_significant_genes

validated = DeseqResultTableSchema.validate(df)
genes = select_significant_genes(validated, padj_thresh=0.1, lfc_thresh=1.0)
```

`select_significant_genes` is the canonical reading: drop rows with no `padj`/`log2_fold_change`/`gene_id`, **then** threshold, **then** collapse transcripts to one row per gene keeping the largest `|log2_fold_change|`.

⚠ **Those steps do not commute.** Collapsing before filtering picks a gene's most-extreme transcript and then tests *that* for significance — a different question, a different gene set, and no error either way. Use the function rather than re-deriving it.

⚠ Ties in `|log2_fold_change|` are **undefined** — the winner depends on row order, and a `+x`/`-x` tie selects opposite signs. This matches the pipelines the function reproduces; impose determinism yourself if you need it.

⚠ The result is the **pre-overlay** gene set. Applying it onto a reference expression table usually yields fewer genes, because genes absent from that reference drop at that later step; a downstream "applied N of M" is expected to be smaller.

## Dependencies

- `pandas`
- `pandera` (added to project dependencies in `pyproject.toml`)
