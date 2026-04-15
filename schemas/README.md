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

## Dependencies

- `pandas`
- `pandera` (added to project dependencies in `pyproject.toml`)
