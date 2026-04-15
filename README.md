# ecoli-sources

Curated, public *E. coli* omics datasets formatted for ingestion into the
[vEcoli](https://github.com/CovertLab/vEcoli) whole-cell model, plus the
schemas, processing code, and analysis utilities that produced them.

## Purpose

`ecoli-sources` is the canonical home for experimental data that feeds into
vEcoli's ParCa (parameter calculator). Moving datasets out of the vEcoli repo
lets the model code stay open while keeping a clean boundary around data
provenance, schema validation, and dataset perturbation.

Current scope:

- **Transcriptomics** — per-condition RNA-seq TPM tables with a samples
  manifest.
- **ParCa input tables** — schemas for the parameter / adjustment / half-life
  / regulation flat files that will migrate out of
  `omics-vEcoli/reconstruction/ecoli/flat/`.
- **Perturbation operators** — deterministic transforms that generate
  derived dataset variants for sensitivity campaigns.

Proteomics and flux data may be added later.

## Repository layout

```
data/                              Per-dataset TSVs + manifest (ingested by vEcoli)
  manifest.tsv                     Index of all RNA-seq datasets
  *.tsv                            TPM tables (gene_id, tpm_mean[, tpm_std])
  perturbations/                   Generated variants (gitignored; regenerated locally)

schemas/                           Pandera validation schemas
  rnaseq.py                        RnaseqTpmTableSchema, RnaseqSamplesManifestSchema
  adjustments.py                   Adjustment / pathway / metabolite schemas
  parameters.py                    Growth-rate-dependent parameters, dry mass
  half_lives.py                    RNA & protein half-life schemas
  translation.py                   Translation efficiency schemas
  regulation.py                    Transcription factors, fold changes, ppGpp
  validate.py                      CLI: validate a TSV against a named schema
  README.md                        Column-level docs for every schema

processing/                        Dataset construction and transformation
  post_processing.py               Dataset variant generators (gene exclusion, etc.)
  perturbations.py                 Perturbation operator library + campaign driver

analysis/                          Cross-dataset analysis utilities
  compare_datasets.py
  meta_analysis.py
  refdata_tpm_to_model_expression.py   TPM ↔ parca-fitted count diagnostics

scripts/                           Repo-level entry points
  validate_all.py                  Validate manifest + every referenced TPM file

datapackage.json                   Frictionless Data package descriptor
CITATION.cff                       GitHub-native citation metadata
.github/workflows/validate.yml     CI: runs schema validation on every PR
```

## Dataset format

Each RNA-seq TSV is tab-separated, one row per gene:

| Column     | Type  | Required | Description |
|------------|-------|----------|-------------|
| `gene_id`  | str   | yes      | EcoCyc gene identifier (e.g. `EG10001`); unique within file |
| `tpm_mean` | float | yes      | Mean TPM across replicates; ≥ 0 |
| `tpm_std`  | float | no       | Std dev of TPM across replicates; ≥ 0 if present |

Extra columns (e.g. gene symbol) are allowed — schemas use `strict="filter"`.

All datasets are registered in `data/manifest.tsv`. See `schemas/README.md`
for the full manifest column spec, including the provenance columns used by
generated perturbation variants (`parent_dataset_id`, `operator`,
`operator_params_json`, `seed`).

## Data sources (public only)

| Source           | Description                                         | Strain |
|------------------|-----------------------------------------------------|--------|
| vEcoli reference | Covert lab in-house reference transcriptome         | MG1655 |
| PRECISE-MG1655   | Public MG1655 RNA-seq compendium across conditions  | MG1655 |

Additional private / proprietary datasets are held in a separate sibling repo
and referenced via a secondary manifest; they are never committed here.

## Validating data

Validate a single file against a named schema:

```bash
uv run python -m schemas.validate --list
uv run python -m schemas.validate RnaseqTpmTableSchema data/vecoli_m9_glucose_plus_aas.tsv
uv run python -m schemas.validate RnaseqSamplesManifestSchema data/manifest.tsv
```

Validate the manifest and every TPM file it references (also run in CI):

```bash
uv run python scripts/validate_all.py
```

## Use from vEcoli

vEcoli's ParCa reads datasets via `rnaseq_manifest_path` +
`rnaseq_basal_dataset_id` in its config. Point it at this repo's
`data/manifest.tsv` (via the `$ECOLI_SOURCES` environment variable or a
direct path).

## Generating perturbation variants

`processing/perturbations.py` provides deterministic operators that map a
source TPM table to a perturbed variant and append a provenance row to the
manifest. Variants live under `data/perturbations/` and are **not** committed
— regenerate them locally from vEcoli's
`runscripts/run_sensitivity_campaign.py`. Pull the canonical `manifest.tsv`
from git before running any campaign, since the campaign meta-runner rewrites
it with generated rows.

## Sensitivity campaigns

The end-to-end workflow (spec → perturbed datasets → multi-parca Nextflow run
→ cross-variant sensitivity plots) is documented in vEcoli's
[`doc/sensitivity_campaigns.rst`](https://github.com/CovertLab/vEcoli/blob/multi-parca-aws/doc/sensitivity_campaigns.rst).
This repo provides the primary data + operator library; vEcoli drives the run
and runs the analyses.

## Environment

```bash
uv sync
```

Requires Python ≥ 3.12.
