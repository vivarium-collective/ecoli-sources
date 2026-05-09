# Bundles

This document is the authoring guide for **bundles** — the data packaging
format ecoli-sources ships to its consumers (primarily the
[vEcoli](https://github.com/vivarium-collective/vEcoli) whole-cell model).

For the *consumer-side* picture (how vEcoli loads bundle data, how to
configure ParCa to use a custom bundle), see vEcoli's
[`doc/data_ingestion.rst`](https://github.com/vivarium-collective/vEcoli/blob/main/doc/data_ingestion.rst).

---

## What a bundle is

A bundle is a TSV manifest of the form `(canonical_key, source_path,
description, schema_name)`, plus the data files it points at. Each row
binds a **canonical key** — an addressable role in the consumer's model
— to a specific file under the bundle's data root.

```
canonical_key             source_path                                          description       schema_name
rnaseq_basal_tpms         rnaseq_experimental/vecoli_m9_glucose_minus_aas.tsv  Basal TPMs        RnaseqTpmTableSchema
rnaseq_experimental_tpms  rnaseq_experimental/vecoli_m9_glucose_minus_aas.tsv  Expt TPMs         RnaseqTpmTableSchema
genes                     flat/genes.tsv                                       Gene metadata
metabolic_reactions       flat/metabolic_reactions.tsv                         Reaction network
metabolite_concentrations flat/metabolite_concentrations.tsv                   Pool sizes
...
```

The default reference bundle ships at
`ecoli_sources/data/reference_bundle.tsv` and is exposed via
`ecoli_sources.BUNDLE_PATH` for consumers that install the package.

The set of canonical keys vEcoli's loader expects is a **contract**
between the supplier (this repo) and the consumer (vEcoli). The
contract is enforced at bundle-load time — any bundle missing a
required key fails validation with an error naming it.

---

## Why this shape

- **Reproducibility**: a bundle's manifest *is* the contract; the
  provenance of a model run is captured as a single hash.
- **Like-for-like substitution**: variant bundles (perturbed
  metabolite concentrations, alternative RNA-seq conditions, swapped
  kinetic parameters) are cheap to express — copy the default
  manifest, override the rows you want, point the consumer at the
  copy.
- **No override / merge semantics at runtime**: the consumer receives
  one complete manifest; there is no diff-against-default logic in the
  loading path. Manifest construction (variant generators, diff
  tools) is a build-time concern outside the runtime contract.
- **Schemas live with the data**: each canonical key can declare a
  Pandera schema name; the validation pipeline applies it.

---

## Manifest schema

The manifest itself is validated by `ReferenceBundleSchema` in
`schemas/reference_bundle.py`. Required columns:

| Column | Type | Description |
|---|---|---|
| `canonical_key` | string | Slot name (snake_case). Primary key; stable across file renames. |
| `source_path` | string | Path to the source file or directory, relative to the bundle root (`ecoli_sources/data/`). |
| `description` | string | One-liner explaining the slot's role. |
| `schema_name` | string | Optional. Pandera schema class name (resolved via `getattr(schemas, schema_name)`); empty for rows without a schema yet. |

The schema enforces uniqueness on `canonical_key` and a DataFrame-level
check that every entry in `REQUIRED_CANONICAL_KEYS` is present in the
manifest. Variants are free to add extra rows beyond the contract list.

---

## Canonical key naming

The convention for keys derived from files under `flat/` is
**path-relative-to-flat with `/` → `__` and the extension dropped**:

| File | Canonical key |
|---|---|
| `flat/genes.tsv` | `genes` |
| `flat/metabolic_reactions.tsv` | `metabolic_reactions` |
| `flat/mass_fractions/glycogen_fractions.tsv` | `mass_fractions__glycogen_fractions` |
| `flat/rrna_options/remove_rrff/genes_removed.tsv` | `rrna_options__remove_rrff__genes_removed` |
| `flat/sequence.fasta` | `sequence` |

Per-file granularity (rather than directory-level keys) gives the
bundle full path-resolution responsibility — option-driven subdirs
(`new_gene_data/<subdir>/`, `rrna_options/<subdir>/`) compose canonical
key strings from runtime options + filename, no directory traversal in
the consumer.

Multi-file groups (e.g., `metabolic_reactions{,_added,_modified,_removed}`)
become separate canonical keys; the consumer's loader joins them.

For RNA-seq and other category-specific data, keys are named by their
*role* in the consumer rather than by file path
(`rnaseq_basal_tpms`, `rnaseq_experimental_tpms`).

---

## Validation pipeline

Three layers, all available via `scripts/validate_all.py` (CI) and
`scripts/validate_bundle.py` (CLI for ad-hoc bundle validation):

1. **Manifest schema** — `ReferenceBundleSchema` checks the TSV's
   columns and types and asserts the canonical-key contract.
2. **Path resolution** — every `source_path` must resolve to an
   existing file under the bundle's data root.
3. **Content schemas** — for rows with `schema_name` set, the
   referenced Pandera schema validates the file's contents.

Layers 1 and 2 are eager (run at consumer load time too, via
`SourceBundle.__init__`); layer 3 runs in CI and on demand.

To validate a bundle from the command line:

```bash
uv run python scripts/validate_bundle.py path/to/reference_bundle.tsv
```

---

## Adding a new canonical key

Adding a new key has three coordinated steps:

1. **Drop the file** under `ecoli_sources/data/` (typically inside
   `flat/` or a category subdirectory).
2. **Add the row** to `ecoli_sources/data/reference_bundle.tsv`:
   `canonical_key`, `source_path` (relative to `data/`),
   `description`, optional `schema_name`.
3. **Update `REQUIRED_CANONICAL_KEYS`** in
   `schemas/reference_bundle.py` if the consumer model now requires
   this key.

Run `scripts/validate_all.py` to confirm the bundle validates.

If the consumer (vEcoli) does not yet know about the new key, the
addition is non-breaking — extra rows are allowed by the contract
schema. The consumer must be updated separately to *use* the key (the
canonical-key reference site in `KnowledgeBaseEcoli` or wherever it
lands).

If a Pandera schema is added for the file, drop the schema class in
`schemas/<category>.py`, export it from `schemas/__init__.py`, and set
`schema_name` on the bundle row to the class name.

---

## Creating a variant bundle

A variant is a complete bundle manifest that differs from the default
in some subset of rows (e.g., a different RNA-seq table, perturbed
metabolite concentrations, swapped kinetic parameters).

Two ways to create one:

- **Copy and edit**: copy `reference_bundle.tsv`, replace `source_path`
  values for the rows you want to change, save the result somewhere
  the consumer can read it. The required-canonical-keys contract
  must still be satisfied.
- **Programmatic generation**: build the variant manifest from a base
  + a `{canonical_key: override}` dict in a script. (A standard
  variant-generator helper is on the roadmap; for now this is hand-
  rolled per use case.)

Validate the variant before handing it to the consumer:

```bash
uv run python scripts/validate_bundle.py path/to/my_variant/reference_bundle.tsv
```

The consumer (vEcoli) takes the variant via its `bundle_manifest_path`
config option.

---

## Reproducibility

- The default bundle's content is determined by the ecoli-sources
  commit; the consumer (vEcoli) pins ecoli-sources by SHA in its
  `pyproject.toml`. A bundle bump in the consumer is a visible PR.
- Variant bundles are normally checked in alongside the campaign or
  experiment that uses them, with a deterministic generator so the
  manifest is reproducible from inputs.

---

## See also

- `schemas/README.md` — column-level documentation for every Pandera
  schema in the repo.
- `schemas/reference_bundle.py` — the manifest schema and the
  authoritative `REQUIRED_CANONICAL_KEYS` list.
- `scripts/validate_all.py`, `scripts/validate_bundle.py` —
  validation entry points.
- vEcoli's [`doc/data_ingestion.rst`](https://github.com/vivarium-collective/vEcoli/blob/main/doc/data_ingestion.rst)
  — consumer-side documentation.
