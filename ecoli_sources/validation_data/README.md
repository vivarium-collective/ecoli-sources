# Validation data

Curated **experimental and theoretical reference values that model *outputs*
are graded against** — distinct from the ParCa reference bundle
(`../data/reference_bundle.tsv`), which ships data fed *into* the model.

The two are deliberately separate concerns:

| | ParCa reference bundle | Validation data (this subsystem) |
|---|---|---|
| Direction | data **into** the model | model output graded **against** data |
| Consumer | ParCa / `KnowledgeBaseEcoli` | report cards / behavioral evaluators |
| Contract | every `REQUIRED_CANONICAL_KEYS` key must resolve, or load fails | no required keys — a missing slot makes a card axis `ungraded` |

## Layout

```
validation_data/
  validation_bundle.tsv          # manifest: <condition>__<observable> -> claim table
  references.bib                 # BibTeX entries for every source
  references/                    # provenance layer (the "mini-wiki")
    <source_id>/
      summary.md                 # what it measured, conditions, caveats, verbatim passages
      raw/                       # extracted figure/table data, as published
  data/                          # curated, card-consumable
    <condition>/
      <observable>.tsv           # one row per source (see Claim tables)
```

**Condition-first.** Data is organized by *condition* (`basal/`, later
`glucose_titration/`, `ko_<gene>/`, …) then by *observable*
(`biomass_yield`, `growth_rate`, …). A report card maps to one condition, so
it reads one directory; and `data/<condition>/` is the natural unit to tie to
a perturbation data bundle later (same condition key on both sides).

## Claim tables — a metric is a *table*, not a value

Each `data/<condition>/<observable>.tsv` holds **one row per source**, validated
by `ScalarClaimSchema` (`schemas/validation_claim.py`). Multiple sources per
observable are the norm. The `kind` column separates the three reference classes
a card treats differently:

- **`measured`** — an experimental observation. The grade target.
- **`theoretical_max`** — a first-principles ceiling (thermodynamic or
  stoichiometric / FBA). A model output that **exceeds** a `theoretical_max` is
  unambiguously broken — a *harder, differentiated failure* than deviating from
  a measured value, requiring no adequacy judgment.
- **`model_predicted`** — a value from another model, for context.

Every `source_id` must resolve to a `references/<source_id>/` page (enforced by
`scripts/validate_all.py`), keeping each curated number tied to its provenance.

### Scalar vs vector claims

A scalar observable (`biomass_yield`, `growth_rate`, …) uses `ScalarClaimSchema`
— one `value` per source. An observable that is a **vector** — a per-reaction
flux map, a per-gene proteome — uses a vector schema instead, where one source
contributes many addressable rows keyed by a within-vector id. Today:

- **`ReactionFluxSchema`** (`data/basal/metabolic_fluxes.tsv`) — ¹³C-MFA flux
  maps, keyed by `reaction_id`, one row per `(source, reaction, timepoint)`.
  Each row carries both `value` (absolute flux) and `value_relative_pct` (flux
  normalized to glucose uptake), so a consumer can compare absolute fluxes or
  scale-invariant carbon routing as appropriate. A source reporting a time
  course carries a `timepoint`; such a course is a trajectory rather than
  replicate measurements — see the per-source `summary.md` before aggregating
  across timepoints.
- **`ProteinAbundanceSchema`** (`data/basal/proteome.tsv`) — proteome maps,
  keyed by `gene`, one row per `(source, gene)`; `value` is absolute protein
  abundance (copies/cell) with an optional replicate `cv`. The proteome analogue
  of the flux vector (typically graded as a log-log R² over shared genes).

### Cultivation-keyed vector observations

The two vector schemas above are *per-source* references: one row per
`(source, entity)`, the vector analogue of `ScalarClaimSchema`.
`VectorObservationSchema` is the cultivation-keyed analogue — the same move
`ScalarObservationSchema` makes for scalars. One row per
`(cultivation_group_id, observable, entity_id, units)`, holding a measurement
this program produced rather than one taken from the literature.

Two things distinguish it from the scalar assay table:

- **It ships pre-aggregated.** Where a scalar table keeps one row per replicate
  and lets the consumer aggregate, a vector observation summarizes an interval
  named by `phase`/`window`, and carries its own dispersion (`sd_log10`) and
  counts (`n`, `n_pos`, `n_total`). Choosing which timepoints represent a
  condition can depend on knowledge absent from the data, so it is a curation
  decision made upstream, not a default a consumer should re-take.
- **It records detection.** A `detection` column separates a genuine measured
  zero from "looked for, absent" — a distinction a bare number cannot carry, and
  one that vanishes if absence is encoded as a missing row or as a null the
  reader coerces.

It is **measurement-only**: no model-side value, no simulation provenance.
Pairing a measurement with a model — resolving the shared entity set,
renormalizing, computing concordance — belongs to whatever grades the
comparison. Baking a model vector in alongside would make the measurement's own
content model-dependent, since the entity set becomes an intersection with that
model's id-space and any normalization over the intersection moves when model
coverage moves.

Vector slots are located, not read, by the loader — a bundle can hold many and
each can hold thousands of rows:

```python
from ecoli_sources.validation import load_vector_observations, read_vector_observations

slots = load_vector_observations()          # {canonical_key: slot}
frame = read_vector_observations(slots["basal__proteome"])
```

## Consuming it

```python
from ecoli_sources import VALIDATION_BUNDLE_PATH
# resolve a (condition, observable) slot -> its claim table, then read + filter
# rows by `kind`. A missing key degrades the card axis to `ungraded`.
```

Because the consumer-side `SourceBundle` is coupled to the ParCa contract, read
this bundle with the ParCa-contract validation **off** (it has none of the
required ParCa keys) — or via a validation-specific loader.

## Adding a source

1. Add `references/<source_id>/summary.md` (+ `raw/` if you extracted data) and
   a `references.bib` entry.
2. Add a row to the relevant `data/<condition>/<observable>.tsv` with the
   `source_id`, `value`, `units`, and `kind`.
3. Add a `validation_bundle.tsv` row if the `(condition, observable)` slot is new.
4. `uv run python scripts/validate_all.py` to confirm.
