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
