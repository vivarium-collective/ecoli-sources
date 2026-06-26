# LaCroix et al. (2015) — wild-type MG1655 aerobic physiology on glucose (ALE study)

- **source_id:** `lacroix_2015`
- **Citation:** LaCroix RA, Sandberg TE, O'Brien EJ, Utrilla J, Ebrahim A,
  Guzman GI, Szubin R, Palsson BO, Feist AM. "Use of adaptive laboratory
  evolution to discover key mutations enabling rapid growth of *Escherichia coli*
  K-12 MG1655 on glucose minimal medium." *Appl. Environ. Microbiol.*
  81(1):17–30, 2015. (BibTeX in `../../references.bib`.)
- **Raw:** `Appl. Environ. Microbiol.-2015-LaCroix-17-30.pdf` (local, gitignored).

## Claims used (aerobic batch, M9 glucose, **wild-type K-12 MG1655**)

| condition · observable | value | ± 95% CI | units | kind |
|---|---|---|---|---|
| basal · growth rate (μ) | 0.69 | 0.02 | 1/h | `measured` |
| basal · glucose uptake (q_glc) | 8.59 | 1.42 | mmol/gDW/h | `measured` |
| basal · biomass yield (Yxs) | 0.44 | 0.07 | gDW/g glucose | `measured` |

## Provenance

**Table 3** ("Phenotypic data from clones isolated from the final flask of each
experiment"), **Wild-type K-12 MG1655** row (the ancestral, pre-evolution
strain): growth rate 0.69 ± 0.02 h⁻¹, glucose uptake 8.59 ± 1.42 mmol gDW⁻¹ h⁻¹,
acetate production 3.91 ± 1.14 mmol gDW⁻¹ h⁻¹, biomass yield 0.44 ± 0.07
gDW g⁻¹ glucose. Mean ± 95% CI over three biological replicates.

Internal consistency: `q_glc = μ / (Yxs · M_glc) = 0.69 / (0.44 · 0.180156)
= 8.7 ≈ 8.59` ✓.

## What this source establishes / caveats

- **Strain-matched** to the v2ecoli model (K-12 **MG1655**), aerobic M9 glucose —
  the closest measured analog to the model's basal condition.
- The value used is the **wild-type ancestor**, *not* an evolved endpoint. The
  evolved strains (rows 3–10) reach μ ≈ 0.89–1.01 h⁻¹ with q_glc ≈ 11–14 and
  similar/slightly-lower yield (0.38–0.49) — available as an `evolved (ALE)`
  condition later if wanted, but distinct from the basal claim.
- Even the wild-type secretes acetate (3.91 mmol/gDW/h) at μ≈0.69 — so the
  measured yield (0.44) already sits below the overflow-free ceiling (0.538).
