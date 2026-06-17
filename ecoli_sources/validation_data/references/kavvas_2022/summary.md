# Kavvas et al. (2022) — WT MG1655 aerobic physiology + 13C fluxome on glucose

- **source_id:** `kavvas_2022`
- **Citation:** Kavvas ES, Long CP, Sastry A, Poudel S, Antoniewicz MR, Ding Y,
  Mohamed ET, Szubin R, Monk JM, Feist AM, Palsson BO. "Experimental evolution
  reveals unifying systems-level adaptations but diversity in driving
  genotypes." *mSystems* 7:e00165-22, 2022. (BibTeX in `../../references.bib`.)
- **Raw (local, gitignored):** the PDF and `…-s0007….csv` (multi-strain 13C-MFA
  flux table + absolute physiology).

## Claims used (aerobic, glucose, **wild-type MG1655**)

| condition · observable | value | ± SE | units | kind |
|---|---|---|---|---|
| basal · growth rate (μ) | 0.81 | 0.03 | 1/h | `measured` |
| basal · glucose uptake (q_glc) | 10.58 | — | mmol/gDW/h | `measured` |
| basal · biomass yield (Yxs) | 0.42 | 0.02 | gDW/g glucose | `measured` |

## Provenance

Supplementary data S7 (the 13C-MFA flux table; strains in columns). The
**MG1655 / WT** column rows: "ABSOLUTE GLUCOSE UPTAKE RATE (mmol/gdw/h)" =
10.58; "Estimated Biomass yield (gdw/g)" = 0.42 ± 0.02; "Estimated Growth rate
(1/h)" = 0.81 ± 0.03.

Internal consistency: `Yxs = μ/(q_glc·M_glc) = 0.81/(10.58·0.180156) = 0.425` ✓.

## Notes

- **Strain-matched** to the model (MG1655). μ here (0.81) is notably higher than
  Long/LaCroix (~0.68); yield (0.42) is consistent. Carry the higher μ as
  honest measured spread.
- **Deferred (available):** the **full multi-strain intracellular 13C fluxome**
  (S7) is the vector-reference data type to collect later — not extracted here.
