# Long et al. (2017) — WT MG1655 aerobic physiology + 13C fluxome on glucose

- **source_id:** `long_2017`
- **Citation:** Long CP, Gonzalez JE, Feist AM, Palsson BO, Antoniewicz MR.
  "Fast growth phenotype of *E. coli* K-12 from adaptive laboratory evolution
  does not require intracellular flux rewiring." *Metabolic Engineering*
  44:100–107, 2017. (BibTeX in `../../references.bib`.)
- **Raw (local, gitignored):** the PDF, `…-mmc4.xlsx` ("Growth Data" — physiology),
  `…-mmc1.xlsx` ("13C-MFA" — intracellular fluxome).

## Claims used (aerobic batch, M9 glucose, **wild-type MG1655**)

| condition · observable | value | ± SE | units | kind |
|---|---|---|---|---|
| basal · growth rate (μ) | 0.677 | 0.002 | 1/h | `measured` |
| basal · glucose uptake (q_glc) | 8.46 | 0.42 | mmol/gDW/h | `measured` |
| basal · biomass yield (Yxs) | 0.444 | 0.022 | gDW/g glucose | `measured` |

## Provenance

Supplement `mmc4.xlsx`, sheet "Growth Data" — *"Physiological data of E. coli
strains grown aerobically in glucose M9 minimal medium."* The **MG1655** row:
growth rate 0.6767 ± 0.0025 h⁻¹, biomass yield 0.4442 ± 0.0218 gDW/g, glucose
uptake 8.464 ± 0.416 mmol/gDW/h, acetate production 5.566 ± 0.274 mmol/gDW/h,
oxygen uptake 11.585 ± 0.599 mmol/gDW/h, acetate yield 0.658 mol/mol.

Internal consistency: `Yxs = μ/(q_glc·M_glc) = 0.677/(8.46·0.180156) = 0.444` ✓.

## Notes

- **Strain-matched** to the model (K-12 MG1655), aerobic M9 glucose, 13C-MFA
  study. Independent of LaCroix 2015 yet near-identical WT values (μ 0.677 vs
  0.69, q_glc 8.46 vs 8.59, Yxs 0.444 vs 0.44) — a cross-validation of the
  basal reference.
- **Deferred (available):** oxygen uptake (11.59) and acetate production (5.57)
  are directly relevant to the model's respiration/overflow behavior and are
  easy future slots. The **full intracellular 13C fluxome** (`mmc1.xlsx`) is the
  vector-reference data type we plan to collect later — not extracted here.
