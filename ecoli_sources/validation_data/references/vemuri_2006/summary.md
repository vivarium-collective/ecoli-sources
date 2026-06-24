# Vemuri et al. (2006) — aerobic acetate-overflow chemostat series (*E. coli* MG1655)

- **source_id:** `vemuri_2006`
- **Citation:** Vemuri GN, Altman E, Sangurdekar DP, Khodursky AB, Eiteman MA.
  "Overflow metabolism in *Escherichia coli* during steady-state growth:
  transcriptional regulation and effect of the redox ratio." *Appl. Environ.
  Microbiol.* 72(5):3653–3661, 2006. DOI 10.1128/AEM.72.5.3653-3661.2006.
  (BibTeX `vemuri_2006` in `../../references.bib`.)
- **Obtained via:** the **Zhuang et al. 2011** Supplementary Information
  (`MOESM2_ESM.xls`, sheet "Growth Data - Fig 2"), which compiled this Vemuri
  chemostat series for its Fig. 2. We hold the Zhuang SI (BibTeX `zhuang_2011`);
  the underlying measurements are Vemuri 2006. **Raw (local, gitignored):** the
  Zhuang SI `.xls` (`MOESM2_ESM.xls`) and the Zhuang paper PDF (`msb.2011.34.pdf`)
  live in this directory; the Vemuri 2006 paper itself is not re-ingested (the data
  is identical, so re-extraction would be redundant).

## What is ingested (`data/perturbation/overflow/acetate_vs_growth.tsv`)

The **7 averaged dilution-rate points** of the WT glucose-limited aerobic
chemostat series (`series = wt_chemostat`, `series_type = chemostat`), one row
per setpoint (each the mean of **n = 3** replicate runs). `ResponseCurveSchema`.

- **`value`** — acetate excretion *q*Ac, **mmol/gDW/h** (native model units, unlike
  Basan's per-OD). **`ci_low`/`ci_high`** = mean ± SD across the 3 reps.
- **`growth_rate_per_h`** — the dilution rate *D* (= μ at steady state), the x-axis.
- **`glucose_uptake`** (GUR), **`o2_uptake`** (*q*O₂), **`co2_evolution`** (*q*CO₂),
  **`biomass_yield`** (Yx/s, g/g) — **co-measured** at the same operating point,
  all mmol/gDW/h (yield g/g).
- **`perturbation` / `level`** — `dilution_rate` / the setpoint *D*.

## Why this is the primary overflow-yield reference

It carries **glucose uptake co-measured with acetate**, so the dimensionless
**acetate-carbon yield** `Y_ac = 2·qAc / (6·GUR)` is computed **directly from each
row** — raw, empirical, no model assumptions and no OD600→gDCW factor (unlike the
Basan-model-derived yield). The series spans the full transition:

| *D* (h⁻¹) | GUR | *q*Ac | **Y_ac** |
|---|---|---|---|
| 0.05–0.30 | 0.68–3.36 | 0 | 0.00 |
| 0.41 | 4.71 | 0.14 | 0.01 |
| 0.49 | 5.96 | 1.82 | 0.10 |
| 0.60 | 7.89 | 6.01 | **0.25** |

The threshold-linear overflow shape, in dimensionless yield. *Y_ac* = 0.25 at
*D* = 0.6 agrees with the independent MG1655 ¹³C-MFA anchors (Long 2017 0.22,
Crown 2015 0.22 at μ≈0.68) — same strain, consistent magnitude.

## Measurement design

- **Strain:** *E. coli* K-12 **MG1655** (WT) — strain-matched to the v2ecoli
  model. (The companion ΔarcA mutant QC2575 in Vemuri 2006 is a separate arm, not
  ingested here.)
- **Culture:** aerobic, **glucose-limited continuous culture (chemostat)** in a
  defined minimal medium (5 g/L glucose). The chemostat **clamps the growth rate**
  (*D*) and lets uptake follow — the inverse causal direction of a batch
  uptake-titration (Basan), but the steady-state (*D*, GUR, *q*Ac) operating
  points are directly comparable as a yield curve.

## Notes

- The overflow **onset** here (*D* ≈ 0.4) is lower than Basan's batch threshold
  (λac ≈ 0.76); onset varies with strain/condition and is treated as a secondary,
  reported feature. The **shape** (0 → linear) and the **yield slope** are the
  graded, strain-robust quantities.
- O₂/CO₂ are retained for a possible respiration-vs-growth axis (a bonus the
  chemostat data supports).
