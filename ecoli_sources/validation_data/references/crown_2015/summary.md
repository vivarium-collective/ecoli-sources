# Crown et al. (2015) — WT *E. coli* MG1655 central-carbon flux map (COMPLETE-MFA, 14 parallel labeling experiments)

- **source_id:** `crown_2015`
- **Citation:** Crown SB, Long CP, Antoniewicz MR. "Integrated ¹³C-metabolic
  flux analysis of 14 parallel labeling experiments in *Escherichia coli*."
  *Metabolic Engineering* 28:151–158, 2015. DOI 10.1016/j.ymben.2015.01.001.
  (BibTeX in `../../references.bib`.)
- **Raw (local, gitignored):** the PDF and `…-mmc1.xlsx` (Supplementary Table S4,
  "Results of ¹³C-metabolic flux analysis").

## What is ingested (`data/basal/metabolic_fluxes.tsv`, `source_id=crown_2015`)

The **COMPLETE-MFA** intracellular net-flux map: **71 reactions**, taken from
the `COMPLETE-MFA` column of Supplementary Table S4. Fluxes are reported
**normalized to a glucose uptake rate of 100**, so for this source `value`
carries the normalized flux (`units = pct_glucose_uptake`) and equals
`value_relative_pct`; `uncertainty` is the flux standard deviation derived from
the reported 95% confidence interval as `(UB95 − LB95) / 3.92`. Reaction
identifiers are `crown_f<N>` (the paper's Flux No.); `reaction_name` carries the
reaction equation. There is no `timepoint` — the values are a single
exponential-steady-state estimate.

The supplement also reports 23 exchange (reversibility) fluxes and per-amino-acid
G-values, and per-tracer flux maps for each of the 14 individual experiments;
only the integrated intracellular net fluxes are ingested here.

## ¹³C-MFA study design (Materials & Methods)

- **Method: COMPLETE-MFA** (complementary parallel labeling experiments) —
  integrated ¹³C-MFA fitting **14 parallel labeling experiments** simultaneously
  to a single flux map (8 experiments from this study + 6 singly-labeled from
  Leighty & Antoniewicz 2013). EMU framework (Metran). 811 mass-isotopomer
  measurements + 2 external flux constraints (glucose, acetate) → 10 net free
  fluxes, 9 exchange fluxes, 129 G-values; 665 redundant measurements; combined
  SSR 690 < 738 (statistically acceptable, 95%). **Accurate 68% and 95%
  confidence intervals computed for every flux.**
- **Parallel labeling, not a single mixed experiment.** Tracers fed to separate
  parallel cultures: [1,2-¹³C], [2,3-¹³C], [4,5,6-¹³C], [2,3,4,5,6-¹³C]glucose,
  mixtures [1-¹³C]+[4,5,6-¹³C] (1:1) and [1-¹³C]+[U-¹³C] (1:1 and 4:1), 20%
  [U-¹³C], plus six singly-labeled [1-]…[6-¹³C]glucose. (Distinct from a single
  culture co-fed a tracer mixture.)
- **Strain:** *E. coli* **K-12 MG1655** (ATCC 700925) — the model's background.
- **Culture:** aerobic batch, **M9 minimal medium**, glucose sole carbon source
  (~2.55 g/L initial), 37 °C, aerated mini-bioreactors (5 mL working volume).
- **Sampling:** during **exponential growth phase**; biomass composition
  confirmed constant over the sampling window (a metabolic/balanced steady
  state). Mass isotopomers of TBDMS-derivatized proteinogenic amino acids by
  GC–MS.
- **Physiology (8 parallel cultures):** specific growth rate 0.72 ± 0.02 h⁻¹,
  biomass yield 0.38 ± 0.02 gDW/g; growth was reproducible across tracers
  (tracer choice did not affect macroscopic growth). Biomass molecular formula
  CH₁.₈O₀.₅N₀.₂ (24.6 gDW/C-mol); OD₆₀₀→0.32 gDW/L.

## Notes

- **Single steady-state vector with rigorous CIs** — unlike a time-resolved
  dataset, there is no sampling-phase trajectory to reconcile; the reported
  uncertainty is the statistical confidence interval of the integrated fit.
- **Strain-matched** to the model (K-12 MG1655), aerobic M9 glucose, exponential
  growth.
- Fluxes are normalized to glucose uptake = 100; the absolute glucose uptake is
  not tabulated per reaction (the integrated fit fixes glucose uptake as the
  normalization basis).
- `reaction_id` (`crown_f<N>`) is this source's native numbering; mapping to
  model / EcoCyc reaction ids (and to other sources' reaction ids) is a
  downstream concern, not applied here.
