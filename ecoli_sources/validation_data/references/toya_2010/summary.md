# Toya et al. (2010) — WT *E. coli* central-carbon ¹³C-MFA flux map (batch glucose)

- **source_id:** `toya_2010`
- **Citation:** Toya Y, Ishii N, Nakahigashi K, Hirasawa T, Soga T, Tomita M,
  Shimizu K. "¹³C-Metabolic flux analysis for batch culture of *Escherichia
  coli* and its *pyk* and *pgi* gene knockout mutants based on mass isotopomer
  distribution of intracellular metabolites." *Biotechnol. Prog.* 26(4):975–992,
  2010. DOI 10.1002/btpr.420. (BibTeX in `../../references.bib`.)
- **Raw (local, gitignored):** the PDF and `btpr_420_sm_suppinfotable1.xls`
  (Supporting Information, Tables S-I … S-X).

## What is ingested (`data/basal/metabolic_fluxes.tsv`, `source_id=toya_2010`)

The **wild-type, glucose-consumption-phase** central-carbon flux map:
**25 reactions × 3 timepoints (5 / 6 / 7 h)** = 75 rows. Reaction identifiers are
the paper's abbreviations (SI Table S-I); `ecocyc_rxn_id` carries an EcoCyc
mapping for 23 of 25 reactions (blank for `Ack`, acetate kinase, and `Mez`,
malic enzyme).

Each row carries both:
- **`value_relative_pct`** — flux normalized to glucose uptake (PTS = 100),
  taken verbatim from SI Table S-IV.
- **`value`** — absolute flux (mmol/gDW/h), computed as
  `value_relative_pct/100 × q_glc(t)`, with the absolute glucose-uptake scale
  `q_glc = {5 h: 11.7, 6 h: 11.1, 7 h: 7.7}` mmol/gDW/h read from Figure 4a.

## ¹³C-MFA study design (Materials & Methods)

- **Nonstationary (INST) ¹³C-MFA** on **intermediate-metabolite mass-isotopomer
  distributions** measured by **CE-TOFMS** (rather than proteinogenic amino-acid
  labeling). MIDs of F16P, DHAP, 3PG, PEP, PYR, Ru5P, R5P, S7P, MAL were used in
  the glucose phase (39 independent data points). Fluxes were optimized by a
  genetic algorithm + sequential quadratic programming and the fit assessed by χ².
- **Tracer:** a mixture — **1.2 g/L (30%) [1-¹³C] + 0.8 g/L (20%) [U-¹³C] +
  2 g/L (50%) naturally-labeled glucose.** Per the paper, [1-¹³C] is released as
  CO₂ in the oxidative PP pathway but conserved in glycolysis (resolving the
  PP/glycolysis split), while [U-¹³C] constrains anaplerotic and oxidative TCA
  fluxes.
- **Culture:** aerobic batch, synthetic M9-type medium, 4 g/L glucose, 37 °C,
  1 L working volume, 1 L/min air, pH 7.0; off-gas O₂/CO₂ monitored.
- **Biomass conversion:** 1 OD₆₀₀ = 0.3 gDCW/L.
- **Reaction model:** SI Table S-I (28 reactions incl. a biomass pseudo-reaction,
  with equations and carbon-atom transitions). The ED pathway and glyoxylate
  shunt were negligible in the wild type on glucose and were excluded from the WT
  model.

## Notes

- **Time-resolved, not a single steady state.** The fluxes are reported at three
  timepoints across the glucose-consumption phase. Most reactions change
  monotonically over 5→7 h as glucose is depleted (TCA-cycle fluxes increase,
  PP-pathway and acetate-secretion fluxes decrease, upper-glycolysis fluxes change
  little), so the spread across the three timepoints reflects this physiological
  change rather than replicate measurement scatter. All three timepoints are
  retained here so the choice of timepoint is left to the consumer.
- **Relation to the 23-reaction subset used in prior model comparisons.** A
  23-reaction subset of this dataset has been distributed with the wcEcoli /
  vivarium-ecoli lineage as a single value per reaction; those values correspond
  to the **mean (± SD) across the 5 / 6 / 7 h timepoints** of the reactions below
  (e.g. glucose uptake 10.17 ± 2.16, citrate synthase 5.46 ± 1.37). This table
  reproduces those means and additionally includes the per-timepoint values and
  the `Ack` and `Mez` reactions, which the subset omitted.
- **Strain:** BW25113 (the Keio-collection parent), distinct from the K-12 MG1655
  background.
- The absolute `value` depends on the Figure-4a glucose-uptake scale;
  `value_relative_pct` is the value reported directly in SI Table S-IV.

## Other data in the same SI (not ingested here)

- Wild-type **acetate-consumption phase** flux map (Table S-V, 8 / 8.5 / 9 h).
- **Pyk** and **Pgi** knockout flux maps (Figures 6–7, Tables S-VI … S-X).
