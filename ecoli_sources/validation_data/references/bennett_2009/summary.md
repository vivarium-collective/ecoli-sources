# Bennett et al. (2009) — absolute intracellular metabolite concentrations (*E. coli*, 3 carbon sources)

- **source_id:** `bennett_2009`
- **Citation:** Bennett BD, Kimball EH, Gao M, Osterhout R, Van Dien SJ,
  Rabinowitz JD. "Absolute metabolite concentrations and implied enzyme active
  site occupancy in *Escherichia coli*." *Nat. Chem. Biol.* 5(8):593–599, 2009.
  DOI 10.1038/nchembio.186. (BibTeX in `../../references.bib`.)
- **Raw (local, gitignored):** `bennett_2009.pdf` (main text) and
  `bennett_2009_SI.pdf` (Supplementary Information).

## What is ingested (`data/basal/metabolite_pools.tsv`, `source_id=bennett_2009`)

The **absolute intracellular concentrations** of **103 metabolites** (mol/L),
across **three carbon sources** — glucose, glycerol, acetate — from
**Supplementary Table 3** (the multi-condition superset of main-text Table 1,
which is glucose-only). One row per `(metabolite, condition)`:

- **`value`** — best-estimate concentration (mol/L).
- **`ci_low` / `ci_high`** — the reported **95% confidence interval** bounds
  (asymmetric; concentrations are roughly log-normal), kept as bounds rather
  than a symmetric ±.
- **`condition`** — `glucose` (103 metabolites) / `glycerol` (70) / `acetate`
  (68); some metabolites were quantified only in glucose.
- **`significance`** — the paper's cross-condition significance flag
  (1 = glucose vs glycerol; 2 = glucose vs acetate; 3 = glycerol vs acetate
  significantly different at FDR 0.05).
- `metabolite` is the name as the paper reports it; `ecocyc_id` carries the
  corresponding EcoCyc / model compound frame id where one exists (**93 of the
  103** metabolites). It is left **blank** for the 10 that have no single id:
  three are footnoted **lumped pools** (`hexose-P` = F6P+G6P+G1P, `pentose-P`
  = R5P+Ru5P+X5P, `isoleucine+leucine` co-eluted) that map to several ids, and
  seven are species not carried as Bennett concentration targets in the
  reconstruction (`flavin mononucleotide`, `glutathione disulfide`,
  `inosine-mono/di/tri-phosphate`, `deoxyribose-5-P`, `PRPP`). The ids were
  assigned by joining the published names against the reconstruction's
  Bennett-curated `metabolite_concentrations.tsv` (name + 2-significant-figure
  value agreement).

The **basal** card uses the **glucose** rows; glycerol/acetate are retained for
future multi-condition work.

## Measurement design (Materials & Methods)

- **Strain:** wild-type K-12 **NCM3722** (prototrophic). *(Strain caveat: the
  v2ecoli model is MG1655; NCM3722 is a different K-12 isolate.)*
- **Culture:** aerobic, exponential growth on **filters atop agarose minimal
  medium** (gentle, fast quenching). Doubling times **77 / 89 / 139 min** in
  glucose / glycerol / acetate. *(The glucose td of 77 min is somewhat slower
  than fast batch growth; filter culture trades growth rate for quench speed.)*
- **Method:** **LC-MS/MS** with an **isotope-ratio internal-standard** approach
  — the unlabeled metabolome was quantified against extracts of cells grown on
  U-¹³C substrate as internal standards; concentrations corrected for incomplete
  ¹³C labeling. Standards limited to single metabolites, freshly prepared.
- Measured species were >90% intracellular. Total pool ≈ **300 mM**, dominated
  on a molar basis by amino acids (49%), nucleotides (15%), central-carbon
  intermediates (15%), redox cofactors + glutathiones (9%). **Glutamate** is the
  single most abundant (~96 mM in glucose).

## Notes

- **Footnoted lumped pools** (per main-text Table 1): `Hexose-P` is the combined
  F6P + G6P + G1P pool; `Pentose-P` the combined R5P + Ru5P + X5P pool; the
  3-phosphoglycerate value may be overestimated by in-prep degradation of
  1,3-bisphosphoglycerate. SI Table 3 lists these under their individual names
  where resolved.
- The **glucose** column of SI Table 3 reproduces main-text Table 1 (with full
  precision and CIs); spot-checks agree (glutamate 9.6×10⁻², acetyl-CoA
  6.06×10⁻⁴, 6-phosphogluconate 3.77×10⁻³).
- **Ingestion fixes:** one name was repaired from a PDF-extraction truncation —
  `phosphate` (3.74×10⁻⁴) → `dihydroxyacetone phosphate` (value + EcoCyc
  `DIHYDROXY-ACETONE-PHOSPHATE` confirm); `isoleucine+leuicine` →
  `isoleucine+leucine`.

## Other data in the same SI (not ingested here)

- Supplementary Table 1/2 — labeling completeness; Supplementary Tables 4–6 —
  intracellular/extracellular classification, ΔG°′ comparisons, computed free
  energies (TMFA). Enzyme-saturation (Kₘ) analysis is in the main text (Fig. 2).
