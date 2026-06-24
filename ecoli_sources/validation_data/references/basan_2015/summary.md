# Basan et al. (2015) — aerobic acetate overflow vs growth rate (*E. coli* NCM3722)

- **source_id:** `basan_2015`
- **Citation:** Basan M, Hui S, Okano H, Zhang Z, Shen Y, Williamson JR, Hwa T.
  "Overflow metabolism in *Escherichia coli* results from efficient proteome
  allocation." *Nature* 528(7580):99–104, 2015. DOI 10.1038/nature15765.
  (BibTeX `basan_2015` in `../../references.bib`.)
- **Raw (local, gitignored):** the main-text PDF + the Fig. 1 source-data
  supplement (`MOESM62`, = Extended Data Table 1).

## What is ingested (`data/perturbation/overflow/acetate_vs_growth.tsv`)

The **Fig. 1 source data** — acetate excretion rate *J*ac vs growth rate *λ*,
for *E. coli* **NCM3722**, as **31 rows** across four **series** (`series_type`).
`ResponseCurveSchema`; one row per `(source, series, level)`, ordered along
`growth_rate_per_h` (the x-axis Basan plots against). Columns:

- **`value`** — acetate excretion *J*ac, in **`mM/OD600/h`** (verbatim header:
  *"Acetate excretion rate Jac (mM/OD600/h)"*). May be ≤ 0 — no detectable
  excretion (measurement noise around zero), below the overflow onset.
- **`growth_rate_per_h`** — growth rate *λ* (1/h); the curve x-axis.
- **`series` / `series_type`** — which curve, and how *λ* was varied:
  - **`uptake_titration`** — carbon uptake set by titrating *transporter
    expression* (not the carbon source): `ptsG_glucose` (strain NQ1243, Pu-*ptsG*,
    glucose uptake titrated by the inducer 3MBA at 0/20/300/800 µM) and
    `lacY_lactose` (strain NQ381, titratable LacY). **The overflow card grades
    the `ptsG_glucose` series** — the apples-to-apples match for a glucose-cache
    model GUR sweep.
  - **`uptake_mutant`** — *glpK* mutants (NQ636/638/640) on glycerol.
  - **`carbon_source`** — WT on 13 carbon sources at their natural growth rates
    (earmarked for a future carbon-source-swap study, not graded here).
  - **`carbon_source_aa`** — WT on carbon sources + 7 amino acids (richer media).
- **`perturbation` / `level`** — the swept variable and its value at each point
  (e.g. `inducer_3MBA_uM`=300; `carbon_source`=glycerol 0.2%).
- **`condition` / `strain`** — media + strain descriptor as reported.

## The claim (for band/onset provenance)

> "acetate overflow is an innate response that depends on the degree of carbon
> influx and not specifically on the nature of carbon sources." (main text)

Operationally: acetate is ~0 below a critical growth rate *λ*ac (≈0.7–0.8 h⁻¹
here) and rises ~linearly with *λ* above it (Basan eq. 1, the "acetate line").
The `curve_response` criterion grades a model's acetate-vs-*λ* response on the
**onset (λac)** and the **slope** above it.

## Measurement design

- **Strain:** wild-type K-12 **NCM3722** and its titratable/mutant derivatives
  (NQ1243 Pu-*ptsG*, NQ381 titratable *lacY*, NQ636/638/640 *glpK*). *(Strain
  caveat: the v2ecoli model is MG1655.)*
- **Culture:** aerobic **batch** growth; acetate accumulation reported **per
  OD600**. Growth rate is **emergent** from the uptake setting (the titration
  sets transporter expression → uptake capacity → *λ*); a chemostat would invert
  this (clamp *λ*, let uptake follow). A model GUR sweep (clamp the uptake bound →
  *λ* emerges) mirrors the titration's causal direction.

## ⚠️ Unit reconciliation (a grade-time assumption, NOT from this paper)

*J*ac is reported in **mM/OD600/h**; the v2ecoli observable is the acetate
exchange flux in **mmol/gDCW/h**. The OD600→gDCW factor is **not stated in Basan
2015** (everything is per-OD). Two faithful options, decided at grade time:
(1) convert the experimental *J*ac with an **externally-sourced** NCM3722 factor
(companion Hwa-lab papers; ~0.4–0.5 gDCW·L⁻¹·OD600⁻¹), flagged with its own
citation; or (2) convert the **model** side to per-OD and grade in the reported
units. **No unsourced factor is baked into the TSV** — values are as-reported.
The overflow **onset** (a growth-rate threshold) is unit-invariant and gradeable
without this factor; the **slope** is not.

## Notes

- 95% CIs are not extracted (the figure carries error bars; `ci_low`/`ci_high`
  are left blank). The onset/slope fit does not require them.
