# Basan et al. (2015) — aerobic acetate overflow vs growth rate (*E. coli* NCM3722)

- **source_id:** `basan_2015`
- **Citation:** Basan M, Hui S, Okano H, Zhang Z, Shen Y, Williamson JR, Hwa T.
  "Overflow metabolism in *Escherichia coli* results from efficient proteome
  allocation." *Nature* 528(7580):99–104, 2015. DOI 10.1038/nature15765.
  (BibTeX `basan_2015` in `../../references.bib`.)
- **Raw (local, gitignored):** the main-text PDF, the Fig. 1 source-data
  supplement (`MOESM62`, = Extended Data Table 1), and the Supplementary
  Information PDF (`MOESM68`, = Supplementary Note 1 — the model derivation +
  fitted parameters used for the carbon-yield conversion below).

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

## Units & grading as a dimensionless carbon yield

*J*ac is reported **per OD600** (mM/OD600/h); the v2ecoli observable is the
acetate exchange flux in **mmol/gDCW/h**. Rather than import an unsourced
OD600→gDCW factor, the overflow card grades a **dimensionless acetate-carbon
yield** `Y_ac = acetate-C / total-carbon-influx` — per-OD vs per-gDCW cancels on
both sides, so no factor is needed. The slot stores *J*ac **as-reported**; the
conversion below is applied **in the consumer** (the card render/grader), not
baked into the data.

**The conversion (from Basan's own Supplementary Note 1 — `MOESM68`).** Basan's
proteome-allocation model gives every parameter as a number:

- `Sac = 1/3` — acetate molecules per carbon atom into fermentation, so
  fermentative carbon in = 3·*J*ac and acetate-carbon = 2·*J*ac (1 C → CO₂).
- `ef = 2`, `er = 4.4` — ATP per carbon for fermentation / respiration.
- `σ ≈ 45.7` mM/OD600 — energy demand, `J_E = σ·λ` (= `J_E,f + J_E,r`).
- `β ≈ 28.5` mM/OD600 — biomass carbon demand, `J_C,BM = β·λ`.
- `λac ≈ 0.76` h⁻¹ — the WT acetate-line threshold.

With the carbon balance `J_C,in = J_C,BM + J_C,f + J_C,r`, `J_E,f = ef·J_C,f`,
`J_C,r = J_E,r/er`, this collapses (per point) to:

```
J_C,in(λ) = (β + σ/er)·λ + (3 − 3·ef/er)·Jac  =  38.9·λ + 1.64·Jac    [mM/OD600/h, carbon]
Y_ac      = 2·Jac / J_C,in                                            [dimensionless]
```

Worked on the graded `ptsG_glucose` series: *Y*ac = 0, 0.002, 0.038, **0.10** at
λ = 0.58 / 0.64 / 0.82 / 0.95 — the threshold-linear overflow shape, unit-free.
Sanity vs an independent ¹³C anchor: NCM3722 ≈0.10 at λ≈1.0 vs MG1655 (Long 2017)
~0.22 at λ≈0.68 — NCM3722 overflows ~2× less, the known strain difference, in the
right direction.

**Caveat:** *Y*ac is **model-derived** (it carries Basan's flux-partition
assumptions — `J_E = σλ`, `J_C,BM = βλ`, `ef`/`er`), not a raw ¹³C flux — but it
is their published, fitted model, fully sourced to `MOESM68`. The overflow
**onset** (a growth-rate threshold) is unit-invariant and gradeable directly from
*J*ac without any of this; the conversion is needed only for the **slope/magnitude**.

## Notes

- 95% CIs are not extracted (the figure carries error bars; `ci_low`/`ci_high`
  are left blank). The onset/slope fit does not require them.
