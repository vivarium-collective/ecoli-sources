# Bremer & Dennis (2008) — macromolecular composition vs growth rate (*E. coli* B/r)

- **source_id:** `bremer_dennis_2008`
- **Citation:** Bremer H, Dennis PP. "Modulation of Chemical Composition and
  Other Parameters of the Cell at Different Exponential Growth Rates." *EcoSal
  Plus* 3(1), 2008. DOI 10.1128/ecosal.5.2.3. (BibTeX in `../../references.bib`.)
- **Raw (local, gitignored):** `bremer_dennis_2008.pdf`.

## What is ingested (`data/basal/macromolecular_composition.tsv`, `source_id=bremer_dennis_2008`)

The **protein / RNA / DNA dry-mass fractions** as a function of growth rate,
derived from **Table 2** (macromolecular composition of exponentially growing
*E. coli* B/r at 37 °C). One row per `(component, doubling_time)`:

- **`component`** — `protein`, `rna` (total RNA: rRNA+tRNA+mRNA), `dna`.
- **`value`** — fraction of total dry weight (0–1).
- **`doubling_time_min`** — 100 / 60 / 40 / 30 / 24 / 20 min
  (= growth rate **0.6 / 1.0 / 1.5 / 2.0 / 2.5 / 3.0 doublings/h**).

3 components × 6 growth rates = **18 rows**. A consumer grades against the row
(or interpolation) matching the model's growth rate.

### Derivation (self-verified)

Dry-mass fractions are **not tabulated directly**; they were computed from
Table 2's per-cell amounts:

- protein μg = `(protein residues/cell) × 108 g/mol / Nₐ` (avg aa residue 108 g/mol)
- RNA μg = `(nucleotide residues/cell) × 324 g/mol / Nₐ` (avg nucleotide 324 g/mol)
- DNA μg = `(genome equivalents/cell) × M_genome / Nₐ` (4700 kbp × ~650 g/mol/bp)
- total dry mass μg = `(mass in OD₄₆₀ units/cell) × 173`

The summed protein+RNA+DNA reproduces Table 2's own `Sum (P+R+D)` row at every
growth rate (e.g. 167.6 vs 167 at td 100 min; 382 vs 383 at td 40 min) — an
internal cross-check that the per-cell rows and molecular weights are read
consistently.

### Resulting fractions

| td (min) | μ (dbl/h) | protein | RNA | DNA |
|---:|---:|---:|---:|---:|
| 100 | 0.6 | 0.606 | 0.103 | 0.036 |
| 60 | 1.0 | 0.561 | 0.115 | 0.027 |
| 40 | 1.5 | 0.531 | 0.136 | 0.023 |
| 30 | 2.0 | 0.495 | 0.164 | 0.020 |
| 24 | 2.5 | 0.469 | 0.195 | 0.020 |
| 20 | 3.0 | 0.416 | 0.209 | 0.020 |

The classic trend: as growth rate rises, **RNA fraction increases** (ribosome
demand) and **protein/DNA fractions fall**.

## Caveats (important — read before using)

- **Strain:** *E. coli* **B/r**, not K-12 MG1655.
- **Not a glucose-minimal dataset.** Table 2 is a **smooth growth-rate curve fit
  through multiple media** — minimal media with different carbon sources
  (succinate, glycerol, glucose; glucose fastest) plus amino-acid supplements at
  higher rates. The single round-number row at **µ = 1.0 dbl/h is specifically
  glycerol-minimal** (per the paper). The values apply to glucose-minimal only
  insofar as **composition depends on growth rate, not medium** — the Maaløe
  thesis, which this paper itself partly qualifies (different media at equal
  growth rate can differ). Use the value at the **matched growth rate**, and read
  it as growth-rate-conditioned, not glucose-specific.
- **protein / RNA / DNA only.** B&D measured "mass, RNA, DNA, and nuclei" — it
  does **not** itemize lipid / LPS / murein / glycogen / soluble pool / ions, so
  the unmeasured remainder (`1 − (protein+RNA+DNA)`, ≈ 0.26–0.30 of dry weight)
  is **not** attributed here. Other sources are needed for those components.
- **RNA is total RNA** (rRNA+tRNA+mRNA; B&D's own parameters: 98% stable RNA,
  of which 14% tRNA → ~84% rRNA / ~14% tRNA / ~2% mRNA). Compare like-for-like
  against a model's total-RNA fraction, not mRNA.

## Relation to the ParCa input

This is an **independent validation reference** (model outputs graded *against*
it), distinct from the ParCa **input** `data/flat/dry_mass_composition.tsv`
(`DryMassCompositionSchema`), whose own provenance is *not* a verified direct
copy of this table — the input's values do not match B&D Table 2's dry-mass
fractions, so the input is a separate (wcEcoli-reconstruction) artifact, not
ingested here.
