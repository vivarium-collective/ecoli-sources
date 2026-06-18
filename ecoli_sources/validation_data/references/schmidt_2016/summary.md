# Schmidt et al. (2016) — *E. coli* MG1655 glucose-minimal proteome (absolute copies/cell)

- **source_id:** `schmidt_2016`
- **Citation:** Schmidt A, Kochanowski K, Vedelaar S, Ahrné E, Volkmer B,
  Callipo L, Knoops K, Bauer M, Aebersold R, Heinemann M. "The quantitative and
  condition-dependent *Escherichia coli* proteome." *Nat. Biotechnol.*
  34(1):104–110, 2016. DOI 10.1038/nbt.3418. (BibTeX in `../../references.bib`.)
  (Published online Dec 2015; the v2ecoli/wcEcoli lineage refers to it as
  "schmidt2015".)
- **Raw (local, gitignored):** the PDF and `41587_2016_BFnbt3418_MOESM18_ESM.xlsx`
  (Supplementary Tables S1–S29).

## What is ingested (`data/basal/proteome.tsv`, `source_id=schmidt_2016`)

The **wild-type MG1655, glucose-minimal-medium** proteome: **2015 proteins**,
absolute **copies/cell**, keyed by gene name (`gene`; `uniprot_id` carried
alongside). `cv` is the reported coefficient of variation across biological
triplicates. Source: Supplementary **Table S9**, column
`Copies/Cell_MG1655.Glucose` (with `cv_MG1655.Glucose`).

**Strain choice.** Schmidt's main 22-condition abundance map (Supplementary
Table S6) is for **BW25113** (the Keio parent). Table S9 additionally quantifies
**MG1655** and NCM3722 at glucose and LB by the same method. We take the
**MG1655 glucose** column because it is **strain-matched to the model**; the
BW25113 glucose column (S6, ~2359 proteins — the basis of the legacy v2ecoli
`schmidt2015_javier_table` reference) and S9's same-table `BW25113.Glucose` are
available as alternatives / a strain cross-check.

## Measurement design (Materials & Methods)

- **Method:** MS-based proteomics. Quantitative proteome maps by LC-MS/MS
  (label-free, OGE/OFFGEL-fractionated and unfractionated), anchored to
  **absolute** quantification of 41 glycolysis/TCA enzymes by stable-isotope
  dilution + selected-reaction monitoring (SRM/SID); whole-proteome absolute
  abundances by iBAQ. Per-cell copies obtained by combining concentrations with
  **flow-cytometry cell counts and cell volumes**. Biological triplicates.
- **Condition:** M9 minimal medium with glucose as carbon source, aerobic,
  37 °C, exponential growth. (Per-condition growth rates, ODs, and harvest
  details are in Supplementary Table S23.)
- Quantified ~2359 proteins overall (~55% of predicted ORFs, >95% of proteome
  mass) across 22 conditions; the MG1655 glucose subset ingested here covers
  ~2015 genes.

## Notes

- **Strain-matched** (MG1655) — the cleanest strain match among the metabolic
  references here (Toya/Crown differ on strain/lineage details).
- Absolute copies/cell — directly comparable to a model's per-cell protein
  counts (typically graded as a log-log R² over shared genes).
- A few proteins report 0 copies/cell in this condition (below detection); kept
  as reported.
- **Sibling proteome source:** Wisniewski et al. 2014 (a "proteomic ruler"
  copy-number estimate) is the other proteome reference in the v2ecoli lineage
  and can be added to this slot as an additional `source_id` the same way.
