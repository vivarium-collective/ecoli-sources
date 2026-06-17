# Varma & Palsson (1994) — measured aerobic physiology of E. coli W3110 on glucose

- **source_id:** `varma_1994`
- **Citation:** Varma A, Palsson BO. "Stoichiometric flux balance models
  quantitatively predict growth and metabolic by-product secretion in wild-type
  *Escherichia coli* W3110." *Appl. Environ. Microbiol.* 60(10):3724–3731, 1994.
  (BibTeX in `../../references.bib`.)
- **Raw:** `aem00027-0254.pdf` in this directory (local, gitignored).

## Claims used (aerobic batch, strain W3110)

| condition · observable | value | units | kind | method |
|---|---|---|---|---|
| basal · growth rate (μ) | 0.68 | 1/h | `measured` | aerobic batch |
| basal · glucose uptake (q_glc) | 10.5 | mmol/gDW/h | `measured` | aerobic batch |
| basal · biomass yield (Yxs) | 0.355 | gDW/g glucose | `measured` | aerobic batch |

## Provenance

**Fig. 2** caption: *"Maximum aerobic glucose utilization rate (10.5 mmol of Glc
per g [dry weight] per h) determined from batch experiments as the growth rate
(0.68 h⁻¹; r²=0.999) divided by the biomass yield (0.064 g [dry weight] per mmol
of Glc; r²=0.995)."*

- μ = **0.68 h⁻¹**, q_glc = **10.5 mmol/gDW/h**, Yxs = **0.064 gDW/mmol glc**.
- Yxs to gDW/g: `0.064 / 0.180156 = 0.355 gDW/g glucose`.
- Internal consistency: `q_glc = μ / Yxs = 0.68 / 0.064 = 10.6 ≈ 10.5` ✓.
- (Fig. 1 also gives max aerobic O₂ uptake **15 mmol/gDW/h** — not ingested here.)

## Caveats

- **Strain W3110** — the v2ecoli model is **MG1655**. Close but not identical;
  noted on every row.
- **Unit-reading caution:** the **0.43 h⁻¹** that appears in this paper is the
  **anaerobic** growth rate (Fig. 3: anaerobic μ=0.43, Yxs=0.023 gDW/mmol,
  q_glc=18.5), *not* an aerobic yield. The aerobic yield is 0.355 gDW/g.
- Aerobic batch near μmax → some acetate overflow is expected, which lowers the
  realized yield below the stoichiometric ceiling (Varma 1993, 0.538 gDW/g).
