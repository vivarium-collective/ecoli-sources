# Varma, Boesch & Palsson (1993) — max theoretical biomass yield on glucose

- **source_id:** `varma_1993`
- **Citation:** Varma A, Boesch BW, Palsson BO. "Biochemical production
  capabilities of *Escherichia coli*." *Biotechnology and Bioengineering*
  42(1):59–73, 1993. (BibTeX in `../../references.bib`.)
- **Raw:** `1993-Biochemical_Production_Capabilities_of_Escherichia_coli.pdf`
  in this directory (local, gitignored).

## Claim used

| condition · observable | value | units | kind | method |
|---|---|---|---|---|
| basal · biomass yield on glucose | 0.538 | gDW / g glucose | `theoretical_max` | FBA (LP flux balance) |

## Provenance

**Table III** — *"Maximum theoretical yield of amino acids and nucleotides on a
glucose substrate in absence of constant maintenance energy requirement"*,
**Biomass** row:

- Maximum yield **0.097**. The column header reads "mol/mol Glc", but **footnote
  *d* states "Biomass yield is in g DW/mmol Glc"** — so the biomass figure is
  **0.097 gDW / mmol glucose** (the mol/mol header applies to the amino-acid /
  nucleotide rows, not biomass).
- CO₂ evolved 1.910 mol/mol Glc; constraint E + S (energy + stoichiometry).

Unit conversion to match the model's `Yxs = μ / (q_glc · M_glc)`:

```
0.097 gDW/mmol ÷ 0.180156 g/mmol  =  0.538 gDW / g glucose
```

## What this source does and does not establish

- A **maximal theoretical** (linear-programming flux-balance) yield of the
  *E. coli* metabolic network — an **upper bound**, not a measured culture yield.
- Computed **without a constant maintenance-energy requirement** (table title)
  and under **fully aerobic** conditions → a deliberately generous stoichiometric
  ceiling. A real *measured* aerobic yield is lower (~0.4–0.5 gDW/g; a separate
  `measured` row to follow).
- Use: the `theoretical_max` rung for basal biomass yield. **A model output above
  0.538 gDW/g exceeds the network's stoichiometric ceiling — a first-principles
  violation**, distinct from merely deviating from a measured value.

**Verbatim** (Results → *Maximal Theoretical Performance* → *Glucose*):
> "The maximal conversion of glucose into amino acids and nucleotides is listed
> in Table III. In addition, the maximal biomass yield under fully aerobic
> conditions is shown. These values, determined without including the constant
> maintenance energy, represent the maximal theoretical stoichiometric
> production capability of the metabolic network."
