# References (provenance layer)

One directory per source, keyed by `source_id` — the same key the curated
claim tables (`../data/<condition>/<observable>.tsv`) reference.

```
references/
  <source_id>/
    summary.md     # what the source measured/derived, conditions, caveats,
                   # and verbatim passage(s) supporting each extracted value
    raw/           # extracted figure/table data as published (csv/xlsx/etc.)
```

`source_id` convention: `firstauthor_year` (e.g. `basan_2015`), matching the
BibTeX key in `../references.bib`.

`summary.md` is the human/agent-facing record that makes a curated number
auditable: it should name the figure/table the value came from, the strain and
growth conditions, any unit conversions applied, and a short verbatim quote of
the claim. Keep it honest about what the source does and does not establish.
