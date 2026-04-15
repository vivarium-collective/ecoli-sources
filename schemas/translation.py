"""
Pandera schema for per-gene translation efficiency measurements.
"""

import pandera.pandas as pa


TranslationEfficiencySchema = pa.DataFrameSchema(
    name="translation_efficiency",
    columns={
        "geneId": pa.Column(
            str, nullable=True,
            description=(
                "Gene id (e.g. 'G7686'). Nullable because the source file "
                "encodes unmatched genes as '#N/A', which pandas coerces to "
                "NaN under default read_csv settings."
            ),
        ),
        "name": pa.Column(
            str, nullable=True, required=False,
            description="Gene symbol (free-form; may be absent).",
        ),
        "translationEfficiency": pa.Column(
            str, nullable=True,
            description=(
                "Relative translation efficiency (ribosome loading / mRNA). "
                "Stored as string because the source file uses 'NA' for "
                "missing values; coerce downstream."
            ),
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Per-gene translation efficiency (Li et al. 2014 Table S4). "
        "'NA' strings are legal and must be handled by the consumer."
    ),
)
