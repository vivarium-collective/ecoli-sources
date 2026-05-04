"""
Pandera schema for the reference bundle manifest.

The bundle manifest is a flat TSV mapping canonical keys (model-side slot
names) to source paths within the data/ package directory. ParCa-time
consumers in vEcoli resolve each canonical key through this bundle to find
the file (or directory) for that data role. The default reference bundle
ships at ``data/reference_bundle.tsv``; alternative bundles (variants)
follow the same shape.
"""

import pandera.pandas as pa


ReferenceBundleSchema = pa.DataFrameSchema(
    name="reference_bundle",
    columns={
        "canonical_key": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description=(
                "Slot name for an addressable data role in the model "
                "(snake_case, primary key). Stable across file renames."
            ),
        ),
        "source_path": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Path to source file or directory, relative to the bundle "
                "root (i.e., the ``data/`` directory of the ecoli-sources "
                "package). Files and directories are both allowed; the "
                "consumer (vEcoli) determines whether to load directly or "
                "traverse based on the canonical key's expected shape."
            ),
        ),
        "description": pa.Column(
            dtype=str,
            nullable=False,
            description="One-liner explaining the slot's role.",
        ),
        "schema_name": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional: Pandera schema name from ecoli-sources/schemas/ "
                "for content validation (e.g. ``RnaseqTpmTableSchema``). "
                "Empty for canonical keys whose source has no associated "
                "schema yet."
            ),
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Bundle manifest mapping canonical keys to source paths. One row "
        "per canonical key. The default reference bundle "
        "(``data/reference_bundle.tsv``) is shipped with the package; "
        "alternative bundles (variants) follow the same shape."
    ),
)
