"""
Post-processing utilities for formatted RNAseq datasets.
"""

import os

import pandas as pd

from schemas import RnaseqTpmTableSchema, RnaseqSamplesManifestSchema


def make_gene_exclusion_variant(
    manifest: pd.DataFrame,
    source_dataset_id: str,
    genes_to_remove: list[str],
    new_dataset_id: str,
    data_dir: str,
    new_description: str | None = None,
    manifest_filename: str = "manifest.tsv",
) -> pd.DataFrame:
    """
    Create a variant of an existing dataset with specified genes removed.

    Loads the source dataset TSV, drops the given gene IDs, writes the new TSV,
    appends a row to the manifest, and writes the manifest. TPM values are not
    renormalized after removal.

    Args:
        manifest: Current manifest DataFrame.
        source_dataset_id: dataset_id of the dataset to derive from.
        genes_to_remove: List of gene_id values to exclude.
        new_dataset_id: dataset_id for the new variant.
        data_dir: Directory containing TSV files and the manifest.
        new_description: Human-readable description for the new dataset.
            Defaults to "<source description> (excluding <N> genes)".
        manifest_filename: Manifest filename within data_dir. Default "manifest.tsv".

    Returns:
        Updated manifest DataFrame (with new row appended).

    Raises:
        ValueError: If source_dataset_id is not found in the manifest, or if
            new_dataset_id already exists in the manifest.
    """
    if source_dataset_id not in manifest["dataset_id"].values:
        raise ValueError(f"source_dataset_id '{source_dataset_id}' not found in manifest")
    if new_dataset_id in manifest["dataset_id"].values:
        raise ValueError(f"new_dataset_id '{new_dataset_id}' already exists in manifest")

    source_row = manifest.loc[manifest["dataset_id"] == source_dataset_id].iloc[0]

    # Load source TPM table
    source_path = os.path.join(data_dir, source_row["file_path"])
    tpm = pd.read_csv(source_path, sep="\t")

    # Remove genes
    genes_to_remove_set = set(genes_to_remove)
    removed = tpm[tpm["gene_id"].isin(genes_to_remove_set)]
    tpm_filtered = tpm[~tpm["gene_id"].isin(genes_to_remove_set)].copy()
    n_removed = len(removed)
    n_not_found = len(genes_to_remove_set - set(removed["gene_id"]))

    if n_not_found:
        missing = genes_to_remove_set - set(removed["gene_id"])
        print(f"Warning: {n_not_found} gene(s) in genes_to_remove not found in dataset: {missing}")

    # Validate and write new TSV
    tpm_validated = RnaseqTpmTableSchema.validate(tpm_filtered)
    out_path = os.path.join(data_dir, f"{new_dataset_id}.tsv")
    tpm_validated.to_csv(out_path, sep="\t", index=False)

    # Build manifest row from source, overriding id/path/description
    if new_description is None:
        new_description = f"{source_row['dataset_description']} (excluding {n_removed} genes)"

    new_row = source_row.to_dict()
    new_row["dataset_id"] = new_dataset_id
    new_row["file_path"] = f"{new_dataset_id}.tsv"
    new_row["dataset_description"] = new_description

    manifest_updated = pd.concat(
        [manifest, pd.DataFrame([new_row])], ignore_index=True
    )
    manifest_validated = RnaseqSamplesManifestSchema.validate(manifest_updated)
    manifest_validated.to_csv(
        os.path.join(data_dir, manifest_filename), sep="\t", index=False
    )

    print(f"Created '{new_dataset_id}': removed {n_removed} genes. Wrote {out_path}")

    return manifest_validated
