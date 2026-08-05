"""
Genotype variant generators and bundle composition.

A **genotype** is a bundle: the manifest a ParCa build reads *is* the genome
and dataset it was built from, so perturbing the model's inputs means writing
a variant manifest rather than passing a flag to the consumer.

This module provides the substrate every genotype generator shares, plus the
first generator (knockdown). The contract:

    generator(...)                  -> GeneratorResult(rows, provenance)
    compose_variant_bundle(results)  -> a complete, validated variant manifest
    genotype_id(manifest)            -> the content hash identifying that genotype

Design notes, because two of these are easy to get wrong:

* **Generators return ROWS, not manifests.** A generator writes its variant
  files and reports which canonical keys they replace; a separate composer
  assembles the complete manifest. This is what makes composition work — a
  knockout combined with a knockdown is two results composed, not two complete
  manifests diffed against each other. The consumer still receives one complete
  manifest and keeps no merge logic, per ``BUNDLES.md``: *"Manifest construction
  (variant generators, diff tools) is a build-time concern outside the runtime
  contract."*

* **Identity is a CONTENT hash, not the manifest's bytes.** A manifest is a list
  of pointers. Editing a variant file in place leaves the manifest byte-identical,
  so a bytes-hash would silently reuse a stale ParCa build — precisely the failure
  the identity exists to prevent.

* **Materialization is cache-like, not archival.** The durable artifact of a
  genotype is its spec and its id; the files are a reproducible derivative that
  may be evicted and rebuilt. A knockout writes ~5.5 MB (the genome sequence is
  4.5 MB of it), so keeping every variant of a thousand-genotype campaign on disk
  is ~5.5 GB for data that is fully determined by (base id, generator, params).
  The gitignored ``rnaseq_experimental/perturbations/`` tree already works this
  way; genotype variants follow it.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

from schemas import ReferenceBundleSchema

__all__ = [
    "GeneratorResult",
    "GENERATOR_VERSION",
    "compose_variant_bundle",
    "genotype_id",
    "knockdown",
    "read_bundle",
]

# Bumped when a generator's output changes for identical inputs. Recorded in
# each genotype's sidecar so a stale materialization is attributable.
GENERATOR_VERSION = "1"

_MANIFEST_NAME = "reference_bundle.tsv"
_SIDECAR_NAME = "genotype.json"


def _default_base_manifest() -> Path:
    from ecoli_sources import BUNDLE_PATH

    return Path(BUNDLE_PATH)


def read_bundle(manifest_path: PathLike = None) -> pd.DataFrame:
    """Read a bundle manifest, with each ``source_path`` resolved to absolute.

    Paths in a manifest are relative to *the manifest's own directory* (that is
    how ``SourceBundle`` resolves them), so a manifest cannot be copied to
    another directory without rewriting them. Resolving on read and re-relativizing
    on write is what makes :func:`compose_variant_bundle` location-independent.
    """
    manifest_path = Path(manifest_path or _default_base_manifest()).resolve()
    root = manifest_path.parent
    df = pd.read_csv(manifest_path, sep="\t", comment="#")
    df["source_path"] = [str((root / p).resolve()) for p in df["source_path"]]
    return df


PathLike = str | os.PathLike | None


@dataclass(frozen=True)
class GeneratorResult:
    """What one generator contributes to a genotype.

    Attributes:
        rows: ``canonical_key -> absolute path`` of the variant file that key
            should resolve to. Generators deal in absolute paths; the composer
            owns relativization.
        provenance: free-form record of how the files were produced — generator
            name, params, and the inputs read. Lands in the genotype sidecar.
    """

    rows: Mapping[str, Path]
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        missing = [k for k, v in self.rows.items() if not Path(v).is_file()]
        if missing:
            raise FileNotFoundError(
                "GeneratorResult references files that do not exist: "
                + ", ".join(f"{k} -> {self.rows[k]}" for k in missing)
            )


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def genotype_id(manifest_path: PathLike = None, *, short: int = 16) -> str:
    """Content-addressed identity of the genotype a manifest describes.

    Hashes ``(canonical_key, sha256(file contents))`` for every row, in sorted
    key order — NOT the manifest's own bytes. Consequences that matter:

    * two genotypes reached by different routes collide correctly, so an already
      built one is reused rather than rebuilt;
    * touching any input invalidates the id automatically;
    * the id is stable across machines and across re-materialization, because it
      depends on no path, mtime, or ordering.

    That last property is what makes eviction safe: a rebuilt variant can be
    *verified* to be the same genotype rather than assumed to be.
    """
    df = read_bundle(manifest_path)
    h = hashlib.sha256()
    for key, path in sorted(zip(df["canonical_key"], df["source_path"])):
        h.update(f"{key}\t{_sha256_file(Path(path))}\n".encode())
    digest = h.hexdigest()
    return digest[:short] if short else digest


def compose_variant_bundle(
    results: Iterable[GeneratorResult],
    out_dir: PathLike,
    *,
    base_manifest: PathLike = None,
    name: str | None = None,
    validate: bool = True,
) -> Path:
    """Assemble one complete variant manifest from a base plus generator results.

    Writes ``reference_bundle.tsv`` and a ``genotype.json`` sidecar into
    ``out_dir`` and returns the manifest path.

    Raises:
        ValueError: if two results write the same canonical key, or if a result
            names a key absent from the base bundle.
    """
    results = list(results)
    base_path = Path(base_manifest or _default_base_manifest()).resolve()
    df = read_bundle(base_path)
    known = set(df["canonical_key"])

    # Collisions are a real situation -- a knockout and an insertion both touch
    # `genes` -- and resolving them by iteration order would be a correctness bug
    # wearing the clothes of a dict update. Make the caller decide instead.
    claimed: dict[str, int] = {}
    for i, result in enumerate(results):
        for key in result.rows:
            if key not in known:
                raise ValueError(
                    f"generator result {i} writes canonical key {key!r}, which is "
                    f"not in the base bundle {base_path}. A variant may only "
                    "override keys the contract already defines."
                )
            if key in claimed:
                raise ValueError(
                    f"canonical key {key!r} is written by two generator results "
                    f"({claimed[key]} and {i}). Composition is not last-writer-wins; "
                    "combine them into one generator or compose them in separate "
                    "bundles."
                )
            claimed[key] = i

    overrides = {k: str(Path(p).resolve()) for r in results for k, p in r.rows.items()}
    df["source_path"] = [
        overrides.get(key, path)
        for key, path in zip(df["canonical_key"], df["source_path"])
    ]

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Re-relativize against the NEW manifest location. Skipping this is the
    # silent failure mode: the manifest validates, and every unoverridden key
    # then resolves to a path that does not exist.
    df["source_path"] = [
        os.path.relpath(p, start=out_dir) for p in df["source_path"]
    ]

    if validate:
        # NB: the schema is strict="filter" -- it DROPS unknown columns rather
        # than rejecting them, so provenance must never ride in the manifest.
        # That is what genotype.json is for.
        df = ReferenceBundleSchema.validate(df)

    manifest_path = out_dir / _MANIFEST_NAME
    df.to_csv(manifest_path, sep="\t", index=False)

    sidecar = {
        "name": name,
        "genotype_id": genotype_id(manifest_path),
        "base_manifest": str(base_path),
        "base_genotype_id": genotype_id(base_path),
        "generator_version": GENERATOR_VERSION,
        "generators": [dict(r.provenance) for r in results],
        "overridden_keys": sorted(overrides),
    }
    (out_dir / _SIDECAR_NAME).write_text(json.dumps(sidecar, indent=2, sort_keys=True, default=str))
    return manifest_path


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------


def knockdown(
    gene_ids: Iterable[str],
    factor: float,
    out_dir: PathLike,
    *,
    base_manifest: PathLike = None,
    canonical_key: str = "rnaseq_experimental_tpms",
    renormalize: bool = False,
) -> GeneratorResult:
    """Knock a gene set down by scaling its measured expression.

    A knockdown is a **single-key** perturbation: it rewrites one expression
    table and nothing else. That is the whole reason it is implemented this way
    rather than as genome surgery — there is no natural "partial chromosome", and
    an expression swap composes trivially with every other rnaseq-derived
    generator.

    Reuses the ``scale_gene_set`` operator, deliberately at the OPERATOR layer
    rather than through ``make_perturbation_variant``. That driver appends a row
    to the shared, committed rnaseq dataset manifest — appropriate for curated DI
    datasets, wrong here: a thousand-genotype campaign would mutate one shared
    file a thousand times, and genotype variants are evictable derivatives that
    should leave no trace in a committed registry. Lineage lives in the genotype
    sidecar instead.

    Args:
        gene_ids: genes to knock down.
        factor: multiplier on ``tpm_mean``. ``0.1`` is a 10x knockdown; ``0.0``
            is expression-level ablation (distinct from a chromosomal knockout,
            which removes the gene entirely).
        out_dir: where to write the variant table.
        base_manifest: bundle to read the source expression table from.
        canonical_key: which expression key to perturb.
        renormalize: rescale the table back to 1e6 TPM after scaling. Off by
            default: renormalizing redistributes the freed capacity across every
            other gene, which is a modelling choice the caller should make
            explicitly rather than inherit.

    Returns:
        A :class:`GeneratorResult` overriding ``canonical_key``.
    """
    from processing.perturbations import scale_gene_set

    gene_ids = list(gene_ids)
    base = read_bundle(base_manifest)
    row = base.loc[base["canonical_key"] == canonical_key]
    if row.empty:
        raise KeyError(
            f"canonical key {canonical_key!r} is not in the base bundle; "
            "cannot knock down an expression table that does not exist."
        )
    source_path = Path(row.iloc[0]["source_path"])

    tpm = pd.read_csv(source_path, sep="\t")
    unknown = sorted(set(gene_ids) - set(tpm["gene_id"]))
    if unknown:
        raise KeyError(
            f"gene ids absent from {source_path.name}: {unknown[:10]}"
            f"{'...' if len(unknown) > 10 else ''}. A knockdown of a gene the "
            "expression table does not measure would silently do nothing."
        )

    out = scale_gene_set(tpm, gene_ids, factor, renormalize=renormalize)

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{canonical_key}__kd.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    return GeneratorResult(
        rows={canonical_key: out_path},
        provenance={
            "generator": "knockdown",
            "operator": "scale_gene_set",
            "params": {
                "gene_ids": gene_ids,
                "factor": factor,
                "renormalize": renormalize,
            },
            "reads": {canonical_key: str(source_path)},
        },
    )
