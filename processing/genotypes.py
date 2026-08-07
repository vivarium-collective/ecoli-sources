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
import warnings
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
    "knockout",
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


# ---------------------------------------------------------------------------
# Flat-file editing, line-oriented
# ---------------------------------------------------------------------------
#
# The ParCa flat files are JSON-typed TSV: every scalar is JSON-encoded
# ("EG10526"), lists are JSON arrays, an absent value is the literal `null` or
# an empty field, and each file opens with `#` provenance comment lines.
#
# Nothing here re-serializes a file. Rows are held as their original text and
# only the FIELDS that change are rewritten, so every untouched row passes
# through byte-identical BY CONSTRUCTION rather than by a correct writer's good
# luck. That property is load-bearing: a variant that silently renormalized
# ~4,700 untouched rows (`null` -> empty, quoting stripped, list columns
# reformatted, `#` headers dropped) would still validate, still produce a
# stable genotype id, and only fail at ParCa build time.
#
# Splitting a data line on tabs is safe because the fields are JSON: a literal
# tab inside a JSON string is encoded `\t`, two characters, so it cannot appear
# raw. Field count is checked against the header anyway.


def _decode(field_text: str) -> Any:
    """Decode one JSON-typed field. An empty field means absent."""
    if field_text == "":
        return None
    return json.loads(field_text)


class _Row:
    """One data line, edited field-wise, rendered from its original text unless
    a field actually changed."""

    __slots__ = ("_raw", "_columns", "_fields", "_dirty", "removed")

    def __init__(self, raw: str, columns: list[str], path: Path, line_no: int):
        self._raw = raw
        self._columns = columns
        self._fields = raw.split("\t")
        if len(self._fields) != len(columns):
            raise ValueError(
                f"{path.name} line {line_no}: {len(self._fields)} fields for "
                f"{len(columns)} columns. Tab-splitting a JSON-typed TSV assumes "
                "no field contains a raw tab; this row breaks that assumption."
            )
        self._dirty = False
        self.removed = False

    def get(self, column: str) -> Any:
        if column not in self._columns:
            return None
        return _decode(self._fields[self._columns.index(column)])

    def set(self, column: str, value: Any) -> None:
        i = self._columns.index(column)
        encoded = json.dumps(value)
        if encoded != self._fields[i]:
            self._fields[i] = encoded
            self._dirty = True

    def render(self) -> str:
        return "\t".join(self._fields) if self._dirty else self._raw


class _FlatFile:
    """A JSON-typed TSV held as lines, with its comment header preserved."""

    def __init__(self, path: PathLike):
        self.path = Path(path)
        text = self.path.read_text(encoding="utf-8")
        self.trailing_newline = text.endswith("\n")
        lines = text.split("\n")
        if self.trailing_newline:
            lines.pop()

        header_i = next(
            (i for i, line in enumerate(lines) if not line.lstrip().startswith("#")),
            None,
        )
        if header_i is None:
            raise ValueError(f"{self.path} has no header row")

        self.preamble = lines[: header_i + 1]
        self.columns = [_decode(f) for f in lines[header_i].split("\t")]
        self.rows = [
            _Row(line, self.columns, self.path, header_i + 1 + n)
            for n, line in enumerate(lines[header_i + 1 :])
            if line != ""
        ]

    def write(self, path: Path) -> Path:
        out = self.preamble + [r.render() for r in self.rows if not r.removed]
        text = "\n".join(out) + ("\n" if self.trailing_newline else "")
        path.write_text(text, encoding="utf-8")
        return path


def _read_fasta(path: Path) -> tuple[str, str, int]:
    """Return (header line, sequence, line width) of a single-record FASTA.

    Hand-rolled rather than biopython: ``SeqIO.write`` re-wraps at 60 columns
    and this file is wrapped at 70, which would rewrite all 66,311 lines of a
    4.5 MB file and destroy the untouched-content property for the largest key.
    """
    lines = path.read_text(encoding="utf-8").split("\n")
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"{path} is not a FASTA (no '>' header)")
    header = lines[0]
    seq_lines = [line for line in lines[1:] if line != ""]
    if any(line.startswith(">") for line in seq_lines):
        raise ValueError(f"{path} holds more than one record; expected a single genome")
    width = len(seq_lines[0]) if seq_lines else 70
    return header, "".join(seq_lines), width


def _write_fasta(path: Path, header: str, sequence: str, width: int) -> Path:
    body = "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))
    path.write_text(f"{header}\n{body}\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The chromosome-deletion transform
# ---------------------------------------------------------------------------
#
# COORDINATE CONVENTION: left_end_pos / right_end_pos are 1-based and
# INCLUSIVE. A feature spanning [L, R] occupies sequence[L - 1 : R] and has
# length R - L + 1. Getting this wrong is SILENT -- feature lengths stay
# correct under an off-by-one while the sequence shifts underneath them.

_SEQUENCE_KEY = "sequence"
_GENES_KEY = "genes"
_DNA_SITES_KEY = "dna_sites"
_TU_KEYS = (
    "transcription_units",
    "transcription_units_added",
    "transcription_units_modified",
)
KNOCKOUT_KEYS = (_SEQUENCE_KEY, _GENES_KEY, _DNA_SITES_KEY) + _TU_KEYS

# The transcription-unit OVERLAYS are why this is six keys and not four.
# v2ecoli runs its deletion after _join_data/_modify_data, so it sees merged
# structures; a generator reads the raw flat files, where the overlays are
# separate canonical keys. `transcription_units_modified` carries real
# coordinates that `_modify_data` applies OVER the base, so omitting it leaves
# those rows in pre-deletion coordinate space and the stale values win.
#
# Not handled, deliberately: the conditional families `rrna_options__*` and
# `new_gene_data__*` apply only when those options are enabled, which is off by
# default and out of scope here.


def _has_coordinates(row: _Row) -> bool:
    """True when a row carries usable left AND right coordinates.

    Tolerant of a half-populated row, which would otherwise raise on
    comparison, and of the two spellings of absent that the shipped data uses:
    `null` in ``genes``, an empty field in ``transcription_units_added``.
    """
    left, right = row.get("left_end_pos"), row.get("right_end_pos")
    return left is not None and right is not None and left != "" and right != ""


def _classify_against_deletion(
    left: int, right: int, del_left: int, del_right: int
) -> str:
    """Classify a feature's position relative to a deletion.

    TOTAL over all (left <= right, del_left <= del_right): every input returns
    exactly one of six labels, so callers need no fallthrough case.

        'before'         entirely upstream; unaffected
        'after'          entirely downstream; shifts left by del_len
        'contained'      wholly inside the deletion; loses its position
        'spans'          starts before and ends after; loses its middle
        'overlaps_left'  starts before, ends inside; truncated at the cut
        'overlaps_right' starts inside, ends after; 5' portion removed
    """
    if right < del_left:
        return "before"
    if left > del_right:
        return "after"
    if left >= del_left and right <= del_right:
        return "contained"
    if left < del_left and right > del_right:
        return "spans"
    if left < del_left:
        return "overlaps_left"
    return "overlaps_right"


def _annotate_removed(row: _Row, del_gene_id: str) -> None:
    """Record that a feature lost content to a deletion.

    Applied only to features truncated or stripped of a member gene -- never to
    ones that merely shifted, which would tag everything downstream of the cut.
    Idempotent per deleted gene.
    """
    marker = f"_removed_{del_gene_id}"
    previous = row.get("common_name") or ""
    if previous.endswith(marker):
        return
    row.set("common_name", previous + marker)


def _update_coordinates(
    rows: list[_Row], kind: str, del_gene_id: str, del_left: int, del_right: int
) -> None:
    """Shift, truncate or retire every feature affected by one deletion.

    Rows are flagged rather than deleted from the list, so nothing is ever
    mutated mid-iteration -- the defect that silently skips the element after
    each removal.
    """
    del_len = del_right - del_left + 1

    for row in rows:
        if row.removed or not _has_coordinates(row):
            continue

        left, right = row.get("left_end_pos"), row.get("right_end_pos")
        if left > right:
            raise ValueError(
                f"{kind} {row.get('id')!r} has left_end_pos {left} > "
                f"right_end_pos {right}"
            )

        case = _classify_against_deletion(left, right, del_left, del_right)

        if case == "before":
            continue

        if case == "contained":
            warnings.warn(
                f"{row.get('id')!r} is contained within the deletion of "
                f"{del_gene_id} and loses its position with it."
            )
            if kind == "gene":
                # D16, extended: any GENE that loses its position is
                # TOMBSTONED, never removed. `rnas` and `proteins` reference
                # genes by id and this transform does not touch them, so
                # dropping the row leaves a dangling cistron -- exactly the
                # failure the tombstone exists to prevent. A nulled gene is
                # already excluded from the transcribed and translated sets:
                # 21 of 4,747 shipped genes ship in that state permanently.
                row.set("left_end_pos", None)
                row.set("right_end_pos", None)
            else:
                # Nothing references a TU or DNA site by id the way `rnas`
                # references genes, so removal creates no dangling reference.
                row.removed = True
            continue

        if case == "after":
            row.set("left_end_pos", left - del_len)
            row.set("right_end_pos", right - del_len)
            continue

        # The remaining cases all LOSE sequence to the deletion.
        if case == "overlaps_left":
            # Starts before, ends inside: truncate at the cut. The kept portion
            # lies entirely upstream, so it does not shift.
            new_left, new_right = left, del_left - 1
        elif case == "overlaps_right":
            # Starts inside, ends after: the surviving 3' portion begins where
            # the deletion used to start.
            new_left, new_right = del_left, right - del_len
        else:  # 'spans'
            # Starts before, ends after: keeps both flanks, loses del_len.
            new_left, new_right = left, right - del_len

        row.set("left_end_pos", new_left)
        row.set("right_end_pos", new_right)
        if kind in ("tu", "dna_site"):
            _annotate_removed(row, del_gene_id)


def _delete_one_gene(
    files: dict[str, _FlatFile], sequence: str, gene_id: str
) -> tuple[str, int]:
    """Splice one gene out of the chromosome and update every coupled key.

    Returns the new sequence and the deleted length.
    """
    gene_row = next((r for r in files[_GENES_KEY].rows if r.get("id") == gene_id), None)
    if gene_row is None:
        raise KeyError(
            f"cannot delete {gene_id!r}: no such gene in the knowledge base. A "
            "knockout of a gene the flat files do not carry would silently "
            "produce a variant identical to its base."
        )
    if not _has_coordinates(gene_row):
        raise KeyError(
            f"cannot delete {gene_id!r}: it has no coordinates on the "
            "chromosome (it may be a phantom gene, or already knocked out)."
        )

    del_left, del_right = gene_row.get("left_end_pos"), gene_row.get("right_end_pos")
    del_len = del_right - del_left + 1

    # 1-based inclusive [L, R] is sequence[L - 1 : R], so the flanks to KEEP
    # are [:L - 1] and [R:]. An off-by-one here preserves every feature LENGTH
    # while corrupting the start codon of the next gene downstream.
    spliced = sequence[: del_left - 1] + sequence[del_right:]
    if len(spliced) != len(sequence) - del_len:
        raise AssertionError(
            f"splicing {gene_id} removed {len(sequence) - len(spliced)} bases, "
            f"expected {del_len}"
        )

    # Detach the gene from every transcription unit that carries it, across the
    # base table and both overlays. A TU left with no genes has no meaningful
    # coordinates.
    for key in _TU_KEYS:
        for tu in files[key].rows:
            members = tu.get("genes") or []
            if gene_id not in members:
                continue
            if len(members) == 1:
                if _has_coordinates(tu):
                    tu.set("left_end_pos", None)
                    tu.set("right_end_pos", None)
            else:
                tu.set("genes", [g for g in members if g != gene_id])
                _annotate_removed(tu, gene_id)

    # Null the target's own coordinates BEFORE the coordinate pass: that is
    # what keeps it out of the `contained` branch below, since it then exits at
    # the no-coordinates guard.
    gene_row.set("left_end_pos", None)
    gene_row.set("right_end_pos", None)

    _update_coordinates(files[_GENES_KEY].rows, "gene", gene_id, del_left, del_right)
    for key in _TU_KEYS:
        _update_coordinates(files[key].rows, "tu", gene_id, del_left, del_right)
    _update_coordinates(
        files[_DNA_SITES_KEY].rows, "dna_site", gene_id, del_left, del_right
    )

    return spliced, del_len


def knockout(
    gene_ids: Iterable[str],
    out_dir: PathLike,
    *,
    base_manifest: PathLike = None,
) -> GeneratorResult:
    """Delete a gene set from the chromosome.

    A knockout is genome surgery, and the hard generator: it writes SIX coupled
    canonical keys that must stay mutually consistent -- the spliced genome,
    the gene table, the DNA sites, and all three transcription-unit tables --
    where a knockdown rewrites one expression table and nothing else.

    What it does NOT write, deliberately: a single gene id is referenced in
    9-10 flat files, so 7-8 dangling references survive every knockout by
    design (`rnas`, `fold_changes`, `translation_efficiency`, the rna_seq
    tables, ...). Pruning them would cascade rnas -> proteins -> complexation
    -> reactions, which carries its own biology decisions; a half-pruned
    knowledge base is a worse failure mode than a stale reference. What those
    dangling references do to a build is what the descriptive studies exist to
    characterize.

    Multi-gene knockouts apply in DESCENDING ``left_end_pos``, so each deletion
    sits above every remaining target and nothing a later step needs has moved.
    That is a correctness requirement, not tidiness: the genotype id is a
    content hash over the produced files, so caller-order-dependent bytes would
    make ``{trpR, tnaA}`` and ``{tnaA, trpR}`` two genotypes -- and two ParCa
    builds -- for one strain.

    Args:
        gene_ids: genes to delete. Order is irrelevant; the canonical order is
            imposed internally.
        out_dir: where to write the variant files (~5.5 MB per genotype).
        base_manifest: bundle to read the base flat files from.

    Returns:
        A :class:`GeneratorResult` overriding the six keys in
        :data:`KNOCKOUT_KEYS`.

    Raises:
        KeyError: a gene is absent from ``genes``, or already has no
            coordinates.
        ValueError: two targets overlap, or a row has left > right.
    """
    gene_ids = list(gene_ids)
    if not gene_ids:
        raise ValueError("knockout() needs at least one gene id")
    duplicates = sorted({g for g in gene_ids if gene_ids.count(g) > 1})
    if duplicates:
        raise ValueError(
            f"repeated gene ids in a single knockout: {duplicates}. The second "
            "deletion of a gene would fail on its nulled coordinates."
        )

    base = read_bundle(base_manifest)
    paths = dict(zip(base["canonical_key"], base["source_path"]))
    missing = [k for k in KNOCKOUT_KEYS if k not in paths]
    if missing:
        raise KeyError(f"base bundle is missing canonical keys {missing}")

    files = {k: _FlatFile(paths[k]) for k in KNOCKOUT_KEYS if k != _SEQUENCE_KEY}
    fasta_path = Path(paths[_SEQUENCE_KEY])
    header, sequence, width = _read_fasta(fasta_path)

    # Resolve every target's coordinates up front, against the UNMODIFIED
    # table, so overlap detection and ordering both see one coordinate space.
    targets = {}
    for gene_id in gene_ids:
        row = next((r for r in files[_GENES_KEY].rows if r.get("id") == gene_id), None)
        if row is None:
            raise KeyError(
                f"cannot delete {gene_id!r}: no such gene in the knowledge base. "
                "A knockout of a gene the flat files do not carry would silently "
                "produce a variant identical to its base."
            )
        if not _has_coordinates(row):
            raise KeyError(
                f"cannot delete {gene_id!r}: it has no coordinates on the "
                "chromosome (it may be a phantom gene, or already knocked out)."
            )
        targets[gene_id] = (row.get("left_end_pos"), row.get("right_end_pos"))

    ordered = sorted(targets, key=lambda g: targets[g][0], reverse=True)

    # Overlapping targets raise. No application order deletes two overlapping
    # genes cleanly -- the second target's coordinates have already been
    # truncated or nulled by the first -- and 926 overlapping gene pairs exist
    # in the shipped data, 60 of them fully nested, so this is reachable rather
    # than theoretical. Interval-merging is the fix if a pairwise screen ever
    # needs it.
    for earlier, later in zip(ordered, ordered[1:]):
        hi_left, hi_right = targets[earlier]
        lo_left, lo_right = targets[later]
        if lo_right >= hi_left and hi_right >= lo_left:
            raise ValueError(
                f"knockout targets {later!r} ({lo_left}-{lo_right}) and "
                f"{earlier!r} ({hi_left}-{hi_right}) overlap. Deleting "
                "overlapping genes in one genotype is not supported: no order "
                "removes both cleanly."
            )

    for gene_id in ordered:
        sequence, _ = _delete_one_gene(files, sequence, gene_id)

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = {
        _SEQUENCE_KEY: _write_fasta(
            out_dir / "sequence__ko.fasta", header, sequence, width
        )
    }
    for key, flat in files.items():
        rows[key] = flat.write(out_dir / f"{key}__ko.tsv")

    return GeneratorResult(
        rows=rows,
        provenance={
            "generator": "knockout",
            "operator": "chromosome_deletion",
            "params": {
                # Recorded in the canonical APPLIED order, not as given, so the
                # sidecar is order-independent alongside the files it describes.
                "gene_ids": ordered,
                "coordinates": {g: list(targets[g]) for g in ordered},
            },
            "reads": {k: str(paths[k]) for k in KNOCKOUT_KEYS},
        },
    )
