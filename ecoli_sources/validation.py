"""
Overlay-union loader for the validation-data subsystem.

Mirrors the ParCa overlay pattern (vEcoli ``wholecell/io/ingestion.py``'s
``ingest_rnaseq_manifest_with_overlays``, unique on ``dataset_id``) for the
*validation* flow: union a primary validation bundle with any overlay bundles,
unique on ``canonical_key``. The reusable loader lives here; overlay repos ship
only their own ``validation_data/`` and are unioned in at load time.

Two claim-table conventions are resolved transparently to one in-memory band:

* **per-observable file** (public ``ScalarClaimSchema``): one file per
  ``(condition, observable)``, e.g. ``data/basal/o2_uptake.tsv``; the measured
  band is the ``kind == measured`` rows.
* **tidy assay table** (private ``ScalarObservationSchema``): one table holding
  many observables across cultivation groups, e.g. ``fermentation.tsv``; the
  bundle row carries a ``(cultivation_group_id, observable)`` filter selecting
  the slot's rows.

Overlay bundle paths come from ``overlay_paths=`` or the
``$ECOLI_SOURCES_VALIDATION_OVERLAYS`` environment variable (``;`` preferred,
``:`` legacy; cloud ``scheme://`` URIs preserved) — the validation analogue of
``$ECOLI_SOURCES_OVERLAYS``.

Note on naming: ``cultivation_group_id`` here identifies the DERIVED group a
card grades against (e.g. strain × condition, aggregating replicate reactors).
It is distinct from a raw, LabKey-atomic per-reactor ``cultivation_id``, which
this loader does not model — that tier lives upstream, in the assay-ingestion
layer that produces this bundle.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from ecoli_sources import VALIDATION_BUNDLE_PATH

# Scalar claim-table schemas this loader resolves to a measured band.
_SCALAR_SCHEMAS = {"ScalarClaimSchema", "ScalarObservationSchema"}

# Vector claim/observation tables. Resolved to a located slot (path + filter)
# rather than a band: a scalar slot is a handful of numbers a consumer grades
# directly, a vector slot is thousands of rows, so one function returning either
# shape would be convenient to call and unpleasant to consume.
_VECTOR_SCHEMAS = {
    "VectorObservationSchema",       # cultivation-keyed, tidy (this subsystem)
    "VectorReplicateSchema",         # per-sample sibling of the above
    "ReactionFluxSchema",            # per-source vector claims
    "ProteinAbundanceSchema",
    "MetaboliteConcentrationSchema",
    "MacromoleculeCompositionSchema",
}

VALIDATION_OVERLAYS_ENV = "ECOLI_SOURCES_VALIDATION_OVERLAYS"


# ---------------------------------------------------------------------------
# Overlay path parsing (mirrors vEcoli ingestion._split_overlay_separator)
# ---------------------------------------------------------------------------

def _split_overlay_separator(raw: str) -> list[str]:
    """Split a separator-delimited overlay list. ``;`` is preferred (URI-safe);
    ``:`` is the legacy separator, with ``scheme://`` fragments reattached so a
    colon-split does not tear ``s3://bucket/...`` apart."""
    raw = raw.strip()
    if not raw:
        return []
    if ";" in raw:
        return [p.strip() for p in raw.split(";") if p.strip()]
    parts = raw.split(":")
    out: list[str] = []
    for part in parts:
        if out and out[-1].endswith("/") and part.startswith("/"):
            # reattach a <scheme>://… that the colon-split tore apart
            out[-1] = out[-1] + ":" + part
        elif out and out[-1] in ("s3", "gs", "gcs", "http", "https", "file"):
            out[-1] = out[-1] + ":" + part
        else:
            out.append(part)
    return [p.strip() for p in out if p.strip()]


def validation_overlay_paths(env: str = VALIDATION_OVERLAYS_ENV) -> list[Path]:
    """Overlay validation-bundle paths from ``$ECOLI_SOURCES_VALIDATION_OVERLAYS``."""
    return [Path(p) for p in _split_overlay_separator(os.environ.get(env, ""))]


# ---------------------------------------------------------------------------
# Band extraction from one bundle row (both conventions)
# ---------------------------------------------------------------------------

def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def _measured_band(df: pd.DataFrame, *, source_col: str) -> dict:
    """From a claim/observation frame (already filtered to one slot), pull the
    measured band + optional theoretical-max ceiling + provenance."""
    meas = df[df["kind"] == "measured"]
    tmax = df[df["kind"] == "theoretical_max"] if "kind" in df else df.iloc[0:0]
    tmax_row = tmax.loc[tmax["value"].idxmin()] if len(tmax) else None
    unc = (meas["uncertainty"] if "uncertainty" in meas
           else pd.Series([None] * len(meas)))
    return {
        "measured": [float(v) for v in meas["value"]],
        "measured_unc": [(float(u) if pd.notna(u) else None) for u in unc],
        "theoretical_max": (float(tmax_row["value"]) if tmax_row is not None else None),
        "theoretical_source": (str(tmax_row[source_col]) if tmax_row is not None
                               and source_col in tmax_row else None),
        "sources": [str(s) for s in meas.get(source_col, [])],
        "units": (str(meas["units"].iloc[0]) if len(meas) and "units" in meas else None),
    }


def _row_band(row: dict, root: Path) -> dict:
    """Resolve one validation-bundle row to a measured band, handling both the
    per-observable-file and tidy-assay-table conventions."""
    table = _read_tsv(root / str(row["source_path"]))
    cultivation_group_id = row.get("cultivation_group_id")
    observable = row.get("observable")
    is_tidy = (pd.notna(cultivation_group_id) if cultivation_group_id is not None else False) and (
        "observable" in table.columns)
    if is_tidy:
        sub = table[table["cultivation_group_id"].astype(str) == str(cultivation_group_id)]
        if pd.notna(observable):
            sub = sub[sub["observable"].astype(str) == str(observable)]
        band = _measured_band(sub, source_col="replicate")
        band["cultivation_group_id"] = str(cultivation_group_id)
        band["observable"] = (str(observable) if pd.notna(observable) else None)
        # derived-tier temporal provenance (constant across replicates of a slot)
        for col in ("phase", "window"):
            vals = ([v for v in sub[col].dropna().unique()] if col in sub.columns else [])
            band[col] = (str(vals[0]) if vals else None)
    else:
        band = _measured_band(table, source_col="source_id")
        band["cultivation_group_id"] = None
        # observable inferred from the canonical_key tail (``<condition>__<obs>``)
        key = str(row["canonical_key"])
        band["observable"] = key.split("__", 1)[1] if "__" in key else key
    band["canonical_key"] = str(row["canonical_key"])
    band["schema_name"] = str(row.get("schema_name") or "")
    band.setdefault("phase", None)    # per-observable convention carries no temporal tier
    band.setdefault("window", None)
    return band


# ---------------------------------------------------------------------------
# Union loader
# ---------------------------------------------------------------------------

def _union_bundle_rows(bundles: list[Path]) -> dict[str, tuple[dict, Path]]:
    """Union every row of every bundle on ``canonical_key``.

    Returns ``{canonical_key: (row, bundle_root)}``. Uniqueness is enforced
    across the WHOLE union — scalar and vector slots share one key namespace, so
    a vector row cannot silently shadow a scalar one (they address the same
    thing: "this observable, for this condition"). Filtering by schema happens
    after the union, so a collision is caught even between rows that no single
    caller would have looked at together."""
    out: dict[str, tuple[dict, Path]] = {}
    for bundle_path in bundles:
        manifest = _read_tsv(bundle_path)
        root = bundle_path.parent
        for _, r in manifest.iterrows():
            row = r.to_dict()
            key = str(row["canonical_key"])
            if key in out:
                raise ValueError(
                    f"duplicate validation canonical_key '{key}' across the "
                    f"primary+overlay union (second seen in {bundle_path})")
            out[key] = (row, root)
    return out


def load_scalar_observations(
    primary_bundle: Path | str | None = None,
    overlay_paths: list[Path | str] | None = None,
    *,
    include_primary: bool = True,
    prefix: str | None = None,
) -> dict[str, dict]:
    """Union scalar validation slots across the primary + overlay bundles.

    Returns ``{canonical_key: band}`` where ``band`` carries ``measured`` (the
    band), ``measured_unc``, ``theoretical_max``/``theoretical_source``,
    ``sources``, ``units``, ``observable``, ``cultivation_group_id``,
    ``schema_name``. Only scalar slots (``ScalarClaimSchema`` /
    ``ScalarObservationSchema``) are returned; vector slots are skipped.
    ``canonical_key`` must be unique across the union — a collision raises
    ``ValueError`` (mirrors the ParCa overlay ``dataset_id`` contract).
    ``prefix`` optionally filters to keys starting with it (e.g. a
    ``"<cultivation_group_id>__"`` prefix)."""
    if primary_bundle is None:
        primary_bundle = VALIDATION_BUNDLE_PATH
    if overlay_paths is None:
        overlay_paths = validation_overlay_paths()

    bundles: list[Path] = []
    if include_primary:
        bundles.append(Path(primary_bundle))
    bundles.extend(Path(p) for p in overlay_paths)

    out: dict[str, dict] = {}
    for key, (row, root) in _union_bundle_rows(bundles).items():
        if str(row.get("schema_name") or "") not in _SCALAR_SCHEMAS:
            continue
        if prefix is not None and not key.startswith(prefix):
            continue
        out[key] = _row_band(row, root)
    return out


def load_vector_observations(
    primary_bundle: Path | str | None = None,
    overlay_paths: list[Path | str] | None = None,
    *,
    include_primary: bool = True,
    prefix: str | None = None,
) -> dict[str, dict]:
    """Locate vector validation slots across the primary + overlay bundles.

    Returns ``{canonical_key: slot}``, where a slot carries ``path`` (the
    resolved table), ``schema_name``, and the filter that selects this slot's
    rows within it (``cultivation_group_id``, ``observable``, ``units``), plus
    ``phase``/``window``/``description``.

    **Located, not read.** Unlike the scalar loader — which resolves a slot all
    the way to a measured band — this stops at the file. A vector table can hold
    thousands of rows, a bundle can hold many slots, and a consumer typically
    wants one; materializing them all to answer "what exists?" is the wrong
    trade. Call :func:`read_vector_observations` on the slot you want.

    Covers both vector conventions: the cultivation-keyed tidy table
    (``VectorObservationSchema``) and the per-source vector claim tables
    (``ReactionFluxSchema``, ``ProteinAbundanceSchema``, …)."""
    if primary_bundle is None:
        primary_bundle = VALIDATION_BUNDLE_PATH
    if overlay_paths is None:
        overlay_paths = validation_overlay_paths()

    bundles: list[Path] = []
    if include_primary:
        bundles.append(Path(primary_bundle))
    bundles.extend(Path(p) for p in overlay_paths)

    out: dict[str, dict] = {}
    for key, (row, root) in _union_bundle_rows(bundles).items():
        schema = str(row.get("schema_name") or "")
        if schema not in _VECTOR_SCHEMAS:
            continue
        if prefix is not None and not key.startswith(prefix):
            continue
        out[key] = {
            "canonical_key": key,
            "schema_name": schema,
            "path": root / str(row["source_path"]),
            "cultivation_group_id": _opt(row, "cultivation_group_id", "cultivation_id"),
            "observable": _opt(row, "observable"),
            "units": _opt(row, "units"),
            "phase": _opt(row, "phase"),
            "window": _opt(row, "window"),
            "description": _opt(row, "description"),
        }
    return out


def _opt(row: dict, *names: str) -> str | None:
    """First present, non-null value among ``names`` in a bundle row.

    Accepts several spellings so a bundle stays readable across the
    ``cultivation_id`` → ``cultivation_group_id`` rename — the group key is the
    one a vector slot filters on."""
    for name in names:
        if name in row:
            val = row[name]
            if val is not None and pd.notna(val):
                return str(val)
    return None


def read_vector_observations(slot: dict, *, validate: bool = True) -> pd.DataFrame:
    """Read one located vector slot, filtered to its rows.

    Applies whichever of ``cultivation_group_id`` / ``observable`` / ``units``
    the slot pins AND the table carries, so the same call works for a tidy
    cultivation-keyed table and for a per-source claim table that pins none of
    them. Validates against whichever schema the slot names, via
    :func:`_schema_for` (``validate=False`` to skip — e.g. at render time, where
    a payload defect should be surfaced by CI rather than break a report).

    ⚠ ``validate`` does not only decide whether a check runs; the vector schemas
    use ``strict="filter"``, so validating also DROPS columns the schema does not
    declare. A caller passing ``validate=False`` can therefore see columns the
    CI path does not."""
    table = _read_tsv(Path(slot["path"]))
    for col in ("cultivation_group_id", "observable", "units"):
        want = slot.get(col)
        if want is not None and col in table.columns:
            table = table[table[col].astype(str) == str(want)]
    # ⚠ Dispatch on the slot's schema, and keep ``_schema_for`` exhaustive over
    # the schemas this repo defines. Registering a name in ``_VECTOR_SCHEMAS``
    # makes a slot RESOLVABLE; it does not make it VALIDATED. A name added to
    # one and not the other yields a caller that passed ``validate=True`` and
    # silently got nothing checked — see test_vector_replicates.py.
    if validate:
        schema = _schema_for(str(slot.get("schema_name") or ""))
        if schema is not None:
            table = schema.validate(table, lazy=True)
    return table.reset_index(drop=True)


def _schema_for(schema_name: str):
    """The pandera schema for a vector slot, or ``None`` if this repo defines none.

    ``None`` is a real answer: several ``_VECTOR_SCHEMAS`` entries are per-source
    claim tables validated elsewhere. What must never happen is a schema that
    EXISTS here going unused because this mapping was not extended alongside
    ``_VECTOR_SCHEMAS``.
    """
    if schema_name == "VectorObservationSchema":
        from schemas.vector_observation import VectorObservationSchema
        return VectorObservationSchema
    if schema_name == "VectorReplicateSchema":
        from schemas.vector_replicate import VectorReplicateSchema
        return VectorReplicateSchema
    return None


#: ⛔ Registered in ``_VECTOR_SCHEMAS`` and NOT dispatched by ``_schema_for``, so
#: a caller passing ``validate=True`` on one of these slots gets no check.
#:
#: This is a KNOWN, PRE-EXISTING gap, enumerated rather than left to an accident.
#: All four are defined in ``schemas/validation_claim.py``; a cold review
#: (2026-08-21) demonstrated the hole by execution — a ``ReactionFluxSchema``
#: slot carrying an invalid ``kind`` was returned clean — and found no other
#: place in this repo that validates them.
#:
#: They are NOT dispatched here because switching four previously-unchecked
#: schemas on is a behaviour change that can red existing payloads, and belongs
#: in its own change with its own verification. **Do not extend this list to
#: silence a test.** A new schema goes in ``_schema_for``.
_UNDISPATCHED_VECTOR_SCHEMAS = frozenset({
    "ReactionFluxSchema",
    "ProteinAbundanceSchema",
    "MetaboliteConcentrationSchema",
    "MacromoleculeCompositionSchema",
})


def observations_for_cultivation_group(
    cultivation_group_id: str,
    primary_bundle: Path | str | None = None,
    overlay_paths: list[Path | str] | None = None,
    *,
    include_primary: bool = True,
) -> dict[str, dict]:
    """Scalar observations for one cultivation group, keyed by ``observable``.

    Convenience view for a ``vs_experiment`` card: returns
    ``{observable: band}`` for the given ``cultivation_group_id`` (the shape
    the basal ``vs_literature`` render consumes, sourced from a cultivation
    group instead of literature)."""
    allk = load_scalar_observations(
        primary_bundle, overlay_paths, include_primary=include_primary)
    return {b["observable"]: b for b in allk.values()
            if b.get("cultivation_group_id") == cultivation_group_id and b.get("observable")}


def load_cultivation_groups(
    primary_registry: Path | str | None = None,
    overlay_registries: list[Path | str] | None = None,
) -> dict[str, dict]:
    """Union the cultivation-group registry across primary + overlays, keyed by
    ``cultivation_group_id`` (unique across the union or ``ValueError``).
    Registry paths default to ``cultivations.tsv`` beside each validation
    bundle."""
    registries: list[Path] = []
    if primary_registry is not None:
        registries.append(Path(primary_registry))
    elif (VALIDATION_BUNDLE_PATH.parent / "cultivations.tsv").exists():
        registries.append(VALIDATION_BUNDLE_PATH.parent / "cultivations.tsv")
    if overlay_registries is None:
        overlay_registries = [p.parent / "cultivations.tsv"
                              for p in validation_overlay_paths()]
    registries.extend(Path(p) for p in overlay_registries if Path(p).exists())

    out: dict[str, dict] = {}
    for reg in registries:
        for _, r in _read_tsv(reg).iterrows():
            cid = str(r["cultivation_group_id"])
            if cid in out:
                raise ValueError(f"duplicate cultivation_group_id '{cid}' across registries ({reg})")
            out[cid] = {k: (None if pd.isna(v) else v) for k, v in r.to_dict().items()}
    return out
