"""Write an immutable snapshot: sequences, raw records table, manifest.

A snapshot directory is never modified after it is written; taking a snapshot into a
directory that already exists is refused.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

from ..errors import ContractError
from ..io.fasta import write_fasta
from ..io.manifest import (
    SourceRecordInfo,
    build_manifest,
    make_snapshot_id,
    write_manifest,
)
from .base import RawRecord, SourceAdapter


def latest_snapshot(raw_root: Path | str) -> Path:
    """Most recent snapshot directory under ``raw_root`` that has a manifest.

    Snapshot ids start with an ISO date, so name order is chronological. Pin
    ``snapshot_id`` in the config for a reproducible rerun instead of relying on this.
    """
    raw_root = Path(raw_root)
    candidates = sorted(p for p in raw_root.glob("*") if (p / "manifest.json").exists())
    if not candidates:
        raise ContractError(f"no snapshots with a manifest under {raw_root}")
    return candidates[-1]


def write_pointer(raw_root: Path | str, scope: str, snapshot_dir: Path | str) -> Path:
    """Record which snapshot a scope's build uses: ``<raw_root>/<scope>.current``.

    The pointer is a tiny text file holding the snapshot id, so Snakemake has a
    fixed output path even though snapshot directories are named dynamically.
    """
    pointer = Path(raw_root) / f"{scope}.current"
    pointer.parent.mkdir(parents=True, exist_ok=True)
    pointer.write_text(Path(snapshot_dir).name + "\n")
    return pointer


def resolve_snapshot(raw_root: Path | str, scope: str, pinned_id: str | None = None) -> Path:
    """Pinned id wins (reproducible rerun), then the scope pointer, then the newest."""
    raw_root = Path(raw_root)
    if pinned_id:
        pinned = raw_root / pinned_id
        if not (pinned / "manifest.json").exists():
            raise ContractError(f"pinned snapshot {pinned_id!r} not found under {raw_root}")
        return pinned
    pointer = raw_root / f"{scope}.current"
    if pointer.exists():
        target = raw_root / pointer.read_text().strip()
        if not (target / "manifest.json").exists():
            raise ContractError(f"{pointer} points at {target.name}, which has no manifest")
        return target
    return latest_snapshot(raw_root)


def load_snapshot_records(snapshot_dir: Path | str) -> list[RawRecord]:
    """Rebuild raw records from a snapshot's FASTA and records table."""
    from ..io.fasta import read_fasta
    from ..io.tables import read_tsv

    snapshot_dir = Path(snapshot_dir)
    sequences = dict(read_fasta(snapshot_dir / "sequences.fasta"))
    records: list[RawRecord] = []
    for row in read_tsv(snapshot_dir / "records.tsv"):
        accession = row.pop("accession")
        source = row.pop("source")
        if accession not in sequences:
            raise ContractError(f"{snapshot_dir.name}: {accession} in records.tsv has no sequence")
        records.append(RawRecord(accession, sequences[accession], source, row))
    return records


def take_snapshot(
    adapters: Sequence[SourceAdapter],
    raw_root: Path | str,
    *,
    date: str | None = None,
    config_hash: str | None = None,
) -> Path:
    """Fetch from every adapter and write ``<raw_root>/<snapshot_id>/``.

    Returns the snapshot directory. Duplicate accessions across sources keep the
    first occurrence; the cross-source dedup that follows happens in curation.
    """
    if not adapters:
        raise ContractError("take_snapshot needs at least one adapter")

    now = datetime.now(UTC)
    date = date or now.strftime("%Y-%m-%d")
    queries = [a.query_description() for a in adapters]
    snapshot_id = make_snapshot_id(date, queries)

    out_dir = Path(raw_root) / snapshot_id
    if (out_dir / "manifest.json").exists():
        raise ContractError(f"snapshot already exists and is immutable: {out_dir}")
    # A folder with no manifest is the leftover of a run that failed before writing anything
    # (nothing is created until every record is fetched), so it is safe to reuse.

    records: list[RawRecord] = []
    seen: set[str] = set()
    infos: list[SourceRecordInfo] = []
    for adapter in adapters:
        retrieved = 0
        for record in adapter.fetch():
            retrieved += 1
            if record.accession in seen:
                continue
            seen.add(record.accession)
            records.append(record)
        infos.append(
            SourceRecordInfo(
                adapter.name, adapter.query_description(), retrieved, adapter.access_tier
            )
        )

    records.sort(key=lambda r: r.accession)

    out_dir.mkdir(parents=True, exist_ok=True)
    seq_path = out_dir / "sequences.fasta"
    write_fasta(seq_path, ((r.accession, r.sequence) for r in records))

    columns = sorted({key for r in records for key in r.fields})
    table_path = out_dir / "records.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "source", *columns])
        for r in records:
            writer.writerow([r.accession, r.source, *(r.fields.get(c, "") for c in columns)])

    manifest = build_manifest(
        snapshot_id=snapshot_id,
        created_utc=now.isoformat(timespec="seconds"),
        sources=infos,
        files={"sequences.fasta": seq_path, "records.tsv": table_path},
        config_hash=config_hash,
    )
    write_manifest(out_dir / "manifest.json", manifest)
    return out_dir
