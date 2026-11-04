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
    if out_dir.exists():
        raise ContractError(f"snapshot already exists and is immutable: {out_dir}")
    out_dir.mkdir(parents=True)

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
