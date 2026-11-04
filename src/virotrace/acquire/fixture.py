"""Offline fixture adapter, used by the dev profile and CI (no network).

Reads a TSV with a header row. Required columns: ``accession``, ``sequence``.
Every other column is passed through as a raw field.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path

from ..errors import ContractError
from .base import RawRecord


class FixtureAdapter:
    name = "fixture"
    access_tier = 1

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)

    def query_description(self) -> str:
        return f"fixture:{self.path.name}"

    def fetch(self) -> Iterator[RawRecord]:
        if not self.path.exists():
            raise ContractError(f"fixture file not found: {self.path}")
        with self.path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = {"accession", "sequence"} - set(reader.fieldnames or [])
            if missing:
                raise ContractError(f"{self.path}: missing columns {sorted(missing)}")
            for row in reader:
                extra = {k: (v or "") for k, v in row.items() if k not in {"accession", "sequence"}}
                yield RawRecord(
                    accession=row["accession"].strip(),
                    sequence=row["sequence"].strip().upper(),
                    source=self.name,
                    fields=extra,
                )
