"""Cross-source deduplication with a decision log.

Levels, in order:

1. **Accession** -- the same record from two sources (version suffix ignored).
2. **Exact sequence identity** -- the same sequence under different accessions, which
   is how one study's deposits appear in both LANL and GenBank.
3. **Same isolate name, country and year, different sequence** -- *flagged, never
   collapsed*: these may be distinct clones or a re-sequencing, and only a person
   can tell. Collapsing them would be an unlogged judgement call.

Every collapsed record is written to the log with the record it lost to, so no
sequence disappears silently.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from ..acquire.base import RawRecord


@dataclass(frozen=True, slots=True)
class DedupDecision:
    accession: str
    kept_instead: str
    level: str  # "accession" | "sequence" | "flag-isolate"
    detail: str = ""


@dataclass(slots=True)
class DedupResult:
    kept: list[RawRecord]
    log: list[DedupDecision]
    n_in: int

    def assert_conserved(self) -> None:
        collapsed = sum(1 for d in self.log if d.level != "flag-isolate")
        if self.n_in != len(self.kept) + collapsed:
            raise AssertionError(
                f"record conservation violated: in={self.n_in} "
                f"kept={len(self.kept)} collapsed={collapsed}"
            )


def base_accession(accession: str) -> str:
    """``AB123456.2`` -> ``AB123456``."""
    return accession.strip().split(".")[0]


def _ambiguity(sequence: str) -> int:
    return sum(1 for b in sequence if b not in "ACGT")


def _rank(record: RawRecord) -> tuple[int, int, int, str]:
    """Lower is better: GenBank first (it carries the study fields LANL lacks), then fewer
    ambiguous bases, then longer, then accession for stability."""
    return (
        0 if record.source == "genbank" else 1,
        _ambiguity(record.sequence),
        -len(record.sequence),
        record.accession,
    )


def deduplicate(records: Iterable[RawRecord]) -> DedupResult:
    records = list(records)
    log: list[DedupDecision] = []

    # Level 1: accession
    by_accession: dict[str, RawRecord] = {}
    for record in sorted(records, key=lambda r: (base_accession(r.accession), _rank(r))):
        key = base_accession(record.accession)
        if key in by_accession:
            log.append(DedupDecision(record.accession, by_accession[key].accession, "accession"))
        else:
            by_accession[key] = record

    # Level 2: exact sequence identity
    by_sequence: dict[str, RawRecord] = {}
    for record in sorted(by_accession.values(), key=_rank):
        fingerprint = record.sequence.upper()
        if fingerprint in by_sequence:
            log.append(
                DedupDecision(record.accession, by_sequence[fingerprint].accession, "sequence")
            )
        else:
            by_sequence[fingerprint] = record

    kept = sorted(by_sequence.values(), key=lambda r: r.accession)

    # Level 3: flag only
    groups: dict[tuple[str, str, str], list[RawRecord]] = defaultdict(list)
    for record in kept:
        isolate = record.fields.get("isolate", "").strip().lower()
        country = record.fields.get("country", "").split(":")[0].strip().lower()
        year = record.fields.get("collection_date", "").strip()[-4:]
        if isolate and country and year:
            groups[(isolate, country, year)].append(record)
    for members in groups.values():
        if len(members) > 1:
            anchor = members[0].accession
            for other in members[1:]:
                log.append(
                    DedupDecision(
                        other.accession,
                        anchor,
                        "flag-isolate",
                        "same isolate, country and year but a different sequence; kept both",
                    )
                )

    result = DedupResult(kept=kept, log=log, n_in=len(records))
    result.assert_conserved()
    return result
