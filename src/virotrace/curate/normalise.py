"""Raw records -> the validated metadata contract (one row per sequence).

Records that cannot satisfy the contract are quarantined with a reason code, never
dropped. If too many quarantine, the stage aborts: silent mass-quarantine is the
failure mode the RecordError/ContractError split would otherwise invite.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from ..acquire.base import RawRecord
from ..errors import ContractError, RecordError
from ..io.fasta import STRAIN_RE
from .dates import parse_collection_date
from .geo import normalise_country
from .study import study_key, treatment_from_text

REQUIRED_FIELDS = (
    "strain",
    "accession",
    "source",
    "access_tier",
    "country",
    "country_code",
    "region",
    "date",
    "date_precision",
    "year",
    "ambiguity_frac",
    "qc_status",
    "snapshot_id",
)
DATE_PRECISIONS = {"day", "month", "year", "unknown"}
SOURCES = {"lanl", "genbank", "fixture"}
_NON_ID = re.compile(r"[^A-Za-z0-9_.\-]")


def canonical_strain(accession: str) -> str:
    """Canonical ID used in FASTA headers and every join."""
    strain = _NON_ID.sub("_", accession.strip())
    if not STRAIN_RE.match(strain):
        raise RecordError(f"cannot form a strain id from {accession!r}", rule="N-1")
    return strain


_CLONE_WORDS = re.compile(r"\b(molecular clone|infectious clone|clone|proviral)\b", re.IGNORECASE)


def is_clone_derived(fields: dict[str, Any]) -> bool:
    """True if the record says it comes from a clone or a provirus.

    A flag, not a verdict: lab infectious clones and within-host PCR clones both set it.
    """
    if fields.get("clone") or fields.get("proviral"):
        return True
    return bool(_CLONE_WORDS.search(fields.get("definition", "") or ""))


def ambiguity_fraction(sequence: str) -> float:
    ungapped = [b for b in sequence if b not in "-. "]
    if not ungapped:
        return 1.0
    return sum(1 for b in ungapped if b.upper() not in "ACGT") / len(ungapped)


def normalise_record(
    record: RawRecord,
    *,
    snapshot_id: str,
    access_tier: int = 1,
    min_year: int = 1985,
    max_year: int = 2026,
) -> dict[str, Any]:
    """One metadata row, or ``RecordError`` naming the rule that rejected it."""
    strain = canonical_strain(record.accession)

    geo = normalise_country(record.fields.get("country"))
    if geo is None:
        raise RecordError(
            f"unmapped country {record.fields.get('country')!r}", record_id=strain, rule="N-2"
        )

    parsed = parse_collection_date(
        record.fields.get("collection_date"), min_year=min_year, max_year=max_year
    )
    if parsed.precision == "unknown":
        # Kept, not quarantined: an undated sequence is still valid for the DRM layer.
        # It is simply ineligible for the time tree, which date_precision records.
        pass

    row: dict[str, Any] = {
        "strain": strain,
        "accession": record.accession.strip(),
        "source": record.source,
        "access_tier": access_tier,
        "country": geo.country,
        "country_code": geo.country_code,
        "region": geo.region,
        "admin1": geo.admin1_raw or None,
        "admin1_resolved": geo.admin1_resolved,
        "date": parsed.date,
        "date_precision": parsed.precision,
        "year": parsed.year,
        "isolate": record.fields.get("isolate", "") or None,
        "subtype_reported": record.fields.get("subtype", "") or None,
        "ambiguity_frac": round(ambiguity_fraction(record.sequence), 6),
        "clone_derived": is_clone_derived(record.fields),
        # Only from an explicit statement in the record's own text; never from mutations.
        "treatment_status": treatment_from_text(record.fields),
        "study_key": study_key(record.fields),
        "study_title": (record.fields.get("study_title", "") or "")[:200] or None,
        "subtype_reported": record.fields.get("subtype", "") or None,
        "ambiguity_frac": round(ambiguity_fraction(record.sequence), 6),
        # Never inferred from mutations: that would be circular.
        "treatment_status": "unknown",
        "qc_status": "pass",
        "qc_reason": None,
        "snapshot_id": snapshot_id,
    }
    validate_row(row)
    return row


def validate_row(row: dict[str, Any]) -> None:
    """Check one row against the metadata contract. Raises ``ContractError``."""
    for name in REQUIRED_FIELDS:
        if name == "year":
            continue  # nullable when the date is unknown
        if row.get(name) in (None, ""):
            if name == "date":
                continue  # empty string is the unknown-date value
            raise ContractError(f"{row.get('strain')}: required field {name!r} is empty")
    if row["date_precision"] not in DATE_PRECISIONS:
        raise ContractError(f"{row['strain']}: bad date_precision {row['date_precision']!r}")
    if row["source"] not in SOURCES:
        raise ContractError(f"{row['strain']}: bad source {row['source']!r}")
    if row["access_tier"] not in (1, 2, 3):
        raise ContractError(f"{row['strain']}: bad access_tier {row['access_tier']!r}")
    if row["date_precision"] != "unknown" and row["year"] is None:
        raise ContractError(f"{row['strain']}: dated record has no year")
    if not 0.0 <= row["ambiguity_frac"] <= 1.0:
        raise ContractError(f"{row['strain']}: ambiguity_frac outside 0..1")


@dataclass(slots=True)
class NormaliseResult:
    rows: list[dict[str, Any]]
    quarantine: list[dict[str, str]]
    n_in: int

    def assert_conserved(self) -> None:
        if self.n_in != len(self.rows) + len(self.quarantine):
            raise AssertionError(
                f"record conservation violated: in={self.n_in} "
                f"rows={len(self.rows)} quarantined={len(self.quarantine)}"
            )


def normalise_records(
    records: Iterable[RawRecord],
    *,
    snapshot_id: str,
    max_record_error_rate: float = 0.25,
    min_year: int = 1985,
    max_year: int = 2026,
) -> NormaliseResult:
    records = list(records)
    rows: list[dict[str, Any]] = []
    quarantine: list[dict[str, str]] = []
    seen: set[str] = set()

    for record in records:
        try:
            row = normalise_record(
                record, snapshot_id=snapshot_id, min_year=min_year, max_year=max_year
            )
            if row["strain"] in seen:
                raise RecordError("duplicate strain id", record_id=row["strain"], rule="N-3")
            seen.add(row["strain"])
            rows.append(row)
        except RecordError as exc:
            quarantine.append(
                {
                    "accession": record.accession,
                    "rule": exc.rule or "N-0",
                    "reason": str(exc),
                }
            )

    result = NormaliseResult(rows=rows, quarantine=quarantine, n_in=len(records))
    result.assert_conserved()
    if records and len(quarantine) / len(records) > max_record_error_rate:
        by_rule = Counter(q["rule"] for q in quarantine)
        examples = "; ".join(
            f"{q['accession']} [{q['rule']}] {q['reason']}" for q in quarantine[:3]
        )
        error = ContractError(
            f"{len(quarantine)}/{len(records)} records quarantined, above "
            f"max_record_error_rate={max_record_error_rate}; refusing to continue. "
            f"By rule: {dict(by_rule)}. Examples: {examples}"
        )
        error.quarantine = quarantine  # type: ignore[attr-defined]
        raise error
        raise ContractError(
            f"{len(quarantine)}/{len(records)} records quarantined, above "
            f"max_record_error_rate={max_record_error_rate}; refusing to continue"
        )
    return result
