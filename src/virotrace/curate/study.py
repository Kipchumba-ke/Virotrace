"""Study identity and treatment context, read from what the record itself says.

Prevalence pooled over GenBank is a mixture of studies: surveillance in ART-naive people
next to treatment-failure cohorts. Two things make that visible without inventing data:

* ``study_key`` groups records that belong to one publication or project, so one study
  cannot dominate a stratum unnoticed.
* ``treatment_from_text`` reads an explicit statement ("treatment failure", "ART-naive")
  from the title, definition or notes. It is never inferred from mutations, which would
  be circular, and a record with no statement, or contradictory ones, stays ``unknown``.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

_EXPERIENCED = re.compile(
    r"treatment[- ]failure|virolog\w*[- ]failure|\bfailing\b|\bfailure\b|"
    r"(?:treatment|art|haart)[- ]experienced|second[- ]line|\bhaart\b|"
    r"art[- ]treated|on (?:art|haart|antiretroviral)|acquired (?:drug )?resistance",
    re.IGNORECASE,
)
_NAIVE = re.compile(
    r"na[iï]ve|untreated|newly diagnosed|transmitted (?:drug )?resistance|"
    r"primary (?:drug )?resistance|pre-?treatment|before (?:starting )?(?:art|treatment)",
    re.IGNORECASE,
)
_TEXT_FIELDS = ("study_title", "definition", "note", "isolation_source")


def treatment_from_text(fields: Mapping[str, Any]) -> str:
    text = " ".join(str(fields.get(name) or "") for name in _TEXT_FIELDS)
    naive, experienced = bool(_NAIVE.search(text)), bool(_EXPERIENCED.search(text))
    if naive and not experienced:
        return "naive"
    if experienced and not naive:
        return "experienced"
    return "unknown"


def study_key(fields: Mapping[str, Any]) -> str:
    """PubMed id, else BioProject, else first author plus title; ``unknown_study`` if none."""
    pubmed = str(fields.get("pubmed") or "").strip()
    if pubmed:
        return f"PMID:{pubmed}"
    bioproject = str(fields.get("bioproject") or "").strip()
    if bioproject:
        return f"BP:{bioproject}"
    authors = str(fields.get("study_authors") or "").strip()
    title = str(fields.get("study_title") or "").strip()
    if authors or title:
        first = authors.split(",")[0].split(" ")[0] if authors else ""
        slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")[:40]
        return f"REF:{first}|{slug}"
    return "unknown_study"


TREATMENT_CONTEXTS = ("naive", "experienced", "mixed", "unknown")


def load_study_annotations(path: Any) -> dict[str, dict[str, str]]:
    """Hand-curated ``study_key -> {treatment_context, note}``; absent file = no annotations."""
    import csv
    from pathlib import Path

    from ..errors import ContractError

    path = Path(path)
    if not path.exists():
        return {}
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        rows = (line for line in handle if not line.lstrip().startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        if not {"study_key", "treatment_context"} <= set(reader.fieldnames or []):
            raise ContractError(f"{path} needs columns study_key and treatment_context")
        for row in reader:
            context = row["treatment_context"].strip().lower()
            if context not in TREATMENT_CONTEXTS:
                raise ContractError(
                    f"{path}: {row['study_key']}: treatment_context {context!r} not in "
                    f"{TREATMENT_CONTEXTS}"
                )
            out[row["study_key"].strip()] = {
                "treatment_context": context,
                "note": (row.get("note") or "").strip(),
            }
    return out
