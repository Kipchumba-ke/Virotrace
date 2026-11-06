"""Clone policy: which clone-derived records count as independent observations.

A GenBank record with a ``/clone`` tag, or "clone"/"proviral" in its definition, is
either a lab infectious clone or one of several PCR clones from a single patient.
Prevalence is a statement about people, so counting every clone would let a few
heavily cloned patients dominate it. The policy is a configuration choice, not a
hard-coded rule:

* ``exclude``   -- clone-derived records leave the prevalence and tree layers. They stay
  in the metadata and the DRM table, so nothing is lost and the choice is reversible.
* ``keep_all``  -- no clone filtering.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from ..errors import ConfigError

POLICIES = ("exclude", "keep_all")


def is_excluded(row: Mapping[str, Any], policy: str) -> bool:
    if policy not in POLICIES:
        raise ConfigError(f"qc.clone_policy must be one of {POLICIES}, got {policy!r}")
    return policy == "exclude" and bool(row.get("clone_derived"))


def split_by_policy(
    rows: Iterable[Mapping[str, Any]], policy: str
) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]]]:
    """``(kept, excluded)``; every input row lands in exactly one."""
    kept: list[Mapping[str, Any]] = []
    excluded: list[Mapping[str, Any]] = []
    for row in rows:
        (excluded if is_excluded(row, policy) else kept).append(row)
    return kept, excluded


INDIVIDUAL_POLICIES = ("keep_all", "one_per_individual")


def individual_id(isolate: str | None, patterns: Sequence[str]) -> str | None:
    """The person an isolate name identifies, or ``None`` when no pattern matches."""
    for pattern in patterns:
        match = re.search(pattern, isolate or "")
        if match:
            return match.group(1)
    return None


def person_key(row: Mapping[str, Any], patterns: Sequence[str]) -> str:
    """The person a record belongs to.

    A matching pattern gives the person id. Otherwise records with the same isolate name in
    the same study are the same sample (repeated clones of one isolate); with no isolate name
    the record is its own person.
    """
    person = individual_id(row.get("isolate"), patterns)
    if person is not None:
        return person
    isolate = str(row.get("isolate") or "").strip()
    if isolate:
        return f"{row.get('study_key') or ''}|{isolate}"
    return str(row["strain"])


def person_map(rows: Iterable[Mapping[str, Any]], patterns: Sequence[str]) -> dict[str, str]:
    """``strain -> person``; a record matching no pattern is its own person (its strain id)."""
    return {r["strain"]: individual_id(r.get("isolate"), patterns) or r["strain"] for r in rows}


def per_person_counts(
    values: Mapping[str, str], person_of: Mapping[str, str]
) -> tuple[dict[str, int], int]:
    """Count people, not sequences, by a per-strain value.

    A person takes their most common value (ties broken alphabetically). Returns
    ``({value: n_people}, n_people_whose_sequences_disagree)``.
    """
    per: dict[str, list[str]] = {}
    for strain, value in values.items():
        per.setdefault(person_of.get(strain, strain), []).append(value)
    counts: dict[str, int] = {}
    mixed = 0
    for vals in per.values():
        tally: dict[str, int] = {}
        for v in vals:
            tally[v] = tally.get(v, 0) + 1
        mixed += len(tally) > 1
        top = sorted(tally.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        counts[top] = counts.get(top, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1])), mixed


def one_per_individual(
    rows: Iterable[Mapping[str, Any]], patterns: Sequence[str]
) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]]]:
    """``(kept, excluded)``: one sequence for each identifiable person.

    Within a person the earliest dated sequence wins, then the least ambiguous, then the
    smallest strain id, so the choice is deterministic and does not depend on input order.
    A record matching no pattern is grouped with other records of the same isolate name in
    the same study (see ``person_key``); with no isolate name it is its own person. Every
    input row lands in exactly one list.
    """
    rows = list(rows)
    groups: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        groups.setdefault(person_key(row, patterns), []).append(row)

    def rank(row: Mapping[str, Any]) -> tuple[str, float, str]:
        return (
            str(row.get("date") or "9999"),
            float(row.get("ambiguity_frac") or 0.0),
            str(row["strain"]),
        )

    winners = {id(min(members, key=rank)) for members in groups.values()}
    kept: list[Mapping[str, Any]] = []
    excluded: list[Mapping[str, Any]] = []
    for row in rows:
        (kept if id(row) in winners else excluded).append(row)
    return kept, excluded


def apply_population_policies(
    rows: Iterable[Mapping[str, Any]],
    *,
    clone_policy: str,
    individual_policy: str,
    patterns: Sequence[str],
) -> tuple[list[Mapping[str, Any]], dict[str, int]]:
    """Clone policy, then individual policy. Returns the kept rows and how many each removed."""
    if individual_policy not in INDIVIDUAL_POLICIES:
        raise ConfigError(
            f"qc.individual_policy must be one of {INDIVIDUAL_POLICIES}, got {individual_policy!r}"
        )
    after_clone, clones = split_by_policy(rows, clone_policy)
    if individual_policy == "one_per_individual":
        kept, repeats = one_per_individual(after_clone, patterns)
    else:
        kept, repeats = list(after_clone), []
    return kept, {"clone_policy": len(clones), "individual_policy": len(repeats)}
