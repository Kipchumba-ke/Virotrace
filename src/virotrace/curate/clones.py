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


def one_per_individual(
    rows: Iterable[Mapping[str, Any]], patterns: Sequence[str]
) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]]]:
    """``(kept, excluded)``: one sequence for each identifiable person.

    Within a person the earliest dated sequence wins, then the least ambiguous, then the
    smallest strain id, so the choice is deterministic and does not depend on input order.
    A record whose isolate matches no pattern is its own person and is always kept. Every
    input row lands in exactly one list.
    """
    rows = list(rows)
    groups: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        person = individual_id(row.get("isolate"), patterns)
        if person is not None:
            groups.setdefault(person, []).append(row)

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
        identified = individual_id(row.get("isolate"), patterns) is not None
        (kept if not identified or id(row) in winners else excluded).append(row)
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
