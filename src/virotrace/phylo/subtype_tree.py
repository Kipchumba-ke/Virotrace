"""A second, tree-based subtype call, independent of the distance placement.

The first call (``curate/subtype.py``) assigns each sequence to the pure subtype whose panel
references are nearest, with site-bootstrap support. This one asks a different question of a
maximum-likelihood tree built from the panel and the tree set together: does the sequence sit
*inside* a clade of references of one subtype? Disagreement between the two is the finding,
not an error to resolve.

Rule. Prune the tree to the panel references plus one query sequence. The query attaches to one
edge of that pruned tree. The call is subtype L when one side of that edge contains only L
references and is a proper part of them (or L has a single reference): the query then nests
inside L. Attaching at the stem of a whole subtype's clade is not a call (a new lineage or a
recombinant sits there). If no sub-subtype qualifies, the same rule is tried at group level
(A1, A2, A3, A6 -> A), because sub-subtype references are often paraphyletic in pol. A support
value on the joining branch below ``min_support`` makes the call unresolved.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Mapping
from typing import Any

import numpy as np

from .newick import Tree

UNRESOLVED = "unresolved"
_NUMBER = re.compile(r"^\s*(\d+(?:\.\d+)?)")


def support_value(label: str) -> float | None:
    """First number of a branch label: ``95`` or ``95.5/100`` -> 95.0; none -> ``None``."""
    match = _NUMBER.match(label or "")
    return float(match.group(1)) if match else None


def _classify(
    tree: Tree, panel: Mapping[str, str], tip_ids: Mapping[str, int], labels: list[str]
) -> dict[int, tuple[str, float | None]]:
    """Per query tip id: (call, support) at one labelling of the panel."""
    index = {lab: i for i, lab in enumerate(sorted(set(panel.values())))}
    n_label = np.zeros(len(index), dtype=int)
    counts = np.zeros((tree.n_nodes, len(index)), dtype=int)
    for name, lab in panel.items():
        node = tip_ids[name]
        counts[node, index[lab]] += 1
        n_label[index[lab]] += 1
    for node in range(tree.n_nodes - 1, 0, -1):  # preorder ids: children come after parents
        counts[tree.parent[node]] += counts[node]
    total = counts.sum(axis=1)
    names = sorted(index, key=index.get)  # type: ignore[arg-type]

    def inside(side: np.ndarray) -> str | None:
        if side.sum() == 0:
            return None
        (present,) = np.nonzero(side)
        if len(present) != 1:
            return None
        k = int(present[0])
        return names[k] if side[k] < n_label[k] or n_label[k] == 1 else None

    out: dict[int, tuple[str, float | None]] = {}
    for query in labels:
        node = tip_ids[query]
        a = tree.parent[node]
        while a != -1 and total[a] == 0:
            node, a = a, tree.parent[a]
        if a == -1:
            out[tip_ids[query]] = (UNRESOLVED, None)
            continue
        bearing = [c for c in tree.children[a] if total[c] > 0]
        if len(bearing) != 1:  # the query joins at a multifurcation: no single edge
            out[tip_ids[query]] = (UNRESOLVED, None)
            continue
        side = counts[bearing[0]]
        call = inside(side) or inside(n_label - side)
        out[tip_ids[query]] = (call or UNRESOLVED, support_value(tree.label[a]))
    return out


def tree_subtype_calls(
    tree: Tree,
    panel: Mapping[str, str],
    groups: Mapping[str, str],
    *,
    min_support: float | None = None,
) -> dict[str, dict[str, Any]]:
    """Call every non-panel tip. ``panel`` maps reference tip -> label, ``groups`` label -> group."""
    tip_ids = {tree.name[t]: t for t in tree.tips}
    missing = sorted(set(panel) - set(tip_ids))
    if missing:
        raise ValueError(f"panel references not in the tree: {missing[:3]}")
    queries = [tree.name[t] for t in tree.tips if tree.name[t] not in panel]
    fine = _classify(tree, panel, tip_ids, queries)
    coarse = _classify(tree, {n: groups[lab] for n, lab in panel.items()}, tip_ids, queries)
    result: dict[str, dict[str, Any]] = {}
    for q in queries:
        call, support = fine[tip_ids[q]]
        level = "subtype"
        if call == UNRESOLVED:
            call, support = coarse[tip_ids[q]]
            level = "group"
        if call == UNRESOLVED:
            level = "none"
        elif min_support is not None and support is not None and support < min_support:
            call, level = UNRESOLVED, "weak_support"
        result[q] = {"call": call, "level": level, "support": support}
    return result


def concordance(
    rows: list[dict[str, str]], agree: Any
) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    """Counts of agree / disagree / unresolved / not_comparable, and a confusion table."""
    summary: Counter[str] = Counter()
    confusion: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        confusion[row["subtype_called"] or "none"][row["subtype_tree"]] += 1
        if row["subtype_tree"] == UNRESOLVED:
            row["concordance"] = UNRESOLVED
        elif row["subtype_called"] in ("", "unassigned"):
            row["concordance"] = "not_comparable"
        else:
            verdict = agree(row["subtype_tree"], row["subtype_called"])
            row["concordance"] = "agree" if verdict == "agree" else (
                "disagree" if verdict == "disagree" else "not_comparable"
            )  # fmt: skip
        summary[row["concordance"]] += 1
    return dict(summary), {k: dict(v) for k, v in confusion.items()}
