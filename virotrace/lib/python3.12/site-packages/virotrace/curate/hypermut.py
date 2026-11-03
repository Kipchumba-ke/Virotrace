"""APOBEC3G/F hypermutation screen (QC-5).

Host APOBEC editing produces G-to-A changes in a sequence-context-dependent way.
Those changes can look exactly like drug-resistance mutations, so hypermutated
sequences are a known source of spurious resistance signal in database surveys.
Omitting this screen would make the headline prevalence number unreliable, which is
why it runs by default and its excluded count is reported.

Method follows the Hypermut-style context-ratio test: compare the rate of G-to-A
substitution in APOBEC context (G followed by G or A) against the rate in control
context (G followed by C or T), using a one-tailed Fisher exact test.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

APOBEC_DOWNSTREAM = {"G", "A"}
CONTROL_DOWNSTREAM = {"C", "T"}


@dataclass(frozen=True, slots=True)
class HypermutResult:
    strain: str
    context_mutated: int
    context_total: int
    control_mutated: int
    control_total: int
    p_value: float
    flagged: bool

    @property
    def context_rate(self) -> float:
        return self.context_mutated / self.context_total if self.context_total else 0.0

    @property
    def control_rate(self) -> float:
        return self.control_mutated / self.control_total if self.control_total else 0.0


def fisher_exact_greater(a: int, b: int, c: int, d: int) -> float:
    """One-tailed Fisher exact p-value for enrichment in the first row.

    Table layout::

        [[a, b],
         [c, d]]

    Returns P(X >= a) under the hypergeometric null.
    """
    if min(a, b, c, d) < 0:
        raise ValueError("counts must be non-negative")
    row1 = a + b
    col1, total = a + c, a + b + c + d
    if total == 0 or row1 == 0 or col1 == 0:
        return 1.0

    denominator = comb(total, row1)
    upper = min(row1, col1)
    tail = sum(comb(col1, k) * comb(total - col1, row1 - k) for k in range(a, upper + 1))
    return min(1.0, tail / denominator)


def screen_hypermutation(
    strain: str, query: str, reference: str, *, p_threshold: float = 0.05
) -> HypermutResult:
    """Screen one aligned query against an aligned reference of equal length."""
    if len(query) != len(reference):
        raise ValueError("query and reference must be aligned to equal length")

    query = query.upper()
    reference = reference.upper()
    counts = {"ctx_mut": 0, "ctx_tot": 0, "ctl_mut": 0, "ctl_tot": 0}

    for i in range(len(reference) - 1):
        ref_base, query_base = reference[i], query[i]
        if ref_base != "G" or query_base not in {"G", "A"}:
            continue

        downstream = next(
            (reference[j] for j in range(i + 1, len(reference)) if reference[j] != "-"), None
        )
        if downstream is None:
            continue

        if downstream in APOBEC_DOWNSTREAM:
            counts["ctx_tot"] += 1
            if query_base == "A":
                counts["ctx_mut"] += 1
        elif downstream in CONTROL_DOWNSTREAM:
            counts["ctl_tot"] += 1
            if query_base == "A":
                counts["ctl_mut"] += 1

    p_value = fisher_exact_greater(
        counts["ctx_mut"],
        counts["ctx_tot"] - counts["ctx_mut"],
        counts["ctl_mut"],
        counts["ctl_tot"] - counts["ctl_mut"],
    )

    return HypermutResult(
        strain=strain,
        context_mutated=counts["ctx_mut"],
        context_total=counts["ctx_tot"],
        control_mutated=counts["ctl_mut"],
        control_total=counts["ctl_tot"],
        p_value=p_value,
        flagged=p_value < p_threshold,
    )
