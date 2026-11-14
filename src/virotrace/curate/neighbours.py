"""Nearest-neighbour identity, to find tips that have no close relative in the tree set.

A tip whose best match over the shared sites is far below what any subtype pair shows
(roughly 0.85 to 0.92 for HIV-1 pol) is either misaligned, from a different virus, or the
only representative of an unrelated lineage. IQ-TREE gives such a tip a very long terminal
branch (up to its ceiling of 10), which stretches the tree axis and drags a midpoint root onto
it. The screen is descriptive: it keeps the tree readable and says what it dropped.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

CHUNK = 256


def nearest_identity(seqs: Sequence[str], *, min_shared: int) -> list[float]:
    """Best identity of each sequence to any other, over sites both have as A, C, G or T.

    A partner needs at least ``min_shared`` such sites. A sequence with no such partner gets
    ``nan`` (not assessable), never 0, so it cannot be dropped for lack of comparisons.
    """
    n = len(seqs)
    if n < 2:
        return [float("nan")] * n
    width = len(seqs[0])
    if any(len(s) != width for s in seqs):
        raise ValueError("sequences must be aligned to one length")
    arr = np.frombuffer("".join(s.upper() for s in seqs).encode(), dtype=np.uint8).reshape(n, width)
    onehot = [(arr == ord(b)).astype(np.float32) for b in "ACGT"]
    called = sum(onehot)  # 1.0 where the site is A, C, G or T
    out: list[float] = []
    for lo in range(0, n, CHUNK):
        hi = min(lo + CHUNK, n)
        shared = called[lo:hi] @ called.T
        same = sum(o[lo:hi] @ o.T for o in onehot)
        shared[np.arange(hi - lo), np.arange(lo, hi)] = 0  # not itself
        with np.errstate(divide="ignore", invalid="ignore"):
            ident = np.where(shared >= min_shared, same / shared, np.nan)
        best = np.full(hi - lo, np.nan)
        has = ~np.all(np.isnan(ident), axis=1)
        best[has] = np.nanmax(ident[has], axis=1)
        out.extend(float(x) for x in best)
    return out
