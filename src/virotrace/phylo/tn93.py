"""Pairwise Tamura-Nei (TN93) distances, vectorised, with pairwise base frequencies.

Only sites where both sequences have an unambiguous A/C/G/T are compared (ambiguity codes
and gaps are skipped, not resolved). A pair needs ``min_overlap`` such sites to be compared
at all; below that it is treated as "unknown", never as "far" or "close".

The pairwise counts come from one-hot matrix products: ``N[a][b][i, j]`` is the number of
shared sites where sequence i has base a and sequence j has base b.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

_CODE = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i
    _CODE[ord(_b.lower())] = _i

A, C, G, T = range(4)


def encode(sequences: Sequence[str]) -> np.ndarray:
    widths = {len(s) for s in sequences}
    if len(widths) != 1:
        raise ValueError(f"sequences differ in length: {sorted(widths)}")
    return np.stack([_CODE[np.frombuffer(s.encode("ascii"), dtype=np.uint8)] for s in sequences])


def _tn93(counts: np.ndarray, shared: np.ndarray) -> np.ndarray:
    """TN93 from the 4x4 count blocks; ``inf`` where the formula saturates."""
    with np.errstate(divide="ignore", invalid="ignore"):
        pi = np.stack(
            [(counts[b].sum(axis=0) + counts[:, b].sum(axis=0)) / (2 * shared) for b in range(4)]
        )
        pA, pC, pG, pT = pi
        pR, pY = pA + pG, pC + pT
        p1 = (counts[A, G] + counts[G, A]) / shared
        p2 = (counts[C, T] + counts[T, C]) / shared
        q = (
            counts[A, C] + counts[C, A] + counts[A, T] + counts[T, A]
            + counts[G, C] + counts[C, G] + counts[G, T] + counts[T, G]
        ) / shared  # fmt: skip
        t1 = 1 - pR * p1 / (2 * pA * pG) - q / (2 * pR)
        t2 = 1 - pY * p2 / (2 * pC * pT) - q / (2 * pY)
        t3 = 1 - q / (2 * pR * pY)
        d = (
            -2 * pA * pG / pR * np.log(t1)
            - 2 * pC * pT / pY * np.log(t2)
            - 2 * (pR * pY - pA * pG * pY / pR - pC * pT * pR / pY) * np.log(t3)
        )
        bad = ~np.isfinite(d) | (t1 <= 0) | (t2 <= 0) | (t3 <= 0)
        d = np.where(bad, np.inf, d)
        # identical sequences: all counts on the diagonal -> exactly zero
        return np.where((p1 == 0) & (p2 == 0) & (q == 0), 0.0, d)


def tn93_block(codes: np.ndarray, rows: slice, min_overlap: int) -> tuple[np.ndarray, np.ndarray]:
    """Distances from ``codes[rows]`` to all sequences. ``nan`` where overlap < min_overlap."""
    onehot = [(codes == b).astype(np.float32) for b in range(4)]
    counts = np.empty((4, 4, rows.stop - rows.start, len(codes)), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            counts[a, b] = onehot[a][rows] @ onehot[b].T
    shared = counts.sum(axis=(0, 1))
    ok = shared >= min_overlap
    d = _tn93(counts, np.where(ok, shared, 1.0))
    return np.where(ok, d, np.nan), shared


def tn93_edges(
    sequences: Sequence[str], threshold: float, min_overlap: int, *, block: int = 512
) -> tuple[list[tuple[int, int, float]], int]:
    """All pairs ``i < j`` with TN93 <= threshold, and how many pairs could be compared."""
    codes = encode(sequences)
    n = len(codes)
    edges: list[tuple[int, int, float]] = []
    compared = 0
    for start in range(0, n, block):
        stop = min(start + block, n)
        d, _ = tn93_block(codes, slice(start, stop), min_overlap)
        for r in range(stop - start):
            i = start + r
            row = d[r, i + 1 :]
            compared += int(np.isfinite(row).sum() + np.isinf(row).sum())
            for off in np.flatnonzero(row <= threshold):
                edges.append((i, i + 1 + int(off), float(row[off])))
    return edges, compared


def components(n: int, edges: Sequence[tuple[int, int, float]]) -> list[list[int]]:
    """Connected components (single linkage), each sorted, largest first."""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j, _ in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    groups: dict[int, list[int]] = {}
    for k in range(n):
        groups.setdefault(find(k), []).append(k)
    return sorted(groups.values(), key=lambda g: (-len(g), g[0]))
