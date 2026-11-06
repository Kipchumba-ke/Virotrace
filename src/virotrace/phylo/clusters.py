"""TN93 clusters as a descriptive population statistic -- never as transmission claims.

Two people are linked when their TN93 distance is at or below ``phylo.clusters.threshold``;
a cluster is a connected component of at least ``min_size`` people (single linkage). The
summary reports cluster sizes, subtype and study composition, how SDRMs are distributed
inside versus outside clusters, and how the pair distances are spread (very close
cross-study pairs can be the same person deposited twice). Membership is written per
sequence only so the Auspice export can colour by it.

Resistance mutations are shared by unrelated people on the same regimen and pull their
sequences together, so by default the whole analysis is repeated with SDRM codons masked
(``sensitivity_mask_sdrm``); the primary result leaves them in.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from ..curate.study import load_study_annotations
from ..drm.prevalence import wilson_ci
from ..drm.sdrm import load_sdrm_list
from ..errors import ContractError
from ..io.fasta import read_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .tn93 import components, tn93_edges

EDGE_BINS = (0.001, 0.0025, 0.005, 0.01)


def sdrm_columns(settings: Settings, repo_root: Path) -> set[int]:
    """0-based alignment columns of every codon at a headline-list SDRM position."""
    genes = settings.reference.genes
    starts = {"PR": genes.pr[0], "RT": genes.rt[0], "IN": genes.in_[0]}
    sdrm = load_sdrm_list(
        repo_root / settings.paths["reference"] / settings.drm.sdrm.list,
        version=settings.drm.sdrm.list_version,
    )
    cols: set[int] = set()
    for gene, position in sdrm.entries:
        first = starts[gene] + 3 * (position - 1) - genes.pol[0]
        cols.update(range(first, first + 3))
    return cols


def _mask(sequence: str, columns: set[int]) -> str:
    return "".join("N" if i in columns else c for i, c in enumerate(sequence))


def _cluster(
    strains: list[str], sequences: list[str], threshold: float, min_overlap: int, min_size: int
) -> dict[str, Any]:
    edges, compared = tn93_edges(sequences, threshold, min_overlap)
    groups = [g for g in components(len(strains), edges) if len(g) >= min_size]
    members = {strains[k]: (f"CL{r:03d}", len(g)) for r, g in enumerate(groups, 1) for k in g}
    return {"edges": edges, "compared": compared, "groups": groups, "members": members}


def _drm_block(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n, k = len(rows), sum(r["has_sdrm"] for r in rows)
    lo, hi = wilson_ci(k, n) if n else (None, None)
    return {
        "n": n, "with_sdrm": k, "proportion": round(k / n, 4) if n else None,
        "ci_low": None if lo is None else round(lo, 4),
        "ci_high": None if hi is None else round(hi, 4),
    }  # fmt: skip


def _by_membership(rows: list[dict[str, Any]], members: dict[str, Any]) -> dict[str, Any]:
    inside = [r for r in rows if r["strain"] in members]
    outside = [r for r in rows if r["strain"] not in members]
    return {"clustered": _drm_block(inside), "unclustered": _drm_block(outside)}


def mantel_haenszel(
    tables: list[tuple[int, int, int, int]], z: float = 1.96
) -> dict[str, float | None]:
    """Common odds ratio over strata (a, b, c, d) = (clustered with/without, unclustered
    with/without), with the Robins-Breslow-Greenland interval."""
    import math

    sr = ss = 0.0
    sum_pr = sum_ps_qr = sum_qs = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        r, s_ = a * d / n, b * c / n
        p, q = (a + d) / n, (b + c) / n
        sr, ss = sr + r, ss + s_
        sum_pr += p * r
        sum_ps_qr += p * s_ + q * r
        sum_qs += q * s_
    if sr == 0 or ss == 0:
        return {"odds_ratio": None, "ci_low": None, "ci_high": None}
    se = math.sqrt(sum_pr / (2 * sr**2) + sum_ps_qr / (2 * sr * ss) + sum_qs / (2 * ss**2))
    estimate = sr / ss
    return {
        "odds_ratio": round(estimate, 3),
        "ci_low": round(math.exp(math.log(estimate) - z * se), 3),
        "ci_high": round(math.exp(math.log(estimate) + z * se), 3),
    }


def run_clusters(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    processed = root / settings.paths["processed"] / settings.scope
    results = root / settings.paths["results"] / settings.scope
    needed = {
        processed / "aln.tree.fasta": "subsample",
        processed / "metadata.final.tsv": "subsample",
        results / "drm" / "drm_table.tsv": "drm call",
    }
    for path, stage in needed.items():
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    cfg = settings.phylo.clusters
    meta = {r["strain"]: r for r in read_tsv(processed / "metadata.final.tsv")}
    aln = list(read_fasta(processed / "aln.tree.fasta"))
    strains = [s for s, _ in aln]
    sequences = [q for _, q in aln]
    if set(strains) != set(meta):
        raise ContractError("aln.tree.fasta and metadata.final.tsv list different strains")

    mutations: dict[str, set[str]] = {}
    for r in read_tsv(results / "drm" / "drm_table.tsv"):
        if r["strain"] in meta and r["sdrm_residues"]:
            for res in r["sdrm_residues"].split(","):
                mutations.setdefault(r["strain"], set()).add(
                    f"{r['gene']}:{r['ref_aa']}{r['position']}{res}"
                )

    primary = _cluster(strains, sequences, cfg.threshold, cfg.min_overlap, cfg.min_size)
    edges, groups, members = primary["edges"], primary["groups"], primary["members"]

    rows = []
    for s in strains:
        cid, size = members.get(s, ("", 1))
        rows.append(
            {
                "strain": s, "cluster_id": cid, "cluster_size": size,
                "subtype_called": meta[s].get("subtype_called", ""),
                "study_key": meta[s].get("study_key", ""),
                "has_sdrm": s in mutations,
            }
        )  # fmt: skip
    (results / "phylo").mkdir(parents=True, exist_ok=True)
    write_tsv(results / "phylo" / "clusters.tsv", rows)

    shared = 0
    for group in groups:
        counts = Counter(m for k in group for m in mutations.get(strains[k], ()))
        shared += any(v >= 2 for v in counts.values())

    study = [meta[s].get("study_key") for s in strains]
    bins: dict[str, dict[str, int]] = {}
    lower = -1.0
    for upper in (*EDGE_BINS, cfg.threshold):
        hit = [(i, j) for i, j, d in edges if lower < d <= upper]
        bins[f"<={upper}"] = {
            "n_pairs": len(hit),
            "cross_study": sum(study[i] != study[j] for i, j in hit),
        }
        lower = upper

    # Is the SDRM enrichment inside clusters just which cohorts are densely sampled?
    ref = root / settings.paths["reference"] / settings.prevalence.study_annotations
    annotations = load_study_annotations(ref) if ref.exists() else {}
    for row in rows:
        row["treatment_context"] = annotations.get(row["study_key"], {}).get(
            "treatment_context", "unannotated"
        )
    tables, by_context = [], {}
    for context in sorted({r["treatment_context"] for r in rows}):
        sub = [r for r in rows if r["treatment_context"] == context]
        by_context[context] = _by_membership(sub, members)
        a = sum(r["has_sdrm"] for r in sub if r["cluster_id"])
        b = sum(not r["has_sdrm"] for r in sub if r["cluster_id"])
        c = sum(r["has_sdrm"] for r in sub if not r["cluster_id"])
        d = sum(not r["has_sdrm"] for r in sub if not r["cluster_id"])
        tables.append((a, b, c, d))

    # near-identical pairs: unrelated people almost never differ at <=0.1% of pol, so these
    # are mostly the same person or sample deposited twice -- a duplicate check, not a finding
    near = [
        {
            "strain_a": strains[i], "strain_b": strains[j], "distance": round(d, 6),
            "isolate_a": meta[strains[i]].get("isolate", ""),
            "isolate_b": meta[strains[j]].get("isolate", ""),
            "study_a": study[i], "study_b": study[j],
        }
        for i, j, d in edges
        if d <= EDGE_BINS[0]
    ]  # fmt: skip
    (results / "qc").mkdir(parents=True, exist_ok=True)
    write_tsv(
        results / "qc" / "near_identical_pairs.tsv",
        near,
        ["strain_a", "strain_b", "distance", "isolate_a", "isolate_b", "study_a", "study_b"],
    )

    sensitivity = None
    if cfg.sensitivity_mask_sdrm:
        cols = sdrm_columns(settings, root)
        alt = _cluster(
            strains,
            [_mask(q, cols) for q in sequences],
            cfg.threshold,
            cfg.min_overlap,
            cfg.min_size,
        )
        sensitivity = {
            "note": "SDRM codons masked before TN93; the primary analysis leaves them in",
            "n_masked_columns": len(cols),
            "n_clusters": len(alt["groups"]),
            "n_clustered": len(alt["members"]),
            "proportion_clustered": round(len(alt["members"]) / len(strains), 4),
            "sdrm": _by_membership(rows, alt["members"]),
        }

    clustered = [r for r in rows if r["cluster_id"]]
    summary = {
        "scope": settings.scope,
        "threshold": cfg.threshold,
        "min_size": cfg.min_size,
        "min_overlap": cfg.min_overlap,
        "n_sequences": len(strains),
        "n_pairs_compared": primary["compared"],
        "n_clusters": len(groups),
        "n_clustered": len(clustered),
        "proportion_clustered": round(len(clustered) / len(strains), 4) if strains else None,
        "size_distribution": dict(sorted(Counter(len(g) for g in groups).items())),
        "largest_cluster": len(groups[0]) if groups else 0,
        "multi_study_clusters": sum(len({study[k] for k in g}) > 1 for g in groups),
        "multi_subtype_clusters": sum(
            len({meta[strains[k]].get("subtype_called") for k in g}) > 1 for g in groups
        ),
        "subtype_clustered": {
            s: {
                "n": sum(r["subtype_called"] == s for r in rows),
                "clustered": sum(r["subtype_called"] == s for r in clustered),
            }
            for s in sorted({r["subtype_called"] for r in rows})
        },
        "sdrm": _by_membership(rows, members),
        "clusters_sharing_an_sdrm": shared,
        "sdrm_by_treatment_context": by_context,
        "sdrm_cluster_odds_ratio_adjusted_for_treatment_context": mantel_haenszel(tables),
        "near_identical_pairs": {
            "n": len(near),
            "within_study": sum(r["study_a"] == r["study_b"] for r in near),
            "by_study": dict(Counter(r["study_a"] for r in near).most_common(8)),
        },
        "edges_by_distance": bins,
        "sensitivity_sdrm_masked": sensitivity,
    }
    (results / "phylo" / "clusters_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
