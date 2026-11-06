"""Choose the tree set: who goes into the phylogeny, and why the rest do not.

Order (each step reports how many it removed, so the waterfall adds up):

1. QC-passed sequences                      (``qc.passed.tsv``)
2. clone policy, then individual policy     (same policies as the prevalence population)
3. recombination policy                     (``recombination.tree_policy``)
4. date requirement                         (``subsample.date_requirement``)
5. stratified subsample                      (``subsample.*``; skipped when disabled)

The resistance layer is unaffected: it reads the full QC-passed set. Subsampling is for
tree inference only.

Outputs: ``data/processed/<scope>/aln.tree.fasta``, ``metadata.final.tsv``,
``results/<scope>/qc/subsample_strata.tsv`` and ``subsample_summary.json``.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from ..errors import ContractError
from ..io.fasta import read_fasta, write_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .clones import apply_population_policies
from .subsample import stratified_subsample

_DATED = {
    "strict": {"day", "month"},
    "year": {"day", "month", "year"},
}


def _truthy(value: str) -> bool:
    return value.strip().lower() == "true"


def run_subsample(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    scope = settings.scope
    interim = root / settings.paths["interim"] / scope
    processed = root / settings.paths["processed"] / scope
    results = root / settings.paths["results"] / scope
    needed = {
        "metadata.valid.tsv": "normalise",
        "qc.passed.tsv": "qc",
        "subtype_calls.tsv": "subtype",
        "recombination.tsv": "recombination",
        "aln.pol.masked.fasta": "mask",
    }
    for name, stage in needed.items():
        if not (interim / name).exists():
            raise ContractError(f"{interim / name} not found; run `virotrace {stage}` first")

    cfg = settings.subsample
    passed = {r["strain"] for r in read_tsv(interim / "qc.passed.tsv")}
    meta = [
        {**r, "clone_derived": _truthy(r.get("clone_derived") or "")}
        for r in read_tsv(interim / "metadata.valid.tsv")
        if r["strain"] in passed
    ]
    if len(meta) != len(passed):
        raise ContractError(
            f"subsample lost records: {len(passed)} passed QC, {len(meta)} in metadata"
        )
    waterfall: dict[str, int] = {"qc_passed": len(meta)}

    kept, removed = apply_population_policies(
        meta,
        clone_policy=settings.qc.clone_policy,
        individual_policy=settings.qc.individual_policy,
        patterns=settings.qc.individual_id_patterns,
    )
    waterfall["after_clone_policy"] = len(meta) - removed["clone_policy"]
    waterfall["after_individual_policy"] = len(kept)

    subtypes = {r["strain"]: r for r in read_tsv(interim / "subtype_calls.tsv")}
    flags = {r["strain"]: r for r in read_tsv(interim / "recombination.tsv")}
    candidates = []
    for row in kept:
        flag = flags.get(row["strain"])
        call = subtypes.get(row["strain"])
        if flag is None or call is None:
            raise ContractError(f"{row['strain']}: missing subtype or recombination row")
        if _truthy(flag["exclude_from_tree"]):
            continue
        candidates.append(
            {
                **row,
                "subtype_called": call["subtype_called"],
                "subtype_support": call["subtype_support"],
                "recomb_status": flag["recomb_status"],
            }
        )
    waterfall["after_recombination_policy"] = len(candidates)

    if cfg.date_requirement != "none":
        allowed = _DATED[cfg.date_requirement]
        candidates = [r for r in candidates if r["date_precision"] in allowed and r["year"]]
    waterfall["after_date_requirement"] = len(candidates)

    for row in candidates:  # the pure function wants typed fields
        row["year"] = int(row["year"]) if row["year"] else 0
        row["gene_coverage"] = 1.0
        row["ambiguity_frac"] = float(row["ambiguity_frac"] or 0.0)

    if cfg.enabled:
        retained, summaries = stratified_subsample(
            candidates,
            strata=cfg.strata,
            max_per_stratum=cfg.max_per_stratum,
            min_per_stratum=cfg.min_per_stratum,
            priority=cfg.priority,
            year_bin_width=cfg.year_bin_width,
        )
    else:
        retained = sorted(r["strain"] for r in candidates)
        summaries = []
    waterfall["retained"] = len(retained)
    keep = set(retained)
    final = sorted((r for r in candidates if r["strain"] in keep), key=lambda r: r["strain"])

    aln = {s: q for s, q in read_fasta(interim / "aln.pol.masked.fasta") if s in keep}
    if len(aln) != len(keep):
        raise ContractError(
            f"{len(keep) - len(aln)} tree-set sequences missing from the masked alignment"
        )
    write_fasta(
        processed / "aln.tree.fasta", [(r["strain"], aln[r["strain"]]) for r in final], wrap=0
    )
    write_tsv(processed / "metadata.final.tsv", final)
    write_tsv(
        results / "qc" / "subsample_strata.tsv",
        [
            {
                "stratum": "|".join(str(p) for p in s.key),
                "n_available": s.n_available,
                "n_retained": s.n_retained,
            }
            for s in summaries
        ],
        ["stratum", "n_available", "n_retained"],
    )

    summary = {
        "scope": scope,
        "date_requirement": cfg.date_requirement,
        "subsample_enabled": cfg.enabled,
        "waterfall": waterfall,
        "by_year": dict(sorted(Counter(r["year"] for r in final).items())),
        "by_subtype": dict(Counter(r["subtype_called"] for r in final).most_common()),
        "by_study": dict(Counter(r["study_key"] for r in final).most_common(10)),
        "year_span": [final[0]["year"], final[-1]["year"]] if final else None,
    }
    if final:
        years = sorted(r["year"] for r in final)
        summary["year_span"] = [years[0], years[-1]]
    (results / "qc").mkdir(parents=True, exist_ok=True)
    (results / "qc" / "subsample_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
