"""The normalise stage: snapshot -> dedup -> scope filter -> validated metadata.

Outputs (per scope):
    interim/<scope>/metadata.valid.tsv      one row per retained sequence
    interim/<scope>/seqs.dedup.fasta        sequences keyed by strain
    results/<scope>/qc/quarantine.tsv       records that failed the contract, with reason
    results/<scope>/qc/dedup_log.tsv        every collapse and flag, with the survivor
    results/<scope>/qc/out_of_scope.tsv     valid records outside the scope (counted, not lost)
    results/<scope>/qc/corpus_counts.json   the P1 gate numbers
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..acquire.snapshot import load_snapshot_records, resolve_snapshot
from ..errors import ContractError
from ..io.fasta import write_fasta
from ..io.tables import write_tsv
from ..settings import Settings
from .clones import is_excluded, split_by_policy
from .corpus import count_corpus
from .dedup import deduplicate
from .normalise import normalise_records


def run_normalise(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    scope_name = settings.scope
    scope = settings.scopes[scope_name]
    codes, regions = scope.get("countries", []), scope.get("regions", [])

    snapshot = resolve_snapshot(root / settings.paths["raw"], scope_name, settings.snapshot_id)
    records = load_snapshot_records(snapshot)
    sequences_by_accession = {r.accession: r.sequence for r in records}

    dedup = deduplicate(records)

    dates = (settings.model_extra or {}).get("curate", {}).get("dates", {})
    try:
        result = normalise_records(
            dedup.kept,
            snapshot_id=snapshot.name,
            max_record_error_rate=settings.qc.max_record_error_rate,
            min_year=int(dates.get("min_year", 1985)),
            max_year=int(dates.get("max_year", 2026)),
        )
    except ContractError as exc:
        # Keep the evidence: an abort must not throw away the reasons for it.
        quarantined = getattr(exc, "quarantine", None)
        if quarantined:
            qc_failed = root / settings.paths["results"] / scope_name / "qc"
            write_tsv(qc_failed / "quarantine.tsv", quarantined, ["accession", "rule", "reason"])
        raise

    clone_policy = settings.qc.clone_policy
    for row in result.rows:
        # Recorded per row so every downstream stage filters on one agreed column
        # instead of re-deriving the policy.
        row["clone_policy_excluded"] = is_excluded(row, clone_policy)

    in_scope = [r for r in result.rows if r["country_code"] in codes or r["region"] in regions]
    out_of_scope = [
        {"strain": r["strain"], "country_code": r["country_code"], "region": r["region"]}
        for r in result.rows
        if not (r["country_code"] in codes or r["region"] in regions)
    ]

    interim = root / settings.paths["interim"] / scope_name
    qc_dir = root / settings.paths["results"] / scope_name / "qc"

    write_tsv(interim / "metadata.valid.tsv", in_scope)
    acc_by_strain = {r["strain"]: r["accession"] for r in in_scope}
    write_fasta(
        interim / "seqs.dedup.fasta",
        ((strain, sequences_by_accession[acc]) for strain, acc in sorted(acc_by_strain.items())),
    )
    write_tsv(qc_dir / "quarantine.tsv", result.quarantine, ["accession", "rule", "reason"])
    write_tsv(
        qc_dir / "dedup_log.tsv",
        [
            {
                "accession": d.accession,
                "kept_instead": d.kept_instead,
                "level": d.level,
                "detail": d.detail,
            }
            for d in dedup.log
        ],
        ["accession", "kept_instead", "level", "detail"],
    )
    write_tsv(qc_dir / "out_of_scope.tsv", out_of_scope, ["strain", "country_code", "region"])

    target_code = codes[0] if codes else "KEN"
    counts = count_corpus(in_scope, country_code=target_code)
    after_policy, _ = split_by_policy(in_scope, clone_policy)
    summary: dict[str, Any] = {
        "snapshot_id": snapshot.name,
        "scope": scope_name,
        "n_snapshot": len(records),
        "n_collapsed_by_dedup": sum(1 for d in dedup.log if d.level != "flag-isolate"),
        "n_flagged_isolates": sum(1 for d in dedup.log if d.level == "flag-isolate"),
        "n_quarantined": len(result.quarantine),
        "n_out_of_scope": len(out_of_scope),
        "n_retained": len(in_scope),
        "clone_policy": clone_policy,
        "n_clone_derived": sum(1 for r in in_scope if r["clone_derived"]),
        "n_after_clone_policy": len(after_policy),
        "corpus": counts,
        "corpus_after_clone_policy": count_corpus(after_policy, country_code=target_code),
    }
    accounted = (
        summary["n_collapsed_by_dedup"]
        + summary["n_quarantined"]
        + summary["n_out_of_scope"]
        + summary["n_retained"]
    )
    if accounted != summary["n_snapshot"]:
        raise ContractError(
            f"record conservation violated across the normalise stage: snapshot has "
            f"{summary['n_snapshot']} records but {accounted} are accounted for"
        )

    qc_dir.mkdir(parents=True, exist_ok=True)
    (qc_dir / "corpus_counts.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
