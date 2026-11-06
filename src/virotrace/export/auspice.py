"""Auspice export: one dataset for the pooled tree and one for each per-subtype tree.

For each tree this runs ``augur ancestral`` (nucleotide mutations on branches) and
``augur export v2``. Tip metadata is assembled here from the stage outputs, so the colour-by
choices (subtype, cluster, study, treatment context, SDRM status per drug class) come from
the same tables the reports use, not from a second copy of the logic.

Clusters are the descriptive TN93 groups from ``phylo clusters``; they are shown as a colouring
and carry no claim about who infected whom. The pooled tree is a midpoint-rooted divergence
tree unless its temporal gate passed; each dataset's title says which it is.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..curate.study import load_study_annotations
from ..drm.classes import DRUG_CLASSES, drug_class_of
from ..errors import ContractError
from ..io.tables import read_tsv, write_tsv
from ..phylo.stage import Runner, augur_date, subprocess_runner
from ..settings import Settings

CLASS_COLUMNS = {c: f"sdrm_{c.lower()}" for c in DRUG_CLASSES}

DESCRIPTION = """\
HIV-1 *pol* sequences from public databases, one sequence per individual after the
clone policy. This is a sample of published sequences, not a population survey.

* **Clusters** are groups of sequences within a TN93 distance threshold, reported as a
  descriptive statistic. They are not transmission chains.
* **SDRM** status uses the WHO 2009 surveillance list. Insertions and deletions on the
  list are not detected.
* **Treatment context** comes from hand-curated study annotations; "unannotated" means the
  study has not been reviewed.
* Dated trees give rough timelines. Estimated rates are low and branch dates carry wide
  uncertainty; a divergence tree has no time axis.
"""


def _paths(settings: Settings, root: Path) -> dict[str, Path]:
    return {
        "processed": root / settings.paths["processed"] / settings.scope,
        "results": root / settings.paths["results"] / settings.scope,
        "auspice": root / settings.paths["auspice"],
        "reference": root / settings.paths["reference"],
    }


def tip_metadata(settings: Settings, root: Path | str = ".") -> list[dict[str, Any]]:
    root = Path(root)
    p = _paths(settings, root)
    needed = {
        p["processed"] / "metadata.final.tsv": "subsample",
        p["results"] / "phylo" / "clusters.tsv": "phylo clusters",
        p["results"] / "drm" / "drm_table.tsv": "drm call",
    }
    for path, stage in needed.items():
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    clusters = {r["strain"]: r for r in read_tsv(p["results"] / "phylo" / "clusters.tsv")}
    annotation_path = p["reference"] / settings.prevalence.study_annotations
    annotations = load_study_annotations(annotation_path) if annotation_path.exists() else {}

    found: dict[str, dict[str, set[str]]] = {}
    for r in read_tsv(p["results"] / "drm" / "drm_table.tsv"):
        if not r.get("sdrm_residues"):
            continue
        klass = drug_class_of(r["gene"], int(r["position"]))
        mutations = {
            f"{r['gene']}:{r['ref_aa']}{r['position']}{x}" for x in r["sdrm_residues"].split(",")
        }
        found.setdefault(r["strain"], {}).setdefault(klass, set()).update(mutations)

    rows = []
    for m in read_tsv(p["processed"] / "metadata.final.tsv"):
        strain = m["strain"]
        cluster = clusters.get(strain, {})
        by_class = found.get(strain, {})
        row: dict[str, Any] = {
            "strain": strain,
            "country": m.get("country") or "unknown",
            "region": m.get("region") or "unknown",
            "year": m.get("year") or "",
            "date": augur_date(m.get("date", ""), m.get("date_precision", "")),
            "subtype": m.get("subtype_called") or "unassigned",
            "study": m.get("study_key") or "unknown",
            "treatment_context": annotations.get(m.get("study_key", ""), {}).get(
                "treatment_context", "unannotated"
            ),
            "cluster": cluster.get("cluster_id") or "unclustered",
            "cluster_size": cluster.get("cluster_size") or "1",
            "sdrm_any": "yes" if by_class else "no",
        }
        for klass, column in CLASS_COLUMNS.items():
            row[column] = "yes" if klass in by_class else "no"
        row["sdrm_mutations"] = ",".join(sorted(x for v in by_class.values() for x in v)) or "none"
        rows.append(row)
    return rows


def _coloring(key: str, title: str, kind: str = "categorical") -> dict[str, str]:
    return {"key": key, "title": title, "type": kind}


def auspice_config(title: str, columns: list[str]) -> dict[str, Any]:
    candidates = [
        _coloring("subtype", "Subtype (called)"),
        _coloring("cluster", "TN93 cluster (descriptive)"),
        _coloring("study", "Study"),
        _coloring("treatment_context", "Treatment context (study-level)"),
        _coloring("sdrm_any", "Any SDRM"),
        *(_coloring(c, f"SDRM: {k}") for k, c in CLASS_COLUMNS.items()),
        _coloring("region", "Region"),
        _coloring("country", "Country"),
        _coloring("year", "Year", "ordinal"),
    ]
    present = [c for c in candidates if c["key"] in columns]
    return {
        "title": title,
        "maintainers": [{"name": "Virotrace"}],
        "colorings": present,
        "display_defaults": {"color_by": "subtype", "distance_measure": "div"},
        "filters": [c["key"] for c in present if c["key"] not in ("year", "cluster")]
        + ["sdrm_mutations"],
        "panels": ["tree", "entropy"],
    }


def _views(settings: Settings, p: dict[str, Path]) -> list[tuple[str | None, Path]]:
    pooled = p["results"] / "phylo"
    return [(None, pooled)] + [(s, pooled / s) for s in settings.phylo.refine.subtype_trees]


def run_export(
    settings: Settings, *, repo_root: Path | str = ".", runner: Runner = subprocess_runner
) -> dict[str, Any]:
    root = Path(repo_root)
    p = _paths(settings, root)
    rows = tip_metadata(settings, root)
    columns = list(rows[0]) if rows else []
    phylo = p["results"] / "phylo"
    metadata = phylo / "auspice_metadata.tsv"
    write_tsv(metadata, rows, columns)
    p["auspice"].mkdir(parents=True, exist_ok=True)
    description = p["auspice"] / "description.md"
    description.write_text(DESCRIPTION)

    outputs: dict[str, Any] = {}
    for subtype, folder in _views(settings, p):
        label = subtype or "pooled"
        tree, node_data = folder / "tree.refined.nwk", folder / "branch_lengths.json"
        for path in (tree, node_data, folder / "temporal_gate.json"):
            if not path.exists():
                step = "phylo refine" if subtype is None else "phylo subtrees"
                raise ContractError(f"{path} not found; run `virotrace {step}` first")
        aln = p["processed"] / (
            "aln.tree.fasta" if subtype is None else f"aln.tree.{subtype}.fasta"
        )
        if not aln.exists():
            raise ContractError(f"{aln} not found; run `virotrace subsample` first")
        gate = json.loads((folder / "temporal_gate.json").read_text())
        kind = "time tree" if gate.get("time_tree") else "divergence tree (midpoint-rooted)"
        scope = settings.scope.capitalize()
        title = f"Virotrace: {scope} HIV-1 pol, " + (
            f"subtype {subtype}, {kind}" if subtype else f"all subtypes, {kind}"
        )

        config_path = folder / "auspice_config.json"
        config_path.write_text(
            json.dumps(auspice_config(title, columns), indent=2, sort_keys=True) + "\n"
        )
        muts = folder / "nt_muts.json"
        name = f"virotrace_{settings.scope}" + (f"_{subtype}" if subtype else "")
        output = p["auspice"] / f"{name}.json"
        for stale in (muts, output):
            stale.unlink(missing_ok=True)

        runner([
            "augur", "ancestral", "--tree", str(tree), "--alignment", str(aln),
            "--inference", "joint", "--output-node-data", str(muts),
        ])  # fmt: skip
        runner([
            "augur", "export", "v2", "--tree", str(tree), "--metadata", str(metadata),
            "--metadata-id-columns", "strain", "--node-data", str(node_data), str(muts),
            "--auspice-config", str(config_path), "--title", title,
            "--description", str(description), "--output", str(output),
        ])  # fmt: skip
        outputs[label] = {"dataset": str(output), "title": title, "time_tree": gate["time_tree"]}
    return {"scope": settings.scope, "datasets": outputs, "n_tips": len(rows)}
