"""Subtype assignment by bootstrap-supported nearest-reference placement.

Two tiers. ``subtype_called`` is placed against the **pure subtypes** only (A1..K): in pol
most CRF references sit very close to a pure subtype (CRF50_A1D next to A1, CRF41_CD next to
D), so letting them compete splits the bootstrap and leaves half of a Kenyan corpus
``unassigned``. ``crf_label`` is a second, annotation-only placement against the full panel,
set only when a CRF is the supported nearest label.

Each QC-passed sequence is compared with a curated subtype reference panel (one or more
sequences per subtype/CRF, all projected onto the HXB2 pol coordinates the DRM layer uses).
Distance is the p-distance over sites where both sequences have an unambiguous base. The
call is the label of the nearest panel sequence; its **support** is the share of site
bootstrap replicates that give the same label, which plays the role a clade bootstrap plays
on a reference tree. A call below ``subtype.min_support``, or with too few comparable
sites, is reported as ``unassigned`` rather than guessed.

This is distance-based placement, not tree placement: cheaper, deterministic and testable,
but it can mislabel a sequence whose true subtype is missing from the panel (it takes the
nearest label with low support, or, for a mosaic, the dominant parent -- see
``recombination``). A second, independent call (an established tool) is still planned; until
then the depositor's own label is compared and disagreement is reported, not resolved.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..align.mafft import Runner, align_chunk, chunked, mafft_runner
from ..align.project import project_to_reference
from ..errors import ContractError
from ..io.fasta import read_fasta, write_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings

UNASSIGNED = "unassigned"
CALL_COLUMNS = [
    "strain",
    "subtype_called",
    "subtype_support",
    "n_sites",
    "nearest_distance",
    "runner_up",
    "runner_up_distance",
    "crf_label",
    "crf_support",
    "subtype_reported",
    "subtype_agreement",
]

_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b)] = _i
    _LUT[ord(_b.lower())] = _i


def encode(sequence: str) -> np.ndarray:
    """A,C,G,T -> 0..3; everything else (gap, N, IUPAC) -> 4, i.e. not comparable."""
    return _LUT[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]


def panel_label(header: str) -> str:
    """Subtype label of a LANL header: ``A1.KE.00.name.ACC``, ``01_AE.TH...`` or ``Ref.A1.AU...``.

    A leading ``Ref`` field (the LANL reference-alignment download) is skipped.
    """
    fields = header.split(".")
    head = (
        fields[1] if fields[0].strip().lower() == "ref" and len(fields) > 1 else fields[0]
    ).strip()
    if not head:
        raise ContractError(f"panel header {header!r} has no subtype label")
    return f"CRF{head}" if re.match(r"^\d{2}_", head) else head


_PURE = re.compile(r"^[A-L][0-9]?$")


def is_pure(label: str) -> bool:
    """A pure subtype or sub-subtype (A1, B, F2, K...), as opposed to a CRF or other mosaic."""
    return bool(_PURE.match(label))


def parent_group(label: str) -> str:
    """Sub-subtypes collapse to their group for the recombination scan: A1..A8 -> A, F1/F2 -> F."""
    return label[0] if is_pure(label) else label


def pure_rows(rows: list[tuple[str, str]]) -> list[tuple[str, str]]:
    return [(name, seq) for name, seq in rows if is_pure(panel_label(name))]


@dataclass(frozen=True)
class Panel:
    names: list[str]
    labels: list[str]  # per sequence, sorted so each label is contiguous
    matrix: np.ndarray  # n x L uint8
    label_names: list[str]  # distinct labels in order
    starts: np.ndarray  # first row of each label, for np.minimum.reduceat

    @property
    def length(self) -> int:
        return int(self.matrix.shape[1])


def build_panel(rows: list[tuple[str, str]]) -> Panel:
    if not rows:
        raise ContractError("subtype reference panel is empty")
    width = {len(seq) for _, seq in rows}
    if len(width) != 1:
        raise ContractError(f"panel rows have different lengths: {sorted(width)}")
    ordered = sorted(rows, key=lambda r: (panel_label(r[0]), r[0]))
    labels = [panel_label(name) for name, _ in ordered]
    label_names = list(dict.fromkeys(labels))
    if len(label_names) < 2:
        raise ContractError("the subtype panel needs at least two distinct labels")
    starts = np.array([labels.index(label) for label in label_names])
    return Panel(
        names=[name for name, _ in ordered],
        labels=labels,
        matrix=np.stack([encode(seq) for _, seq in ordered]),
        label_names=label_names,
        starts=starts,
    )


@dataclass(frozen=True)
class Placement:
    label: str
    support: float
    n_sites: int
    distance: float
    runner_up: str
    runner_up_distance: float


def _label_distances(num: np.ndarray, den: np.ndarray, panel: Panel, min_sites: int) -> np.ndarray:
    """Per-label minimum distance; panel members with too few shared sites are ignored."""
    with np.errstate(divide="ignore", invalid="ignore"):
        dist = np.where(den >= min_sites, num / np.maximum(den, 1), np.inf)
    return np.minimum.reduceat(dist, panel.starts, axis=-1)


def place(
    query: np.ndarray,
    panel: Panel,
    *,
    n_bootstrap: int,
    min_sites: int,
    rng: np.random.Generator,
    columns: np.ndarray | None = None,
    exclude_row: int | None = None,
) -> Placement | None:
    """Nearest-label placement of one encoded query; ``None`` if too few comparable sites.

    ``columns`` restricts the comparison to those alignment columns (used for segments).
    """
    usable = query < 4
    if columns is not None:
        mask = np.zeros_like(usable)
        mask[columns] = True
        usable &= mask
    idx = np.flatnonzero(usable)
    if idx.size < min_sites:
        return None
    sub = panel.matrix[:, idx]
    shared = sub < 4
    if exclude_row is not None:  # leave-one-out: the query must not find itself
        shared = shared.copy()
        shared[exclude_row] = False
    mism = ((sub != query[idx]) & shared).astype(np.float32)
    shared_f = shared.astype(np.float32)

    num, den = mism.sum(axis=1), shared_f.sum(axis=1)
    label_dist = _label_distances(num, den, panel, min_sites)
    if not np.isfinite(label_dist).any():
        return None
    order = np.argsort(label_dist)
    best = int(order[0])
    second = int(order[1]) if len(order) > 1 and np.isfinite(label_dist[order[1]]) else best

    if n_bootstrap > 0:
        weights = rng.multinomial(idx.size, np.full(idx.size, 1.0 / idx.size), size=n_bootstrap)
        weights = weights.astype(np.float32)
        boot = _label_distances(weights @ mism.T, weights @ shared_f.T, panel, min_sites)
        winners = np.argmin(boot, axis=1)
        support = float(np.mean(winners == best))
    else:
        support = 1.0
    return Placement(
        label=panel.label_names[best],
        support=support,
        n_sites=int(idx.size),
        distance=float(label_dist[best]),
        runner_up=panel.label_names[second] if second != best else "",
        runner_up_distance=float(label_dist[second]) if second != best else float("nan"),
    )


_REPORTED = re.compile(r"^(?:subtype[_\- ]*)?([A-Za-z0-9_]+)$")


def agreement(called: str, reported: str) -> str:
    """Compare the call with the depositor's label: agree | disagree | unreported | not_comparable.

    A depositor label of ``A`` agrees with ``A1`` and its sub-subtypes (and vice versa for the
    call ``A``): many depositors wrote the group, not the sub-subtype. CRF labels are compared
    after dropping a leading ``CRF`` and underscores.
    """
    reported = (reported or "").strip()
    if not reported:
        return "unreported"
    if called == UNASSIGNED:
        return "not_comparable"

    def norm(label: str) -> str:
        label = label.upper().replace("CRF", "").replace("_", "").replace("-", "")
        return label

    a, b = norm(called), norm(reported)
    if a == b:
        return "agree"
    pure = re.compile(r"^[A-K]$")
    if pure.match(a) and b.startswith(a) and b[1:].isdigit():
        return "agree"
    if pure.match(b) and a.startswith(b) and a[1:].isdigit():
        return "agree"
    return "disagree"


def select_panel(rows: list[tuple[str, str]], cfg: Any) -> list[tuple[str, str]]:
    """Drop labels that are not targets: non-M groups, and CRFs outside ``include_crfs``.

    ``include_crfs`` of ``None`` keeps every CRF. Pure subtypes are always kept.
    """
    excluded = {label.upper() for label in cfg.exclude_labels}
    only = (
        None
        if cfg.include_crfs is None
        else {c.upper().replace("CRF", "") for c in cfg.include_crfs}
    )
    kept = []
    for name, seq in rows:
        label = panel_label(name)
        bare = label.upper().replace("CRF", "")
        if label.upper() in excluded or bare in excluded:
            continue
        if label.startswith("CRF") and only is not None and bare not in only:
            continue
        kept.append((name, seq))
    return kept


def project_panel(
    settings: Settings, root: Path, runner: Runner = mafft_runner
) -> list[tuple[str, str]]:
    """Project the panel onto HXB2 pol coordinates exactly as the query sequences are."""
    from ..drm.reference import load_reference_sequence

    pol_start, pol_end = settings.reference.genes.pol
    reference_dir = root / settings.paths["reference"]
    panel_path = reference_dir / settings.subtype.reference_panel
    if not panel_path.exists() or panel_path.stat().st_size == 0:
        raise ContractError(
            f"subtype reference panel not found or empty: {panel_path}. "
            "Download the LANL HIV-1 subtype reference sequences (pol region or full genomes) "
            "and save them there with headers like 'A1.KE.00.name.ACCESSION'"
        )
    reference_pol = load_reference_sequence(reference_dir / settings.reference.record)[
        pol_start - 1 : pol_end
    ]
    reference = (f"{settings.reference.accession}_pol", reference_pol)
    out: list[tuple[str, str]] = []
    wanted = select_panel(list(read_fasta(panel_path)), settings.subtype)
    if not wanted:
        raise ContractError(f"{panel_path}: no panel sequences left after the label filters")
    for chunk in chunked(wanted, 100):
        aligned_ref, rows = align_chunk(reference, chunk, threads=settings.threads, runner=runner)
        for name, _ in chunk:
            out.append(
                (
                    name,
                    project_to_reference(
                        aligned_ref, rows[name], start=1, end=len(reference_pol)
                    ).projected,
                )
            )
    return out


def panel_self_check(
    panel: Panel, *, n_bootstrap: int, min_sites: int, min_support: float, rng: np.random.Generator
) -> dict[str, Any]:
    """Each reference sequence, placed against the rest of the panel, should find its label.

    The methods contract: "reference sequences of known subtype must recover their subtype".
    Labels with a single member are skipped (nothing left to recover them). A failure is a
    reported finding about the panel (near-identical labels, a mislabelled entry), not an error.
    """
    counts = Counter(panel.labels)
    tested = recovered = 0
    failures: list[dict[str, Any]] = []
    for row, (name, label) in enumerate(zip(panel.names, panel.labels, strict=True)):
        if counts[label] < 2:
            continue
        tested += 1
        placed = place(
            panel.matrix[row],
            panel,
            n_bootstrap=n_bootstrap,
            min_sites=min_sites,
            rng=rng,
            exclude_row=row,
        )
        ok = placed is not None and placed.label == label and placed.support >= min_support
        recovered += ok
        if not ok:
            failures.append(
                {
                    "name": name,
                    "expected": label,
                    "got": None if placed is None else placed.label,
                    "support": None if placed is None else round(placed.support, 2),
                }
            )
    return {
        "n_tested": tested,
        "n_recovered": recovered,
        "recovery_rate": round(recovered / tested, 4) if tested else None,
        "failures": failures[:40],
        "n_failures": len(failures),
    }


def run_subtype(
    settings: Settings, *, repo_root: Path | str = ".", runner: Runner = mafft_runner
) -> dict[str, Any]:
    root = Path(repo_root)
    scope = settings.scope
    interim = root / settings.paths["interim"] / scope
    results = root / settings.paths["results"] / scope
    aln_path, passed_path = interim / "aln.pol.unmasked.fasta", interim / "qc.passed.tsv"
    for path, stage in ((aln_path, "align"), (passed_path, "qc")):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    cfg = settings.subtype
    panel_rows = project_panel(settings, root, runner)
    write_fasta(interim / "panel.pol.fasta", panel_rows, wrap=0)
    panel = build_panel(pure_rows(panel_rows))  # tier 1: pure subtypes decide `subtype_called`
    crf_panel = build_panel(panel_rows)  # tier 2: with CRFs, only to annotate `crf_label`

    passed = {r["strain"] for r in read_tsv(passed_path)}
    reported = {
        r["strain"]: r.get("subtype_reported", "") for r in read_tsv(interim / "metadata.valid.tsv")
    }
    rng = np.random.default_rng(settings.stage_seed("subtype"))

    rows: list[dict[str, Any]] = []
    for strain, sequence in read_fasta(aln_path):
        if strain not in passed:
            continue
        if len(sequence) != panel.length:
            raise ContractError(
                f"{strain}: alignment width {len(sequence)} != panel {panel.length}"
            )
        placed = place(
            encode(sequence),
            panel,
            n_bootstrap=cfg.n_bootstrap,
            min_sites=cfg.min_overlap_sites,
            rng=rng,
        )
        called = UNASSIGNED if placed is None or placed.support < cfg.min_support else placed.label
        crf_label, crf_support = "", ""
        if cfg.annotate_crfs and placed is not None:
            with_crfs = place(
                encode(sequence),
                crf_panel,
                n_bootstrap=cfg.n_bootstrap,
                min_sites=cfg.min_overlap_sites,
                rng=rng,
            )
            if with_crfs is not None:
                crf_support = round(with_crfs.support, 3)
                if not is_pure(with_crfs.label) and with_crfs.support >= cfg.min_support:
                    crf_label = with_crfs.label
        reported_label = reported.get(strain, "") or ""
        agree = agreement(called, reported_label)
        if agree != "agree" and crf_label and agreement(crf_label, reported_label) == "agree":
            agree = "agree"
        rows.append(
            {
                "strain": strain,
                "subtype_called": called,
                "subtype_support": "" if placed is None else round(placed.support, 3),
                "n_sites": 0 if placed is None else placed.n_sites,
                "nearest_distance": "" if placed is None else round(placed.distance, 4),
                "runner_up": "" if placed is None else placed.runner_up,
                "runner_up_distance": (
                    ""
                    if placed is None or not placed.runner_up
                    else round(placed.runner_up_distance, 4)
                ),
                "crf_label": crf_label,
                "crf_support": crf_support,
                "subtype_reported": reported_label,
                "subtype_agreement": agree,
            }
        )
    if len(rows) != len(passed):
        raise ContractError(f"subtype lost records: {len(passed)} passed QC, {len(rows)} called")
    write_tsv(interim / "subtype_calls.tsv", rows, CALL_COLUMNS)

    summary = {
        "scope": scope,
        "n_sequences": len(rows),
        "panel": {
            "n_sequences": len(panel_rows),
            "n_pure": len(panel.names),
            "pure_labels": dict(Counter(panel.labels)),
            "n_crf_labels": len(set(crf_panel.label_names) - set(panel.label_names)),
        },
        "crf_annotated": dict(
            Counter(r["crf_label"] for r in rows if r["crf_label"]).most_common()
        ),
        "panel_self_check_with_crfs": panel_self_check(
            crf_panel,
            n_bootstrap=cfg.n_bootstrap,
            min_sites=cfg.min_overlap_sites,
            min_support=cfg.min_support,
            rng=np.random.default_rng(settings.stage_seed("subtype_self_check_crf")),
        ),
        "panel_self_check": panel_self_check(
            panel,
            n_bootstrap=cfg.n_bootstrap,
            min_sites=cfg.min_overlap_sites,
            min_support=cfg.min_support,
            rng=np.random.default_rng(settings.stage_seed("subtype_self_check")),
        ),
        "min_support": cfg.min_support,
        "n_bootstrap": cfg.n_bootstrap,
        "called": dict(Counter(r["subtype_called"] for r in rows).most_common()),
        "agreement": dict(Counter(r["subtype_agreement"] for r in rows)),
        # which pairs drive the disagreements: "called|reported" -> count
        "disagreement_pairs": dict(
            Counter(
                f"{r['subtype_called']}|{r['subtype_reported']}"
                for r in rows
                if r["subtype_agreement"] == "disagree"
            ).most_common(25)
        ),
        "unassigned_by_n_sites": {
            "under_600": sum(
                1 for r in rows if r["subtype_called"] == UNASSIGNED and int(r["n_sites"]) < 600
            ),
            "600_or_more": sum(
                1 for r in rows if r["subtype_called"] == UNASSIGNED and int(r["n_sites"]) >= 600
            ),
        },
    }
    (results / "qc").mkdir(parents=True, exist_ok=True)
    (results / "qc" / "subtype_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
