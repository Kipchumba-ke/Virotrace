"""Smoke test: the entire dev-profile DAG, end to end.

Asserts the reproducibility claim rather than leaving it in the README: a rerun must
produce identical checksums for deterministic outputs.

Skipped automatically until the pipeline stages are implemented and snakemake is on
PATH, so a fresh clone still gets a green fast suite.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = [
    pytest.mark.smoke,
    pytest.mark.needs_tools,
    pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed"),
]

DETERMINISTIC_OUTPUTS = (
    "results/kenya/drm/drm_table.tsv",
    "results/kenya/drm/prevalence.tsv",
    "results/kenya/qc/excluded.tsv",
    "results/kenya/qc/subsample_strata.tsv",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_dev_build(repo_root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            "config/config.yaml",
            "--cores",
            "2",
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.skip(reason="enable once the pipeline stages are implemented (Phase 2)")
def test_dev_build_completes_and_is_reproducible(repo_root: Path) -> None:
    first = _run_dev_build(repo_root)
    assert first.returncode == 0, first.stderr[-2000:]

    checksums = {}
    for relative in DETERMINISTIC_OUTPUTS:
        path = repo_root / relative
        assert path.exists() and path.stat().st_size > 0, f"missing or empty: {relative}"
        assert Path(f"{path}.provenance.json").exists(), f"no provenance sidecar for {relative}"
        checksums[relative] = _sha256(path)

    report = repo_root / "reports/kenya/index.html"
    assert report.exists()
    # report.include_limitations cannot be configured off; this proves the code honours it.
    assert "Limitations" in report.read_text()

    second = _run_dev_build(repo_root)
    assert second.returncode == 0, second.stderr[-2000:]
    for relative, digest in checksums.items():
        assert _sha256(repo_root / relative) == digest, f"not reproducible: {relative}"


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_dag_is_valid_for_every_profile(repo_root: Path) -> None:
    """Dry-run only: checks the DAG wires up, executes nothing."""
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            "config/config.yaml",
            "--dry-run",
            "--cores",
            "1",
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-3000:]
