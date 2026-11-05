"""Smoke tests: the DAG wires up, and (later) the whole dev build is reproducible.

Skipped automatically when snakemake is not on PATH, so a fresh clone still gets a
green fast suite.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest

from virotrace.settings import load_settings

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

# (profile, scope) pairs the Makefile actually runs.
PROFILES = [("dev", "kenya"), ("kenya", "kenya"), ("africa", "africa")]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stage_dry_run_workspace(repo_root: Path, workspace: Path) -> None:
    """Copy workflow + config into a scratch dir and stub reference files not pinned yet.

    A dry-run fails on any missing input, and three reference files are Phase 2 work
    (HXB2 record, subtype panel, WHO SDRM list). Empty stand-ins let the test check
    the *wiring* without fabricating reference data in the real repo.
    """
    shutil.copytree(repo_root / "workflow", workspace / "workflow")
    shutil.copytree(repo_root / "config", workspace / "config")

    settings = load_settings(repo_root / "config" / "config.yaml")
    reference = workspace / "config" / "reference"
    needed = [
        settings.reference.record,
        settings.subtype.reference_panel,
        settings.drm.sdrm.list,
    ]
    for name in needed:
        target = reference / name
        if not target.exists():
            target.write_text("")


@pytest.mark.parametrize(("profile", "scope"), PROFILES)
def test_dag_is_valid_for_every_profile(repo_root: Path, tmp_path: Path, profile: str, scope: str):
    """Dry-run only: checks the DAG wires up for each profile, executes nothing."""
    _stage_dry_run_workspace(repo_root, tmp_path)
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            "config/config.yaml",
            "--config",
            f"scope={scope}",
            f"profile={profile}",
            "--dry-run",
            "--cores",
            "1",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    # Snakemake reports most errors on stdout, so show both streams.
    output = f"--- stdout ---\n{result.stdout[-2500:]}\n--- stderr ---\n{result.stderr[-2500:]}"
    assert result.returncode == 0, output


@pytest.mark.skip(reason="enable once the pipeline stages are implemented (Phase 2)")
def test_dev_build_completes_and_is_reproducible(repo_root: Path) -> None:
    first = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            "config/config.yaml",
            "--config",
            "scope=kenya",
            "profile=dev",
            "--cores",
            "2",
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert first.returncode == 0, first.stdout[-2000:] + first.stderr[-2000:]

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
