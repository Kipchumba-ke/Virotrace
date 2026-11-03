"""Provenance sidecars.

Every terminal output gets a ``<name>.provenance.json`` recording the config hash,
snapshot, tool versions, seeds and timestamp. Reproducibility is defined as
*same snapshot + same config -> same result*, and this file is what makes that
claim checkable after the fact. Provenance writing is not configurable off.
"""

from __future__ import annotations

import json
import platform
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from . import SCHEMA_VERSION, __version__

EXTERNAL_TOOLS = ("mafft", "iqtree2", "iqtree", "treetime", "augur", "snakemake")


def tool_versions(tools: tuple[str, ...] = EXTERNAL_TOOLS) -> dict[str, str]:
    """Resolved versions of external binaries, or ``'absent'`` when not installed."""
    versions: dict[str, str] = {}
    for tool in tools:
        path = shutil.which(tool)
        if path is None:
            versions[tool] = "absent"
            continue
        for flag in ("--version", "-version", "-v"):
            try:
                proc = subprocess.run(
                    [tool, flag], capture_output=True, text=True, timeout=20, check=False
                )
            except (OSError, subprocess.TimeoutExpired):
                continue
            output = (proc.stdout or proc.stderr).strip().splitlines()
            if output:
                versions[tool] = output[0][:120]
                break
        else:
            versions[tool] = "unknown"
    return versions


def build_provenance(
    *,
    stage: str,
    config_hash: str,
    snapshot_id: str | None,
    seed: int,
    extra: dict[str, Any] | None = None,
    include_tool_versions: bool = True,
) -> dict[str, Any]:
    return {
        "virotrace_version": __version__,
        "schema_version": SCHEMA_VERSION,
        "stage": stage,
        "config_hash": config_hash,
        "snapshot_id": snapshot_id,
        "seed": seed,
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "platform": f"{platform.system()} {platform.machine()} py{platform.python_version()}",
        "tool_versions": tool_versions() if include_tool_versions else {},
        **(extra or {}),
    }


def write_provenance(output_path: Path | str, record: dict[str, Any]) -> Path:
    """Write ``<output_path>.provenance.json`` beside the output it describes."""
    sidecar = Path(f"{output_path}.provenance.json")
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return sidecar
