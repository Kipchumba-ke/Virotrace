"""Snapshot manifests.

A manifest records exactly which queries produced a corpus and the SHA256 of every
file, so "same snapshot + same config -> same result" is checkable. Public databases
grow, so a rerun months later will not reproduce the same corpus; the manifest is
what makes that visible rather than silent.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .. import __version__


def sha256_file(path: Path | str, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def make_snapshot_id(date: str, queries: list[str]) -> str:
    """``YYYY-MM-DD_<6 hex>``: the date plus a short hash of the source queries.

    Sorted, so the id does not depend on adapter order.
    """
    payload = "\n".join(sorted(queries)).encode()
    return f"{date}_{hashlib.sha256(payload).hexdigest()[:6]}"


@dataclass(frozen=True, slots=True)
class SourceRecordInfo:
    name: str
    query: str
    retrieved: int
    access_tier: int = 1


def build_manifest(
    *,
    snapshot_id: str,
    created_utc: str,
    sources: list[SourceRecordInfo],
    files: dict[str, Path],
    config_hash: str | None = None,
) -> dict[str, Any]:
    return {
        "snapshot_id": snapshot_id,
        "created_utc": created_utc,
        "virotrace_version": __version__,
        "config_hash": config_hash,
        "sources": [
            {
                "name": s.name,
                "query": s.query,
                "retrieved": s.retrieved,
                "access_tier": s.access_tier,
            }
            for s in sources
        ],
        "files": [
            {"path": name, "sha256": sha256_file(path)} for name, path in sorted(files.items())
        ],
    }


def write_manifest(path: Path | str, manifest: dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return path


def verify_manifest(manifest_path: Path | str) -> list[str]:
    """Return a list of problems (empty = every listed file matches its checksum)."""
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    problems: list[str] = []
    for entry in manifest["files"]:
        target = manifest_path.parent / entry["path"]
        if not target.exists():
            problems.append(f"missing: {entry['path']}")
        elif sha256_file(target) != entry["sha256"]:
            problems.append(f"checksum mismatch: {entry['path']}")
    return problems
