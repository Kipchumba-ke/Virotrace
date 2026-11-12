"""The acquire stage: settings -> adapters -> an immutable snapshot + a scope pointer."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..curate.geo import countries_for_scope
from ..errors import ConfigError
from ..settings import Settings
from .base import SourceAdapter
from .fixture import FixtureAdapter
from .genbank import GenBankAdapter, build_term
from .lanl import LanlAdapter
from .snapshot import take_snapshot, write_pointer


def _section(settings: Settings, name: str) -> dict[str, Any]:
    return dict((settings.model_extra or {}).get(name, {}))


def build_adapters(
    settings: Settings, *, max_records: int | None = None, repo_root: Path | str = "."
) -> list[SourceAdapter]:
    acquire = _section(settings, "acquire")
    scope = settings.scopes[settings.scope]
    adapters: list[SourceAdapter] = []

    for name in acquire.get("sources", []):
        if name == "fixture":
            adapters.append(
                FixtureAdapter(
                    Path(repo_root) / acquire.get("fixture_path", "tests/fixtures/records_toy.tsv")
                )
            )
        elif name == "genbank":
            countries = countries_for_scope(scope.get("countries", []), scope.get("regions", []))
            if not countries:
                raise ConfigError(f"scope {settings.scope!r} selects no countries")
            term = build_term(
                country=countries,
                min_length=int(acquire.get("min_length", 1000)),
                max_length=int(acquire.get("max_length", 10000)),
                extra=acquire.get("extra_term"),
            )
            adapters.append(
                GenBankAdapter(
                    term,
                    retries=int(acquire.get("retries", 3)),
                    max_records=max_records,
                    cache_dir=(
                        Path(repo_root) / "data" / "raw" / ".fetch_cache"
                        if acquire.get("cache", True)
                        else None
                    ),
                )
            )
        elif name == "lanl":
            adapters.append(
                LanlAdapter(Path(repo_root) / acquire.get("lanl_path", "data/external/lanl"))
            )
        else:
            raise ConfigError(f"unknown acquire source {name!r}")

    if not adapters:
        raise ConfigError("acquire.sources is empty")
    return adapters


def run_acquire(
    settings: Settings, *, max_records: int | None = None, repo_root: Path | str = "."
) -> Path:
    """Take a snapshot and point the scope at it. Returns the snapshot directory."""
    raw_root = Path(repo_root) / settings.paths["raw"]
    snapshot = take_snapshot(
        build_adapters(settings, max_records=max_records, repo_root=repo_root),
        raw_root,
        config_hash=settings.config_hash(),
    )
    write_pointer(raw_root, settings.scope, snapshot)
    return snapshot
