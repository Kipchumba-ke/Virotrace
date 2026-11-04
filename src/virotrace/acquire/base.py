"""Source adapter contract.

Adapters talk to the outside world and nothing else: they yield raw records and
describe their query. They never interpret or filter biology, which is why a
controlled-access adapter could be added later without touching anything downstream.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Protocol


@dataclass(frozen=True, slots=True)
class RawRecord:
    """One record as the source supplied it. Nothing here is normalised."""

    accession: str
    sequence: str
    source: str
    fields: dict[str, str] = field(default_factory=dict)


class SourceAdapter(Protocol):
    name: str
    access_tier: int

    def query_description(self) -> str:
        """The exact query, recorded verbatim in the snapshot manifest."""
        ...

    def fetch(self) -> Iterator[RawRecord]: ...
