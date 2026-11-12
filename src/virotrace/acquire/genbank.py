"""NCBI GenBank adapter (Entrez esearch + efetch), stdlib only.

* The query string is recorded verbatim in the snapshot manifest, because re-running
  it later returns a superset and the query is the reproducibility anchor.
* NCBI requires a contact email for programmatic access. It is read from the
  environment (``VIROTRACE_NCBI_EMAIL``), never from a committed file.
* Transport is injectable so tests run offline and the first live run is a one-line
  swap, not a rewrite.
"""

from __future__ import annotations

import hashlib
import os
import time
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

from ..errors import ConfigError, ContractError
from .base import RawRecord
from .genbank_parse import parse_genbank

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
Transport = Callable[[str, dict[str, str]], str]
ID_CAP = 100000  # most ids one esearch call returns

# Kenya HIV-1 sequences long enough to plausibly cover pol (>= 1000 nt).
DEFAULT_TERM = "txid11676[Organism] AND Kenya[Country] AND 1000:10000[SLEN]"


def build_term(
    *,
    country: str | list[str] = "Kenya",
    min_length: int = 1000,
    max_length: int = 10000,
    taxid: int = 11676,
    extra: str | None = None,
) -> str:
    """Entrez term. A list of countries becomes an OR group (the Africa scope).

    ``extra`` is an extra Entrez clause ANDed on, for narrowing a very large query (for example
    to records that name pol); it is part of the recorded query, so it changes the snapshot id.
    """
    countries = [country] if isinstance(country, str) else list(country)
    if not countries:
        raise ConfigError("build_term needs at least one country")
    if len(countries) == 1:
        place = f"{countries[0]}[Country]"
    else:
        place = "(" + " OR ".join(f"{c}[Country]" for c in countries) + ")"
    term = f"txid{taxid}[Organism] AND {place} AND {min_length}:{max_length}[SLEN]"
    if extra and extra.strip():
        term += f" AND ({extra.strip()})"
    return term


def _evenly_spaced(ids: list[str], n: int) -> list[str]:
    """Deterministic spread of ``n`` ids across the whole ordered list.

    A dev cap taken from the front would be the lowest accessions, i.e. the oldest and
    least-annotated submissions, which is not a sample of anything useful.
    """
    if n <= 0:
        return []
    if n >= len(ids):
        return ids
    return [ids[(i * len(ids)) // n] for i in range(n)]


def default_transport(url: str, params: dict[str, str]) -> str:
    """Plain HTTPS GET. Raises ``ContractError`` on failure after the caller's retries."""
    full = f"{url}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(full, headers={"User-Agent": "virotrace/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


class GenBankAdapter:
    name = "genbank"
    access_tier = 1

    def __init__(
        self,
        term: str = DEFAULT_TERM,
        *,
        email: str | None = None,
        api_key: str | None = None,
        batch_size: int = 200,
        retries: int = 3,
        transport: Transport = default_transport,
        sleep: Callable[[float], None] = time.sleep,
        max_records: int | None = None,
        cache_dir: Path | str | None = None,
    ) -> None:
        self.term = term
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.email = email or os.environ.get("VIROTRACE_NCBI_EMAIL")
        self.api_key = api_key or os.environ.get("VIROTRACE_NCBI_API_KEY")
        self.batch_size = batch_size
        self.retries = retries
        self.transport = transport
        self.sleep = sleep
        self.max_records = max_records
        # NCBI allows 3 requests/s without a key and 10/s with one.
        self.min_interval = 0.11 if self.api_key else 0.34

    def query_description(self) -> str:
        # A capped trial run must not share a snapshot id with the full query,
        # or the immutable-snapshot rule would block the real run afterwards.
        cap = f" max_records={self.max_records}" if self.max_records is not None else ""
        return f"db=nuccore term={self.term}{cap}"

    # ------------------------------------------------------------------ internals
    def _common(self) -> dict[str, str]:
        if not self.email:
            raise ConfigError(
                "NCBI requires a contact email: set VIROTRACE_NCBI_EMAIL (see .env.example)"
            )
        params = {"tool": "virotrace", "email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        return params

    def _call(self, endpoint: str, params: dict[str, str]) -> str:
        last: Exception | None = None
        for attempt in range(self.retries):
            self.sleep(self.min_interval)
            try:
                return self.transport(f"{EUTILS}/{endpoint}", params)
            except Exception as exc:  # network errors are the expected failure here
                last = exc
                self.sleep(2.0**attempt)
        raise ContractError(f"Entrez {endpoint} failed after {self.retries} attempts: {last}")

    def _esearch(self, term: str, retmax: int) -> tuple[list[str], int]:
        import json

        params = {
            **self._common(),
            "db": "nuccore",
            "term": term,
            "retmode": "json",
            "retmax": str(retmax),
            "sort": "accession",
        }
        payload: dict[str, Any] = json.loads(self._call("esearch.fcgi", params))
        result = payload.get("esearchresult", {})
        ids = [str(i) for i in result.get("idlist", [])]
        return ids, int(result.get("count", len(ids)))

    def _search(self) -> tuple[list[str], int]:
        """UIDs for the query, in order, plus NCBI's reported total.

        One esearch returns at most ``ID_CAP`` ids. A larger query (the Africa scope) is cut into
        disjoint publication-date ranges that each fit, and the ranges are joined in date order.
        """
        ids, total = self._esearch(self.term, ID_CAP)
        if total <= ID_CAP:
            return ids, total
        joined = self._search_by_date(1900, 2100)
        if len(joined) != total:
            raise ContractError(
                f"date-split search found {len(joined)} ids but NCBI reports {total} for the "
                "whole query (records may have been added during the run); retry"
            )
        return joined, total

    def _search_by_date(self, first_year: int, last_year: int) -> list[str]:
        term = f"({self.term}) AND {first_year}/01/01:{last_year}/12/31[PDAT]"
        ids, count = self._esearch(term, ID_CAP)
        if count == 0:
            return []
        if count <= ID_CAP:
            return ids
        if first_year == last_year:
            raise ContractError(
                f"{count} records were published in {first_year} alone; narrow the query"
            )
        middle = (first_year + last_year) // 2
        return self._search_by_date(first_year, middle) + self._search_by_date(
            middle + 1, last_year
        )

    # ---------------------------------------------------------------------- public
    def fetch(self) -> Iterator[RawRecord]:
        ids, total = self._search()
        if total > len(ids):
            raise ContractError(
                f"query matched {total} records but only {len(ids)} ids were returned; "
                "narrow the query rather than silently truncating the corpus"
            )
        if self.max_records is not None:
            ids = _evenly_spaced(ids, self.max_records)

        for start in range(0, len(ids), self.batch_size):
            chunk = ids[start : start + self.batch_size]
            yield from parse_genbank(self._fetch_batch(chunk), source=self.name)

    def _fetch_batch(self, chunk: list[str]) -> str:
        """One efetch batch, kept on disk so an interrupted download resumes where it stopped."""
        cached = None
        if self.cache_dir is not None:
            digest = hashlib.sha1(("|".join(chunk)).encode()).hexdigest()[:20]
            cached = self.cache_dir / f"{digest}.gb"
            if cached.exists():
                return cached.read_text()
        text = self._call(
            "efetch.fcgi",
            {
                **self._common(),
                "db": "nuccore",
                "id": ",".join(chunk),
                "rettype": "gb",
                "retmode": "text",
            },
        )
        if cached is not None:
            cached.parent.mkdir(parents=True, exist_ok=True)
            partial = cached.with_suffix(".part")
            partial.write_text(text)
            partial.replace(cached)  # a half-written batch is never mistaken for a finished one
        return text
