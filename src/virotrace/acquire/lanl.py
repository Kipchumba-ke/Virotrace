"""LANL HIV Sequence Database adapter, reading search downloads from disk.

LANL has no query API, so the person runs the search in the web interface, downloads the
result as FASTA with the default header (``Subtype.Country.Year.Name.Accession``), and puts
the files in ``acquire.lanl_path``. The adapter reads them and nothing else; it never reaches
the network. The snapshot manifest records the file names and their SHA256, which is the
reproducibility anchor for a source that returns a superset when re-queried.

What it does to the data is format translation only: two-letter country codes become country
names, a two-digit header year becomes a year, and alignment gaps are removed. A LANL header
carries no title, authors or PubMed id, so records seen only in LANL have no study key;
a record that is also in GenBank is resolved to the GenBank copy by deduplication.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterator
from pathlib import Path

from ..errors import ContractError
from .base import RawRecord

SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}

ISO2_TO_ISO3 = {
    "AO": "AGO", "BI": "BDI", "BJ": "BEN", "BF": "BFA", "BW": "BWA", "CF": "CAF",
    "CI": "CIV", "CM": "CMR", "CD": "COD", "CG": "COG", "KM": "COM", "CV": "CPV",
    "DJ": "DJI", "DZ": "DZA", "EG": "EGY", "ER": "ERI", "ET": "ETH", "GA": "GAB",
    "GH": "GHA", "GN": "GIN", "GM": "GMB", "GW": "GNB", "GQ": "GNQ", "KE": "KEN",
    "LR": "LBR", "LY": "LBY", "LS": "LSO", "MA": "MAR", "MG": "MDG", "ML": "MLI",
    "MZ": "MOZ", "MR": "MRT", "MU": "MUS", "MW": "MWI", "NA": "NAM", "NE": "NER",
    "NG": "NGA", "RW": "RWA", "SD": "SDN", "SN": "SEN", "SL": "SLE", "SO": "SOM",
    "SS": "SSD", "ST": "STP", "SZ": "SWZ", "SC": "SYC", "TD": "TCD", "TG": "TGO",
    "TN": "TUN", "TZ": "TZA", "UG": "UGA", "ZA": "ZAF", "ZM": "ZMB", "ZW": "ZWE",
}  # fmt: skip
_ACCESSION = re.compile(r"^(?:[A-Z]{1,2}\d{5,6}|[A-Z]{4}\d{8,10})(?:\.\d+)?$")
_YEAR_PIVOT = 40  # 00-40 -> 2000s, 83-99 -> 1900s, 41-82 -> unknown (no HIV sequence is that old)
_UNSAFE = re.compile(r"[^A-Za-z0-9_.\-]")


def parse_year(text: str) -> str:
    text = text.strip()
    if re.fullmatch(r"(19|20)\d{2}", text):
        return text
    if re.fullmatch(r"\d{2}", text):
        yy = int(text)
        if yy <= _YEAR_PIVOT:
            return str(2000 + yy)
        if yy >= 83:
            return str(1900 + yy)
    return ""


def parse_header(header: str) -> dict[str, str]:
    """``A1.KE.00.name.ACC`` -> subtype, country name, year, isolate name, accession."""
    parts = header.strip().split(".")
    if len(parts) < 5:
        raise ContractError(
            f"LANL header {header!r} is not Subtype.Country.Year.Name.Accession; download the "
            "FASTA with the default header format"
        )
    subtype, code, year, accession = parts[0], parts[1], parts[2], parts[-1]
    name = ".".join(parts[3:-1])
    iso3 = ISO2_TO_ISO3.get(code.upper())
    from ..curate.geo import _AFRICA  # country names are defined once, in the geo table

    country = _AFRICA[iso3][0] if iso3 else code  # unmapped codes reach normalise and are logged
    has_accession = bool(_ACCESSION.match(accession))
    return {
        "subtype": subtype,
        "country": country,
        "collection_date": parse_year(year),
        "isolate": name,
        "accession": accession if has_accession else "",
    }


class LanlAdapter:
    name = "lanl"
    access_tier = 1

    def __init__(self, directory: Path | str) -> None:
        self.directory = Path(directory)

    def _files(self) -> list[Path]:
        if not self.directory.is_dir():
            raise ContractError(
                f"LANL directory {self.directory} not found; run the search at hiv.lanl.gov, "
                "download the FASTA, and put it there (docs/03-data-acquisition.md)"
            )
        files = sorted(p for p in self.directory.iterdir() if p.suffix.lower() in SUFFIXES)
        if not files:
            raise ContractError(f"no FASTA files in {self.directory}")
        return files

    def query_description(self) -> str:
        parts = [
            f"{p.name}:{hashlib.sha256(p.read_bytes()).hexdigest()[:16]}" for p in self._files()
        ]
        return "lanl-download:" + ",".join(parts)

    def fetch(self) -> Iterator[RawRecord]:
        for path in self._files():
            header, chunks = None, []
            for line in [*path.read_text().splitlines(), ">"]:
                if line.startswith(">"):
                    if header is not None:
                        yield self._record(header, "".join(chunks))
                    header, chunks = line[1:].strip(), []
                else:
                    chunks.append(line.strip())

    def _record(self, header: str, raw: str) -> RawRecord:
        fields = parse_header(header)
        accession = fields.pop("accession") or "LANL_" + _UNSAFE.sub("_", fields["isolate"])
        sequence = re.sub(r"[-.\s]", "", raw).upper()
        return RawRecord(accession=accession, sequence=sequence, source=self.name, fields=fields)
