"""Minimal GenBank flat-file parser (stdlib only).

Extracts what Virotrace needs from each record: accession, the nucleotide sequence,
and the source-feature qualifiers (/country, /collection_date, /note, /isolate,
/strain). Deliberately small: it does not interpret biology, and anything it cannot
read is left empty for normalisation to quarantine or mark unknown.
"""

from __future__ import annotations

import re
from collections.abc import Iterator

from .base import RawRecord

_QUALIFIER = re.compile(r"^\s{21}/(?P<key>[A-Za-z_]+)(?:=(?P<value>.*))?$")
_WANTED = {
    "country",
    "geo_loc_name",
    "collection_date",
    "note",
    "isolate",
    "strain",
    "organism",
    "host",
    "clone",
    "proviral",
    "isolation_source",
}


def _clean(value: str) -> str:
    return " ".join(value.strip().strip('"').split())


def _parse_source_qualifiers(lines: list[str]) -> dict[str, str]:
    """Qualifiers of the first ``source`` feature, with wrapped values rejoined."""
    out: dict[str, str] = {}
    in_features = False
    in_source = False
    key: str | None = None
    buf: list[str] = []

    def flush() -> None:
        nonlocal key, buf
        if key in _WANTED and key not in out:
            out[key] = _clean(" ".join(buf))
        key, buf = None, []

    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN"):
            break
        if re.match(r"^\s{5}\S", line):  # a feature key line, e.g. "     source   1..1300"
            flush()
            if in_source:
                break  # only the first source feature
            in_source = line.split()[0] == "source"
            continue
        if not in_source:
            continue
        match = _QUALIFIER.match(line)
        if match:
            flush()
            key = match["key"]
            buf = [match["value"] or ""]
        elif key is not None and line.startswith(" " * 21):
            buf.append(line.strip())
    flush()
    return out


def _parse_definition(lines: list[str]) -> str:
    """The DEFINITION line with its wrapped continuation lines rejoined."""
    parts: list[str] = []
    for line in lines:
        if line.startswith("DEFINITION"):
            parts.append(line[len("DEFINITION") :].strip())
        elif parts:
            if line.startswith(" " * 12):
                parts.append(line.strip())
            else:
                break
    return " ".join(parts)


_REF_FIELD = re.compile(r"^\s{2,3}(AUTHORS|TITLE|JOURNAL|PUBMED)\s+(.*)$")


def _parse_first_reference(lines: list[str]) -> dict[str, str]:
    """AUTHORS, TITLE, JOURNAL and PUBMED of REFERENCE 1, wrapped lines rejoined."""
    out: dict[str, list[str]] = {}
    in_ref = False
    key: str | None = None
    for line in lines:
        if line.startswith("REFERENCE"):
            if in_ref:
                break
            in_ref = True
            continue
        if not in_ref:
            continue
        if line.startswith(("FEATURES", "ORIGIN", "COMMENT")):
            break
        match = _REF_FIELD.match(line)
        if match:
            key = match[1]
            out[key] = [match[2].strip()]
        elif key is not None and line.startswith(" " * 12):
            out[key].append(line.strip())
    return {k: " ".join(v) for k, v in out.items()}


def _parse_bioproject(lines: list[str]) -> str:
    for line in lines:
        if line.startswith("DBLINK") or line.startswith(" " * 12):
            match = re.search(r"BioProject:\s*(PRJ[A-Z]+\d+)", line)
            if match:
                return match[1]
    return ""


def _parse_sequence(lines: list[str]) -> str:
    seq: list[str] = []
    in_origin = False
    for line in lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            seq.append("".join(ch for ch in line if ch.isalpha()))
    return "".join(seq).upper()


def parse_genbank(text: str, *, source: str = "genbank") -> Iterator[RawRecord]:
    """Yield one :class:`RawRecord` per ``//``-terminated GenBank record."""
    block: list[str] = []
    for line in text.splitlines():
        if line.startswith("//"):
            if block:
                record = _record_from_block(block, source)
                if record is not None:
                    yield record
            block = []
        elif line.strip() or block:
            block.append(line)
    if block and any(line.startswith("LOCUS") for line in block):
        record = _record_from_block(block, source)
        if record is not None:
            yield record


def _record_from_block(lines: list[str], source: str) -> RawRecord | None:
    accession = ""
    version = ""
    for line in lines:
        if line.startswith("ACCESSION") and not accession:
            parts = line.split()
            accession = parts[1] if len(parts) > 1 else ""
        elif line.startswith("VERSION") and not version:
            parts = line.split()
            version = parts[1] if len(parts) > 1 else ""
    # Prefer the versioned id: a changed sequence under one accession must not
    # silently masquerade as the same record in a later snapshot.
    identifier = version or accession
    if not identifier:
        return None

    sequence = _parse_sequence(lines)
    if not sequence:
        return None

    qualifiers = _parse_source_qualifiers(lines)
    reference = _parse_first_reference(lines)
    fields = {
        # INSDC renamed /country to /geo_loc_name; newer records only carry the latter.
        "country": qualifiers.get("country", "") or qualifiers.get("geo_loc_name", ""),
        "collection_date": qualifiers.get("collection_date", ""),
        "subtype": _subtype_from_note(qualifiers),
        "isolate": qualifiers.get("isolate", "") or qualifiers.get("strain", ""),
        "host": qualifiers.get("host", ""),
        "definition": _parse_definition(lines),
        # Recorded as stated. Whether a clone is a lab strain or a within-host PCR clone
        # is not decidable here, so the policy lives in config, not in the parser.
        "clone": qualifiers.get("clone", ""),
        "proviral": "yes" if "proviral" in qualifiers else "",
        # Which study a record belongs to, and what it says about the people sampled.
        "note": qualifiers.get("note", ""),
        "isolation_source": qualifiers.get("isolation_source", ""),
        "study_authors": reference.get("AUTHORS", ""),
        "study_title": reference.get("TITLE", ""),
        "study_journal": reference.get("JOURNAL", ""),
        "pubmed": reference.get("PUBMED", ""),
        "bioproject": _parse_bioproject(lines),
    }
    return RawRecord(accession=identifier, sequence=sequence, source=source, fields=fields)


_SUBTYPE = re.compile(r"\bsubtype[:\s]+([A-Za-z0-9_]+)", re.IGNORECASE)


def _subtype_from_note(qualifiers: dict[str, str]) -> str:
    """Subtype as *reported* by the depositor; Virotrace makes its own call later."""
    note = qualifiers.get("note", "")
    match = _SUBTYPE.search(note)
    return match.group(1) if match else ""
