"""Country normalisation and region assignment (UN M49 sub-regions for Africa).

GenBank ``/country`` looks like ``Kenya: Nairobi``. Only the country is trusted to
resolve; the part after the colon is kept as a *raw* admin1 hint and is never
promoted to a resolved county. The map legend must state the country-only fraction,
so an unmapped or vague location stays unresolved rather than being guessed.
"""

from __future__ import annotations

from dataclasses import dataclass

# alpha3 -> (canonical name, region)
_AFRICA: dict[str, tuple[str, str]] = {
    # eastern
    "BDI": ("Burundi", "eastern_africa"),
    "COM": ("Comoros", "eastern_africa"),
    "DJI": ("Djibouti", "eastern_africa"),
    "ERI": ("Eritrea", "eastern_africa"),
    "ETH": ("Ethiopia", "eastern_africa"),
    "KEN": ("Kenya", "eastern_africa"),
    "MDG": ("Madagascar", "eastern_africa"),
    "MWI": ("Malawi", "eastern_africa"),
    "MUS": ("Mauritius", "eastern_africa"),
    "MOZ": ("Mozambique", "eastern_africa"),
    "RWA": ("Rwanda", "eastern_africa"),
    "SYC": ("Seychelles", "eastern_africa"),
    "SOM": ("Somalia", "eastern_africa"),
    "SSD": ("South Sudan", "eastern_africa"),
    "TZA": ("Tanzania", "eastern_africa"),
    "UGA": ("Uganda", "eastern_africa"),
    "ZMB": ("Zambia", "eastern_africa"),
    "ZWE": ("Zimbabwe", "eastern_africa"),
    # southern
    "BWA": ("Botswana", "southern_africa"),
    "SWZ": ("Eswatini", "southern_africa"),
    "LSO": ("Lesotho", "southern_africa"),
    "NAM": ("Namibia", "southern_africa"),
    "ZAF": ("South Africa", "southern_africa"),
    # western
    "BEN": ("Benin", "western_africa"),
    "BFA": ("Burkina Faso", "western_africa"),
    "CPV": ("Cabo Verde", "western_africa"),
    "CIV": ("Cote d'Ivoire", "western_africa"),
    "GMB": ("Gambia", "western_africa"),
    "GHA": ("Ghana", "western_africa"),
    "GIN": ("Guinea", "western_africa"),
    "GNB": ("Guinea-Bissau", "western_africa"),
    "LBR": ("Liberia", "western_africa"),
    "MLI": ("Mali", "western_africa"),
    "MRT": ("Mauritania", "western_africa"),
    "NER": ("Niger", "western_africa"),
    "NGA": ("Nigeria", "western_africa"),
    "SEN": ("Senegal", "western_africa"),
    "SLE": ("Sierra Leone", "western_africa"),
    "TGO": ("Togo", "western_africa"),
    # middle
    "AGO": ("Angola", "middle_africa"),
    "CMR": ("Cameroon", "middle_africa"),
    "CAF": ("Central African Republic", "middle_africa"),
    "TCD": ("Chad", "middle_africa"),
    "COG": ("Republic of the Congo", "middle_africa"),
    "COD": ("Democratic Republic of the Congo", "middle_africa"),
    "GNQ": ("Equatorial Guinea", "middle_africa"),
    "GAB": ("Gabon", "middle_africa"),
    "STP": ("Sao Tome and Principe", "middle_africa"),
    # northern
    "DZA": ("Algeria", "northern_africa"),
    "EGY": ("Egypt", "northern_africa"),
    "LBY": ("Libya", "northern_africa"),
    "MAR": ("Morocco", "northern_africa"),
    "SDN": ("Sudan", "northern_africa"),
    "TUN": ("Tunisia", "northern_africa"),
}

# A few common non-African sources, so they classify as non_africa rather than unmapped.
_NON_AFRICA: dict[str, str] = {
    "USA": "United States",
    "GBR": "United Kingdom",
    "FRA": "France",
    "BRA": "Brazil",
    "IND": "India",
    "CHN": "China",
    "THA": "Thailand",
    "ESP": "Spain",
    "DEU": "Germany",
    "BEL": "Belgium",
    "NLD": "Netherlands",
    "AUS": "Australia",
    "CAN": "Canada",
}

_ALIASES: dict[str, str] = {
    "ivory coast": "CIV",
    "cote d'ivoire": "CIV",
    "côte d'ivoire": "CIV",
    "united republic of tanzania": "TZA",
    "swaziland": "SWZ",
    "the gambia": "GMB",
    "cape verde": "CPV",
    "congo": "COG",
    "drc": "COD",
    "dr congo": "COD",
    "democratic republic of congo": "COD",
    "republic of congo": "COG",
    "usa": "USA",
    "united states of america": "USA",
    "uk": "GBR",
    "great britain": "GBR",
    "south sudan": "SSD",
}


def _build_index() -> dict[str, str]:
    index: dict[str, str] = {}
    for code, (name, _) in _AFRICA.items():
        index[name.lower()] = code
    for code, name in _NON_AFRICA.items():
        index[name.lower()] = code
    index.update(_ALIASES)
    return index


_INDEX = _build_index()


@dataclass(frozen=True, slots=True)
class GeoResult:
    country: str
    country_code: str
    region: str
    admin1_raw: str  # whatever followed the colon, unresolved
    admin1_resolved: bool


def normalise_country(raw: str | None) -> GeoResult | None:
    """Map a source country string, or return ``None`` if it is not recognised.

    The caller quarantines ``None`` with a reason and logs the unmapped value; an
    unrecognised string is never silently assigned a country.
    """
    text = (raw or "").strip()
    if not text:
        return None
    country_part, _, rest = text.partition(":")
    code = _INDEX.get(country_part.strip().lower())
    if code is None:
        return None
    if code in _AFRICA:
        name, region = _AFRICA[code]
    else:
        name, region = _NON_AFRICA[code], "non_africa"
    return GeoResult(
        country=name,
        country_code=code,
        region=region,
        admin1_raw=rest.strip(),
        admin1_resolved=False,
    )
