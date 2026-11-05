"""The clone policy: configuration decides, and no record is ever lost."""

from __future__ import annotations

import pytest

from virotrace.curate.clones import is_excluded, split_by_policy
from virotrace.errors import ConfigError

ROWS = [
    {"strain": "A", "clone_derived": True},
    {"strain": "B", "clone_derived": False},
    {"strain": "C", "clone_derived": True},
    {"strain": "D"},  # a row from before the flag existed
]


def test_exclude_policy_removes_only_clone_derived_rows():
    kept, excluded = split_by_policy(ROWS, "exclude")
    assert [r["strain"] for r in kept] == ["B", "D"]
    assert [r["strain"] for r in excluded] == ["A", "C"]


def test_keep_all_policy_excludes_nothing():
    kept, excluded = split_by_policy(ROWS, "keep_all")
    assert len(kept) == 4 and excluded == []


@pytest.mark.parametrize("policy", ["exclude", "keep_all"])
def test_every_row_lands_in_exactly_one_side(policy):
    kept, excluded = split_by_policy(ROWS, policy)
    assert len(kept) + len(excluded) == len(ROWS)


def test_unknown_policy_is_a_config_error_not_a_silent_default():
    with pytest.raises(ConfigError, match="clone_policy"):
        is_excluded({"clone_derived": True}, "one_per_isolate")


# --- individual policy ---------------------------------------------------------------

from virotrace.curate.clones import (  # noqa: E402
    apply_population_policies,
    individual_id,
    one_per_individual,
)

RV = [r"^(RV217_\d+)_"]


def _r(strain, isolate, date="2010-05", amb=0.0, clone=False):
    return {"strain": strain, "isolate": isolate, "date": date, "ambiguity_frac": amb,
            "clone_derived": clone}  # fmt: skip


def test_individual_id_uses_the_first_matching_pattern():
    assert individual_id("RV217_20225_v01_GP_1_90", RV) == "RV217_20225"
    assert individual_id("Ke-123", RV) is None
    assert individual_id(None, RV) is None


def test_one_per_individual_keeps_earliest_then_least_ambiguous_then_smallest_id():
    rows = [
        _r("c", "RV217_1_v02_GP_1", date="2010-06"),
        _r("b", "RV217_1_v01_GP_2", date="2010-05", amb=0.02),
        _r("a", "RV217_1_v01_GP_3", date="2010-05", amb=0.01),
        _r("z", "RV217_1_v01_GP_4", date="2010-05", amb=0.01),
        _r("solo", "unrelated"),
        _r("other", "RV217_2_v01_GP_1"),
    ]
    kept, dropped = one_per_individual(rows, RV)
    assert {r["strain"] for r in kept} == {"a", "solo", "other"}
    assert len(kept) + len(dropped) == len(rows)
    again, _ = one_per_individual(list(reversed(rows)), RV)
    assert {r["strain"] for r in again} == {"a", "solo", "other"}  # order independent


def test_unidentifiable_records_are_never_collapsed():
    rows = [_r("x1", ""), _r("x2", ""), _r("x3", None)]
    kept, dropped = one_per_individual(rows, RV)
    assert len(kept) == 3 and not dropped


def test_policies_apply_clone_then_individual_and_report_counts():
    rows = [
        _r("a", "RV217_1_v01", date="2010-05"),
        _r("b", "RV217_1_v02", date="2010-06"),
        _r("c", "plain", clone=True),
        _r("d", "plain2"),
    ]
    kept, removed = apply_population_policies(
        rows, clone_policy="exclude", individual_policy="one_per_individual", patterns=RV
    )
    assert {r["strain"] for r in kept} == {"a", "d"}
    assert removed == {"clone_policy": 1, "individual_policy": 1}
    kept_all, removed_all = apply_population_policies(
        rows, clone_policy="keep_all", individual_policy="keep_all", patterns=RV
    )
    assert len(kept_all) == 4 and removed_all == {"clone_policy": 0, "individual_policy": 0}
    with pytest.raises(ConfigError):
        apply_population_policies(rows, clone_policy="exclude", individual_policy="x", patterns=RV)
