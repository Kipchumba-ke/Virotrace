"""The tree-based subtype call: nesting inside a reference clade, stems, groups, support."""

from __future__ import annotations

import pytest

from virotrace.curate.subtype import agreement, parent_group
from virotrace.phylo.newick import parse_newick
from virotrace.phylo.subtype_tree import (
    UNRESOLVED,
    concordance,
    support_value,
    tree_subtype_calls,
)

PANEL = {
    "A1a": "A1", "A1b": "A1", "A6a": "A6", "A6b": "A6",
    "Da": "D", "Db": "D", "Ba": "B", "Bb": "B",
}  # fmt: skip
GROUPS = {lab: parent_group(lab) for lab in set(PANEL.values())}
NEWICK = (
    "(((((A1a:1,q_in:1)90:1,A1b:1)95:1,q_stem:1)80:1,(A6a:1,A6b:1)99:1)98:1,"
    "((Da:1,Db:1)100:1,q_sister_D:1)97:1,"
    "(Ba:1,Bb:1,q_poly:1)100:1);"
)


def _calls(newick=NEWICK, **kw):
    return tree_subtype_calls(parse_newick(newick), PANEL, GROUPS, **kw)


def test_a_query_nested_inside_a_reference_clade_takes_its_subtype():
    got = _calls()["q_in"]
    assert (got["call"], got["level"], got["support"]) == ("A1", "subtype", 90.0)


def test_a_query_on_a_subtype_stem_is_not_called_that_subtype():
    got = _calls()
    # q_stem joins the whole A1 clade; at group level A1+A6 are not a pure side either way
    assert got["q_stem"]["call"] in {"A", UNRESOLVED}
    assert got["q_sister_D"]["call"] == UNRESOLVED and got["q_sister_D"]["level"] == "none"


def test_group_level_rescues_paraphyletic_sub_subtypes():
    newick = "((A1a:1,A6a:1)90:1,q:1)95:1,(A1b:1,A6b:1)99:1,(Da:1,Db:1,Ba:1,Bb:1)100:1);"
    newick = "(" + newick
    got = _calls(newick)["q"]
    assert (got["call"], got["level"]) == ("A", "group")


def test_a_multifurcation_gives_no_call_and_weak_support_is_unresolved():
    assert _calls()["q_poly"]["call"] == UNRESOLVED
    weak = _calls(min_support=95.0)["q_in"]
    assert (weak["call"], weak["level"]) == (UNRESOLVED, "weak_support")
    assert _calls(min_support=80.0)["q_in"]["call"] == "A1"


def test_a_missing_reference_is_an_error_and_support_labels_parse():
    with pytest.raises(ValueError, match="not in the tree"):
        tree_subtype_calls(parse_newick("(A:1,B:1,C:1);"), PANEL, GROUPS)
    assert (support_value("95.5/100"), support_value("98"), support_value("")) == (95.5, 98.0, None)


def test_concordance_uses_the_group_aware_agreement_rule():
    rows = [
        {"strain": "1", "subtype_tree": "A", "subtype_called": "A1"},
        {"strain": "2", "subtype_tree": "D", "subtype_called": "A1"},
        {"strain": "3", "subtype_tree": UNRESOLVED, "subtype_called": "C"},
        {"strain": "4", "subtype_tree": "C", "subtype_called": "unassigned"},
    ]
    summary, confusion = concordance(rows, agreement)
    assert summary == {"agree": 1, "disagree": 1, "unresolved": 1, "not_comparable": 1}
    assert confusion["A1"] == {"A": 1, "D": 1}


def test_disagreements_are_counted_by_study_country_and_label_pair():
    from virotrace.curate.subtype import disagreements_by_study

    def row(strain, called, reported, agree, support="0.9"):
        return {
            "strain": strain, "subtype_called": called, "subtype_reported": reported,
            "subtype_agreement": agree, "subtype_support": support,
            "nearest_distance": "0.05", "n_sites": "900",
        }  # fmt: skip

    rows = [row("a", "C", "A1", "disagree"), row("b", "C", "A1", "disagree", "0.7"),
            row("c", "C", "A1", "disagree"), row("d", "A1", "A1", "agree")]  # fmt: skip
    meta = {
        "a": {"study_key": "S1", "country": "KEN"}, "b": {"study_key": "S1", "country": "KEN"},
        "c": {"study_key": "S2", "country": "KEN"}, "d": {"study_key": "S1", "country": "KEN"},
    }  # fmt: skip
    out = disagreements_by_study(rows, meta)
    assert [(r["study_key"], r["n"]) for r in out] == [("S1", 2), ("S2", 1)]
    assert out[0]["median_support"] == 0.8 and out[0]["subtype_called"] == "C"
