import math

import pytest

from virotrace.curate.neighbours import nearest_identity


def test_each_tip_gets_its_best_match_and_the_odd_one_out_scores_low():
    base = "ACGT" * 150
    near = base[:-1] + "A"  # 1 difference in 600
    far = "".join({"A": "C", "C": "G", "G": "T", "T": "A"}[c] for c in base)  # every site differs
    ident = nearest_identity([base, near, far], min_shared=500)
    assert ident[0] == pytest.approx(599 / 600)
    assert ident[1] == pytest.approx(599 / 600)
    assert ident[2] < 0.01


def test_a_tip_with_no_partner_reaching_min_shared_is_not_assessable_not_zero():
    a = "ACGT" * 150
    b = "-" * 400 + "ACGT" * 50  # only 200 called sites
    ident = nearest_identity([a, b], min_shared=500)
    assert all(math.isnan(x) for x in ident)


def test_ambiguity_codes_and_gaps_do_not_count_as_shared_sites():
    a = "ACGT" * 150
    b = "N" * 100 + a[100:]
    assert nearest_identity([a, b], min_shared=500)[0] == pytest.approx(1.0)


def test_single_sequence_and_ragged_input():
    assert math.isnan(nearest_identity(["ACGT"], min_shared=1)[0])
    with pytest.raises(ValueError):
        nearest_identity(["ACGT", "ACG"], min_shared=1)
