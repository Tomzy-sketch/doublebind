import math

import pytest

from doublebind.csi.categories import COLLATERAL, CROSS, INCONCLUSIVE, SPARING, label_pair
from doublebind.csi.csi import (csi_from_ddg_bind, ddg_bind_for_csi, ddg_from_log10_ic50_delta,
                                ic50_fold_change, mean_needed_for_call, prob_csi_at_least, rt)


def test_rt_at_body_temperature():
    assert rt(310.15) == pytest.approx(0.616, abs=0.001)


def test_rt_at_room_temperature():
    assert rt(298.15) == pytest.approx(0.593, abs=0.001)


def test_four_fold_susceptibility_needs_minus_0_85():
    # PRD Section 1 item 5 and Section 19
    assert ddg_bind_for_csi(2.0) == pytest.approx(-0.854, abs=0.002)


def test_one_kcal_is_about_five_fold():
    assert ic50_fold_change(1.0) == pytest.approx(5.07, abs=0.05)


def test_csi_round_trip():
    for csi in (-2.0, -0.5, 0.0, 1.3, 2.0, 3.7):
        assert csi_from_ddg_bind(ddg_bind_for_csi(csi)) == pytest.approx(csi)


def test_negative_ddg_means_more_susceptible():
    assert csi_from_ddg_bind(-1.0) > 0 > csi_from_ddg_bind(1.0)


def test_worked_example_sigma_1_2():
    # PRD Section 19: sigma = 1.2 needs a mean near -1.9 kcal/mol, a 20-fold predicted gain
    mu = mean_needed_for_call(1.2, csi_min=2.0, confidence=0.8)
    assert mu == pytest.approx(-1.86, abs=0.02)
    assert 1 / ic50_fold_change(mu) == pytest.approx(20, abs=1.5)
    assert prob_csi_at_least(mu, 1.2) == pytest.approx(0.8, abs=1e-6)


def test_probability_limits():
    assert prob_csi_at_least(-5.0, 0.0) == 1.0
    assert prob_csi_at_least(0.0, 0.0) == 0.0
    assert prob_csi_at_least(-0.854, 1.0) == pytest.approx(0.5, abs=0.01)


def test_log10_ic50_delta_conversion():
    # a 10-fold rise in IC50 is RT*ln(10) kcal/mol
    assert ddg_from_log10_ic50_delta(1.0) == pytest.approx(rt() * math.log(10))
    assert ic50_fold_change(ddg_from_log10_ic50_delta(1.0)) == pytest.approx(10.0)


@pytest.mark.parametrize("csi,p,expected", [
    (2.5, 0.9, COLLATERAL),
    (2.5, 0.5, INCONCLUSIVE),   # large but not confident
    (1.5, 0.99, INCONCLUSIVE),  # weak sensitization
    (0.4, 0.1, SPARING),
    (-0.4, 0.0, SPARING),
    (-1.0, 0.0, CROSS),
    (-3.0, 0.0, CROSS),
])
def test_label_pair(csi, p, expected):
    assert label_pair(csi, p) == expected
