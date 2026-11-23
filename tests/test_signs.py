import pytest

from doublebind.signs import (SOURCES, Convention, check_directions, register_source, to_doublebind)


def test_doublebind_convention_is_unchanged():
    assert to_doublebind(1.5, "doublebind") == 1.5


def test_stabilizing_positive_sources_are_flipped():
    assert to_doublebind(2.0, "megascale") == -2.0
    assert to_doublebind(-0.7, "megascale") == 0.7


def test_destabilizing_positive_sources_are_kept():
    assert to_doublebind(1.0, "thermompnn") == 1.0


def test_unknown_source_is_refused():
    with pytest.raises(KeyError):
        to_doublebind(1.0, "some-new-tool")


def test_register_new_source():
    register_source("mytool", Convention.STABILIZING_POSITIVE, "checked docs 2026-11-21", confirmed=True)
    assert to_doublebind(3.0, "mytool") == -3.0
    del SOURCES["mytool"]


def test_unconfirmed_sources_are_flagged_for_checking():
    assert SOURCES["doublebind"].confirmed
    assert not SOURCES["megascale"].confirmed  # until you read the dataset README


def test_direction_check_counts_agreement():
    expected = {"G238S": "destabilizing", "M182T": "stabilizing", "A42G": "destabilizing"}
    preds = {"G238S": 1.2, "M182T": -0.5}  # A42G missing
    report = check_directions(preds, expected)
    assert report.n_agree == 2 and report.n_total == 3
    assert [r[3] for r in report.rows] == ["agree", "agree", "missing"]


def test_direction_check_flags_a_flipped_model():
    expected = {"G238S": "destabilizing", "M182T": "stabilizing"}
    flipped = {"G238S": -1.2, "M182T": 0.5}
    assert check_directions(flipped, expected).n_agree == 0


def test_direction_check_rejects_bad_labels():
    with pytest.raises(ValueError):
        check_directions({"X1Y": 1.0}, {"X1Y": "bad"})
