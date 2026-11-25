import random

import pytest
from helpers import build_cds, random_protein

from doublebind.baselines import average_ranks, find_esm_file, format_report, run_baselines
from doublebind.data.codons import aggregate_missense, enumerate_variants
from doublebind.tables import COLUMNS, build_variant_table, write_table


def test_average_ranks_with_ties():
    assert average_ranks({"a": 3, "b": 3, "c": 1}) == {"a": 1.5, "b": 1.5, "c": 3.0}
    assert average_ranks({"a": 3, "b": 1}, higher_is_better=False) == {"b": 1.0, "a": 2.0}


def test_find_esm_file_prefers_the_biggest_model(tmp_path):
    for name in ("gyra_esm_esm2_t12_35M_UR50D.csv", "gyra_esm_esm2_t33_650M_UR50D.csv",
                 "gyra_esm_esm2_t33_650M_UR50D_positions.csv", "tem1_esm_esm2_t6_8M_UR50D.csv"):
        (tmp_path / name).write_text("x")
    assert find_esm_file(tmp_path, "gyra").name == "gyra_esm_esm2_t33_650M_UR50D.csv"
    assert find_esm_file(tmp_path, "tem1").name == "tem1_esm_esm2_t6_8M_UR50D.csv"
    assert find_esm_file(tmp_path, "none") is None
    assert find_esm_file(tmp_path, "gyra", explicit="custom.csv").name == "custom.csv"


def _synthetic():
    protein = random_protein(120, seed=11, fixed={83: "S", 87: "D"})
    cds = build_cds(protein, overrides={83: "TCG", 87: "GAC"}, seed=11)
    aa = aggregate_missense(enumerate_variants(cds))
    rng = random.Random(0)
    esm = [{"pos": v.pos, "mut_aa": v.mut_aa, "llr": rng.uniform(-6, 1)} for v in aa]
    for row in esm:                                    # the two known changes look very plausible
        if (row["pos"], row["mut_aa"]) in {(83, "L"), (87, "N")}:
            row["llr"] = 3.0
    pocket = [{"pos": p, "min_dist": 3.0 + abs(p - 85) * 1.0 if abs(p - 85) <= 10 else 30.0} for p in range(1, 121)]
    known = [{"mutation": "S83L"}, {"mutation": "D87N"}]
    return protein, aa, esm, pocket, known


def test_run_baselines_ranks_known_positions_sensibly():
    protein, aa, esm, pocket, known = _synthetic()
    res = run_baselines(protein, aa, known, esm, pocket, cutoff=8.0)
    assert res["n_positions"] == 119 and res["numbering_offset"] == 0
    assert all(k["reachable"] for k in res["known"])
    assert res["random"]["expected_rank"] == pytest.approx(60.0)
    prox = res["position_level"]["proximity"]["known_ranks"]
    assert prox[83] <= 5 and prox[87] <= 5                        # both sit within 2-4 residues of the pocket centre
    inside = res["position_level"]["proximity_then_esm"]
    assert inside["n_ranked"] == 11 and inside["known_ranks"][83] <= 2 and inside["known_ranks"][87] <= 2
    v = res["variant_level"]["known"]["S83L"]
    assert v["llr"] == 3.0 and v["rank_all"] <= 2 and v["rank_in_pocket"] <= 2
    text = format_report("gyra", res)
    assert "pos 83" in text and "proximity only" in text and "random (expected rank)" in text


def test_run_baselines_without_optional_inputs():
    protein, aa, _esm, _pocket, known = _synthetic()
    res = run_baselines(protein, aa, known)
    assert res["position_level"] == {} and res["variant_level"] == {}
    assert "random" in format_report("gyra", res)


def test_run_baselines_rejects_labels_that_do_not_match_the_sequence():
    protein, aa, *_ = _synthetic()
    with pytest.raises(ValueError):
        run_baselines(protein, aa, [{"mutation": "W83L"}, {"mutation": "W87N"}])


def test_variant_table_merges_everything(tmp_path):
    _protein, aa, esm, pocket, _known = _synthetic()
    rows = build_variant_table(aa, esm[:50], pocket, cutoff=8.0)
    assert len(rows) == len(aa) and list(rows[0]) == COLUMNS
    assert sum(r["esm_llr"] != "" for r in rows) == 50
    near = next(r for r in rows if r["pos"] == 85)
    far = next(r for r in rows if r["pos"] == 2)
    assert near["in_pocket"] == 1 and far["in_pocket"] == 0 and near["dist_min"] == 3.0
    empty = build_variant_table(aa)
    assert empty[0]["esm_llr"] == "" and empty[0]["dist_min"] == "" and empty[0]["in_pocket"] == ""
    out = tmp_path / "t.csv"
    write_table(out, rows)
    assert out.read_text().splitlines()[0] == ",".join(COLUMNS)
