import math

import pytest
from helpers import build_cds, random_protein

from doublebind.data.codons import (ALT_START_CODONS, CODON_TABLE, aggregate_missense, cds_to_protein, check_cds,
                                    effect_counts, enumerate_variants, find_numbering_offsets, parse_mutation,
                                    translate, ts_tv_probability)


def test_genetic_code_basics():
    assert len(CODON_TABLE) == 64
    assert sum(aa == "*" for aa in CODON_TABLE.values()) == 3
    assert CODON_TABLE["ATG"] == "M" and CODON_TABLE["TGG"] == "W" and CODON_TABLE["TAA"] == "*"
    assert CODON_TABLE["TCG"] == "S" and CODON_TABLE["TTG"] == "L"


def test_translate_and_alternative_start():
    assert translate("ATGAAATAA") == "MK*"
    assert cds_to_protein("ATGAAATAA") == "MK"
    assert "GTG" in ALT_START_CODONS and cds_to_protein("GTGAAATAA") == "MK"


def test_bad_dna_is_rejected():
    with pytest.raises(ValueError):
        translate("ATGA")
    with pytest.raises(ValueError):
        translate("ATGNNN")


def test_check_cds_exact_and_mismatch():
    protein = random_protein(150, seed=3)
    cds = build_cds(protein, seed=3)
    assert check_cds(cds, protein).ok
    changed = protein[:40] + ("A" if protein[40] != "A" else "C") + protein[41:]
    chk = check_cds(cds, changed)
    assert not chk.ok and len(chk.differences) == 1 and chk.differences[0][0] == 41
    assert chk.identity == pytest.approx(149 / 150)


def test_ts_tv_probabilities_sum_to_one():
    for kappa in (0.5, 1.0, 2.0, 8.0):
        for ref in "ACGT":
            total = sum(ts_tv_probability(ref, alt, kappa) for alt in "ACGT" if alt != ref)
            assert total == pytest.approx(1.0)
    assert ts_tv_probability("C", "T", 2.0) == pytest.approx(0.5)      # transition
    assert ts_tv_probability("C", "A", 2.0) == pytest.approx(0.25)     # transversion


def test_nine_changes_per_codon_and_phe_example():
    variants = enumerate_variants("ATG" + "TTT" + "AAA" + "TAA")
    phe = [v for v in variants if v.pos == 2]
    assert len(phe) == 9
    counts = effect_counts(phe)
    assert counts["synonymous"] == 1 and counts["missense"] == 8 and counts["nonsense"] == 0
    aa = [a for a in aggregate_missense(variants) if a.pos == 2]
    assert {a.mut_aa for a in aa} == set("LSYCIV")
    assert next(a for a in aa if a.mut_aa == "L").n_routes == 3


def test_lys_codon_has_one_stop_neighbour():
    variants = [v for v in enumerate_variants("ATG" + "AAA" + "TAA") if v.pos == 2]
    assert effect_counts(variants)["nonsense"] == 1


def test_start_codon_and_terminal_stop_are_skipped_by_default():
    variants = enumerate_variants("ATG" + "TTT" + "TAA")
    assert {v.pos for v in variants} == {2}
    with_start = enumerate_variants("ATG" + "TTT" + "TAA", include_start=True)
    assert {v.pos for v in with_start} == {1, 2}


def test_weights_follow_kappa():
    v = next(x for x in enumerate_variants("ATG" + "TTT" + "TAA", kappa=4.0)
             if x.pos == 2 and x.codon_site == 2 and x.alt_nt == "C")  # T->C transition
    assert v.transition and v.weight == pytest.approx(4 / 6)


def test_standard_code_global_counts():
    """Over all 61 sense codons: 549 single-nucleotide changes split 134 / 392 / 23."""
    sense = [c for c, aa in CODON_TABLE.items() if aa != "*"]
    cds = "ATG" + "".join(sense) + "TAA"
    counts = effect_counts(enumerate_variants(cds))
    assert sum(counts.values()) == 61 * 9
    assert (counts["synonymous"], counts["missense"], counts["nonsense"]) == (134, 392, 23)


def test_aggregate_weights_sum_to_missense_weight():
    protein = random_protein(80, seed=5)
    variants = enumerate_variants(build_cds(protein, seed=5))
    total_missense = sum(v.weight for v in variants if v.effect == "missense")
    assert sum(a.weight for a in aggregate_missense(variants)) == pytest.approx(total_missense)


def test_parse_mutation():
    assert parse_mutation("S83L") == ("S", 83, "L")
    with pytest.raises(ValueError):
        parse_mutation("83L")


def test_find_numbering_offsets():
    protein = "MKTAYIAKQRQISFVKSHFSRQ"  # T3, Y5, K8 in sequence numbering
    literature = ["T5A", "Y7F", "K10R"]  # literature numbering is shifted by +2
    assert find_numbering_offsets(protein, literature)[0] == -2
    assert find_numbering_offsets(protein, ["T3A", "Y5F"])[0] == 0
    assert find_numbering_offsets(protein, ["T3A", "W5F"]) == []
