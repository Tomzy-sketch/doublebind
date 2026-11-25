import numpy as np
import pytest

from doublebind.scoring.esm_llr import AA20, MAX_RESIDUES, llr_table, position_logits

V, CLS, EOS, MASK = 33, 0, 2, 32
AA_IDS = {aa: 4 + i for i, aa in enumerate(AA20)}
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQ"


def tokens(protein):
    return [CLS] + [AA_IDS[a] for a in protein] + [EOS]


def favoured(pos):
    return AA20[(pos * 3) % 20]


def fake_model(calls):
    """Logits are all zero except at a masked position, where the 'favoured' amino acid gets +5."""
    def fn(batch):
        calls.append(batch.shape)
        out = np.zeros((batch.shape[0], batch.shape[1], V))
        for b in range(batch.shape[0]):
            for t in range(batch.shape[1]):
                if batch[b, t] == MASK:
                    out[b, t, AA_IDS[favoured(t)]] = 5.0
        return out
    return fn


def test_masked_mode_reads_the_masked_position():
    calls = []
    rows = position_logits(fake_model(calls), tokens(PROTEIN), AA_IDS, mode="masked", mask_id=MASK, batch_size=4)
    assert sorted(rows) == list(range(1, len(PROTEIN) + 1))
    for pos, row in rows.items():
        assert int(np.argmax(row)) == AA20.index(favoured(pos)), f"position {pos} read from the wrong index"
    assert len(calls) == -(-len(PROTEIN) // 4)  # ceil(L / batch)
    assert all(shape[1] == len(PROTEIN) + 2 for shape in calls)


def test_wt_mode_needs_a_single_pass_and_sees_no_mask():
    calls = []
    rows = position_logits(fake_model(calls), tokens(PROTEIN), AA_IDS, mode="wt")
    assert len(calls) == 1 and calls[0][0] == 1
    assert all(np.allclose(r, 0) for r in rows.values())


def test_subset_of_positions_and_progress_callback():
    seen = []
    rows = position_logits(fake_model([]), tokens(PROTEIN), AA_IDS, mode="masked", mask_id=MASK,
                           batch_size=2, positions=[3, 5, 9], progress=lambda d, t: seen.append((d, t)))
    assert sorted(rows) == [3, 5, 9] and seen[-1] == (3, 3)


def test_errors():
    with pytest.raises(ValueError):
        position_logits(fake_model([]), tokens(PROTEIN), AA_IDS, mode="masked")  # mask_id missing
    with pytest.raises(ValueError):
        position_logits(fake_model([]), tokens(PROTEIN), AA_IDS, mode="other", mask_id=MASK)
    with pytest.raises(ValueError):
        position_logits(fake_model([]), [CLS] + [AA_IDS["A"]] * (MAX_RESIDUES + 1) + [EOS], AA_IDS, mode="wt")


def test_llr_table_values():
    protein = "ACD"
    rows = {1: np.zeros(20), 2: np.zeros(20), 3: np.zeros(20)}
    rows[2][AA20.index("L")] = 5.0     # position 2 (wt C): L is favoured
    rows[3][AA20.index("D")] = 5.0     # position 3 (wt D): the wild type itself is favoured
    subs, stats = llr_table(protein, rows)
    assert len(subs) == 3 * 19
    lookup = {(p, m): v for p, _w, m, v in subs}
    assert lookup[(1, "G")] == pytest.approx(0.0)
    assert lookup[(2, "L")] == pytest.approx(5.0)
    assert lookup[(3, "A")] == pytest.approx(-5.0)
    by_pos = {p: (lp, ent) for p, _w, lp, ent in stats}
    assert by_pos[1][1] == pytest.approx(np.log(20))            # flat distribution: maximum entropy
    assert by_pos[3][0] > by_pos[2][0]                          # wild type is more likely at position 3


def test_unknown_residues_are_skipped():
    subs, stats = llr_table("AXC", {1: np.zeros(20), 2: np.zeros(20), 3: np.zeros(20)})
    assert {p for p, *_ in subs} == {1, 3} and {p for p, *_ in stats} == {1, 3}
