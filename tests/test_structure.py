import numpy as np
import pytest
from helpers import random_protein
from pdb_helpers import build_pdb

pytest.importorskip("Bio")
from doublebind.structure import proximity as prox  # noqa: E402


def _model(tmp_path, seq, **kw):
    text, lig, res = build_pdb(seq, **kw)
    path = tmp_path / "t.pdb"
    path.write_text(text)
    return prox.first_model(prox.load_structure(path)), lig, res


def test_chains_and_hetero_groups(tmp_path):
    seq = random_protein(40, seed=1)
    model, _, _ = _model(tmp_path, seq, first_resseq=101, ligand_index=20)
    chains = prox.protein_chains(model)
    assert len(chains) == 1 and chains[0].sequence == seq and chains[0].residues[0].id[1] == 101
    groups = prox.hetero_groups(model)
    assert {g.name for g in groups} == {"LIG", "MG"}          # water is hidden
    assert [g.name for g in prox.candidate_ligands(groups)] == ["LIG"]  # the ion is not a drug


def test_small_hetero_groups_are_not_candidates(tmp_path):
    text, _, _ = build_pdb(random_protein(30, seed=2), ligand_index=10)
    (tmp_path / "s.pdb").write_text(text)
    groups = prox.hetero_groups(prox.first_model(prox.load_structure(tmp_path / "s.pdb")))
    assert prox.candidate_ligands(groups, min_heavy=50) == []


def test_mapping_handles_numbering_offset_and_flanks(tmp_path):
    seq = random_protein(80, seed=3)
    model, _, _ = _model(tmp_path, seq, first_resseq=101)
    target = "WWWWW" + seq + "WWWWW"
    m = prox.map_chain_to_target(prox.protein_chains(model)[0], target)
    assert m.identity == 1.0 and m.aligned == 80
    assert m.mapping[6].id[1] == 101 and m.mapping[85].id[1] == 180   # target 6 -> first residue
    assert 1 not in m.mapping and 86 not in m.mapping


def test_mapping_reports_a_mismatch(tmp_path):
    seq = random_protein(80, seed=4)
    changed = seq[:40] + ("W" if seq[40] != "W" else "A") + seq[41:]
    model, _, _ = _model(tmp_path, seq)
    m = prox.map_chain_to_target(prox.protein_chains(model)[0], changed)
    assert m.aligned == 80 and m.identity == pytest.approx(79 / 80)


def test_distances_pick_the_residue_next_to_the_ligand(tmp_path):
    seq = random_protein(60, seed=5)
    model, lig, res = _model(tmp_path, seq, first_resseq=10, ligand_index=30)
    ligands = prox.candidate_ligands(prox.hetero_groups(model))
    dists, maps = prox.ligand_distances(model, seq, ligands, min_aligned=20)
    assert len(maps) == 1 and len(dists) == 60
    nearest = min(dists, key=lambda d: d.min_dist)
    assert nearest.pos == 31 and nearest.resseq == 40          # index 30 -> position 31, residue number 40
    lig_xyz = np.array(lig)
    expected = min(np.linalg.norm(np.array(xyz) - lig_xyz, axis=1).min() for xyz in res[40])
    assert nearest.min_dist == pytest.approx(expected, abs=1e-3)
    far = max(dists, key=lambda d: d.min_dist)
    assert far.min_dist > nearest.min_dist + 20


def test_distances_respect_the_chain_filter_and_identity_floor(tmp_path):
    seq = random_protein(60, seed=6)
    model, _, _ = _model(tmp_path, seq, ligand_index=30)
    ligands = prox.candidate_ligands(prox.hetero_groups(model))
    with pytest.raises(ValueError):
        prox.ligand_distances(model, seq, [], min_aligned=20)
    dists, maps = prox.ligand_distances(model, seq, ligands, chains={"Z"}, min_aligned=20)
    assert dists == [] and maps == []
    unrelated = random_protein(60, seed=99)
    dists, maps = prox.ligand_distances(model, unrelated, ligands, min_identity=0.9, min_aligned=20)
    assert dists == []


def test_mmcif_path_matches_the_pdb_path(tmp_path):
    """RCSB serves .cif files, so the .cif reader has to give the same answers as the .pdb reader."""
    from pdb_helpers import build_cif

    from doublebind.data.fetch import pdb_title

    seq = random_protein(60, seed=7)
    cif = tmp_path / "t.cif"
    cif.write_text(build_cif(seq, first_resseq=10, ligand_index=30, title="A TITLE WITH 5'-D(AT)-3' IN IT"))
    model = prox.first_model(prox.load_structure(cif))
    chains = prox.protein_chains(model)
    assert len(chains) == 1 and chains[0].sequence == seq and chains[0].residues[0].id[1] == 10
    groups = prox.hetero_groups(model)
    assert {g.name for g in groups} == {"LIG", "MG"}
    ligands = prox.candidate_ligands(groups)
    dists, _ = prox.ligand_distances(model, seq, ligands, min_aligned=20)
    nearest = min(dists, key=lambda d: d.min_dist)
    assert nearest.pos == 31 and nearest.resseq == 40
    assert pdb_title(cif) == "A TITLE WITH 5'-D(AT)-3' IN IT"
