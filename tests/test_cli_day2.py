import csv
import json
import random

import pytest
from helpers import build_cds, random_protein
from pdb_helpers import build_pdb

pytest.importorskip("Bio")
from doublebind.cli import main  # noqa: E402


@pytest.fixture()
def workspace(tmp_path):
    protein = random_protein(120, seed=11, fixed={83: "S", 87: "D"})
    cds = build_cds(protein, overrides={83: "TCG", 87: "GAC"}, seed=11)
    raw = tmp_path / "raw"
    (raw / "cds").mkdir(parents=True)
    (raw / "pdb").mkdir()
    (raw / "cds" / "gyra_cds.fasta").write_text(f">gyra\n{cds}\n")
    # structure covers residues 20..110 of the target, numbered from 1000, drug 5 A above residue 85
    text, _, _ = build_pdb(protein[19:110], first_resseq=1000, ligand_index=85 - 20)
    (raw / "pdb" / "9TST.pdb").write_text(text)
    known = tmp_path / "known.csv"
    known.write_text("target,mutation,numbering,note,evidence_status\ngyra,S83L,uniprot,t,t\ngyra,D87N,uniprot,t,t\n")
    return {"raw": str(raw), "out": str(tmp_path / "processed"), "known": str(known), "tmp": tmp_path, "protein": protein}


def _args(ws, *extra):
    return ["--raw-dir", ws["raw"], "--out-dir", ws["out"], "--known", ws["known"], *extra]


def test_inspect_lists_chains_and_ligand(workspace, capsys):
    assert main(["inspect", "9TST", "--target", "gyra", "--raw-dir", workspace["raw"]]) == 0
    text = capsys.readouterr().out
    assert "protein chains" in text and "identity 1.00" in text
    assert "LIG" in text and "candidate ligand" in text and "MG" in text


def test_pocket_baselines_and_table_flow(workspace, capsys):
    assert main(["pocket", "gyra", "--pdb", "9TST", *_args(workspace)]) == 0
    text = capsys.readouterr().out
    assert "ligand(s) LIG" in text and "known S83L" in text and "known D87N" in text
    pocket = list(csv.DictReader(open(f"{workspace['out']}/gyra_pocket.csv")))
    assert len(pocket) == 91 and {r["pos"] for r in pocket} >= {"83", "87"}
    near = {int(r["pos"]) for r in pocket if r["in_pocket"] == "1"}
    assert {83, 84, 85, 86, 87} <= near

    # fabricate ESM scores so the rest of the flow can run without model weights
    rng = random.Random(1)
    esm_path = workspace["tmp"] / "processed" / "gyra_esm_esm2_t12_35M_UR50D.csv"
    with open(esm_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pos", "wt_aa", "mut_aa", "llr", "reachable", "mode"])
        for pos in range(2, 121):
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                llr = 3.0 if (pos, aa) in {(83, "L"), (87, "N")} else rng.uniform(-6, 1)
                w.writerow([pos, workspace["protein"][pos - 1], aa, f"{llr:.3f}", 1, "masked"])
    assert main(["baselines", "gyra", *_args(workspace)]) == 0
    out = capsys.readouterr().out
    assert "proximity only" in out and "ESM tolerated" in out and "proximity <= 8 A, then ESM tolerated" in out
    result = json.load(open(f"{workspace['out']}/gyra_baselines.json"))
    assert result["position_level"]["proximity_then_esm"]["known_ranks"]["83"] <= 2

    assert main(["table", "gyra", *_args(workspace)]) == 0
    assert "with an ESM score" in capsys.readouterr().out
    rows = list(csv.DictReader(open(f"{workspace['out']}/gyra_variant_table.csv")))
    assert rows and {"esm_llr", "dist_min", "in_pocket"} <= set(rows[0])


def test_pocket_without_a_ligand_lists_what_is_there(workspace, capsys):
    text, _, _ = build_pdb(workspace["protein"][19:110], first_resseq=1000, ligand_index=None)
    (workspace["tmp"] / "raw" / "pdb" / "8NOL.pdb").write_text(text)
    assert main(["pocket", "gyra", "--pdb", "8NOL", *_args(workspace)]) == 1
    assert "no ligand found" in capsys.readouterr().out
