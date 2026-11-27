import csv
import json

import pytest
from helpers import build_cds, random_protein
from pdb_helpers import build_pdb

pytest.importorskip("Bio")
from doublebind.cli import main  # noqa: E402
from doublebind.stability import AA20  # noqa: E402


@pytest.fixture()
def ws(tmp_path):
    protein = random_protein(120, seed=31)
    raw = tmp_path / "raw"
    (raw / "cds").mkdir(parents=True)
    (raw / "alphafold").mkdir()
    (raw / "cds" / "gyra_cds.fasta").write_text(f">gyra\n{build_cds(protein, seed=31)}\n")
    text, _, _ = build_pdb(protein, first_resseq=1, with_extras=False)
    (raw / "alphafold" / "AF-P0AES4-F1.pdb").write_text(text)
    directions = tmp_path / "directions.csv"
    directions.write_text("target,mutation,direction,numbering,evidence_status\n")
    return {"raw": str(raw), "out": str(tmp_path / "processed"), "dir": str(directions), "tmp": tmp_path, "protein": protein}


def _common(ws):
    return ["--raw-dir", ws["raw"], "--out-dir", ws["out"]]


def _write_thermompnn_csv(path, seq, sign=1.0):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "Model", "Dataset", "ddG_pred", "position", "wildtype", "mutation", "pdb", "chain"])
        i = 0
        for p, aa in enumerate(seq):
            for m in AA20:
                w.writerow([i, "ThermoMPNN", "x", sign * ((2.0 if m == "P" else 0.3) + 0.001 * p), p, aa, m, "x", "A"])
                i += 1


def test_prep_import_table_flow(ws, capsys):
    assert main(["stability-prep", "gyra", *_common(ws)]) == 0
    text = capsys.readouterr().out
    assert "identity 1.00" in text and "covers 100%" in text and "ThermoMPNN_inference_gyra_AF-P0AES4-F1_A.csv" in text
    stab = ws["tmp"] / "processed" / "stability"
    assert (stab / "gyra_AF-P0AES4-F1_A.pdb").exists() and (stab / "gyra_AF-P0AES4-F1_A.map.json").exists()

    csv_path = ws["tmp"] / "ThermoMPNN_inference_gyra_AF-P0AES4-F1_A.csv"
    _write_thermompnn_csv(csv_path, ws["protein"])
    assert main(["stability-import", "gyra", "--csv", str(csv_path), "--directions", ws["dir"], *_common(ws)]) == 0
    text = capsys.readouterr().out
    assert "coverage of reachable variants" in text and "100.0%" in text and "PASS (positive = destabilizing)" in text
    rows = list(csv.DictReader(open(f"{ws['out']}/gyra_ddg_fold.csv")))
    assert rows and {"pos", "wt_aa", "mut_aa", "ddg_fold", "source"} <= set(rows[0])
    report = json.load(open(f"{ws['out']}/gyra_stability_report.json"))
    assert report["canary"]["passed"] is True and report["coverage"]["all"][0] == report["coverage"]["all"][1]

    assert main(["table", "gyra", *_common(ws), "--known", ws["dir"]]) == 0
    assert "with a stability score" in capsys.readouterr().out
    table = list(csv.DictReader(open(f"{ws['out']}/gyra_variant_table.csv")))
    assert all(r["ddg_fold"] != "" for r in table)


def test_import_exits_with_an_error_when_the_sign_looks_wrong(ws, capsys):
    assert main(["stability-prep", "gyra", *_common(ws)]) == 0
    capsys.readouterr()
    csv_path = ws["tmp"] / "flipped.csv"
    _write_thermompnn_csv(csv_path, ws["protein"], sign=-1.0)
    assert main(["stability-import", "gyra", "--csv", str(csv_path), "--directions", ws["dir"], *_common(ws)]) == 1
    assert "FAIL: check the sign" in capsys.readouterr().out


def test_prep_without_the_alphafold_file_says_what_to_do(ws, capsys):
    (ws["tmp"] / "raw" / "alphafold" / "AF-P0AES4-F1.pdb").unlink()
    assert main(["stability-prep", "gyra", *_common(ws)]) == 1
    assert "data --alphafold" in capsys.readouterr().out


def test_prep_accepts_a_structure_file(ws, capsys):
    text, _, _ = build_pdb(ws["protein"], first_resseq=7, with_extras=False)
    mine = ws["tmp"] / "mine.pdb"
    mine.write_text(text)
    assert main(["stability-prep", "gyra", "--file", str(mine), *_common(ws)]) == 0
    assert (ws["tmp"] / "processed" / "stability" / "gyra_mine_A.pdb").exists()
