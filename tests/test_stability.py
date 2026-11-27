import csv
import json

import pytest
from helpers import build_cds, random_protein
from pdb_helpers import build_pdb

pytest.importorskip("Bio")
from doublebind.data import fetch  # noqa: E402
from doublebind.data.codons import aggregate_missense, enumerate_variants  # noqa: E402
from doublebind.stability import (AA20, coverage, direction_report, import_thermompnn, proline_canary,  # noqa: E402
                                  read_thermompnn_csv)
from doublebind.structure import proximity as prox  # noqa: E402
from doublebind.structure.prep import write_clean_pdb  # noqa: E402

REAL_SAMPLE = """,Model,Dataset,ddG_pred,position,wildtype,mutation,pdb,chain
0,ThermoMPNN,2OCJ,-0.04087644815444946,0,S,A,2OCJ,A
1,ThermoMPNN,2OCJ,-0.3005809187889099,0,S,C,2OCJ,A
2,ThermoMPNN,2OCJ,0.009091079235076904,0,S,D,2OCJ,A
3,ThermoMPNN,2OCJ,-0.051891446113586426,0,S,E,2OCJ,A
"""


def fake_rows(seq, value=lambda p, m: 0.01 * p + 0.1 * AA20.index(m) + (2.0 if m == "P" else 0.0)):
    """Rows shaped like real ThermoMPNN output: all 20 letters at every position, including wild type -> wild type."""
    return [{"": str(i), "Model": "ThermoMPNN", "Dataset": "x", "ddG_pred": repr(value(p, m)), "position": str(p),
             "wildtype": aa, "mutation": m, "pdb": "x", "chain": "A"}
            for p, aa in enumerate(seq) for i, m in enumerate(AA20)]


def prepared(tmp_path, seq, first_resseq=40, **kw):
    text, _, _ = build_pdb(seq, first_resseq=first_resseq, ligand_index=20, **kw)
    src = tmp_path / "in.pdb"
    src.write_text(text)
    return write_clean_pdb(src, tmp_path / "out" / "s_A.pdb")


def test_real_thermompnn_header_parses(tmp_path):
    f = tmp_path / "real.csv"
    f.write_text(REAL_SAMPLE)
    rows = read_thermompnn_csv(f)
    assert len(rows) == 4 and rows[1]["wildtype"] == "S" and rows[1]["mutation"] == "C"
    assert float(rows[1]["ddG_pred"]) == pytest.approx(-0.3005809187889099)
    bad = tmp_path / "bad.csv"
    bad.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError):
        read_thermompnn_csv(bad)


def test_clean_pdb_renumbers_and_keeps_a_way_back(tmp_path):
    seq = random_protein(50, seed=21)
    prep = prepared(tmp_path, seq, first_resseq=101)
    assert prep.sequence == seq and prep.n_skipped == 0
    assert prep.index_map[0]["resseq"] == 101 and prep.index_map[49]["resseq"] == 150
    model = prox.first_model(prox.load_structure(prep.pdb_path))
    chain = prox.protein_chains(model)[0]
    assert chain.chain_id == "A" and [r.id[1] for r in chain.residues] == list(range(1, 51))
    assert {a.get_id() for a in chain.residues[0]} == {"N", "CA", "C", "O"}
    assert not prox.hetero_groups(model)                     # the ligand is not written
    assert json.loads(prep.map_path.read_text())["sequence"] == seq


def test_residues_without_a_full_backbone_are_skipped(tmp_path):
    seq = random_protein(50, seed=22)
    text, _, _ = build_pdb(seq, first_resseq=1, with_extras=False)
    kept = [ln for ln in text.splitlines() if not (ln.startswith("ATOM") and ln[12:16] == " O  " and int(ln[22:26]) == 10)]
    src = tmp_path / "gap.pdb"
    src.write_text("\n".join(kept) + "\n")
    prep = write_clean_pdb(src, tmp_path / "gap_A.pdb")
    assert prep.n_skipped == 1 and len(prep.index_map) == 49
    assert prep.sequence == seq[:9] + seq[10:] and prep.index_map[9]["resseq"] == 11


def test_no_protein_chain_is_an_error(tmp_path):
    src = tmp_path / "empty.pdb"
    src.write_text("END\n")
    with pytest.raises(ValueError):
        write_clean_pdb(src, tmp_path / "x_A.pdb")


def test_import_maps_structure_index_to_target_numbering(tmp_path):
    seq = random_protein(70, seed=23)
    prep = prepared(tmp_path, seq)
    target = "WWWWW" + seq + "WW"                     # structure residue i is target position i + 6
    imp = import_thermompnn(fake_rows(seq), prep.index_map, target)
    assert imp.identity == 1.0 and imp.n_structure_mismatch == 0 and imp.n_unmapped == 0
    assert len(imp.rows) == 70 * 19                    # wild type -> wild type rows are dropped
    lookup = {(r["pos"], r["mut_aa"]): r["ddg_fold"] for r in imp.rows}
    assert lookup[(6 + 10, "P")] == pytest.approx(0.01 * 10 + 0.1 * AA20.index("P") + 2.0)
    assert all(r["wt_aa"] == target[r["pos"] - 1] for r in imp.rows)
    assert min(r["pos"] for r in imp.rows) == 6 and max(r["pos"] for r in imp.rows) == 75


def test_sign_is_converted_for_a_stabilizing_positive_source(tmp_path):
    seq = random_protein(40, seed=24)
    prep = prepared(tmp_path, seq)
    plain = {(r["pos"], r["mut_aa"]): r["ddg_fold"] for r in import_thermompnn(fake_rows(seq), prep.index_map, seq).rows}
    flipped = {(r["pos"], r["mut_aa"]): r["ddg_fold"]
               for r in import_thermompnn(fake_rows(seq), prep.index_map, seq, source="megascale").rows}
    assert all(flipped[k] == pytest.approx(-v) for k, v in plain.items())


def test_bad_rows_are_counted_not_used(tmp_path):
    seq = random_protein(40, seed=25)
    prep = prepared(tmp_path, seq)
    rows = fake_rows(seq)
    rows[1]["wildtype"] = "W" if seq[0] != "W" else "A"          # disagrees with the structure
    rows[2]["position"] = "999"                                  # outside the chain
    rows[3]["mutation"] = "X"                                    # ThermoMPNN's unknown-residue letter
    imp = import_thermompnn(rows, prep.index_map, seq)
    assert imp.n_structure_mismatch == 2 and len(imp.rows) == 40 * 19 - 3     # two bad rows plus the 'X' row
    homolog = seq[:20] + ("A" if seq[20] != "A" else "C") + seq[21:]   # the target differs at one position
    imp2 = import_thermompnn(fake_rows(seq), prep.index_map, homolog)
    assert imp2.n_target_mismatch == 19 and all(r["pos"] != 21 for r in imp2.rows)
    imp3 = import_thermompnn(fake_rows(seq), prep.index_map, homolog, allow_mismatch=True)
    assert any(r["pos"] == 21 for r in imp3.rows)


def test_coverage_with_and_without_a_folded_region(tmp_path):
    protein = random_protein(100, seed=26)
    aa = aggregate_missense(enumerate_variants(build_cds(protein, seed=26)))
    rows = [{"pos": v.pos, "mut_aa": v.mut_aa, "ddg_fold": 0.0} for v in aa if v.pos >= 11]   # positions 2-10 uncovered
    cov = coverage(rows, aa, region=(11, 100))
    assert cov["region"][0] == cov["region"][1]                  # everything inside the region is covered
    assert cov["all"][0] < cov["all"][1]                         # but not the whole protein
    assert "region" not in coverage(rows, aa)


def test_proline_canary_passes_and_fails():
    good = [{"mut_aa": m, "ddg_fold": (2.0 if m == "P" else 0.5)} for m in AA20 for _ in range(5)]
    assert proline_canary(good)["passed"] is True
    flipped = [{"mut_aa": m, "ddg_fold": (-2.0 if m == "P" else -0.5)} for m in AA20 for _ in range(5)]
    assert proline_canary(flipped)["passed"] is False
    assert proline_canary([])["passed"] is None


def test_direction_report_uses_the_numbering_offset():
    protein = random_protein(60, seed=27, fixed={38: "G", 28: "M"})   # literature numbers are sequence + 2
    rows = [{"pos": 38, "mut_aa": "S", "ddg_fold": 1.2}, {"pos": 28, "mut_aa": "T", "ddg_fold": -0.6}]
    expected = [{"mutation": "G40S", "direction": "destabilizing"}, {"mutation": "M30T", "direction": "stabilizing"}]
    rep = direction_report(rows, protein, expected)
    assert rep["checked"] == 2 and rep["agree"] == 2 and rep["offset"] == -2
    flipped = [dict(r, ddg_fold=-r["ddg_fold"]) for r in rows]
    assert direction_report(flipped, protein, expected)["agree"] == 0
    assert direction_report(rows, protein, [])["checked"] == 0
    impossible = [{"mutation": f"W{n}A", "direction": "stabilizing"} for n in (5, 6, 7, 8)]   # no run of four W
    assert "problem" in direction_report(rows, protein, impossible)


class _Resp:
    def __init__(self, content=b"", data=None):
        self.content, self._data = content, data

    def json(self):
        return self._data


def _target():
    from doublebind.data.targets import Target
    return Target(key="demo", protein="Demo", gene="demo", organism="Test", uniprot="P0AES4", expected_length=10,
                  genome_accession="TEST.1", pdb_ids=())


def test_alphafold_download_uses_the_api_url(monkeypatch, tmp_path):
    from doublebind.data.manifest import Manifest
    calls = []

    def fake(url, params=None, **kw):
        calls.append(url)
        if url.endswith("/api/prediction/P0AES4"):
            return _Resp(data=[{"pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P0AES4-F1-model_v9.pdb"}])
        if url.endswith("model_v9.pdb"):
            return _Resp(content=b"ATOM\n")
        raise RuntimeError("unexpected " + url)

    monkeypatch.setattr(fetch, "http_get", fake)
    manifest = Manifest(tmp_path / "MANIFEST.json")
    path = fetch.fetch_alphafold(_target(), tmp_path, manifest)
    assert path.name == "AF-P0AES4-F1.pdb" and path.read_bytes() == b"ATOM\n" and len(calls) == 2
    assert all(status == "ok" for _, status in manifest.verify())


def test_alphafold_download_falls_back_and_explains_failure(monkeypatch, tmp_path):
    from doublebind.data.manifest import Manifest

    def api_down(url, params=None, **kw):
        if "/api/" in url:
            raise RuntimeError("api down")
        if url.endswith("model_v6.pdb"):
            return _Resp(content=b"ATOM v6\n")
        raise RuntimeError("missing")

    monkeypatch.setattr(fetch, "http_get", api_down)
    assert fetch.fetch_alphafold(_target(), tmp_path, Manifest(tmp_path / "M.json")).read_bytes() == b"ATOM v6\n"

    def everything_down(url, params=None, **kw):
        raise RuntimeError("offline")

    monkeypatch.setattr(fetch, "http_get", everything_down)
    with pytest.raises(RuntimeError, match="Download it by hand"):
        fetch.fetch_alphafold(_target(), tmp_path / "other", Manifest(tmp_path / "other" / "M.json"))


def test_pdb_ca_sequence_and_verify_check(tmp_path):
    seq = random_protein(40, seed=28)
    text, _, _ = build_pdb(seq, with_extras=False)
    af = tmp_path / "raw" / "alphafold"
    af.mkdir(parents=True)
    (af / "AF-P0AES4-F1.pdb").write_text(text)
    assert fetch.pdb_ca_sequence(af / "AF-P0AES4-F1.pdb") == seq
    (tmp_path / "raw" / "uniprot").mkdir()
    (tmp_path / "raw" / "uniprot" / "P0AES4.fasta").write_text(f">x\n{seq}\n")
    (tmp_path / "raw" / "uniprot" / "P0AES4.json").write_text(json.dumps({"sequence": {"length": 40}}))
    from doublebind.data import targets as targets_mod
    targets_mod.TARGETS["demo"] = _target()
    try:
        check = {c.name: c for c in fetch.verify_all(["demo"], tmp_path / "raw")}["demo: AlphaFold model"]
        assert check.status == "PASS" and "matches UniProt" in check.detail
        (tmp_path / "raw" / "uniprot" / "P0AES4.fasta").write_text(f">x\n{seq[:-1]}A\n")
        check = {c.name: c for c in fetch.verify_all(["demo"], tmp_path / "raw")}["demo: AlphaFold model"]
        assert check.status == "WARN"
        (af / "AF-P0AES4-F1.pdb").unlink()
        check = {c.name: c for c in fetch.verify_all(["demo"], tmp_path / "raw")}["demo: AlphaFold model"]
        assert check.status == "WARN" and "data --alphafold" in check.detail
    finally:
        del targets_mod.TARGETS["demo"]
