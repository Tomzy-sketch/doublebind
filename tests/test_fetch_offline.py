import json

import pytest
from helpers import build_cds, random_protein

from doublebind.data import fetch
from doublebind.data.manifest import Manifest
from doublebind.data.targets import Target


def test_parse_fasta():
    text = ">a one\nATG\nAAA\n\n>b two\nCCC\n"
    assert fetch.parse_fasta(text) == [("a one", "ATGAAA"), ("b two", "CCC")]


def test_best_cds_match_picks_the_right_gene():
    protein = random_protein(200, seed=1)
    decoys = [(f"decoy{i}", build_cds(random_protein(200, seed=10 + i), seed=i)) for i in range(5)]
    records = decoys[:2] + [("target gene", build_cds(protein, seed=99))] + decoys[2:]
    m = fetch.best_cds_match(records, protein)
    assert m.header == "target gene" and m.identity == 1.0 and m.differences == ()


def test_best_cds_match_reports_a_single_difference():
    protein = random_protein(300, seed=2)
    cds_protein = protein[:99] + ("A" if protein[99] != "A" else "C") + protein[100:]
    m = fetch.best_cds_match([("variant gene", build_cds(cds_protein, seed=4))], protein)
    assert len(m.differences) == 1 and m.differences[0][0] == 100
    assert m.identity == pytest.approx(299 / 300)


def test_best_cds_match_returns_none_when_nothing_is_close():
    protein = random_protein(200, seed=1)
    short = ("short", build_cds(random_protein(50, seed=7), seed=7))
    assert fetch.best_cds_match([short], protein) is None


class _Resp:
    def __init__(self, data):
        self._d = data

    def json(self):
        return self._d


def test_chembl_pages_follow_next(monkeypatch):
    base = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    pages = {
        base: {"activities": [1, 2], "page_meta": {"next": "/chembl/api/data/activity.json?limit=2&offset=2"}},
        base + "?limit=2&offset=2": {"activities": [3], "page_meta": {"next": None}},
    }
    monkeypatch.setattr(fetch, "http_get", lambda url, params=None, **kw: _Resp(pages[url]))
    got = [a for page in fetch.chembl_pages("activity") for a in page["activities"]]
    assert got == [1, 2, 3]


def _fake_target(key="demo"):
    return Target(key=key, protein="Demo", gene="demo", organism="Test", uniprot="P00000",
                  expected_length=120, genome_accession="TEST.1", pdb_ids=())


def test_fetch_chembl_writes_files_and_counts(monkeypatch, tmp_path):
    def fake_pages(path, params=None):
        if path == "target":
            yield {"targets": [{"target_chembl_id": "CHEMBL1"}]}
        elif path == "activity":
            yield {"activities": [{"molecule_chembl_id": "CHEMBLA"}, {"molecule_chembl_id": "CHEMBLB"}]}
        elif path == "molecule":
            yield {"molecules": [
                {"molecule_chembl_id": "CHEMBLA", "max_phase": "4.0", "pref_name": "CIPROFLOXACIN",
                 "molecule_structures": {"canonical_smiles": "CCO"}},
                {"molecule_chembl_id": "CHEMBLB", "max_phase": None, "pref_name": None, "molecule_structures": None},
            ]}

    monkeypatch.setattr(fetch, "chembl_pages", fake_pages)
    raw = tmp_path / "raw"
    manifest = Manifest(raw / "MANIFEST.json")
    summary = fetch.fetch_chembl(_fake_target(), raw, manifest)
    assert summary == {"targets": 1, "activities": 2, "molecules": 2, "approved": 1}
    assert all(status == "ok" for _, status in manifest.verify())
    again = fetch.fetch_chembl(_fake_target(), raw, manifest)  # second call uses the cache
    assert again == {"cached": True}


def test_extract_and_verify_end_to_end(tmp_path):
    protein = random_protein(120, seed=8)
    target = _fake_target()
    raw = tmp_path / "raw"
    (raw / "uniprot").mkdir(parents=True)
    (raw / "ncbi").mkdir()
    (raw / "uniprot" / "P00000.fasta").write_text(f">sp|P00000|DEMO\n{protein}\n")
    (raw / "uniprot" / "P00000.json").write_text(json.dumps({
        "sequence": {"length": 120}, "genes": [{"geneName": {"value": "demo"}}],
        "organism": {"scientificName": "Test organism"},
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Demo protein"}}}}))
    genome = [(f"decoy{i}", build_cds(random_protein(120, seed=20 + i), seed=i)) for i in range(3)]
    genome.append(("[gene=demo]", build_cds(protein, seed=5)))
    (raw / "ncbi" / "TEST.1_cds.fasta").write_text("".join(f">{h}\n{s}\n" for h, s in genome))

    manifest = Manifest(raw / "MANIFEST.json")
    for p in (raw / "uniprot" / "P00000.fasta", raw / "uniprot" / "P00000.json"):
        manifest.record(p, "test")
    match = fetch.extract_target_cds(target, raw, manifest)
    assert match.identity == 1.0
    assert (raw / "cds" / "demo_cds.fasta").exists()

    from doublebind.data import targets as targets_mod
    targets_mod.TARGETS["demo"] = target
    try:
        checks = {c.name: c for c in fetch.verify_all(["demo"], raw)}
    finally:
        del targets_mod.TARGETS["demo"]
    assert checks["manifest"].status == "PASS"
    assert checks["demo: coding sequence"].status == "PASS"
    assert checks["demo: UniProt P00000"].status == "PASS"
    assert checks["demo: ChEMBL activities"].status == "WARN"  # not downloaded in this test


def test_fetch_all_collects_errors_instead_of_crashing(monkeypatch, tmp_path):
    def broken(url, params=None, **kw):
        raise RuntimeError("network blocked")

    monkeypatch.setattr(fetch, "http_get", broken)
    messages = []
    errors = fetch.fetch_all(["tem1"], tmp_path / "raw", chembl=True, log=messages.append)
    steps = [step for step, _ in errors]
    assert any(s.startswith("uniprot") for s in steps) and any(s.startswith("ncbi") for s in steps)
    assert any("finished: " in m and "failed" in m for m in messages)


def test_pdb_title_parsing(tmp_path):
    quoted = tmp_path / "a.cif"
    quoted.write_text("data_X\n#\n_struct.entry_id   X\n_struct.title   'BETA-LACTAMASE TEM-1'\n")
    block = tmp_path / "b.cif"
    block.write_text("data_Y\n_struct.title\n;Crystal structure of a ternary\ncomplex with ciprofloxacin\n;\n")
    assert fetch._pdb_title(quoted) == "BETA-LACTAMASE TEM-1"
    assert fetch._pdb_title(block) == "Crystal structure of a ternary complex with ciprofloxacin"
