from helpers import build_cds, random_protein

from doublebind.cli import main


def test_numbers_command(capsys):
    assert main(["numbers"]) == 0
    out = capsys.readouterr().out
    assert "0.616" in out and "-0.854" in out and "20.6x more susceptible" in out


def test_variants_command_on_a_synthetic_gyra(tmp_path, capsys):
    # S83 (TCG) and D87 (GAC) can reach L and N by one nucleotide change
    protein = random_protein(120, seed=11, fixed={83: "S", 87: "D"})
    cds = build_cds(protein, overrides={83: "TCG", 87: "GAC"}, seed=11)
    fasta = tmp_path / "gyra_cds.fasta"
    fasta.write_text(f">gyra\n{cds}\n")
    known = tmp_path / "known.csv"
    known.write_text("target,mutation,numbering,note,evidence_status\n"
                     "gyra,S83L,uniprot,test,test\ngyra,D87N,uniprot,test,test\ngyra,S83W,uniprot,test,test\ngyra,S83K,uniprot,test,test\n")
    out = tmp_path / "out"
    rc = main(["variants", "gyra", "--cds", str(fasta), "--out-dir", str(out), "--known", str(known)])
    text = capsys.readouterr().out
    assert rc == 0
    assert (out / "gyra_variants_nt.csv").exists() and (out / "gyra_variants_aa.csv").exists()
    assert "numbering offset +0" in text
    assert "S83L = residue 83: reachable" in text
    assert "D87N = residue 87: reachable" in text
    assert "S83W = residue 83: reachable" in text      # TCG -> TGG is one change (C>G)
    assert "S83K = residue 83: NOT reachable" in text  # TCG -> AAG needs two changes
