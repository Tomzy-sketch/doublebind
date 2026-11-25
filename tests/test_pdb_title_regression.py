from doublebind.data.fetch import pdb_title


def test_title_far_down_the_file_is_found(tmp_path):
    """Day 1 bug: verify printed '(title not parsed)' because only the first 20,000 characters were read."""
    filler = "".join(f"_entity_poly.seq{i} {'A' * 60}\n" for i in range(600))
    path = tmp_path / "big.cif"
    path.write_text("data_BIG\n" + filler + "_struct.title   'COMPLEX WITH 5'-D(AT)-3' DNA'\n#\n")
    assert path.stat().st_size > 30000
    assert pdb_title(path) == "COMPLEX WITH 5'-D(AT)-3' DNA"


def test_missing_title(tmp_path):
    path = tmp_path / "none.cif"
    path.write_text("data_X\n_entry.id X\n")
    assert pdb_title(path) == "(title not found)"
