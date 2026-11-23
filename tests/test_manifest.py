from doublebind.data.manifest import Manifest, sha256_file


def test_record_and_verify(tmp_path):
    raw = tmp_path / "raw"
    (raw / "sub").mkdir(parents=True)
    f = raw / "sub" / "a.txt"
    f.write_text("hello")
    m = Manifest(raw / "MANIFEST.json")
    m.record(f, "test")
    assert m.verify() == [("sub/a.txt", "ok")]
    assert Manifest(raw / "MANIFEST.json").entries["sub/a.txt"]["sha256"] == sha256_file(f)


def test_changed_and_missing_files_are_caught(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    a, b = raw / "a.txt", raw / "b.txt"
    a.write_text("one")
    b.write_text("two")
    m = Manifest(raw / "MANIFEST.json")
    m.record(a, "t")
    m.record(b, "t")
    a.write_text("ONE")
    b.unlink()
    assert dict(m.verify()) == {"a.txt": "changed", "b.txt": "missing"}
