"""Download raw data for the MVP targets (Day 1).

All network access goes through http_get, so everything else here is testable offline.
The endpoints below were written from memory of the public documentation and could not be
tested from the build sandbox. If a call fails, the error prints the exact URL.

  UniProt REST   https://rest.uniprot.org/uniprotkb/{accession}.json | .fasta
  NCBI E-utils   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi  (rettype=fasta_cds_na)
  RCSB files     https://files.rcsb.org/download/{PDB_ID}.cif
  ChEMBL API     https://www.ebi.ac.uk/chembl/api/data/...
  CARD download  https://card.mcmaster.ca/latest/data   (academic use only)
"""
from __future__ import annotations

import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path

from doublebind.data.codons import check_cds
from doublebind.data.manifest import Manifest
from doublebind.data.targets import Target, get_targets

USER_AGENT = "doublebind-starter/0.0.1 (research prototype)"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RCSB = "https://files.rcsb.org/download"
CHEMBL_HOST = "https://www.ebi.ac.uk"
CHEMBL_API = CHEMBL_HOST + "/chembl/api/data"
CARD_URL = "https://card.mcmaster.ca/latest/data"


def http_get(url, params=None, timeout=180, retries=4, backoff=3.0):
    import requests  # imported here so offline tests do not need it at import time

    last = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, params=params, timeout=timeout, headers={"User-Agent": USER_AGENT})
            if r.status_code == 429 or r.status_code >= 500:
                last = RuntimeError(f"HTTP {r.status_code}")
            else:
                r.raise_for_status()
                return r
        except requests.HTTPError as e:  # other 4xx errors: retrying will not help
            raise RuntimeError(f"GET {url} failed: {e}") from e
        except requests.RequestException as e:
            last = e
        time.sleep(backoff * attempt)
    raise RuntimeError(f"GET {url} failed after {retries} tries: {last}")


# ---------- small helpers ----------

def parse_fasta(text: str) -> list[tuple[str, str]]:
    records, header, chunks = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def _ensure_recorded(manifest: Manifest, dest: Path, source: str) -> None:
    if dest.relative_to(manifest.root).as_posix() not in manifest.entries:
        manifest.record(dest, source)


@dataclass
class CdsMatch:
    header: str
    cds: str
    identity: float
    differences: tuple
    length_protein: int
    length_cds: int


def best_cds_match(records, protein: str, length_tolerance: float = 0.1):
    """Find the coding sequence whose translation is closest to `protein` (None if nothing is close)."""
    best = None
    slack = 3 * max(10, int(length_tolerance * len(protein)))
    for header, seq in records:
        if abs(len(seq) - 3 * len(protein)) > slack + 3:
            continue
        try:
            chk = check_cds(seq, protein)
        except ValueError:
            continue  # not a clean coding sequence
        if best is None or chk.identity > best.identity:
            best = CdsMatch(header, seq, chk.identity, chk.differences, chk.length_protein, chk.length_cds)
    return best


# ---------- downloads ----------

def fetch_uniprot(target: Target, raw_dir: Path, manifest: Manifest, force=False) -> None:
    folder = raw_dir / "uniprot"
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ("json", "fasta"):
        dest = folder / f"{target.uniprot}.{ext}"
        url = f"{UNIPROT}/{target.uniprot}.{ext}"
        if not dest.exists() or force:
            dest.write_bytes(http_get(url).content)
        _ensure_recorded(manifest, dest, url)


def fetch_genome_cds(accession: str, raw_dir: Path, manifest: Manifest, force=False) -> Path:
    dest = raw_dir / "ncbi" / f"{accession}_cds.fasta"
    dest.parent.mkdir(parents=True, exist_ok=True)
    params = {"db": "nuccore", "id": accession, "rettype": "fasta_cds_na", "retmode": "text"}
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    if os.environ.get("NCBI_EMAIL"):
        params["email"] = os.environ["NCBI_EMAIL"]
        params["tool"] = "doublebind"
    source = f"{EFETCH}?db=nuccore&id={accession}&rettype=fasta_cds_na"
    if not dest.exists() or force:
        text = http_get(EFETCH, params=params).text
        if not text.lstrip().startswith(">"):
            raise RuntimeError(f"NCBI did not return FASTA for {accession}: {text[:200]!r}")
        dest.write_text(text, encoding="utf-8")
    _ensure_recorded(manifest, dest, source)
    return dest


def read_uniprot_sequence(raw_dir: Path, target: Target) -> str:
    path = raw_dir / "uniprot" / f"{target.uniprot}.fasta"
    records = parse_fasta(path.read_text(encoding="utf-8"))
    if not records:
        raise RuntimeError(f"empty FASTA: {path}")
    return records[0][1]


def extract_target_cds(target: Target, raw_dir: Path, manifest: Manifest) -> CdsMatch:
    """Pick the coding sequence that encodes the UniProt protein (matched by translation, not gene tag)."""
    genome = raw_dir / "ncbi" / f"{target.genome_accession}_cds.fasta"
    protein = read_uniprot_sequence(raw_dir, target)
    match = best_cds_match(parse_fasta(genome.read_text(encoding="utf-8")), protein)
    if match is None or match.identity < 0.9:
        raise RuntimeError(f"[{target.key}] no coding sequence in {target.genome_accession} resembles "
                           f"{target.uniprot}; check genome_accession in targets.py")
    folder = raw_dir / "cds"
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / f"{target.key}_cds.fasta"
    fasta.write_text(f">{target.key} {match.header}\n{match.cds}\n", encoding="utf-8")
    info = {
        "target": target.key, "uniprot": target.uniprot, "genome_accession": target.genome_accession,
        "header": match.header, "identity": match.identity,
        "differences": [list(d) for d in match.differences],
        "length_protein": match.length_protein, "length_cds_protein": match.length_cds,
    }
    info_path = folder / f"{target.key}_cds_match.json"
    info_path.write_text(json.dumps(info, indent=2), encoding="utf-8")
    for p in (fasta, info_path):
        manifest.record(p, f"derived from {target.genome_accession} by best translation match to {target.uniprot}")
    return match


def fetch_pdb(pdb_id: str, raw_dir: Path, manifest: Manifest, force=False) -> Path:
    pdb_id = pdb_id.upper()
    dest = raw_dir / "pdb" / f"{pdb_id}.cif"
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"{RCSB}/{pdb_id}.cif"
    if not dest.exists() or force:
        dest.write_bytes(http_get(url).content)
    _ensure_recorded(manifest, dest, url)
    return dest


def chembl_pages(path: str, params=None):
    """Yield each JSON page of a ChEMBL list endpoint, following page_meta.next."""
    url = f"{CHEMBL_API}/{path}.json"
    query = {"limit": 1000}
    query.update(params or {})
    while url:
        data = http_get(url, params=query).json()
        yield data
        nxt = (data.get("page_meta") or {}).get("next")
        url = None if not nxt else (nxt if nxt.startswith("http") else CHEMBL_HOST + nxt)
        query = None  # the next URL already carries its query string


def _is_phase4(value) -> bool:
    try:
        return float(value) == 4.0
    except (TypeError, ValueError):
        return False


def fetch_chembl(target: Target, raw_dir: Path, manifest: Manifest, force=False) -> dict:
    """Targets, activities and molecule phases for one UniProt accession. The Drug B pool is built on Day 4."""
    folder = raw_dir / "chembl"
    folder.mkdir(parents=True, exist_ok=True)
    f_targets = folder / f"{target.key}_targets.json"
    f_acts = folder / f"{target.key}_activities.jsonl"
    f_mols = folder / f"{target.key}_molecules.jsonl"
    source = f"{CHEMBL_API} (target_components__accession={target.uniprot})"
    if f_targets.exists() and f_acts.exists() and f_mols.exists() and not force:
        for p in (f_targets, f_acts, f_mols):
            _ensure_recorded(manifest, p, source)
        return {"cached": True}

    targets = []
    for page in chembl_pages("target", {"target_components__accession": target.uniprot}):
        targets.extend(page.get("targets", []))
    f_targets.write_text(json.dumps(targets, indent=2), encoding="utf-8")

    n_act, mol_ids = 0, set()
    with open(f_acts, "w", encoding="utf-8") as f:
        for t in targets:
            for page in chembl_pages("activity", {"target_chembl_id": t["target_chembl_id"]}):
                for a in page.get("activities", []):
                    f.write(json.dumps(a) + "\n")
                    n_act += 1
                    if a.get("molecule_chembl_id"):
                        mol_ids.add(a["molecule_chembl_id"])

    n_mol = n_approved = 0
    ids = sorted(mol_ids)
    with open(f_mols, "w", encoding="utf-8") as f:
        for i in range(0, len(ids), 50):
            for page in chembl_pages("molecule", {"molecule_chembl_id__in": ",".join(ids[i:i + 50])}):
                for m in page.get("molecules", []):
                    row = {
                        "molecule_chembl_id": m.get("molecule_chembl_id"),
                        "pref_name": m.get("pref_name"),
                        "max_phase": m.get("max_phase"),
                        "canonical_smiles": (m.get("molecule_structures") or {}).get("canonical_smiles"),
                    }
                    f.write(json.dumps(row) + "\n")
                    n_mol += 1
                    n_approved += int(_is_phase4(row["max_phase"]))
    for p in (f_targets, f_acts, f_mols):
        manifest.record(p, source)
    return {"targets": len(targets), "activities": n_act, "molecules": n_mol, "approved": n_approved}


def fetch_card(raw_dir: Path, manifest: Manifest, force=False, log=print) -> Path:
    log("CARD data is free for academic, government and non-profit use only; commercial use needs a license.")
    dest = raw_dir / "card" / "card-data.tar.bz2"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or force:
        dest.write_bytes(http_get(CARD_URL, timeout=600).content)
    _ensure_recorded(manifest, dest, CARD_URL)
    return dest


def fetch_all(keys, raw_dir, chembl=True, card=False, force=False, log=print) -> list[tuple[str, str]]:
    """Run every download step; a failing step is reported at the end and does not stop the others."""
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest = Manifest(raw_dir / "MANIFEST.json")
    targets = get_targets(keys)
    errors: list[tuple[str, str]] = []

    def run(step, fn):
        try:
            return fn()
        except Exception as e:  # keep going; report at the end
            errors.append((step, str(e)))
            log(f"  FAILED {step}: {e}")
            return None

    for t in targets:
        log(f"[{t.key}] UniProt {t.uniprot}")
        run(f"uniprot {t.uniprot}", lambda t=t: fetch_uniprot(t, raw_dir, manifest, force))
    for acc in sorted({t.genome_accession for t in targets}):
        log(f"[ncbi] coding sequences of {acc} (a few MB, can take a minute)")
        run(f"ncbi {acc}", lambda acc=acc: fetch_genome_cds(acc, raw_dir, manifest, force))
    for t in targets:
        log(f"[{t.key}] matching the coding sequence to the protein")
        m = run(f"cds {t.key}", lambda t=t: extract_target_cds(t, raw_dir, manifest))
        if m:
            log(f"    identity {m.identity:.4f}, {len(m.differences)} difference(s)")
        for pdb in t.pdb_ids:
            log(f"[{t.key}] structure {pdb}")
            run(f"pdb {pdb}", lambda pdb=pdb: fetch_pdb(pdb, raw_dir, manifest, force))
        if chembl:
            log(f"[{t.key}] ChEMBL activities (can take a few minutes)")
            s = run(f"chembl {t.key}", lambda t=t: fetch_chembl(t, raw_dir, manifest, force))
            if s:
                log(f"    {s}")
    if card:
        run("card", lambda: fetch_card(raw_dir, manifest, force, log))
    log("finished: " + (f"{len(errors)} step(s) failed" if errors else "all steps ok")
        + "; next run: python -m doublebind verify")
    return errors


# ---------- verification ----------

@dataclass
class Check:
    name: str
    status: str  # PASS | WARN | FAIL
    detail: str


def _pdb_title(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")[:20000]
    m = re.search(r"_struct\.title\s+(?:'([^']*)'|\"([^\"]*)\"|;\s*(.*?)\s*\n;|(\S+))", text, re.DOTALL)
    if not m:
        return "(title not parsed)"
    return " ".join(next(g for g in m.groups() if g).split())


def verify_all(keys, raw_dir) -> list[Check]:
    raw_dir = Path(raw_dir)
    checks: list[Check] = []
    manifest = Manifest(raw_dir / "MANIFEST.json")
    if not manifest.entries:
        checks.append(Check("manifest", "FAIL", "no MANIFEST.json yet: run the data command first"))
    else:
        results = manifest.verify()
        bad = [(p, s) for p, s in results if s != "ok"]
        checks.append(Check("manifest", "FAIL" if bad else "PASS",
                            f"problems: {bad}" if bad else f"{len(results)} files match their checksums"))
    for t in get_targets(keys):
        up = raw_dir / "uniprot" / f"{t.uniprot}.json"
        if not up.exists():
            checks.append(Check(f"{t.key}: UniProt", "FAIL", f"{up.name} missing"))
        else:
            d = json.loads(up.read_text(encoding="utf-8"))
            length = (d.get("sequence") or {}).get("length")
            gene = ((d.get("genes") or [{}])[0].get("geneName") or {}).get("value")
            organism = (d.get("organism") or {}).get("scientificName")
            name = (((d.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}).get("value")
            status = "PASS" if length == t.expected_length else "WARN"
            extra = "" if status == "PASS" else f" (expected {t.expected_length}: update targets.py if the record is right)"
            checks.append(Check(f"{t.key}: UniProt {t.uniprot}", status, f"{name} | gene {gene} | {organism} | {length} aa{extra}"))
        mj = raw_dir / "cds" / f"{t.key}_cds_match.json"
        if not mj.exists():
            checks.append(Check(f"{t.key}: coding sequence", "FAIL", "not extracted yet"))
        else:
            info = json.loads(mj.read_text(encoding="utf-8"))
            ident, diffs = info["identity"], info["differences"]
            if ident == 1.0:
                checks.append(Check(f"{t.key}: coding sequence", "PASS", f"translates exactly to {t.uniprot} ({info['header'][:70]})"))
            elif ident >= 0.99:
                checks.append(Check(f"{t.key}: coding sequence", "WARN", f"identity {ident:.4f}; differences (pos, UniProt, DNA): {diffs}"))
            else:
                checks.append(Check(f"{t.key}: coding sequence", "FAIL", f"identity only {ident:.4f}"))
        for pdb in t.pdb_ids:
            p = raw_dir / "pdb" / f"{pdb.upper()}.cif"
            checks.append(Check(f"{t.key}: structure {pdb}", "PASS" if p.exists() else "FAIL",
                                _pdb_title(p) if p.exists() else "missing: check the title matches what you expect"))
        acts = raw_dir / "chembl" / f"{t.key}_activities.jsonl"
        if acts.exists():
            n = sum(1 for _ in open(acts, encoding="utf-8"))
            checks.append(Check(f"{t.key}: ChEMBL activities", "PASS" if n else "WARN", f"{n} activity records"))
        else:
            checks.append(Check(f"{t.key}: ChEMBL activities", "WARN", "not downloaded (needed on Day 4)"))
    return checks
