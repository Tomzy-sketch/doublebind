"""MVP targets. Every value below is a starting point that `python -m doublebind verify` checks.

expected_length and the PDB choices were written from memory: a mismatch is reported as a
warning, so you can correct this file instead of trusting it.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Target:
    key: str
    protein: str
    gene: str
    organism: str
    uniprot: str
    expected_length: int
    genome_accession: str   # nucleotide record that contains the coding sequence
    pdb_ids: tuple
    notes: str = ""


TARGETS: dict[str, Target] = {
    "gyra": Target(
        key="gyra", protein="DNA gyrase subunit A", gene="gyrA",
        organism="Escherichia coli K-12", uniprot="P0AES4", expected_length=875,
        genome_accession="U00096.3", pdb_ids=("2XCT",),
        notes="2XCT is believed to be an S. aureus gyrase-DNA-ciprofloxacin complex: a homolog, "
              "not E. coli. Search RCSB for an E. coli complex and replace it here if you find one. "
              "If you keep the homolog, residue numbers need an alignment to E. coli numbering."),
    "gyrb": Target(
        key="gyrb", protein="DNA gyrase subunit B", gene="gyrB",
        organism="Escherichia coli K-12", uniprot="P0AES6", expected_length=804,
        genome_accession="U00096.3", pdb_ids=(), notes="No structure needed for the MVP."),
    "tem1": Target(
        key="tem1", protein="Beta-lactamase TEM (TEM-1)", gene="bla",
        organism="Escherichia coli plasmid (pBR322 / Tn3 lineage)", uniprot="P62593", expected_length=286,
        genome_accession="J01749.1", pdb_ids=("1BTL",),
        notes="Literature numbering for TEM-1 is Ambler numbering, which differs from UniProt numbering "
              "by a constant offset for most residues. The variants command finds that offset for you."),
}


def get_targets(keys):
    unknown = [k for k in keys if k not in TARGETS]
    if unknown:
        raise KeyError(f"unknown target(s) {unknown}; choose from {sorted(TARGETS)}")
    return [TARGETS[k] for k in keys]
