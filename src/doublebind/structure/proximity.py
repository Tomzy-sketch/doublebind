"""Ligand proximity from an experimental structure (Day 2). Needs Biopython and numpy.

Definition fixed before looking at any results: the distance of a residue is the smallest
distance between any heavy atom of that residue and any heavy atom of the chosen ligand(s).
A structure from a homolog is mapped onto the target protein by a local alignment, and a
target position takes the smallest distance over all chains that map to it.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
}
WATER = {"HOH", "WAT", "DOD", "H2O"}
# Ions, solvents, buffers and other common crystallization additives are never "the drug".
NOT_A_DRUG = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "MPD", "ACT", "ACY", "FMT", "DMS", "TRS", "MES",
    "EPE", "BME", "DTT", "IMD", "CIT", "TAR", "MLI", "BCT", "EOH", "MOH", "IPA", "NO3", "SCN", "AZI",
    "CL", "NA", "K", "MG", "MN", "ZN", "CA", "CD", "IOD", "BR", "NI", "CO", "FE", "FE2", "CU", "HG",
    "SR", "CS", "RB", "LI", "BA", "F",
}


def _is_heavy(atom) -> bool:
    element = (atom.element or "").strip().upper()
    if not element:
        element = atom.get_id().strip()[:1].upper()
    return element not in ("H", "D")


def load_structure(path):
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.PDBParser import PDBParser

    p = Path(path)
    parser = MMCIFParser(QUIET=True) if p.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    return parser.get_structure(p.stem, str(p))


def first_model(structure):
    for model in structure.get_models():
        return model
    raise ValueError("the structure file contains no atoms")


@dataclass
class ChainInfo:
    chain_id: str
    sequence: str
    residues: list


@dataclass
class HeteroGroup:
    chain_id: str
    name: str
    resseq: int
    n_heavy: int
    residue: object


@dataclass
class ChainMapping:
    chain_id: str
    identity: float    # identical residues / aligned residues
    coverage: float    # aligned residues / target length
    aligned: int
    mapping: dict      # target position (1-based) -> Biopython residue


@dataclass
class PositionDistance:
    pos: int
    wt_aa: str         # residue in the target protein
    struct_aa: str     # residue found at the aligned place in the structure
    chain_id: str
    resseq: int
    min_dist: float


def protein_chains(model, min_length: int = 10) -> list[ChainInfo]:
    chains = []
    for chain in model:
        seq, residues = [], []
        for r in chain:
            name = r.get_resname().strip().upper()
            if name in AA3_TO_1 and (r.id[0] == " " or name == "MSE"):
                seq.append(AA3_TO_1[name])
                residues.append(r)
        if len(seq) >= min_length:
            chains.append(ChainInfo(chain.id, "".join(seq), residues))
    return chains


def hetero_groups(model) -> list[HeteroGroup]:
    groups = []
    for chain in model:
        for r in chain:
            if r.id[0] == " ":
                continue
            name = r.get_resname().strip().upper()
            if name in WATER or name in AA3_TO_1:
                continue
            groups.append(HeteroGroup(chain.id, name, r.id[1], sum(_is_heavy(a) for a in r), r))
    return groups


def candidate_ligands(groups: list[HeteroGroup], min_heavy: int = 12) -> list[HeteroGroup]:
    return [g for g in groups if g.name not in NOT_A_DRUG and g.n_heavy >= min_heavy]


def map_chain_to_target(chain: ChainInfo, target: str) -> ChainMapping:
    """Local alignment of the target protein to the chain's observed residues."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = next(iter(aligner.align(target, chain.sequence)))
    target_blocks, chain_blocks = alignment.aligned
    mapping, same = {}, 0
    for (t0, t1), (c0, _c1) in zip(target_blocks, chain_blocks):
        for k in range(int(t1 - t0)):
            t_i, c_i = int(t0) + k, int(c0) + k
            mapping[t_i + 1] = chain.residues[c_i]
            same += int(target[t_i] == chain.sequence[c_i])
    aligned = len(mapping)
    return ChainMapping(chain.chain_id, same / aligned if aligned else 0.0, aligned / len(target), aligned, mapping)


def ligand_distances(model, target: str, ligands: list[HeteroGroup], chains=None,
                     min_identity: float = 0.3, min_aligned: int = 50):
    """Return (list of PositionDistance sorted by position, list of ChainMapping that were used)."""
    ligand_xyz = np.array([a.coord for g in ligands for a in g.residue if _is_heavy(a)], dtype=float)
    if ligand_xyz.size == 0:
        raise ValueError("the chosen ligand has no heavy atoms")
    mappings = []
    for chain in protein_chains(model):
        if chains is not None and chain.chain_id not in chains:
            continue
        m = map_chain_to_target(chain, target)
        if m.identity >= min_identity and m.aligned >= min_aligned:
            mappings.append(m)
    best: dict[int, PositionDistance] = {}
    for m in mappings:
        for pos, res in m.mapping.items():
            xyz = np.array([a.coord for a in res if _is_heavy(a)], dtype=float)
            if xyz.size == 0:
                continue
            d = float(np.min(np.linalg.norm(xyz[:, None, :] - ligand_xyz[None, :, :], axis=2)))
            if pos not in best or d < best[pos].min_dist:
                struct_aa = AA3_TO_1.get(res.get_resname().strip().upper(), "X")
                best[pos] = PositionDistance(pos, target[pos - 1], struct_aa, m.chain_id, int(res.id[1]), d)
    return [best[p] for p in sorted(best)], mappings
