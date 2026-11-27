"""Write a clean, ThermoMPNN-ready PDB file from any structure (Day 3). Needs Biopython.

ThermoMPNN labels every prediction with 'position', the 0-based index of the residue in the chain it was given.
Here the chain is written with residues renumbered 1..N, so index = line order, and a sidecar file maps each
index back to the original residue. That removes every numbering pitfall (insertion codes, gaps, Ambler numbering).
Only the backbone (N, CA, C, O) is written, which is all ProteinMPNN-based models read.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from doublebind.structure.proximity import AA3_TO_1, first_model, load_structure, protein_chains

BACKBONE = ("N", "CA", "C", "O")
AA1_TO_3 = {v: k for k, v in AA3_TO_1.items() if k != "MSE"}


@dataclass
class PreparedChain:
    chain_id: str
    sequence: str
    index_map: list     # one dict per written residue: index, chain, resseq, icode, resname, aa
    n_skipped: int      # residues left out because a backbone atom was missing
    pdb_path: Path
    map_path: Path


def _atom_line(serial: int, name: str, resname: str, resseq: int, xyz, element: str) -> str:
    x, y, z = xyz
    return (f"ATOM  {serial:>5}  {name:<3} {resname:>3} A{resseq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}")


def write_clean_pdb(structure_path, out_pdb, chain_id=None, min_length: int = 30) -> PreparedChain:
    model = first_model(load_structure(structure_path))
    chains = [c for c in protein_chains(model) if len(c.sequence) >= min_length]
    if chain_id:
        chains = [c for c in chains if c.chain_id == chain_id]
    if not chains:
        raise ValueError(f"no protein chain of at least {min_length} residues"
                         + (f" with id {chain_id}" if chain_id else "") + f" in {structure_path}")
    chain = chains[0]
    lines, index_map, serial, skipped = [], [], 1, 0
    for res in chain.residues:
        atoms = {}
        for a in res:
            name = a.get_id().strip()
            if name in BACKBONE and name not in atoms:
                atoms[name] = a
        if any(k not in atoms for k in BACKBONE):
            skipped += 1
            continue
        original = res.get_resname().strip().upper()
        aa = AA3_TO_1[original]
        new_number = len(index_map) + 1
        for name in BACKBONE:
            lines.append(_atom_line(serial, name, AA1_TO_3[aa], new_number, atoms[name].coord, name[0]))
            serial += 1
        index_map.append({"index": new_number - 1, "chain": chain.chain_id, "resseq": int(res.id[1]),
                          "icode": str(res.id[2]).strip(), "resname": original, "aa": aa})
    lines += ["TER", "END"]
    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")
    map_path = out_pdb.with_suffix(".map.json")
    sequence = "".join(e["aa"] for e in index_map)
    map_path.write_text(json.dumps({"source": str(structure_path), "chain": chain.chain_id,
                                    "sequence": sequence, "index_map": index_map}, indent=1), encoding="utf-8")
    return PreparedChain(chain.chain_id, sequence, index_map, skipped, out_pdb, map_path)
