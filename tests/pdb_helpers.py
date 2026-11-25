"""Build tiny synthetic PDB files so structure code can be tested without downloads."""
from doublebind.structure.proximity import AA3_TO_1

ONE_TO_3 = {v: k for k, v in AA3_TO_1.items() if k != "MSE"}


def atom_line(record, serial, name, resname, chain, resseq, xyz, element):
    x, y, z = xyz
    return (f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}")


def residue_atoms(i, aa):
    """Backbone plus CB for residue index i, laid out along the x axis (3.8 A apart)."""
    x = 3.8 * i
    atoms = [(" N  ", (x - 1.2, 0.5, 0.0), "N"), (" CA ", (x, 0.0, 0.0), "C"),
             (" C  ", (x + 1.2, 0.5, 0.0), "C"), (" O  ", (x + 1.2, 1.7, 0.0), "O")]
    if aa != "G":
        atoms.append((" CB ", (x, -1.0, 1.0), "C"))
    return atoms


def ligand_atoms(center_index, n=12, name="LIG"):
    """A compact cluster of heavy atoms 5 A above residue `center_index`."""
    xl = 3.8 * center_index
    coords = [(xl + 0.3 * (k % 3 - 1), 5.0 + 0.3 * (k // 3), 0.0) for k in range(n)]
    return [(f" C{k + 1:<2}", c, "C") for k, c in enumerate(coords)], name


def build_pdb(seq, first_resseq=1, chain="A", ligand_index=None, ligand_name="LIG", with_extras=True):
    """Return (pdb text, list of ligand coordinates, dict resseq -> list of atom coordinates)."""
    lines, serial = [], 1
    residue_xyz = {}
    for i, aa in enumerate(seq):
        resseq = first_resseq + i
        residue_xyz[resseq] = []
        for name, xyz, el in residue_atoms(i, aa):
            lines.append(atom_line("ATOM", serial, name, ONE_TO_3[aa], chain, resseq, xyz, el))
            residue_xyz[resseq].append(xyz)
            serial += 1
    lig_xyz = []
    if ligand_index is not None:
        atoms, _ = ligand_atoms(ligand_index)
        for name, xyz, el in atoms:
            lines.append(atom_line("HETATM", serial, name, ligand_name, chain, 900, xyz, el))
            lig_xyz.append(xyz)
            serial += 1
    if with_extras:
        lines.append(atom_line("HETATM", serial, " O  ", "HOH", chain, 901, (0.0, 10.0, 0.0), "O"))
        lines.append(atom_line("HETATM", serial + 1, "MG  ", "MG", chain, 902, (50.0, 50.0, 50.0), "MG"))
    lines.append("END")
    return "\n".join(lines) + "\n", lig_xyz, residue_xyz


def build_cif(seq, first_resseq=1, chain="A", ligand_index=None, ligand_name="LIG", title="SYNTHETIC TEST STRUCTURE"):
    """Same geometry as build_pdb, written as a minimal mmCIF file (the format RCSB serves)."""
    head = [
        "data_TEST", "#", f"_struct.title   '{title}'", "#", "loop_",
        "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
        "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
        "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv", "_atom_site.auth_seq_id", "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id", "_atom_site.auth_atom_id", "_atom_site.pdbx_PDB_model_num",
    ]
    rows, serial = [], 1

    def add(group, atom, resname, resseq, xyz, element, label_seq):
        nonlocal serial
        x, y, z = xyz
        rows.append(f"{group} {serial} {element} {atom} . {resname} {chain} 1 {label_seq} ? "
                    f"{x:.3f} {y:.3f} {z:.3f} 1.00 0.00 {resseq} {resname} {chain} {atom} 1")
        serial += 1

    for i, aa in enumerate(seq):
        for name, xyz, el in residue_atoms(i, aa):
            add("ATOM", name.strip(), ONE_TO_3[aa], first_resseq + i, xyz, el, i + 1)
    if ligand_index is not None:
        atoms, _ = ligand_atoms(ligand_index)
        for name, xyz, el in atoms:
            add("HETATM", name.strip(), ligand_name, 900, xyz, el, ".")
    add("HETATM", "O", "HOH", 901, (0.0, 10.0, 0.0), "O", ".")
    add("HETATM", "MG", "MG", 902, (50.0, 50.0, 50.0), "MG", ".")
    return "\n".join(head + rows + ["#"]) + "\n"
