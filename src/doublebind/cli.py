"""Command line entry point: python -m doublebind <command>.

Output is plain ASCII on purpose, so it prints correctly in any Windows console.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

from doublebind.csi.csi import ddg_bind_for_csi, ic50_fold_change, mean_needed_for_call, rt
from doublebind.data.codons import (aggregate_missense, cds_to_protein, effect_counts, enumerate_variants,
                                    find_numbering_offsets, parse_mutation)
from doublebind.data.fetch import fetch_all, parse_fasta, pdb_title, verify_all
from doublebind.data.targets import TARGETS


# ---------- shared helpers ----------

def _load_cds(key, raw_dir, cds=None) -> str:
    path = Path(cds) if cds else Path(raw_dir) / "cds" / f"{key}_cds.fasta"
    records = parse_fasta(path.read_text(encoding="utf-8"))
    if not records:
        raise SystemExit(f"no sequence found in {path}")
    return records[0][1]


def _load_known(path: Path, key: str) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return [row for row in csv.DictReader(f) if row["target"] == key]


def _structure_path(raw_dir, pdb: str) -> Path:
    folder = Path(raw_dir) / "pdb"
    for ext in (".cif", ".pdb"):
        p = folder / f"{pdb.upper()}{ext}"
        if p.exists():
            return p
    return folder / f"{pdb.upper()}.cif"


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------- Day 1 ----------

def cmd_data(args) -> int:
    keys = [k.strip() for k in args.targets.split(",") if k.strip()]
    errors = fetch_all(keys, args.raw_dir, chembl=not args.no_chembl, card=args.card, force=args.force)
    return 1 if errors else 0


def cmd_verify(args) -> int:
    keys = [k.strip() for k in args.targets.split(",") if k.strip()]
    checks = verify_all(keys, args.raw_dir)
    for c in checks:
        print(f"{c.status:<5} {c.name:<34} {c.detail}")
    fails = sum(c.status == "FAIL" for c in checks)
    warns = sum(c.status == "WARN" for c in checks)
    print(f"\n{len(checks) - fails - warns} passed, {warns} warning(s), {fails} failed")
    return 1 if fails else 0


def cmd_variants(args) -> int:
    key = args.target
    cds = _load_cds(key, args.raw_dir, args.cds)
    protein = cds_to_protein(cds)
    variants = enumerate_variants(cds, kappa=args.kappa)
    aa = aggregate_missense(variants)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    nt_csv, aa_csv = out / f"{key}_variants_nt.csv", out / f"{key}_variants_aa.csv"
    with open(nt_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pos", "wt_codon", "mut_codon", "codon_site", "ref_nt", "alt_nt", "wt_aa", "mut_aa",
                    "effect", "transition", "weight"])
        for v in variants:
            w.writerow([v.pos, v.wt_codon, v.mut_codon, v.codon_site, v.ref_nt, v.alt_nt, v.wt_aa, v.mut_aa,
                        v.effect, int(v.transition), f"{v.weight:.6f}"])
    with open(aa_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pos", "wt_aa", "mut_aa", "label", "weight", "n_routes"])
        for v in aa:
            w.writerow([v.pos, v.wt_aa, v.mut_aa, v.label, f"{v.weight:.6f}", v.n_routes])

    counts = effect_counts(variants)
    print(f"{key}: {len(protein)} residues, kappa={args.kappa}")
    print(f"  single-nucleotide changes: {len(variants)} "
          f"(missense {counts['missense']}, synonymous {counts['synonymous']}, nonsense {counts['nonsense']})")
    print(f"  distinct amino-acid variants reachable by one nucleotide change: {len(aa)}")
    print(f"  wrote {nt_csv} and {aa_csv}")

    known = _load_known(Path(args.known), key)
    if known:
        labels = [row["mutation"] for row in known]
        offsets = find_numbering_offsets(protein, labels)
        print("  known resistance mutations (data/known_mutations.csv):")
        if not offsets:
            print("    no single numbering offset makes every wild-type residue line up: check the sequence "
                  "or the labels before trusting any of these")
        else:
            off = offsets[0]
            print(f"    numbering offset {off:+d} (sequence position = literature position {off:+d})")
            index = {(v.pos, v.mut_aa): v for v in aa}
            for row in known:
                wt, pos, mut = parse_mutation(row["mutation"])
                v = index.get((pos + off, mut))
                if v:
                    print(f"    {row['mutation']:>6} = residue {pos + off}: reachable by one nucleotide change "
                          f"({v.n_routes} route(s), supply weight {v.weight:.3f})")
                else:
                    print(f"    {row['mutation']:>6} = residue {pos + off}: NOT reachable by a single nucleotide change")
    return 0


def cmd_numbers(args) -> int:
    t = args.temp
    print(f"RT at {t} K = {rt(t):.3f} kcal/mol")
    print(f"CSI of 2.0 (4x more susceptible) needs ddG_bind <= {ddg_bind_for_csi(2.0, t):.3f} kcal/mol")
    print(f"ddG_bind of +1.0 kcal/mol = {ic50_fold_change(1.0, t):.1f}x higher IC50")
    print("To call collateral sensitivity with 80% confidence, the predicted mean ddG_bind must be at or below:")
    for sigma in (0.6, 1.0, 1.2, 1.5):
        mu = mean_needed_for_call(sigma, 2.0, 0.8, t)
        fold = ic50_fold_change(mu, t)
        print(f"  model error sigma = {sigma:.1f} kcal/mol -> {mu:+.2f} kcal/mol ({1 / fold:.1f}x more susceptible)")
    return 0


# ---------- Day 2 ----------

def cmd_inspect(args) -> int:
    from doublebind.structure import proximity as prox

    path = _structure_path(args.raw_dir, args.pdb)
    if not path.exists():
        print(f"{path} not found: run the data command first")
        return 1
    model = prox.first_model(prox.load_structure(path))
    print(f"{args.pdb.upper()}: {pdb_title(path)}")
    target = cds_to_protein(_load_cds(args.target, args.raw_dir)) if args.target else None
    print("protein chains:")
    for chain in prox.protein_chains(model):
        line = (f"  {chain.chain_id:<3}{len(chain.sequence):>5} residues "
                f"(numbers {chain.residues[0].id[1]} to {chain.residues[-1].id[1]})")
        if target:
            m = prox.map_chain_to_target(chain, target)
            line += (f"   vs {args.target}: identity {m.identity:.2f}, covers {100 * m.coverage:.0f}% "
                     f"of the target ({m.aligned} residues aligned)")
        print(line)
    groups = prox.hetero_groups(model)
    candidates = {id(g) for g in prox.candidate_ligands(groups)}
    print("other molecules (water hidden):")
    for g in groups:
        flag = "<- candidate ligand" if id(g) in candidates else "(ion, solvent or small)"
        print(f"  {g.name:<5} chain {g.chain_id} number {g.resseq:<6}{g.n_heavy:>4} heavy atoms  {flag}")
    return 0


def cmd_pocket(args) -> int:
    from doublebind.structure import proximity as prox

    key = args.target
    protein = cds_to_protein(_load_cds(key, args.raw_dir, args.cds))
    pdb = args.pdb or (TARGETS[key].pdb_ids[0] if TARGETS[key].pdb_ids else None)
    if not pdb:
        print(f"no structure configured for {key}: pass --pdb")
        return 1
    path = _structure_path(args.raw_dir, pdb)
    if not path.exists():
        print(f"{path} not found: run the data command first")
        return 1
    model = prox.first_model(prox.load_structure(path))
    groups = prox.hetero_groups(model)
    if args.ligand:
        names = {n.strip().upper() for n in args.ligand.split(",")}
        ligands = [g for g in groups if g.name in names]
    else:
        ligands = prox.candidate_ligands(groups)
    if not ligands:
        print("no ligand found. Molecules in this structure:")
        for g in groups:
            print(f"  {g.name} chain {g.chain_id} number {g.resseq} ({g.n_heavy} heavy atoms)")
        print("Pass the right one with --ligand NAME")
        return 1
    chains = {c.strip() for c in args.chains.split(",")} if args.chains else None
    dists, maps = prox.ligand_distances(model, protein, ligands, chains=chains)
    if not dists:
        print("no chain of this structure aligns to the target protein well enough (identity >= 0.3, 50+ residues)")
        return 1

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    out_csv = out / f"{key}_pocket.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pos", "wt_aa", "struct_aa", "chain", "resseq", "min_dist", "in_pocket"])
        for d in dists:
            w.writerow([d.pos, d.wt_aa, d.struct_aa, d.chain_id, d.resseq, f"{d.min_dist:.2f}", int(d.min_dist <= args.cutoff)])

    used = ", ".join(sorted({f"{g.name} (chain {g.chain_id})" for g in ligands}))
    print(f"{key} pocket from {pdb.upper()}: ligand(s) {used}")
    for m in maps:
        print(f"  chain {m.chain_id}: identity {m.identity:.2f}, covers {100 * m.coverage:.0f}% of the target")
    inside = [d for d in dists if d.min_dist <= args.cutoff]
    print(f"  {len(dists)} of {len(protein)} positions have a structure; {len(inside)} are within {args.cutoff:g} A of the ligand")
    print("  nearest positions (target residue, distance, residue found in the structure):")
    for d in sorted(dists, key=lambda x: x.min_dist)[:15]:
        print(f"    {d.wt_aa}{d.pos:<5} {d.min_dist:5.2f} A   structure {d.struct_aa}{d.resseq} chain {d.chain_id}")

    known = _load_known(Path(args.known), key)
    if known:
        offsets = find_numbering_offsets(protein, [k["mutation"] for k in known])
        if offsets:
            order = {d.pos: i + 1 for i, d in enumerate(sorted(dists, key=lambda x: x.min_dist))}
            by_pos = {d.pos: d for d in dists}
            for k in known:
                pos = parse_mutation(k["mutation"])[1] + offsets[0]
                d = by_pos.get(pos)
                if d:
                    print(f"  known {k['mutation']}: {d.min_dist:.2f} A, closest-first rank {order[pos]} of {len(dists)}")
                else:
                    print(f"  known {k['mutation']}: no structure coverage at position {pos}")
    print(f"  wrote {out_csv}")
    return 0


def cmd_esm(args) -> int:
    from doublebind.scoring.esm_llr import score_protein

    key = args.target
    cds = _load_cds(key, args.raw_dir, args.cds)
    protein = cds_to_protein(cds)
    reachable = {(v.pos, v.mut_aa) for v in aggregate_missense(enumerate_variants(cds))}
    tag = args.model.split("/")[-1]
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    def progress(done, total):
        print(f"\r  scored {done}/{total} positions", end="", flush=True)

    print(f"{key}: {len(protein)} residues, model {args.model}, mode {args.mode}")
    start = time.time()
    subs, stats, device = score_protein(protein, args.model, mode=args.mode, device=args.device,
                                        batch_size=args.batch_size, progress=progress)
    print(f"\n  device {device}, {time.time() - start:.0f} s")

    llr_csv = out / f"{key}_esm_{tag}.csv"
    with open(llr_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pos", "wt_aa", "mut_aa", "llr", "reachable", "mode"])
        for pos, wt, mut, llr in subs:
            w.writerow([pos, wt, mut, f"{llr:.4f}", int((pos, mut) in reachable), args.mode])
    pos_csv = out / f"{key}_esm_{tag}_positions.csv"
    with open(pos_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pos", "wt_aa", "wt_logp", "entropy"])
        for pos, wt, logp, entropy in stats:
            w.writerow([pos, wt, f"{logp:.4f}", f"{entropy:.4f}"])

    reach_vals = [llr for pos, _wt, mut, llr in subs if (pos, mut) in reachable]
    if reach_vals:
        print(f"  {len(subs)} substitutions scored, {len(reach_vals)} reachable by one nucleotide change "
              f"(mean LLR {sum(reach_vals) / len(reach_vals):+.2f})")
    known = _load_known(Path(args.known), key)
    if known:
        offsets = find_numbering_offsets(protein, [k["mutation"] for k in known])
        if offsets:
            lookup = {(pos, mut): llr for pos, _wt, mut, llr in subs}
            for k in known:
                _wt, p, m = parse_mutation(k["mutation"])
                val = lookup.get((p + offsets[0], m))
                print(f"  known {k['mutation']}: LLR {'n/a' if val is None else format(val, '+.2f')}")
    print(f"  wrote {llr_csv} and {pos_csv}")
    return 0


def cmd_baselines(args) -> int:
    from doublebind.baselines import find_esm_file, format_report, read_rows, run_baselines

    key = args.target
    cds = _load_cds(key, args.raw_dir, args.cds)
    protein = cds_to_protein(cds)
    aa = aggregate_missense(enumerate_variants(cds))
    processed = Path(args.out_dir)
    esm_file = find_esm_file(processed, key, args.esm_file)
    pocket_file = Path(args.pocket_file) if args.pocket_file else processed / f"{key}_pocket.csv"
    esm_rows = read_rows(esm_file) if esm_file and Path(esm_file).exists() else None
    pocket_rows = read_rows(pocket_file) if pocket_file.exists() else None
    known = _load_known(Path(args.known), key)
    print(f"ESM scores: {esm_file if esm_rows else 'none found (run the esm command)'}")
    print(f"pocket distances: {pocket_file if pocket_rows else 'none found (run the pocket command)'}\n")
    result = run_baselines(protein, aa, known, esm_rows, pocket_rows, cutoff=args.cutoff)
    print(format_report(key, result))
    out_json = processed / f"{key}_baselines.json"
    out_json.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote {out_json}")
    return 0


def cmd_table(args) -> int:
    from doublebind.baselines import find_esm_file, read_rows
    from doublebind.tables import build_variant_table, write_table

    key = args.target
    cds = _load_cds(key, args.raw_dir, args.cds)
    aa = aggregate_missense(enumerate_variants(cds))
    processed = Path(args.out_dir)
    esm_file = find_esm_file(processed, key, args.esm_file)
    pocket_file = processed / f"{key}_pocket.csv"
    esm_rows = read_rows(esm_file) if esm_file and Path(esm_file).exists() else None
    pocket_rows = read_rows(pocket_file) if pocket_file.exists() else None
    rows = build_variant_table(aa, esm_rows, pocket_rows, cutoff=args.cutoff)
    out = processed / f"{key}_variant_table.csv"
    write_table(out, rows)
    with_esm = sum(1 for r in rows if r["esm_llr"] != "")
    with_pocket = sum(1 for r in rows if r["dist_min"] != "")
    print(f"{key}: {len(rows)} variants; {with_esm} with an ESM score, {with_pocket} with a pocket distance")
    print(f"wrote {out}")
    return 0


# ---------- parser ----------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="doublebind", description="DoubleBind research prototype")
    sub = p.add_subparsers(dest="cmd", required=True)
    default_targets = ",".join(TARGETS)

    def common(sp, cds=True):
        sp.add_argument("--raw-dir", default="data/raw")
        sp.add_argument("--out-dir", default="data/processed")
        sp.add_argument("--known", default="data/known_mutations.csv")
        if cds:
            sp.add_argument("--cds", default=None, help="path to a CDS FASTA (default: data/raw/cds/<target>_cds.fasta)")

    d = sub.add_parser("data", help="download raw data (Day 1)")
    d.add_argument("--targets", default=default_targets)
    d.add_argument("--raw-dir", default="data/raw")
    d.add_argument("--no-chembl", action="store_true", help="skip ChEMBL (needed on Day 4, can wait)")
    d.add_argument("--card", action="store_true", help="also download CARD (academic use only)")
    d.add_argument("--force", action="store_true", help="download again even if files exist")
    d.set_defaults(func=cmd_data)

    v = sub.add_parser("verify", help="check checksums, accessions and coding sequences")
    v.add_argument("--targets", default=default_targets)
    v.add_argument("--raw-dir", default="data/raw")
    v.set_defaults(func=cmd_verify)

    m = sub.add_parser("variants", help="list every amino-acid change one nucleotide can reach")
    m.add_argument("target", choices=sorted(TARGETS))
    m.add_argument("--kappa", type=float, default=2.0, help="transition/transversion rate ratio (placeholder 2.0)")
    common(m)
    m.set_defaults(func=cmd_variants)

    n = sub.add_parser("numbers", help="print the PRD energy-to-susceptibility numbers")
    n.add_argument("--temp", type=float, default=310.15, help="temperature in kelvin")
    n.set_defaults(func=cmd_numbers)

    i = sub.add_parser("inspect", help="show the chains and ligands of a downloaded structure (Day 2)")
    i.add_argument("pdb")
    i.add_argument("--target", choices=sorted(TARGETS), default=None, help="also align each chain to this target")
    i.add_argument("--raw-dir", default="data/raw")
    i.set_defaults(func=cmd_inspect)

    pk = sub.add_parser("pocket", help="distance of every residue to the bound drug (Day 2)")
    pk.add_argument("target", choices=sorted(TARGETS))
    pk.add_argument("--pdb", default=None, help="structure ID (default: the first one in targets.py)")
    pk.add_argument("--ligand", default=None, help="ligand name(s), comma separated (default: auto-detect)")
    pk.add_argument("--chains", default=None, help="chains to use, comma separated (default: all that align)")
    pk.add_argument("--cutoff", type=float, default=8.0, help="pocket radius in angstrom")
    common(pk)
    pk.set_defaults(func=cmd_pocket)

    e = sub.add_parser("esm", help="ESM-2 plausibility score for every substitution (Day 2)")
    e.add_argument("target", choices=sorted(TARGETS))
    e.add_argument("--model", default="facebook/esm2_t30_150M_UR50D")
    e.add_argument("--mode", choices=["masked", "wt"], default="masked")
    e.add_argument("--device", default="auto", help="auto, cpu or cuda")
    e.add_argument("--batch-size", type=int, default=None)
    common(e)
    e.set_defaults(func=cmd_esm)

    b = sub.add_parser("baselines", help="rank the known resistance mutations under each baseline (Day 2)")
    b.add_argument("target", choices=sorted(TARGETS))
    b.add_argument("--esm-file", default=None)
    b.add_argument("--pocket-file", default=None)
    b.add_argument("--cutoff", type=float, default=8.0)
    common(b)
    b.set_defaults(func=cmd_baselines)

    t = sub.add_parser("table", help="merge variants, ESM scores and pocket distances into one table (Day 2)")
    t.add_argument("target", choices=sorted(TARGETS))
    t.add_argument("--esm-file", default=None)
    t.add_argument("--cutoff", type=float, default=8.0)
    common(t)
    t.set_defaults(func=cmd_table)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)
