"""Command line entry point: python -m doublebind <command>."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from doublebind.csi.csi import (ddg_bind_for_csi, ic50_fold_change, mean_needed_for_call, rt)
from doublebind.data.codons import (aggregate_missense, cds_to_protein, effect_counts, enumerate_variants,
                                    find_numbering_offsets, parse_mutation)
from doublebind.data.fetch import fetch_all, parse_fasta, verify_all
from doublebind.data.targets import TARGETS


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


def _load_known(path: Path, key: str) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return [row for row in csv.DictReader(f) if row["target"] == key]


def cmd_variants(args) -> int:
    key = args.target
    raw = Path(args.raw_dir)
    cds_path = Path(args.cds) if args.cds else raw / "cds" / f"{key}_cds.fasta"
    records = parse_fasta(cds_path.read_text(encoding="utf-8"))
    if not records:
        print(f"no sequence found in {cds_path}")
        return 1
    cds = records[0][1]
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


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="doublebind", description="DoubleBind research prototype")
    sub = p.add_subparsers(dest="cmd", required=True)
    default_targets = ",".join(TARGETS)

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

    m = sub.add_parser("variants", help="list every amino-acid change one nucleotide can reach (Day 2)")
    m.add_argument("target", choices=sorted(TARGETS))
    m.add_argument("--kappa", type=float, default=2.0, help="transition/transversion rate ratio (placeholder 2.0)")
    m.add_argument("--raw-dir", default="data/raw")
    m.add_argument("--cds", default=None, help="path to a CDS FASTA (default: data/raw/cds/<target>_cds.fasta)")
    m.add_argument("--out-dir", default="data/processed")
    m.add_argument("--known", default="data/known_mutations.csv")
    m.set_defaults(func=cmd_variants)

    n = sub.add_parser("numbers", help="print the PRD energy-to-susceptibility numbers")
    n.add_argument("--temp", type=float, default=310.15, help="temperature in kelvin")
    n.set_defaults(func=cmd_numbers)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)
