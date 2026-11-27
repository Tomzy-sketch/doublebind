"""The Day 2 variant table: one row per reachable amino-acid variant (grows a column per day)."""
from __future__ import annotations

import csv

COLUMNS = ["pos", "wt_aa", "mut_aa", "label", "weight", "n_routes", "esm_llr", "ddg_fold", "dist_min", "in_pocket"]


def build_variant_table(aa_variants, esm_rows=None, pocket_rows=None, cutoff: float = 8.0, ddg_rows=None) -> list[dict]:
    ddg = {(int(r["pos"]), r["mut_aa"]): float(r["ddg_fold"]) for r in (ddg_rows or [])}
    llr = {(int(r["pos"]), r["mut_aa"]): float(r["llr"]) for r in (esm_rows or [])}
    dist = {int(r["pos"]): float(r["min_dist"]) for r in (pocket_rows or []) if r.get("min_dist") not in (None, "")}
    rows = []
    for v in aa_variants:
        d = dist.get(v.pos)
        rows.append({
            "pos": v.pos, "wt_aa": v.wt_aa, "mut_aa": v.mut_aa, "label": v.label,
            "weight": round(v.weight, 6), "n_routes": v.n_routes,
            "esm_llr": round(llr[(v.pos, v.mut_aa)], 4) if (v.pos, v.mut_aa) in llr else "",
            "ddg_fold": round(ddg[(v.pos, v.mut_aa)], 4) if (v.pos, v.mut_aa) in ddg else "",
            "dist_min": "" if d is None else round(d, 2),
            "in_pocket": "" if d is None else int(d <= cutoff),
        })
    return rows


def write_table(path, rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
        