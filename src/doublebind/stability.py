"""ThermoMPNN stability scores (Day 3): read, map onto the target protein, convert the sign, check.

Convention: ddg_fold > 0 means destabilizing (see signs.py). ThermoMPNN already uses it, so the conversion is
the identity, but it still goes through to_doublebind() so a different tool cannot slip in unconverted.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from statistics import mean

from doublebind.data.codons import find_numbering_offsets, parse_mutation
from doublebind.signs import check_directions, to_doublebind
from doublebind.structure.proximity import ChainInfo, map_chain_to_target

AA20 = "ACDEFGHIKLMNPQRSTVWY"


def read_thermompnn_csv(path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    needed = {"ddG_pred", "position", "wildtype", "mutation"}
    if not rows or not needed <= set(rows[0]):
        found = sorted(rows[0]) if rows else "an empty file"
        raise ValueError(f"{path} does not look like ThermoMPNN output (needs {sorted(needed)}, found {found})")
    return rows


@dataclass
class StabilityImport:
    rows: list = field(default_factory=list)   # dicts: pos, wt_aa, mut_aa, ddg_fold
    n_input: int = 0
    n_structure_mismatch: int = 0              # CSV wild type disagrees with the structure's residue: should be 0
    n_unmapped: int = 0                        # structure residue has no counterpart in the target
    n_target_mismatch: int = 0                 # structure residue differs from the target residue (homolog)
    identity: float = 0.0
    coverage_of_target: float = 0.0


def import_thermompnn(csv_rows, index_map, target_protein, source="thermompnn", allow_mismatch=False) -> StabilityImport:
    sequence = "".join(e["aa"] for e in index_map)
    mapping = map_chain_to_target(ChainInfo("A", sequence, list(range(len(sequence)))), target_protein)
    index_to_target = {idx: tpos for tpos, idx in mapping.mapping.items()}
    out = StabilityImport(n_input=len(csv_rows), identity=mapping.identity, coverage_of_target=mapping.coverage)
    seen = {}
    for r in csv_rows:
        mut = r["mutation"].strip()
        wt = r["wildtype"].strip()
        if mut not in AA20 or mut == wt:
            continue  # 'X', or the wild type itself
        p = int(float(r["position"]))
        if not 0 <= p < len(index_map) or index_map[p]["aa"] != wt:
            out.n_structure_mismatch += 1
            continue
        tpos = index_to_target.get(p)
        if tpos is None:
            out.n_unmapped += 1
            continue
        target_wt = target_protein[tpos - 1]
        if target_wt != wt:
            out.n_target_mismatch += 1
            if not allow_mismatch:
                continue
        seen[(tpos, mut)] = {"pos": tpos, "wt_aa": target_wt, "mut_aa": mut,
                             "ddg_fold": float(to_doublebind(float(r["ddG_pred"]), source))}
    out.rows = [seen[k] for k in sorted(seen)]
    return out


def coverage(ddg_rows, aa_variants, region=None) -> dict:
    """Share of reachable variants that have a ddG, overall and inside region=(first, last) of the target."""
    have = {(r["pos"], r["mut_aa"]) for r in ddg_rows}
    reach = {(v.pos, v.mut_aa) for v in aa_variants}
    res = {"all": [len(reach & have), len(reach)]}
    if region:
        lo, hi = region
        inside = {k for k in reach if lo <= k[0] <= hi}
        res["region"] = [len(inside & have), len(inside)]
        res["region_bounds"] = [lo, hi]
    return res


def proline_canary(ddg_rows, margin: float = 0.2) -> dict:
    """Substitutions to proline are destabilizing on average; if they do not score higher than the rest,
    the sign (or the mapping) is probably wrong."""
    to_p = [r["ddg_fold"] for r in ddg_rows if r["mut_aa"] == "P"]
    rest = [r["ddg_fold"] for r in ddg_rows if r["mut_aa"] != "P"]
    if not to_p or not rest:
        return {"passed": None, "difference": None, "note": "not enough substitutions to test"}
    diff = mean(to_p) - mean(rest)
    return {"passed": diff > margin, "difference": diff, "mean_to_proline": mean(to_p), "mean_other": mean(rest)}


def direction_report(ddg_rows, protein, expected_rows) -> dict:
    """Compare predicted signs with literature directions (rows of data/known_directions.csv for one target)."""
    if not expected_rows:
        return {"checked": 0}
    labels = [r["mutation"] for r in expected_rows]
    offsets = find_numbering_offsets(protein, labels)
    if not offsets:
        return {"checked": 0, "problem": "no numbering offset lines the labels up with the sequence"}
    off = offsets[0]
    lookup = {(r["pos"], r["mut_aa"]): r["ddg_fold"] for r in ddg_rows}
    preds, expected = {}, {}
    for r in expected_rows:
        wt, pos, mut = parse_mutation(r["mutation"])
        value = lookup.get((pos + off, mut))
        if value is not None:
            preds[r["mutation"]] = value
        expected[r["mutation"]] = r["direction"]
    rep = check_directions(preds, expected)
    return {"checked": rep.n_total, "agree": rep.n_agree, "offset": off, "rows": [list(x) for x in rep.rows]}


def format_summary(key, csv_name, imp: StabilityImport, cov: dict, canary: dict, directions: dict) -> str:
    lines = [f"{key}: imported {len(imp.rows)} substitutions from {csv_name} ({imp.n_input} rows read)",
             f"  structure aligned to the target: identity {imp.identity:.2f}, covers {100 * imp.coverage_of_target:.0f}% of the target",
             f"  skipped: {imp.n_structure_mismatch} rows whose wild type disagrees with the structure (should be 0), "
             f"{imp.n_unmapped} with no target position, {imp.n_target_mismatch} where the structure residue differs from the target"]
    a, b = cov["all"]
    lines.append(f"  coverage of reachable variants: {a}/{b} ({100 * a / b:.1f}%)" if b else "  no reachable variants")
    if "region" in cov:
        a, b = cov["region"]
        lo, hi = cov["region_bounds"]
        pct = 100 * a / b if b else 0.0
        lines.append(f"  coverage inside the folded region {lo}-{hi}: {a}/{b} ({pct:.1f}%)  -> "
                     + ("meets the 95% goal" if pct >= 95 else "BELOW the 95% goal"))
    if canary["passed"] is None:
        lines.append(f"  sign canary: {canary['note']}")
    else:
        lines.append(f"  sign canary: substitutions to proline score {canary['difference']:+.2f} kcal/mol above the rest -> "
                     + ("PASS (positive = destabilizing)" if canary["passed"] else "FAIL: check the sign and the mapping"))
    if directions.get("checked"):
        lines.append(f"  literature directions: {directions['agree']} of {directions['checked']} agree")
        for mut, expected, value, status in directions["rows"]:
            shown = "n/a" if value is None else f"{value:+.2f}"
            lines.append(f"    {mut}: predicted {shown} kcal/mol, literature says {expected}: {status}")
    elif directions.get("problem"):
        lines.append(f"  literature directions: {directions['problem']}")
    return "\n".join(lines)


def save_report(path, key, imp, cov, canary, directions) -> None:
    payload = {"target": key, "n_substitutions": len(imp.rows), "identity": imp.identity,
               "coverage_of_target": imp.coverage_of_target, "n_structure_mismatch": imp.n_structure_mismatch,
               "n_unmapped": imp.n_unmapped, "n_target_mismatch": imp.n_target_mismatch,
               "coverage": cov, "canary": canary, "directions": directions}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=str)
