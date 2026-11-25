"""Day 2 baselines for the resistance-mutation recall metric (PRD Section 9).

Known mutations are few (two for GyrA), so every number here is an anecdote and not a statistic:
use it to check that the pipeline behaves sensibly, not to claim performance.

Scorers, all decided before looking at results:
  random                    expected rank (N+1)/2, chance of the top 10 = 10/N
  proximity only            closest to the ligand first
  ESM tolerated             positions whose most plausible reachable substitution has the highest LLR
  ESM constrained           positions whose reachable substitutions have the lowest mean LLR (diagnostic)
  proximity then ESM        inside the pocket cutoff, ranked like 'ESM tolerated' (first real candidate)
"""
from __future__ import annotations

import csv
from pathlib import Path

from doublebind.data.codons import find_numbering_offsets, parse_mutation

TOP_N = 10
ESM_PREFERENCE = ("t33_650M", "t30_150M", "t12_35M", "t6_8M")


def read_rows(path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def find_esm_file(processed_dir, key, explicit=None):
    """Pick an ESM score file: explicit path, else the largest model available."""
    if explicit:
        return Path(explicit)
    files = [p for p in sorted(Path(processed_dir).glob(f"{key}_esm_*.csv")) if not p.name.endswith("_positions.csv")]
    for tag in ESM_PREFERENCE:
        for p in files:
            if tag in p.name:
                return p
    return files[0] if files else None


def average_ranks(scores: dict, higher_is_better: bool = True) -> dict:
    """Rank 1 = best. Ties share the average of the ranks they span."""
    items = sorted(scores.items(), key=lambda kv: kv[1], reverse=higher_is_better)
    ranks, i = {}, 0
    while i < len(items):
        j = i
        while j + 1 < len(items) and items[j + 1][1] == items[i][1]:
            j += 1
        for k in range(i, j + 1):
            ranks[items[k][0]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def _summarize(label, scores, higher_is_better, wt, known_pos, top_n):
    ranks = average_ranks(scores, higher_is_better)
    order = sorted(scores, key=lambda p: ranks[p])
    return {
        "label": label,
        "n_ranked": len(scores),
        "known_ranks": {p: ranks.get(p) for p in known_pos},
        "top": [{"pos": p, "wt_aa": wt.get(p, "?"), "score": scores[p], "rank": ranks[p]} for p in order[:top_n]],
    }


def run_baselines(protein, aa_variants, known_rows, esm_rows=None, pocket_rows=None, cutoff=8.0, top_n=TOP_N) -> dict:
    """Pure function: no file access. aa_variants are AaVariant objects from codons.aggregate_missense."""
    wt = {v.pos: v.wt_aa for v in aa_variants}
    positions = sorted(wt)
    n_pos = len(positions)
    reachable = {(v.pos, v.mut_aa) for v in aa_variants}

    labels = [r["mutation"] for r in known_rows]
    offset = 0
    if labels:
        offsets = find_numbering_offsets(protein, labels)
        if not offsets:
            raise ValueError("no single numbering offset lines up every known mutation with the protein sequence")
        offset = offsets[0]
    known = []
    for r in known_rows:
        w, p, m = parse_mutation(r["mutation"])
        pos = p + offset
        known.append({"mutation": r["mutation"], "pos": pos, "wt_aa": w, "mut_aa": m,
                      "reachable": (pos, m) in reachable})
    known_pos = sorted({k["pos"] for k in known})

    llr = {}
    for r in esm_rows or []:
        key = (int(r["pos"]), r["mut_aa"])
        if key in reachable:
            llr[key] = float(r["llr"])
    by_pos: dict[int, list] = {}
    for (pos, _mut), val in llr.items():
        by_pos.setdefault(pos, []).append(val)
    esm_max = {p: max(v) for p, v in by_pos.items()}
    esm_mean = {p: sum(v) / len(v) for p, v in by_pos.items()}

    dist = {int(r["pos"]): float(r["min_dist"]) for r in (pocket_rows or []) if r.get("min_dist") not in (None, "")}
    pocket_positions = {p for p, d in dist.items() if d <= cutoff}

    position_level = {}
    if dist:
        position_level["proximity"] = _summarize("proximity only", {p: -d for p, d in dist.items()}, True, wt, known_pos, top_n)
    if esm_max:
        position_level["esm_tolerated"] = _summarize("ESM tolerated (max LLR)", esm_max, True, wt, known_pos, top_n)
        position_level["esm_constrained"] = _summarize("ESM constrained (mean LLR, low first)", esm_mean, False, wt, known_pos, top_n)
    if esm_max and dist:
        inside = {p: esm_max[p] for p in pocket_positions if p in esm_max}
        if inside:
            position_level["proximity_then_esm"] = _summarize(
                f"proximity <= {cutoff:g} A, then ESM tolerated", inside, True, wt, known_pos, top_n)

    variant_level = {}
    if llr:
        ranks_all = average_ranks(llr, True)
        in_pocket = {k: v for k, v in llr.items() if k[0] in pocket_positions}
        ranks_pocket = average_ranks(in_pocket, True) if in_pocket else {}
        variant_level = {
            "n_variants": len(llr),
            "n_pocket_variants": len(in_pocket),
            "known": {k["mutation"]: {
                "reachable": k["reachable"],
                "llr": llr.get((k["pos"], k["mut_aa"])),
                "rank_all": ranks_all.get((k["pos"], k["mut_aa"])),
                "rank_in_pocket": ranks_pocket.get((k["pos"], k["mut_aa"])),
            } for k in known},
        }

    return {
        "n_positions": n_pos, "numbering_offset": offset, "cutoff": cutoff, "top_n": top_n,
        "known": known,
        "random": {"expected_rank": (n_pos + 1) / 2, "p_top_n_each": min(1.0, top_n / n_pos) if n_pos else 0.0},
        "n_pocket_positions": len(pocket_positions),
        "position_level": position_level,
        "variant_level": variant_level,
    }


def format_report(key: str, res: dict) -> str:
    known_pos = sorted({k["pos"] for k in res["known"]})
    lines = [f"{key} baselines: {res['n_positions']} positions, numbering offset {res['numbering_offset']:+d}, "
             f"pocket cutoff {res['cutoff']:g} A ({res['n_pocket_positions']} positions inside)"]
    for k in res["known"]:
        lines.append(f"  known {k['mutation']} = position {k['pos']} ({'reachable' if k['reachable'] else 'NOT reachable'} by one nucleotide change)")
    lines += ["", f"Position-level rank of the known positions (1 = best; the top {res['top_n']} is what counts)"]
    head = f"  {'method':<42}{'ranked':>7}" + "".join(f"{('pos ' + str(p)):>10}" for p in known_pos)
    lines.append(head)
    rnd = res["random"]
    lines.append(f"  {'random (expected rank)':<42}{res['n_positions']:>7}" + "".join(f"{rnd['expected_rank']:>10.1f}" for _ in known_pos)
                 + f"   chance of top {res['top_n']}: {100 * rnd['p_top_n_each']:.1f}% each")
    for s in res["position_level"].values():
        cells = "".join(f"{('n/a' if s['known_ranks'][p] is None else format(s['known_ranks'][p], '.1f')):>10}" for p in known_pos)
        lines.append(f"  {s['label']:<42}{s['n_ranked']:>7}{cells}")
    for s in res["position_level"].values():
        top = ", ".join(f"{t['wt_aa']}{t['pos']}" for t in s["top"])
        lines.append(f"  top {res['top_n']} by {s['label']}: {top}")
    v = res["variant_level"]
    if v:
        lines += ["", f"Variant-level rank by ESM LLR (higher = more plausible): {v['n_variants']} reachable variants, "
                      f"{v['n_pocket_variants']} inside the pocket"]
        for name, info in v["known"].items():
            llr = "n/a" if info["llr"] is None else f"{info['llr']:+.2f}"
            ra = "n/a" if info["rank_all"] is None else f"{info['rank_all']:.1f}"
            rp = "n/a" if info["rank_in_pocket"] is None else f"{info['rank_in_pocket']:.1f}"
            lines.append(f"  {name:>6}: LLR {llr}, rank among all {ra}, rank among pocket variants {rp}")
    return "\n".join(lines)
