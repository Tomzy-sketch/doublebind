"""Genetic code and single-nucleotide mutation grammar.

PRD refs: FR-A1 (the gene DNA must translate to the UniProt protein) and FR-A2
(enumerate every single-nucleotide change per codon, keep missense variants, weight
them by a transition/transversion ratio).

Mutation supply is not selection. This module only says which amino-acid changes one
nucleotide can reach and how likely the DNA is to produce them. Fitness is handled
elsewhere.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

BASES = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}
STOP = "*"
_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
# Bacterial genes can start with these codons; the ribosome still inserts Met.
ALT_START_CODONS = {"GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}


def clean_dna(seq: str) -> str:
    s = "".join(seq.split()).upper().replace("U", "T")
    bad = set(s) - set("ACGT")
    if bad:
        raise ValueError(f"unexpected characters in DNA: {sorted(bad)}")
    return s


def translate(dna: str) -> str:
    s = clean_dna(dna)
    if len(s) % 3:
        raise ValueError(f"length {len(s)} is not a multiple of 3")
    return "".join(CODON_TABLE[s[i:i + 3]] for i in range(0, len(s), 3))


def cds_to_protein(dna: str, alt_start: bool = True, strip_stop: bool = True) -> str:
    """Translate a coding sequence the way the cell does (alternative start codons give Met)."""
    s = clean_dna(dna)
    prot = translate(s)
    if alt_start and prot and s[:3] in ALT_START_CODONS:
        prot = "M" + prot[1:]
    if strip_stop and prot.endswith(STOP):
        prot = prot[:-1]
    return prot


@dataclass(frozen=True)
class CdsCheck:
    ok: bool
    identity: float
    differences: tuple  # (position, protein_aa, cds_aa), 1-based
    length_protein: int
    length_cds: int


def check_cds(dna: str, protein: str) -> CdsCheck:
    """FR-A1: does this DNA encode exactly this protein?"""
    translated = cds_to_protein(dna)
    n = min(len(protein), len(translated))
    diffs = tuple((i + 1, protein[i], translated[i]) for i in range(n) if protein[i] != translated[i])
    longest = max(len(protein), len(translated))
    identity = (n - len(diffs)) / longest if longest else 1.0
    ok = len(protein) == len(translated) and not diffs
    return CdsCheck(ok, identity, diffs, len(protein), len(translated))


def ts_tv_probability(ref: str, alt: str, kappa: float = 2.0) -> float:
    """Chance that a mutation at base `ref` becomes `alt`, with kappa = transition rate / transversion rate.

    Each base has 1 transition and 2 transversions available, so the three probabilities
    are kappa/(kappa+2), 1/(kappa+2), 1/(kappa+2) and always sum to 1. kappa=2.0 is a
    placeholder: estimate it from the organism's mutation spectrum when you can.
    """
    if ref == alt:
        raise ValueError("ref and alt must differ")
    return kappa / (kappa + 2.0) if (ref, alt) in _TRANSITIONS else 1.0 / (kappa + 2.0)


@dataclass(frozen=True)
class NtVariant:
    pos: int            # residue number, 1-based, in the numbering of the supplied sequence
    wt_codon: str
    mut_codon: str
    codon_site: int     # 0, 1 or 2
    ref_nt: str
    alt_nt: str
    wt_aa: str
    mut_aa: str
    transition: bool
    weight: float       # relative supply of this exact change (see ts_tv_probability)

    @property
    def effect(self) -> str:
        if self.mut_aa == STOP:
            return "nonsense"
        return "synonymous" if self.mut_aa == self.wt_aa else "missense"

    @property
    def label(self) -> str:
        return f"{self.wt_aa}{self.pos}{self.mut_aa}"


def enumerate_variants(cds: str, kappa: float = 2.0, include_start: bool = False) -> list[NtVariant]:
    """All single-nucleotide substitutions of every sense codon (9 per codon).

    The terminal stop codon is skipped. The start codon is skipped unless include_start,
    because changing it is a start-loss event and not a missense change.
    """
    s = clean_dna(cds)
    if len(s) % 3:
        raise ValueError(f"CDS length {len(s)} is not a multiple of 3")
    n_codons = len(s) // 3
    if n_codons and CODON_TABLE[s[-3:]] == STOP:
        n_codons -= 1
    out: list[NtVariant] = []
    for idx in range(n_codons):
        pos = idx + 1
        if pos == 1 and not include_start:
            continue
        codon = s[3 * idx:3 * idx + 3]
        wt_aa = CODON_TABLE[codon]
        if wt_aa == STOP:
            continue
        if pos == 1 and codon in ALT_START_CODONS:
            wt_aa = "M"
        for site in range(3):
            ref = codon[site]
            for alt in "ACGT":
                if alt == ref:
                    continue
                mut = codon[:site] + alt + codon[site + 1:]
                out.append(NtVariant(pos, codon, mut, site, ref, alt, wt_aa, CODON_TABLE[mut],
                                     (ref, alt) in _TRANSITIONS, ts_tv_probability(ref, alt, kappa)))
    return out


@dataclass(frozen=True)
class AaVariant:
    pos: int
    wt_aa: str
    mut_aa: str
    weight: float   # summed supply over all nucleotide routes to this amino-acid change
    n_routes: int

    @property
    def label(self) -> str:
        return f"{self.wt_aa}{self.pos}{self.mut_aa}"


def aggregate_missense(variants: list[NtVariant]) -> list[AaVariant]:
    """Merge nucleotide routes into amino-acid variants (missense only)."""
    acc: dict[tuple[int, str, str], tuple[float, int]] = {}
    for v in variants:
        if v.effect != "missense":
            continue
        key = (v.pos, v.wt_aa, v.mut_aa)
        w, n = acc.get(key, (0.0, 0))
        acc[key] = (w + v.weight, n + 1)
    return [AaVariant(pos, wt, mut, w, n) for (pos, wt, mut), (w, n) in sorted(acc.items())]


def effect_counts(variants: list[NtVariant]) -> Counter:
    return Counter(v.effect for v in variants)


_MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse_mutation(label: str) -> tuple[str, int, str]:
    m = _MUT_RE.match(label.strip())
    if not m:
        raise ValueError(f"cannot parse mutation '{label}' (expected a form like S83L)")
    return m.group(1), int(m.group(2)), m.group(3)


def find_numbering_offsets(protein: str, labels: list[str], search: int = 40) -> list[int]:
    """Offsets d such that sequence_position = literature_position + d for ALL labels.

    Literature often uses a different numbering from UniProt (for example Ambler numbering
    for TEM-1). An offset is returned only if every wild-type residue lines up. Smallest
    absolute offsets come first. An empty list means no consistent numbering was found.
    """
    parsed = [parse_mutation(label) for label in labels]
    found = []
    for off in sorted(range(-search, search + 1), key=abs):
        if all(1 <= pos + off <= len(protein) and protein[pos + off - 1] == wt for wt, pos, _ in parsed):
            found.append(off)
    return found
