"""Small builders for synthetic genes, so tests never need the internet."""
import random

from doublebind.data.codons import CODON_TABLE

_REVERSE: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    _REVERSE.setdefault(_aa, []).append(_codon)

AMINO = "ACDEFGHIKLMNPQRSTVWY"


def random_protein(n, seed=0, fixed=None):
    rng = random.Random(seed)
    prot = ["M"] + [rng.choice(AMINO) for _ in range(n - 1)]
    for pos, aa in (fixed or {}).items():
        prot[pos - 1] = aa
    return "".join(prot)


def build_cds(protein, overrides=None, seed=0, stop="TAA"):
    """Encode `protein` with random synonymous codons; overrides maps residue number -> codon."""
    rng = random.Random(seed)
    overrides = overrides or {}
    codons = []
    for i, aa in enumerate(protein, start=1):
        if i in overrides:
            codons.append(overrides[i])
        elif i == 1:
            codons.append("ATG")
        else:
            codons.append(rng.choice(_REVERSE[aa]))
    return "".join(codons) + stop
