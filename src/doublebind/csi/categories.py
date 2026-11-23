"""Pair labels (PRD Section 19). Evaluated on terminal states already resistant to Drug A."""
from __future__ import annotations

COLLATERAL = "collateral sensitivity"
SPARING = "resistance-sparing"
CROSS = "cross-resistance"
INCONCLUSIVE = "inconclusive"


def label_pair(csi: float, p_csi_ge_2: float, p_min: float = 0.8, csi_min: float = 2.0,
               sparing_band: float = 1.0) -> str:
    """Return exactly one label for a (Drug A, Drug B) pair.

    csi          : CSI of Drug B at the terminal state (median across trajectories is fine)
    p_csi_ge_2   : probability that CSI >= csi_min under surrogate uncertainty
    """
    if csi >= csi_min and p_csi_ge_2 >= p_min:
        return COLLATERAL
    if csi <= -sparing_band:
        return CROSS
    if -sparing_band < csi < sparing_band:
        return SPARING
    return INCONCLUSIVE
