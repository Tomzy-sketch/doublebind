"""Energy to susceptibility maths (PRD Section 19).

Convention: ddG_bind > 0 means weaker binding, which means more resistance.
IC50 ratio ~ exp(ddG_bind / RT); CSI = log2(IC50_wt / IC50_mut) = -ddG_bind / (RT ln 2).
A positive CSI means the mutant is MORE susceptible to the drug.
"""
from __future__ import annotations

import math
from statistics import NormalDist

from doublebind.constants import DEFAULT_T_K, R_KCAL_PER_MOL_K

LN2 = math.log(2.0)


def rt(temp_k: float = DEFAULT_T_K) -> float:
    """RT in kcal/mol."""
    return R_KCAL_PER_MOL_K * temp_k


def ic50_fold_change(ddg_bind: float, temp_k: float = DEFAULT_T_K) -> float:
    """IC50(mutant) / IC50(wild type), treating IC50 as proportional to Kd."""
    return math.exp(ddg_bind / rt(temp_k))


def csi_from_ddg_bind(ddg_bind: float, temp_k: float = DEFAULT_T_K) -> float:
    return -ddg_bind / (rt(temp_k) * LN2)


def ddg_bind_for_csi(csi: float, temp_k: float = DEFAULT_T_K) -> float:
    """The ddG_bind that corresponds to a given CSI (inverse of csi_from_ddg_bind)."""
    return -csi * rt(temp_k) * LN2


def ddg_from_log10_ic50_delta(delta_log10_ic50: float, temp_k: float = DEFAULT_T_K) -> float:
    """Convert a change in log10(IC50) (for example from Boltz-2) to ddG_bind."""
    return rt(temp_k) * math.log(10.0) * delta_log10_ic50


def prob_csi_at_least(mu: float, sigma: float, csi_min: float = 2.0,
                      temp_k: float = DEFAULT_T_K) -> float:
    """P(CSI >= csi_min) when ddG_bind ~ Normal(mu, sigma)."""
    threshold = ddg_bind_for_csi(csi_min, temp_k)
    if sigma <= 0:
        return 1.0 if mu <= threshold else 0.0
    return NormalDist().cdf((threshold - mu) / sigma)


def mean_needed_for_call(sigma: float, csi_min: float = 2.0, confidence: float = 0.8,
                         temp_k: float = DEFAULT_T_K) -> float:
    """Highest predicted mean ddG_bind that still gives P(CSI >= csi_min) >= confidence."""
    z = NormalDist().inv_cdf(confidence)
    return ddg_bind_for_csi(csi_min, temp_k) - z * sigma
