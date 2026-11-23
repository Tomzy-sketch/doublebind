"""Sign conventions for energy changes (PRD Section 19).

DoubleBind convention: ddG > 0 means destabilizing (fold) or weaker binding (bind).
Convert every outside source ONCE, at ingestion, with to_doublebind().
Sources marked confirmed=False were written from memory or a third-party page:
open the source's own documentation and flip confirmed to True when you have checked.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Mapping


class Convention(str, Enum):
    DESTABILIZING_POSITIVE = "destabilizing_positive"
    STABILIZING_POSITIVE = "stabilizing_positive"


@dataclass(frozen=True)
class Source:
    convention: Convention
    confirmed: bool
    note: str = ""


SOURCES: dict[str, Source] = {
    "doublebind": Source(Convention.DESTABILIZING_POSITIVE, True, "PRD v2.0 Section 19"),
    "thermompnn": Source(Convention.DESTABILIZING_POSITIVE, False,
                         "A third-party page documents 'negative = stabilizing'; confirm in the repository"),
    "megascale": Source(Convention.STABILIZING_POSITIVE, False,
                        "Believed to report unfolding free energy (positive = stabilizing); confirm in the dataset README"),
    "foldx": Source(Convention.DESTABILIZING_POSITIVE, False, "Commonly reported this way; confirm"),
    "rosetta": Source(Convention.DESTABILIZING_POSITIVE, False, "Commonly reported this way; confirm"),
}


def register_source(name: str, convention: Convention, note: str = "", confirmed: bool = False) -> None:
    SOURCES[name.lower()] = Source(convention, confirmed, note)


def to_doublebind(value, source: str):
    """Return value in the DoubleBind convention. Works on floats, numpy arrays and pandas columns."""
    key = source.lower()
    if key not in SOURCES:
        raise KeyError(f"Unknown source '{source}'. Register it with register_source() after "
                       f"reading its documentation. Known: {sorted(SOURCES)}")
    return -value if SOURCES[key].convention is Convention.STABILIZING_POSITIVE else value


@dataclass(frozen=True)
class DirectionReport:
    n_agree: int
    n_total: int
    rows: tuple  # (mutation, expected, predicted, status)

    @property
    def fraction(self) -> float:
        return self.n_agree / self.n_total if self.n_total else float("nan")


def check_directions(predictions: Mapping[str, float], expected: Mapping[str, str]) -> DirectionReport:
    """Compare predicted ddG_fold (DoubleBind convention) with expected directions.

    expected values are 'destabilizing' or 'stabilizing'. Mutations missing from
    predictions are reported as 'missing' and do not count as agreement.
    """
    rows = []
    agree = 0
    for mutation, direction in expected.items():
        if direction not in ("destabilizing", "stabilizing"):
            raise ValueError(f"direction for {mutation} must be 'destabilizing' or 'stabilizing'")
        if mutation not in predictions:
            rows.append((mutation, direction, None, "missing"))
            continue
        value = predictions[mutation]
        predicted = "destabilizing" if value > 0 else "stabilizing"
        ok = predicted == direction
        agree += int(ok)
        rows.append((mutation, direction, value, "agree" if ok else "DISAGREE"))
    return DirectionReport(agree, len(expected), tuple(rows))
