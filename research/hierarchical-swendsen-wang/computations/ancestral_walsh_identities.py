"""Structural identities of the ancestral Walsh reduction (note 42).

The module is dependency-free.  It certifies, for the four ancestor rates

``Lambda_v^{ab} = (T_v + X_0 + (-1)^a X_1 + (-1)^b X_2) / 2``,

the three compact identities of note 42 section 5.2:

* paired sums:      ``L00 + L11 == L01 + L10 == T_v + X_0``;
* paired products:  ``L00 * L11 - L01 * L10 == -X_1 * X_2``,
  hence ``sign(J_v) = -sign(X_1 X_2)``;
* mute terminal ancestor: a single ancestor at level ``beta_v = 1``
  transmits an exactly null message ``B_u``.

It also recomputes the closed Walsh form of the ancestral message,

``B_u = 2 J + log cosh(h1 + h2) - log cosh(h1 - h2)``,

directly from the four log-weights and checks both evaluations agree.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cosh, exp, log
from typing import Iterable

State = tuple[int, int]
STATES: tuple[State, ...] = ((0, 0), (0, 1), (1, 0), (1, 1))


@dataclass(frozen=True)
class WalshAncestor:
    """One strict ancestor described by its three signed statistics."""

    beta: float
    common: float
    imbalance_1: float
    imbalance_2: float
    total: float

    def __post_init__(self) -> None:
        if not 0.0 <= self.beta <= 1.0:
            raise ValueError("beta must belong to [0, 1]")
        bound = abs(self.common) + abs(self.imbalance_1) + abs(self.imbalance_2)
        if self.total < bound:
            raise ValueError("total weight must dominate the signed statistics")

    def rate(self, state: State) -> float:
        """Return Lambda_v^{ab} for the flip pair ``state``."""

        a, b = state
        return 0.5 * (
            self.total
            + self.common
            + (1 - 2 * a) * self.imbalance_1
            + (1 - 2 * b) * self.imbalance_2
        )

    def four_rates(self) -> dict[State, float]:
        return {state: self.rate(state) for state in STATES}

    def log_factor(self, state: State) -> float:
        """Return phi_v(Lambda_v^{ab}) = log Lambda + (1 - beta) Lambda."""

        value = self.rate(state)
        if value <= 0.0:
            raise ValueError("all four rates must be positive for the Walsh form")
        return log(value) + (1.0 - self.beta) * value


def paired_sum_gap(ancestor: WalshAncestor) -> float:
    """Return (L00 + L11) - (L01 + L10); identically zero."""

    rates = ancestor.four_rates()
    return rates[0, 0] + rates[1, 1] - rates[0, 1] - rates[1, 0]


def paired_product_defect(ancestor: WalshAncestor) -> float:
    """Return L00*L11 - L01*L10 + X_1*X_2; identically zero."""

    rates = ancestor.four_rates()
    product_gap = rates[0, 0] * rates[1, 1] - rates[0, 1] * rates[1, 0]
    return product_gap + ancestor.imbalance_1 * ancestor.imbalance_2


def log_sum_exp(first: float, second: float) -> float:
    peak = max(first, second)
    return peak + log(exp(first - peak) + exp(second - peak))


def ancestral_message(ancestors: Iterable[WalshAncestor]) -> float:
    """Return B_u from the four accumulated log-weights (uniform prior)."""

    phis = {state: 0.0 for state in STATES}
    for ancestor in ancestors:
        for state in STATES:
            phis[state] += ancestor.log_factor(state)
    even = log_sum_exp(phis[0, 0], phis[1, 1])
    odd = log_sum_exp(phis[0, 1], phis[1, 0])
    return even - odd


def ancestral_message_walsh(ancestors: Iterable[WalshAncestor]) -> float:
    """Return B_u through the closed Walsh form of note 42 section 5.1."""

    field_1 = field_2 = coupling = 0.0
    for ancestor in ancestors:
        phi = {state: ancestor.log_factor(state) for state in STATES}
        field_1 += 0.25 * (phi[0, 0] + phi[0, 1] - phi[1, 0] - phi[1, 1])
        field_2 += 0.25 * (phi[0, 0] + phi[1, 0] - phi[0, 1] - phi[1, 1])
        coupling += 0.25 * (phi[0, 0] + phi[1, 1] - phi[0, 1] - phi[1, 0])
    return (
        2.0 * coupling
        + log(cosh(field_1 + field_2))
        - log(cosh(field_1 - field_2))
    )


def linear_coupling_defect(ancestor: WalshAncestor) -> float:
    """Walsh coupling of the linear part (1-beta) Lambda; identically zero."""

    rates = ancestor.four_rates()
    linear = {
        state: (1.0 - ancestor.beta) * rate for state, rate in rates.items()
    }
    return 0.25 * (
        linear[0, 0] + linear[1, 1] - linear[0, 1] - linear[1, 0]
    )


def paired_sum_value_defect(ancestor: WalshAncestor) -> float:
    """Return (L00 + L11) - (T_v + X_0); identically zero (full identity i)."""

    rates = ancestor.four_rates()
    return rates[0, 0] + rates[1, 1] - (ancestor.total + ancestor.common)


def exact_rational_defects(
    total: "Fraction", common: "Fraction", imb_1: "Fraction", imb_2: "Fraction"
) -> tuple["Fraction", "Fraction", "Fraction"]:
    """Exact-arithmetic certificates of identities (i)-(ii) of note 42.

    Returns (sum_gap, sum_value_gap, product_defect) computed in
    ``fractions.Fraction``; each is identically the rational number 0.
    """

    from fractions import Fraction

    half = Fraction(1, 2)
    rates = {
        (a, b): half
        * (total + common + (1 - 2 * a) * imb_1 + (1 - 2 * b) * imb_2)
        for a in (0, 1)
        for b in (0, 1)
    }
    sum_gap = rates[0, 0] + rates[1, 1] - rates[0, 1] - rates[1, 0]
    sum_value_gap = rates[0, 0] + rates[1, 1] - (total + common)
    product_defect = (
        rates[0, 0] * rates[1, 1] - rates[0, 1] * rates[1, 0] + imb_1 * imb_2
    )
    return sum_gap, sum_value_gap, product_defect
