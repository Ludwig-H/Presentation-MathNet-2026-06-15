"""Exact stability benchmark of the beta = 0 (Glauber) dynamics on the SBM.

The module is dependency-free and Monte-Carlo-free.  For the binary
symmetric SBM with intra/inter probabilities ``p = a/n`` and ``q = b/n``,
it computes exactly, by binomial enumeration, the genie experiment of note
SBM/07: all labels but one are set to the truth and the single-site heat
bath (the beta = 0 limit of the hierarchical dynamics) resamples the last
label.

Quantities:

* ``affinity(n, a, b)``          -- the exact Bhattacharyya affinity
  ``rho_n ** (n - 1)`` of the local experiment;
* ``exp_half_llr(n, a, b)``      -- exact ``E[exp(-Delta_v / 2)]``; the
  identity of Proposition 3.1 states it EQUALS ``affinity``;
* ``flip_probability(n, a, b)``  -- exact heat-bath flip probability
  ``E[1 / (1 + exp(Delta_v))]`` at the truth;
* ``expected_sweep_flips(n, a, b)`` -- ``n`` times the above: the expected
  number of flips of one full sweep started from the truth;
* ``log_regime_crossing(A, B, n)``  -- the sign of
  ``log(n * flip) / log(n)``, negative above the exact-recovery
  threshold ``(sqrt(A) - sqrt(B))**2 = 2`` and positive below.

The enumeration splits the ``n - 1`` potential incident edges into
``n_same`` same-class slots (edge probability ``p``) and ``n_diff``
cross-class slots (probability ``q``); conditionally on the labels the
slots are independent, so the LLR only depends on the two binomial counts.
"""

from __future__ import annotations

from math import comb, exp, log, sqrt


def _binomial_pmf(count: int, probability: float) -> list[float]:
    """Return the full Binomial(count, probability) pmf."""

    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must belong to [0, 1]")
    return [
        comb(count, k) * probability**k * (1.0 - probability) ** (count - k)
        for k in range(count + 1)
    ]


def _slot_counts(n: int) -> tuple[int, int]:
    """Split the n - 1 other vertices into same/different class slots.

    With an i.i.d. uniform prior the class sizes fluctuate; the benchmark
    uses the balanced split, which matches the planted bisection and the
    leading exponents of note SBM/05.
    """

    others = n - 1
    n_same = others // 2
    return n_same, others - n_same


def _llr_terms(n: int, a: float, b: float) -> tuple[float, float, float]:
    """Return (log p/q, log (1-p)/(1-q), constant offset of the LLR).

    Delta_v = S_same * log(p/q) + (n_same - S_same) * log((1-p)/(1-q))
            - S_diff * log(p/q) - (n_diff - S_diff) * log((1-p)/(1-q)).
    """

    p, q = a / n, b / n
    if not 0.0 < q < p < 1.0:
        raise ValueError("need 0 < b/n < a/n < 1")
    return log(p / q), log((1.0 - p) / (1.0 - q)), 0.0


def _iterate_llr(n: int, a: float, b: float):
    """Yield (probability, Delta_v) over the exact binomial enumeration."""

    log_edge, log_gap, _ = _llr_terms(n, a, b)
    n_same, n_diff = _slot_counts(n)
    p, q = a / n, b / n
    pmf_same = _binomial_pmf(n_same, p)
    pmf_diff = _binomial_pmf(n_diff, q)
    for s_same, w_same in enumerate(pmf_same):
        base = s_same * log_edge + (n_same - s_same) * log_gap
        for s_diff, w_diff in enumerate(pmf_diff):
            delta = base - s_diff * log_edge - (n_diff - s_diff) * log_gap
            yield w_same * w_diff, delta


def affinity(n: int, a: float, b: float) -> float:
    """Return rho_n ** (n - 1), the exact local Bhattacharyya affinity."""

    p, q = a / n, b / n
    rho = sqrt(p * q) + sqrt((1.0 - p) * (1.0 - q))
    return rho ** (n - 1)


def exp_half_llr(n: int, a: float, b: float) -> float:
    """Return E[exp(-Delta_v / 2)] by exact enumeration."""

    return sum(weight * exp(-delta / 2.0) for weight, delta in _iterate_llr(n, a, b))


def _sigmoid_of_minus(delta: float) -> float:
    """Return 1 / (1 + exp(delta)) without overflow."""

    if delta >= 0.0:
        z = exp(-delta)
        return z / (1.0 + z)
    return 1.0 / (1.0 + exp(delta))


def flip_probability(n: int, a: float, b: float) -> float:
    """Return the exact heat-bath flip probability at the truth."""

    return sum(
        weight * _sigmoid_of_minus(delta) for weight, delta in _iterate_llr(n, a, b)
    )


def expected_sweep_flips(n: int, a: float, b: float) -> float:
    """Expected number of flips of one full sweep started from the truth."""

    return n * flip_probability(n, a, b)


def log_regime_parameters(n: int, big_a: float, big_b: float) -> tuple[float, float]:
    """Return (a_n, b_n) = (A log n, B log n) for the exact-recovery regime."""

    return big_a * log(n), big_b * log(n)


def log_regime_crossing(big_a: float, big_b: float, n: int) -> float:
    """Return log(n * flip) / log n in the logarithmic regime.

    Negative values mean the truth is stable under a full sweep (exact
    recovery side); positive values mean the expected flip count grows.
    The first-order theory predicts the sign of
    ``1 - (sqrt(A) - sqrt(B)) ** 2 / 2``.
    """

    a_n, b_n = log_regime_parameters(n, big_a, big_b)
    return log(expected_sweep_flips(n, a_n, b_n)) / log(n)
