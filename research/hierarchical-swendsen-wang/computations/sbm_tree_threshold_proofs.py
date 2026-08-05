"""Exact certificates for the from-scratch proofs of note SBM/08.

The module is dependency-free and Monte-Carlo-free.  It certifies, by
exact enumeration on small deterministic trees, the ingredients of the
two tree theorems of SBM/08:

* ``psi(m)``                -- the binary divergence
  ``((1+m)/2) ln(1+m) + ((1-m)/2) ln(1-m)`` and its power series
  ``sum_k m^{2k} / (2k(2k-1))`` (Lemme A);
* the SDPI inequality ``psi(theta m) <= theta^2 psi(m)`` and the
  quadratic minoration ``psi(m) >= m^2 / 2`` (Lemmes A-B);
* on the complete d-ary tree of depth t with broadcast correlation
  ``theta``: the exact reconstruction weight ``q_t``, the exact mutual
  information ``I(sigma_rho; sigma_L)``, the information bound
  ``I <= ln 2 * d^t theta^{2t}`` (Théorème I.18, whose inductive bound
  is a fixed-tree statement, genuinely instantiated here) and the
  bound ``q_t >= ell_t(d theta^2)`` — on the REGULAR tree this is a
  strictly weakened corollary of Théorème I.17 (the regular variance
  E[N(N-1)] = d(d-1) is smaller than the Poisson d^2, so the
  inequality holds a fortiori and does not certify the Poisson
  recursion itself);
* the Poisson second-moment recursion of Théorème I.17, certified in
  EXACT rational arithmetic: v_t = d v_{t-1} + d^2 theta^2 (d theta)^{2t-2}
  coincides with the closed form (d theta)^{2t} sum_{s<=t} lambda^{-s};
* the scalar identities of Faits I.9 (frozen probability (a-b)/n on
  both pair types), I.12 (residual channel), I.13 (shared-cut
  inflation) and Proposition II.5 (pair affinity value).
"""

from __future__ import annotations

from math import exp, log, sqrt


def psi(m: float) -> float:
    """Binary divergence psi(m) in nats; psi(+-1) = ln 2."""

    if not -1.0 <= m <= 1.0:
        raise ValueError("m must belong to [-1, 1]")
    total = 0.0
    for share, arg in (((1.0 + m) / 2.0, 1.0 + m), ((1.0 - m) / 2.0, 1.0 - m)):
        if arg > 0.0:
            total += share * log(arg)
    return total


def psi_series(m: float, terms: int = 4000) -> float:
    """Partial sum of psi's power series sum_k m^{2k}/(2k(2k-1))."""

    return sum(m ** (2 * k) / (2 * k * (2 * k - 1)) for k in range(1, terms + 1))


def leaf_distributions(depth: int, arity: int, theta: float):
    """Return exact leaf-vector laws given root spin +1 and -1.

    Output: two dicts mapping leaf tuples (in {-1,+1}^{arity**depth}) to
    probabilities, for root spin +1 and -1 respectively.
    """

    p = (1.0 + theta) / 2.0

    def dist(k: int, spin: int):
        if k == 0:
            return {(spin,): 1.0}
        sub = {+1: dist(k - 1, +1), -1: dist(k - 1, -1)}
        base = {(): 1.0}
        for _ in range(arity):
            new = {}
            for vec, weight in base.items():
                for child_spin in (+1, -1):
                    w_child = p if child_spin == spin else 1.0 - p
                    for child_vec, child_weight in sub[child_spin].items():
                        key = vec + child_vec
                        new[key] = new.get(key, 0.0) + weight * w_child * child_weight
            base = new
        return base

    return dist(depth, +1), dist(depth, -1)


def tree_reconstruction_quantities(depth: int, arity: int, theta: float):
    """Return (q_t, information, info_bound, ell_t) for the d-ary tree.

    ``q_t`` is the exact reconstruction weight E[E[sigma_rho|leaves]^2],
    ``information`` the exact I(sigma_rho; leaves) in nats,
    ``info_bound`` the fixed-tree bound ln2 * (arity**depth) * theta^(2 depth),
    ``ell_t`` the second-moment lower bound with lambda = arity*theta^2.
    """

    law_plus, law_minus = leaf_distributions(depth, arity, theta)
    q_t = 0.0
    information = 0.0
    for vec, weight_plus in law_plus.items():
        weight_minus = law_minus.get(vec, 0.0)
        mass = (weight_plus + weight_minus) / 2.0
        if mass > 0.0:
            magnet = (weight_plus - weight_minus) / (weight_plus + weight_minus)
            q_t += mass * magnet * magnet
            information += mass * psi(magnet)
    leaves = arity**depth
    info_bound = log(2.0) * leaves * theta ** (2 * depth)
    lam = arity * theta * theta
    ell_t = 1.0 / sum(lam ** (-s) for s in range(depth + 1))
    return q_t, information, info_bound, ell_t


def frozen_probability_same_class(n: int, a: float, b: float) -> float:
    """Exact freeze probability of an intra-class pair (Fait I.9)."""

    return (a / n) * (1.0 - b / a)


def frozen_probability_cross_class(n: int, a: float, b: float) -> float:
    """Exact freeze probability of a cross-class pair (Fait I.9)."""

    return (1.0 - b / n) * (1.0 - (1.0 - a / n) / (1.0 - b / n))


def residual_channel_identity(theta: float, q: float) -> tuple[float, float]:
    """Return (marginalisation defect, shared-cut inflation defect).

    First value: q*1 + (1-q)*c0 - theta (identically zero, Fait I.12).
    Second: q + (1-q)*c0^2 - theta^2 - q(1-theta)^2/(1-q) (zero, I.13).
    """

    c0 = (theta - q) / (1.0 - q)
    marg = q + (1.0 - q) * c0 - theta
    shared = q + (1.0 - q) * c0 * c0 - theta * theta - q * (1.0 - theta) ** 2 / (1.0 - q)
    return marg, shared


def pair_affinity(n: int, a: float, b: float) -> float:
    """Exact Bhattacharyya factor of one incident pair (Prop. II.5)."""

    return sqrt(a * b) / n + sqrt((1.0 - a / n) * (1.0 - b / n))


def pair_affinity_from_channel(n: int, a: float, b: float) -> float:
    """Same value computed from the channel definition sum_w sqrt(P P)."""

    p_edge_plus, p_edge_minus = a / n, b / n
    return sqrt(p_edge_plus * p_edge_minus) + sqrt(
        (1.0 - p_edge_plus) * (1.0 - p_edge_minus)
    )

def poisson_second_moment_defect(
    d: "Fraction", theta: "Fraction", depth: int
) -> "Fraction":
    """Exact-arithmetic certificate of the recursion of Théorème I.17.

    Returns v_depth(recursion) - v_depth(closed form) in
    ``fractions.Fraction``; identically the rational number 0, where
    v_t = d v_{t-1} + d^2 theta^2 (d theta)^{2(t-1)} with v_0 = 1 and
    the closed form is (d theta)^{2t} sum_{s=0}^{t} (d theta^2)^{-s}.
    """

    from fractions import Fraction

    v = Fraction(1)
    for t in range(1, depth + 1):
        v = d * v + d * d * theta * theta * (d * theta) ** (2 * (t - 1))
    lam = d * theta * theta
    closed = (d * theta) ** (2 * depth) * sum(
        lam ** (-sst) for sst in range(depth + 1)
    )
    return v - closed
