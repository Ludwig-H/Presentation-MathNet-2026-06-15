"""Exact four-state calculus for the ancestors of a Kruskal merge.

The module is dependency-free.  It separates two tasks:

* evaluate all four values of ``Lambda_v`` after flipping either child of a
  target merge ``u``;
* integrate the homogeneous BSC marks conditionally on an edge-unlabelled
  Kruskal skeleton for small, exactly enumerable examples.

Group 0 of an ancestor bucket is unaffected by the two flips.  Groups 1 and 2
contain the edges whose endpoint in the child leading to ``u`` belongs to the
first or second child of ``u`` respectively.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb, exp, inf, isfinite, log, log1p, tanh
from typing import Iterable, Mapping


State = tuple[int, int]
STATES: tuple[State, ...] = ((0, 0), (0, 1), (1, 0), (1, 1))
EVEN_STATES: tuple[State, ...] = ((0, 0), (1, 1))
ODD_STATES: tuple[State, ...] = ((0, 1), (1, 0))


def coupling(p: float) -> float:
    """Return the Nishimori BSC coupling log(p / (1-p))."""

    if not 0.5 < p < 1.0:
        raise ValueError("p must satisfy 0.5 < p < 1")
    return log(p / (1.0 - p))


def closed_satisfaction_probability(p: float, beta: float) -> float:
    """Probability that an edge is satisfied given that it is closed at beta."""

    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must belong to [0, 1]")
    u_p = coupling(p)
    numerator = p * exp(-u_p * beta)
    return numerator / (1.0 - p + numerator)


@dataclass(frozen=True)
class AncestorBucket:
    """Sufficient statistics of one strict ancestor of the target merge."""

    beta: float
    total: tuple[float, float, float]
    satisfied: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not 0.0 <= self.beta <= 1.0:
            raise ValueError("beta must belong to [0, 1]")
        if len(self.total) != 3 or len(self.satisfied) != 3:
            raise ValueError("total and satisfied must have three groups")
        for total, satisfied in zip(self.total, self.satisfied):
            if total < 0.0 or not 0.0 <= satisfied <= total:
                raise ValueError("each group must satisfy 0 <= satisfied <= total")

    def lambda_after(self, state: State) -> float:
        """Return Lambda_v after the two indicated child flips."""

        a, b = state
        if a not in (0, 1) or b not in (0, 1):
            raise ValueError("a flip state must belong to {0, 1}^2")
        _, total_1, total_2 = self.total
        sat_0, sat_1, sat_2 = self.satisfied
        contribution_1 = total_1 - sat_1 if a else sat_1
        contribution_2 = total_2 - sat_2 if b else sat_2
        return sat_0 + contribution_1 + contribution_2

    def four_lambdas(self) -> dict[State, float]:
        """Return the exact four-state vector of ancestor rates."""

        return {state: self.lambda_after(state) for state in STATES}


def log_factor(rate: float, beta: float) -> float:
    """Return log(x exp((1-beta)x)), with log(0) represented by -inf."""

    if rate < 0.0:
        raise ValueError("a rate cannot be negative")
    if rate == 0.0:
        return -inf
    return log(rate) + (1.0 - beta) * rate


def _logsumexp(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return -inf
    maximum = max(values)
    if maximum == -inf:
        return -inf
    return maximum + log(sum(exp(value - maximum) for value in values))


def _log_ratio(log_numerator: float, log_denominator: float) -> float:
    if log_numerator == -inf and log_denominator == -inf:
        raise ValueError("both parity classes have zero conditional mass")
    if log_denominator == -inf:
        return inf
    if log_numerator == -inf:
        return -inf
    return log_numerator - log_denominator


def ancestral_log_weights(
    ancestors: Iterable[AncestorBucket],
    log_prior: Mapping[State, float] | None = None,
) -> dict[State, float]:
    """Accumulate the four outside log-weights above the target merge."""

    weights = {
        state: 0.0 if log_prior is None else log_prior[state] for state in STATES
    }
    for bucket in ancestors:
        rates = bucket.four_lambdas()
        for state in STATES:
            weights[state] += log_factor(rates[state], bucket.beta)
    return weights


def ancestral_message(
    ancestors: Iterable[AncestorBucket],
    log_prior: Mapping[State, float] | None = None,
) -> float:
    """Return B_u using a stable four-state log-sum-exp computation."""

    weights = ancestral_log_weights(ancestors, log_prior)
    log_even = _logsumexp(weights[state] for state in EVEN_STATES)
    log_odd = _logsumexp(weights[state] for state in ODD_STATES)
    return _log_ratio(log_even, log_odd)


@dataclass(frozen=True)
class LocalResult:
    """Exact parity message and persistence at the target merge."""

    log_odds: float
    persistence: float
    reliability: float


@dataclass(frozen=True)
class WalshCoefficients:
    """Four-state log-potential in the basis 1, (-1)^a, (-1)^b, (-1)^(a+b)."""

    constant: float
    field_1: float
    field_2: float
    coupling: float


def walsh_coefficients(log_weights: Mapping[State, float]) -> WalshCoefficients:
    """Collapse four finite log-weights to their exact Walsh coefficients."""

    values = tuple(log_weights[state] for state in STATES)
    if not all(isfinite(value) for value in values):
        raise ValueError("Walsh reduction requires four finite log-weights")
    w00 = log_weights[(0, 0)]
    w01 = log_weights[(0, 1)]
    w10 = log_weights[(1, 0)]
    w11 = log_weights[(1, 1)]
    return WalshCoefficients(
        constant=(w00 + w01 + w10 + w11) / 4.0,
        field_1=(w00 + w01 - w10 - w11) / 4.0,
        field_2=(w00 + w10 - w01 - w11) / 4.0,
        coupling=(w00 + w11 - w01 - w10) / 4.0,
    )


def _log_cosh(value: float) -> float:
    absolute = abs(value)
    return absolute + log1p(exp(-2.0 * absolute)) - log(2.0)


def message_from_coefficients(coefficients: WalshCoefficients) -> float:
    """Return B_u from the exact three nonconstant Walsh coefficients."""

    return (
        2.0 * coefficients.coupling
        + _log_cosh(coefficients.field_1 + coefficients.field_2)
        - _log_cosh(coefficients.field_1 - coefficients.field_2)
    )


def local_result(
    local_rate: float,
    local_total: float,
    local_beta: float,
    ancestors: Iterable[AncestorBucket] = (),
    log_prior: Mapping[State, float] | None = None,
) -> LocalResult:
    """Return L_u, tanh(L_u/2), and eta_u without indeterminate infinities."""

    if not 0.0 <= local_rate <= local_total:
        raise ValueError("local_rate must belong to [0, local_total]")
    if not 0.0 <= local_beta <= 1.0:
        raise ValueError("local_beta must belong to [0, 1]")

    weights = ancestral_log_weights(ancestors, log_prior)
    even_factor = log_factor(local_rate, local_beta)
    odd_factor = log_factor(local_total - local_rate, local_beta)
    log_even = _logsumexp(weights[state] + even_factor for state in EVEN_STATES)
    log_odd = _logsumexp(weights[state] + odd_factor for state in ODD_STATES)
    log_odds = _log_ratio(log_even, log_odd)
    if log_odds == inf:
        persistence = 1.0
    elif log_odds == -inf:
        persistence = -1.0
    else:
        persistence = tanh(log_odds / 2.0)
    return LocalResult(log_odds, persistence, persistence * persistence)


def binomial_pmf(size: int, probability: float) -> dict[int, float]:
    """Return a small exact-enumeration binomial PMF in floating point."""

    if size < 0:
        raise ValueError("size cannot be negative")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must belong to [0, 1]")
    return {
        value: comb(size, value)
        * probability**value
        * (1.0 - probability) ** (size - value)
        for value in range(size + 1)
    }


def conditional_total_count_pmf(
    size: int, p: float, beta: float
) -> dict[int, float]:
    """Law of 1 + Bin(size-1, s_p(beta)) at an unlabelled merge."""

    if size <= 0:
        raise ValueError("a merging bucket must contain at least one edge")
    probability = closed_satisfaction_probability(p, beta)
    return {
        value + 1: mass
        for value, mass in binomial_pmf(size - 1, probability).items()
    }


def grouped_count_pmf(
    group_sizes: tuple[int, int, int], p: float, beta: float
) -> dict[tuple[int, int, int], float]:
    """Exact grouped satisfied-count law conditional on an unlabelled merge.

    A latent winning edge is uniform in the bucket.  It is satisfied, while
    every other edge is independently satisfied with probability s_p(beta)
    conditionally on remaining closed at beta.
    """

    if len(group_sizes) != 3 or any(size < 0 for size in group_sizes):
        raise ValueError("group_sizes must contain three nonnegative integers")
    total_size = sum(group_sizes)
    if total_size == 0:
        raise ValueError("a merging bucket must contain at least one edge")

    probability = closed_satisfaction_probability(p, beta)
    result: dict[tuple[int, int, int], float] = {}
    for winner_group, winner_size in enumerate(group_sizes):
        if winner_size == 0:
            continue
        winner_mass = winner_size / total_size
        group_pmfs: list[dict[int, float]] = []
        for group, size in enumerate(group_sizes):
            winner = int(group == winner_group)
            pmf = {
                value + winner: mass
                for value, mass in binomial_pmf(size - winner, probability).items()
            }
            group_pmfs.append(pmf)
        for counts in product(*(pmf.keys() for pmf in group_pmfs)):
            mass = winner_mass
            for group, count in enumerate(counts):
                mass *= group_pmfs[group][count]
            result[counts] = result.get(counts, 0.0) + mass
    return result


@dataclass(frozen=True)
class AncestorGeometry:
    """Unmarked geometry needed for exact annealed enumeration."""

    beta: float
    group_sizes: tuple[int, int, int]


def exact_reliability_given_skeleton(
    local_size: int,
    local_beta: float,
    ancestor_geometry: Iterable[AncestorGeometry],
    p: float,
) -> float:
    """Enumerate E[eta_u | unmarked skeleton] for a small homogeneous chain."""

    u_p = coupling(p)
    local_pmf = conditional_total_count_pmf(local_size, p, local_beta)
    geometries = tuple(ancestor_geometry)
    ancestor_pmfs = tuple(
        grouped_count_pmf(geometry.group_sizes, p, geometry.beta)
        for geometry in geometries
    )

    expectation = 0.0
    count_spaces = tuple(pmf.items() for pmf in ancestor_pmfs)
    for local_count, local_mass in local_pmf.items():
        for ancestors_with_mass in product(*count_spaces):
            mass = local_mass
            buckets: list[AncestorBucket] = []
            for geometry, (counts, ancestor_mass) in zip(
                geometries, ancestors_with_mass
            ):
                mass *= ancestor_mass
                buckets.append(
                    AncestorBucket(
                        beta=geometry.beta,
                        total=tuple(u_p * size for size in geometry.group_sizes),
                        satisfied=tuple(u_p * count for count in counts),
                    )
                )
            expectation += mass * local_result(
                local_rate=u_p * local_count,
                local_total=u_p * local_size,
                local_beta=local_beta,
                ancestors=buckets,
            ).reliability
    return expectation


def main() -> None:
    p = 0.8
    local_size = 4
    local_beta = 1.0
    no_ancestor = exact_reliability_given_skeleton(
        local_size=local_size,
        local_beta=local_beta,
        ancestor_geometry=(),
        p=p,
    )
    with_ancestors = exact_reliability_given_skeleton(
        local_size=local_size,
        local_beta=0.45,
        ancestor_geometry=(
            AncestorGeometry(beta=0.62, group_sizes=(3, 2, 2)),
            AncestorGeometry(beta=0.81, group_sizes=(5, 1, 2)),
        ),
        p=p,
    )
    print(f"terminal check E[eta | m={local_size}, beta=1] = {no_ancestor:.12f}")
    print(f"expected value 1/m                         = {1/local_size:.12f}")
    print(f"two-ancestor example E[eta | skeleton]    = {with_ancestors:.12f}")


if __name__ == "__main__":
    main()
