"""Conditional estimates for every Lambda above a Kruskal merge.

The unmarked Kruskal skeleton fixes, for each ancestor bucket, its merge time
and its three edge groups.  Conditional on that skeleton, the winning edge is
latent.  With heterogeneous rates it is not uniform: its probability is
proportional to its conditional hazard ``weight * satisfaction_probability``.

This module implements:

* the exact heterogeneous winner mixture;
* exact conditional means and covariances of the three satisfied weights;
* their affine transport to the four rates after the two child flips;
* deterministic bounds on the three useful Walsh contrasts;
* certified transport of an ancestral-message error to the LCA reliability.

No external dependency is required.  Exact enumeration of all residual marks
is intentionally left to small counter-audits; the formulas below scale
linearly in the bucket size.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log, log1p, sqrt, tanh
from typing import Iterable

from ancestral_lambda_chain import AncestorBucket, STATES, State, coupling


RELIABILITY_LIPSCHITZ_CONSTANT = 2.0 / (3.0 * sqrt(3.0))


def lca_reliability(log_odds: float) -> float:
    """Return the exact LCA persistence ``tanh(log_odds / 2)^2``."""

    if not isfinite(log_odds):
        raise ValueError("log_odds must be finite")
    return tanh(0.5 * log_odds) ** 2


def reliability_tail_bound(message_error: float) -> float:
    """Transport a deterministic log-odds error to a reliability error.

    The derivative of ``tanh(x / 2)^2`` is
    ``tanh(x / 2) * sech(x / 2)^2``.  Its sharp global supremum is
    ``2 / (3 sqrt(3))``.
    """

    if not isfinite(message_error) or message_error < 0.0:
        raise ValueError("message_error must be finite and nonnegative")
    return min(1.0, RELIABILITY_LIPSCHITZ_CONSTANT * message_error)


def _logistic(value: float) -> float:
    if value >= 0.0:
        return 1.0 / (1.0 + exp(-value))
    exponential = exp(value)
    return exponential / (1.0 + exponential)


def _log_sigmoid(value: float) -> float:
    if value >= 0.0:
        return -log1p(exp(-value))
    return value - log1p(exp(value))


@dataclass(frozen=True)
class WeightedEdge:
    """One edge in an ancestor bucket.

    ``prior_satisfied`` is the annealed probability of a satisfied mark before
    conditioning its exponential clock to stay above the bucket time.  In a
    Nishimori edge channel it equals ``logistic(weight)``.
    """

    weight: float
    group: int
    prior_satisfied: float

    def __post_init__(self) -> None:
        if not isfinite(self.weight) or self.weight <= 0.0:
            raise ValueError("weight must be finite and strictly positive")
        if self.group not in (0, 1, 2):
            raise ValueError("group must belong to {0, 1, 2}")
        if not 0.0 < self.prior_satisfied < 1.0:
            raise ValueError("prior_satisfied must belong to (0, 1)")

    @classmethod
    def nishimori(cls, weight: float, group: int) -> "WeightedEdge":
        """Construct an edge whose prior log-odds equals its clock rate."""

        return cls(weight, group, _logistic(weight))


def closed_satisfaction_probability(edge: WeightedEdge, beta: float) -> float:
    """Return P(Y_e=1 | T_e > beta) for one marked exponential clock."""

    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must belong to [0, 1]")
    prior_log_odds = log(edge.prior_satisfied / (1.0 - edge.prior_satisfied))
    return _logistic(prior_log_odds - edge.weight * beta)


def winner_probabilities(
    edges: Iterable[WeightedEdge], beta: float
) -> tuple[float, ...]:
    """Exact latent-winner law conditional on an unmarked minimum at beta.

    The probability of edge e is proportional to w_e s_e(beta), not to w_e
    alone and not, in general, to one.
    """

    edges = tuple(edges)
    if not edges:
        raise ValueError("a merging bucket must contain at least one edge")
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must belong to [0, 1]")

    log_hazards: list[float] = []
    for edge in edges:
        prior_log_odds = log(
            edge.prior_satisfied / (1.0 - edge.prior_satisfied)
        )
        posterior_log_odds = prior_log_odds - edge.weight * beta
        log_hazards.append(log(edge.weight) + _log_sigmoid(posterior_log_odds))
    maximum = max(log_hazards)
    hazards = tuple(exp(value - maximum) for value in log_hazards)
    normalizer = sum(hazards)
    return tuple(value / normalizer for value in hazards)


@dataclass(frozen=True)
class GroupMoments:
    """Moments of (lambda_0, lambda_1, lambda_2) given the skeleton."""

    mean: tuple[float, float, float]
    covariance: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    winner_probability: tuple[float, ...]
    closed_probability: tuple[float, ...]


def conditional_group_moments(
    edges: Iterable[WeightedEdge], beta: float
) -> GroupMoments:
    """Return exact first and second moments after hiding the winning edge."""

    edges = tuple(edges)
    probabilities = tuple(
        closed_satisfaction_probability(edge, beta) for edge in edges
    )
    winner = winner_probabilities(edges, beta)

    baseline = [0.0, 0.0, 0.0]
    winner_correction = [0.0, 0.0, 0.0]
    for edge, probability, winner_mass in zip(edges, probabilities, winner):
        baseline[edge.group] += edge.weight * probability
        winner_correction[edge.group] += (
            winner_mass * edge.weight * (1.0 - probability)
        )

    mean = tuple(
        baseline[group] + winner_correction[group] for group in range(3)
    )
    covariance = [[0.0 for _ in range(3)] for _ in range(3)]

    # Expected conditional covariance given the latent winner.
    for edge, probability, winner_mass in zip(edges, probabilities, winner):
        covariance[edge.group][edge.group] += (
            (1.0 - winner_mass)
            * edge.weight**2
            * probability
            * (1.0 - probability)
        )

    # Covariance of the winner correction d_G.  Its vector has exactly one
    # nonzero coordinate, hence the negative off-diagonal terms after mixing.
    for edge, probability, winner_mass in zip(edges, probabilities, winner):
        correction = edge.weight * (1.0 - probability)
        covariance[edge.group][edge.group] += winner_mass * correction**2
    for first in range(3):
        for second in range(3):
            covariance[first][second] -= (
                winner_correction[first] * winner_correction[second]
            )

    covariance_tuple = (
        (covariance[0][0], covariance[0][1], covariance[0][2]),
        (covariance[1][0], covariance[1][1], covariance[1][2]),
        (covariance[2][0], covariance[2][1], covariance[2][2]),
    )
    return GroupMoments(
        mean=(mean[0], mean[1], mean[2]),
        covariance=covariance_tuple,
        winner_probability=winner,
        closed_probability=probabilities,
    )


@dataclass(frozen=True)
class FourRateMoments:
    """Moments of the ordered vector (Lambda^00, Lambda^01, Lambda^10, Lambda^11)."""

    states: tuple[State, ...]
    mean: tuple[float, ...]
    covariance: tuple[tuple[float, ...], ...]


def four_rate_moments(
    totals: tuple[float, float, float], moments: GroupMoments
) -> FourRateMoments:
    """Transport group moments through the exact complementation map."""

    if len(totals) != 3 or any(
        not isfinite(total) or total < 0.0 for total in totals
    ):
        raise ValueError("totals must contain three finite nonnegative values")

    coefficients: list[tuple[float, float, float]] = []
    constants: list[float] = []
    for a, b in STATES:
        coefficients.append((1.0, 1.0 - 2.0 * a, 1.0 - 2.0 * b))
        constants.append(a * totals[1] + b * totals[2])

    means = tuple(
        constant
        + sum(coefficient[group] * moments.mean[group] for group in range(3))
        for constant, coefficient in zip(constants, coefficients)
    )
    covariance: list[list[float]] = []
    for first in coefficients:
        row: list[float] = []
        for second in coefficients:
            row.append(
                sum(
                    first[r] * moments.covariance[r][s] * second[s]
                    for r in range(3)
                    for s in range(3)
                )
            )
        covariance.append(row)
    return FourRateMoments(
        states=STATES,
        mean=means,
        covariance=tuple(tuple(row) for row in covariance),
    )


def homogeneous_four_rate_moments(
    group_sizes: tuple[int, int, int], p: float, beta: float
) -> FourRateMoments:
    """Closed-form four-rate moments for an equal-weight Nishimori bucket."""

    if len(group_sizes) != 3 or any(
        not isinstance(size, int) or size < 0 for size in group_sizes
    ):
        raise ValueError("group_sizes must contain three nonnegative integers")
    if sum(group_sizes) == 0:
        raise ValueError("a merging bucket must contain at least one edge")
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must belong to [0, 1]")
    size = sum(group_sizes)
    weight = coupling(p)
    probability = _logistic(weight * (1.0 - beta))
    imbalance = 2.0 * probability - 1.0
    alpha = imbalance + (1.0 - imbalance) / size
    proportions = tuple(group_size / size for group_size in group_sizes)
    coefficients = tuple(
        (1.0, 1.0 - 2.0 * a, 1.0 - 2.0 * b) for a, b in STATES
    )
    signed_sizes = tuple(
        sum(group_sizes[group] * coefficient[group] for group in range(3))
        for coefficient in coefficients
    )
    means = tuple(
        0.5 * weight * (size + alpha * signed_size)
        for signed_size in signed_sizes
    )
    covariance: list[list[float]] = []
    for first in coefficients:
        row: list[float] = []
        for second in coefficients:
            signed_product = sum(
                group_sizes[group] * first[group] * second[group]
                for group in range(3)
            )
            winner_product = sum(
                proportions[group] * first[group] * second[group]
                for group in range(3)
            )
            winner_first = sum(
                proportions[group] * first[group] for group in range(3)
            )
            winner_second = sum(
                proportions[group] * second[group] for group in range(3)
            )
            row.append(
                weight**2
                * (
                    probability
                    * (1.0 - probability)
                    * (1.0 - 1.0 / size)
                    * signed_product
                    + (1.0 - probability) ** 2
                    * (winner_product - winner_first * winner_second)
                )
            )
        covariance.append(row)
    return FourRateMoments(
        states=STATES,
        mean=means,
        covariance=tuple(tuple(row) for row in covariance),
    )


@dataclass(frozen=True)
class ContrastEnvelope:
    """Deterministic upper bounds for one ancestor's useful contrasts."""

    field_1: float
    field_2: float
    coupling: float
    message: float
    minimum_rate: float


def contrast_envelope(bucket: AncestorBucket) -> ContrastEnvelope:
    """Bound |h_1|, |h_2|, |J| and this ancestor's message impact.

    The bounds follow from exact oriented-integral formulas for phi' and the
    mixed second increment of log.  Four strictly positive corner rates are
    required; zero-rate cases must use the exact four-state calculation.
    """

    rates = bucket.four_lambdas()
    minimum = min(rates.values())
    if minimum <= 0.0:
        raise ValueError("contrast bounds require four strictly positive rates")
    delta_1 = bucket.total[1] - 2.0 * bucket.satisfied[1]
    delta_2 = bucket.total[2] - 2.0 * bucket.satisfied[2]
    derivative_bound = (1.0 - bucket.beta) + 1.0 / minimum
    field_1 = 0.5 * abs(delta_1) * derivative_bound
    field_2 = 0.5 * abs(delta_2) * derivative_bound
    interaction = abs(delta_1 * delta_2) / (4.0 * minimum**2)
    return ContrastEnvelope(
        field_1=field_1,
        field_2=field_2,
        coupling=interaction,
        message=2.0 * (field_1 + field_2 + interaction),
        minimum_rate=minimum,
    )


def ancestral_tail_bound(ancestors: Iterable[AncestorBucket]) -> float:
    """Return a certified bound on the change in B after omitting ancestors."""

    return sum(contrast_envelope(bucket).message for bucket in ancestors)


def main() -> None:
    group_sizes = (5, 3, 2)
    moments = homogeneous_four_rate_moments(group_sizes, p=0.8, beta=0.7)
    print("conditional means of (Lambda^00, Lambda^01, Lambda^10, Lambda^11)")
    print(" ".join(f"{value:.9f}" for value in moments.mean))


if __name__ == "__main__":
    main()
