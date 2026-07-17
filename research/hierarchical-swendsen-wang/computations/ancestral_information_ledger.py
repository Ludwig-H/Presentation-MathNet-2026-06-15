"""Exact information ledger for hierarchical LCA side information.

The generic ledger decomposes the increase of a binary posterior
magnetization as side-information coordinates are revealed.  The triangular
functions give a closed audit of the favorable Kruskal event on one face.
They are diagnostics for the hierarchical route, not a weak-recovery proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, expm1, inf, log, log1p, log2, sqrt, tanh
from typing import Hashable, Mapping, Sequence

from critical_band_thresholds import Q_CRITICAL, beta_critical, coupling
from nishimori_hierarchical_entropy import (
    binary_entropy_bits,
    conjectured_nishimori_root,
    face_syndrome_probability,
)


Outcome = tuple[Hashable, ...]


@dataclass(frozen=True)
class BinaryFiltrationLedger:
    """Scale-by-scale posterior statistics for a finite binary experiment."""

    second_moments: tuple[float, ...]
    conditional_entropies_bits: tuple[float, ...]
    l2_gains: tuple[float, ...]
    information_gains_bits: tuple[float, ...]


@dataclass(frozen=True)
class LogitIncrementAudit:
    """Exact quantities and universal upper bounds for one logit update."""

    kl_nats: float
    squared_magnetization_change: float
    kl_quadratic_bound: float
    magnetization_quadratic_bound: float
    magnetization_pinsker_bound: float


def _logistic(value: float) -> float:
    if value >= 0.0:
        return 1.0 / (1.0 + exp(-value))
    exponential = exp(value)
    return exponential / (1.0 + exponential)


def binary_filtration_ledger(
    joint_weights: Mapping[Outcome, float],
) -> BinaryFiltrationLedger:
    """Reveal the non-target coordinates from left to right.

    Each key is ``(x, y_1, ..., y_h)`` with ``x`` in ``{-1,+1}``.
    Values may be normalized probabilities or arbitrary non-negative weights.
    Stage zero reveals no ``y`` coordinate; stage ``k`` reveals the first
    ``k`` coordinates.
    """

    if not joint_weights:
        raise ValueError("joint_weights must not be empty")
    lengths = {len(outcome) for outcome in joint_weights}
    if len(lengths) != 1 or next(iter(lengths)) < 1:
        raise ValueError("all outcomes must have the same positive length")
    for outcome, weight in joint_weights.items():
        if outcome[0] not in (-1, 1):
            raise ValueError("the first coordinate must belong to {-1,+1}")
        if weight < 0.0:
            raise ValueError("joint weights must be non-negative")
    total = sum(joint_weights.values())
    if total <= 0.0:
        raise ValueError("at least one joint weight must be positive")

    normalized = tuple(
        (outcome, weight / total)
        for outcome, weight in joint_weights.items()
        if weight > 0.0
    )
    side_coordinates = next(iter(lengths)) - 1
    second_moments: list[float] = []
    entropies: list[float] = []

    for revealed in range(side_coordinates + 1):
        cells: dict[tuple[Hashable, ...], list[float]] = {}
        for outcome, probability in normalized:
            prefix = outcome[1 : 1 + revealed]
            mass_and_signed = cells.setdefault(prefix, [0.0, 0.0])
            mass_and_signed[0] += probability
            mass_and_signed[1] += probability * float(outcome[0])
        second_moment = 0.0
        entropy = 0.0
        for mass, signed_mass in cells.values():
            magnetization = signed_mass / mass
            second_moment += mass * magnetization * magnetization
            entropy += mass * binary_entropy_bits((1.0 + magnetization) / 2.0)
        second_moments.append(second_moment)
        entropies.append(entropy)

    l2_gains = tuple(
        second_moments[index + 1] - second_moments[index]
        for index in range(side_coordinates)
    )
    information_gains = tuple(
        entropies[index] - entropies[index + 1]
        for index in range(side_coordinates)
    )
    return BinaryFiltrationLedger(
        second_moments=tuple(second_moments),
        conditional_entropies_bits=tuple(entropies),
        l2_gains=l2_gains,
        information_gains_bits=information_gains,
    )


def logit_increment_audit(
    prior_logit: float, posterior_logit: float
) -> LogitIncrementAudit:
    """Audit the KL and magnetization change between two finite logits."""

    prior_probability = _logistic(prior_logit)
    posterior_probability = _logistic(posterior_logit)
    kl_nats = (
        posterior_probability * log(posterior_probability / prior_probability)
        + (1.0 - posterior_probability)
        * log((1.0 - posterior_probability) / (1.0 - prior_probability))
    )
    prior_magnetization = tanh(prior_logit / 2.0)
    posterior_magnetization = tanh(posterior_logit / 2.0)
    squared_change = (posterior_magnetization - prior_magnetization) ** 2
    squared_logit_change = (posterior_logit - prior_logit) ** 2
    return LogitIncrementAudit(
        kl_nats=kl_nats,
        squared_magnetization_change=squared_change,
        kl_quadratic_bound=squared_logit_change / 8.0,
        magnetization_quadratic_bound=squared_logit_change / 4.0,
        magnetization_pinsker_bound=2.0 * kl_nats,
    )


def censored_exponential_kl_nats(
    rate: float, reference_rate: float, horizon: float
) -> float:
    """KL between exponential clocks observed up to a censoring horizon.

    The output is the ringing time when it is at most ``horizon`` and a
    censoring symbol otherwise.  Zero rates are handled by their exact
    one-sided limits.
    """

    if rate < 0.0 or reference_rate < 0.0:
        raise ValueError("rates must be non-negative")
    if horizon <= 0.0:
        raise ValueError("horizon must be positive")
    if rate == reference_rate:
        return 0.0
    if rate == 0.0:
        return reference_rate * horizon
    if reference_rate == 0.0:
        return inf
    ringing_probability = -expm1(-rate * horizon)
    relative_difference = reference_rate / rate - 1.0
    if abs(relative_difference) < 1e-4:
        value = relative_difference
        contrast = value * value * (
            0.5
            - value / 3.0
            + value * value / 4.0
            - value**3 / 5.0
            + value**4 / 6.0
        )
    else:
        contrast = relative_difference - log1p(relative_difference)
    return ringing_probability * contrast


def binary_censored_exponential_information_bound_nats(
    prior_plus: float, rate_plus: float, rate_minus: float, horizon: float
) -> float:
    """Jeffreys upper bound for a binary censored-exponential channel."""

    if not 0.0 <= prior_plus <= 1.0:
        raise ValueError("prior_plus must belong to [0,1]")
    jeffreys = censored_exponential_kl_nats(
        rate_plus, rate_minus, horizon
    ) + censored_exponential_kl_nats(rate_minus, rate_plus, horizon)
    if prior_plus in (0.0, 1.0):
        return 0.0
    return prior_plus * (1.0 - prior_plus) * jeffreys


def _validate_two_point_mixture(
    rates: Sequence[float], weights: Sequence[float]
) -> None:
    if len(rates) != 2 or len(weights) != 2:
        raise ValueError("rates and weights must each have length two")
    if any(rate < 0.0 for rate in rates):
        raise ValueError("rates must be non-negative")
    if any(weight < 0.0 for weight in weights):
        raise ValueError("weights must be non-negative")
    if abs(sum(weights) - 1.0) > 1e-12:
        raise ValueError("weights must sum to one")


def _optimal_two_by_two_transport_cost(
    left_weights: Sequence[float],
    right_weights: Sequence[float],
    costs: Sequence[Sequence[float]],
) -> float:
    """Exact optimal transport cost between two two-point laws."""

    left_zero = left_weights[0]
    right_zero = right_weights[0]
    lower = max(0.0, left_zero + right_zero - 1.0)
    upper = min(left_zero, right_zero)

    def objective(gamma_zero_zero: float) -> float:
        coupling = (
            (gamma_zero_zero, left_zero - gamma_zero_zero),
            (
                right_zero - gamma_zero_zero,
                1.0 - left_zero - right_zero + gamma_zero_zero,
            ),
        )
        total = 0.0
        for row in range(2):
            for column in range(2):
                mass = coupling[row][column]
                if abs(mass) < 1e-15:
                    continue
                if mass < 0.0:
                    raise ArithmeticError("negative transport mass")
                total += mass * costs[row][column]
        return total

    return min(objective(lower), objective(upper))


def four_rate_censored_information_bound_nats(
    prior_even: float,
    even_rates: Sequence[float],
    even_weights: Sequence[float],
    odd_rates: Sequence[float],
    odd_weights: Sequence[float],
    horizon: float,
) -> float:
    """Transport-KL bound for two parity mixtures with four clock rates."""

    if not 0.0 <= prior_even <= 1.0:
        raise ValueError("prior_even must belong to [0,1]")
    if horizon <= 0.0:
        raise ValueError("horizon must be positive")
    _validate_two_point_mixture(even_rates, even_weights)
    _validate_two_point_mixture(odd_rates, odd_weights)
    if prior_even in (0.0, 1.0):
        return 0.0

    forward_costs = tuple(
        tuple(
            censored_exponential_kl_nats(even_rate, odd_rate, horizon)
            for odd_rate in odd_rates
        )
        for even_rate in even_rates
    )
    reverse_costs = tuple(
        tuple(
            censored_exponential_kl_nats(odd_rate, even_rate, horizon)
            for even_rate in even_rates
        )
        for odd_rate in odd_rates
    )
    forward = _optimal_two_by_two_transport_cost(
        even_weights, odd_weights, forward_costs
    )
    reverse = _optimal_two_by_two_transport_cost(
        odd_weights, even_weights, reverse_costs
    )
    return prior_even * (1.0 - prior_even) * (forward + reverse)


def satisfied_edge_open_probability(p: float, time: float) -> float:
    """Probability that an Exp(log(p/(1-p))) clock rings by ``time``."""

    if time < 0.0:
        raise ValueError("time must be non-negative")
    return 1.0 - exp(-time * coupling(p))


def triangle_full_connection_probability(
    satisfied_edges: int, p: float, time: float
) -> float:
    """Full-tree probability on a triangle with a fixed satisfied-edge count."""

    if satisfied_edges not in (0, 1, 2, 3):
        raise ValueError("satisfied_edges must belong to {0,1,2,3}")
    opened = satisfied_edge_open_probability(p, time)
    if satisfied_edges < 2:
        return 0.0
    if satisfied_edges == 2:
        return opened * opened
    return 3.0 * opened * opened - 2.0 * opened**3


def _high_count_masses(p: float) -> tuple[float, float, float]:
    coupling(p)
    syndrome_plus = face_syndrome_probability(p)
    q = 1.0 - p
    all_three_given_plus = p**3 / syndrome_plus
    exactly_two_given_minus = 3.0 * p * p * q / (1.0 - syndrome_plus)
    return syndrome_plus, all_three_given_plus, exactly_two_given_minus


def triangle_connection_information_bits(p: float, time: float) -> float:
    """Return I(Z; C_t | S), where C_t is full connection by ``time``."""

    syndrome_plus, mass_plus, mass_minus = _high_count_masses(p)
    connection_three = triangle_full_connection_probability(3, p, time)
    connection_two = triangle_full_connection_probability(2, p, time)
    information_plus = (
        binary_entropy_bits(mass_plus * connection_three)
        - mass_plus * binary_entropy_bits(connection_three)
    )
    information_minus = (
        binary_entropy_bits(mass_minus * connection_two)
        - mass_minus * binary_entropy_bits(connection_two)
    )
    return (
        syndrome_plus * information_plus
        + (1.0 - syndrome_plus) * information_minus
    )


def _triangle_edge_information_bits(
    p: float, time: float, connected_entropy_minus: float
) -> float:
    syndrome_plus, mass_plus, mass_minus = _high_count_masses(p)
    connection_three = triangle_full_connection_probability(3, p, time)
    connection_two = triangle_full_connection_probability(2, p, time)

    target_plus_given_plus = (1.0 + 2.0 * mass_plus) / 3.0
    null_plus_mass = 1.0 - mass_plus * connection_three
    target_plus_and_null = (
        mass_plus * (1.0 - connection_three) + (1.0 - mass_plus) / 3.0
    )
    information_plus = binary_entropy_bits(target_plus_given_plus) - (
        null_plus_mass
        * binary_entropy_bits(target_plus_and_null / null_plus_mass)
    )

    target_plus_given_minus = 2.0 * mass_minus / 3.0
    null_minus_mass = 1.0 - mass_minus * connection_two
    target_plus_and_null = (
        2.0 * mass_minus * (1.0 - connection_two) / 3.0
    )
    information_minus = (
        binary_entropy_bits(target_plus_given_minus)
        - mass_minus * connection_two * connected_entropy_minus
        - null_minus_mass
        * binary_entropy_bits(target_plus_and_null / null_minus_mass)
    )
    return (
        syndrome_plus * information_plus
        + (1.0 - syndrome_plus) * information_minus
    )


def triangle_edge_connection_information_bits(p: float, time: float) -> float:
    """Return I(Z_1; C_t | S) for one fixed edge of the triangle."""

    return _triangle_edge_information_bits(
        p, time, binary_entropy_bits(2.0 / 3.0)
    )


def triangle_connected_first_merge_information_bits(
    p: float, time: float
) -> float:
    """Return I(Z; Y_t | S) for a coarse function of the non-marked dendrogram.

    ``Y_t`` is null unless the triangle is fully connected by ``time``.  On
    full connection it records the first merged singleton pair.  Partial
    forests are deliberately collapsed, so this is a lower bound on the
    information carried by the complete dendrogram.
    """

    syndrome_plus, mass_plus, mass_minus = _high_count_masses(p)
    connection_three = triangle_full_connection_probability(3, p, time)
    connection_two = triangle_full_connection_probability(2, p, time)
    information_plus = (
        binary_entropy_bits(mass_plus * connection_three)
        - mass_plus * binary_entropy_bits(connection_three)
    )
    information_minus = (
        binary_entropy_bits(mass_minus * connection_two)
        + mass_minus * connection_two * log2(3.0)
        - mass_minus
        * (binary_entropy_bits(connection_two) + connection_two)
    )
    return (
        syndrome_plus * information_plus
        + (1.0 - syndrome_plus) * information_minus
    )


def triangle_edge_connected_first_merge_information_bits(
    p: float, time: float
) -> float:
    """Return I(Z_1; Y_t | S) for one fixed edge of the triangle."""

    return _triangle_edge_information_bits(p, time, 2.0 / 3.0)


def _triangle_edge_l2_gain(
    p: float, time: float, connected_second_moment_minus: float
) -> float:
    syndrome_plus, mass_plus, mass_minus = _high_count_masses(p)
    connection_three = triangle_full_connection_probability(3, p, time)
    connection_two = triangle_full_connection_probability(2, p, time)

    target_plus_given_plus = (1.0 + 2.0 * mass_plus) / 3.0
    null_plus_mass = 1.0 - mass_plus * connection_three
    target_plus_and_null = (
        mass_plus * (1.0 - connection_three) + (1.0 - mass_plus) / 3.0
    )
    null_plus_magnetization = (
        2.0 * target_plus_and_null / null_plus_mass - 1.0
    )

    target_plus_given_minus = 2.0 * mass_minus / 3.0
    null_minus_mass = 1.0 - mass_minus * connection_two
    target_plus_and_null = (
        2.0 * mass_minus * (1.0 - connection_two) / 3.0
    )
    null_minus_magnetization = (
        2.0 * target_plus_and_null / null_minus_mass - 1.0
    )

    baseline = (
        syndrome_plus * (2.0 * target_plus_given_plus - 1.0) ** 2
        + (1.0 - syndrome_plus)
        * (2.0 * target_plus_given_minus - 1.0) ** 2
    )
    terminal_plus = (
        mass_plus * connection_three
        + null_plus_mass * null_plus_magnetization**2
    )
    terminal_minus = (
        mass_minus * connection_two * connected_second_moment_minus
        + null_minus_mass * null_minus_magnetization**2
    )
    terminal = (
        syndrome_plus * terminal_plus
        + (1.0 - syndrome_plus) * terminal_minus
    )
    return terminal - baseline


def triangle_edge_connection_l2_gain(p: float, time: float) -> float:
    """L2 posterior gain when only full connection is revealed."""

    return _triangle_edge_l2_gain(p, time, 1.0 / 9.0)


def triangle_edge_connected_first_merge_l2_gain(
    p: float, time: float
) -> float:
    """L2 posterior gain when connection and the first merge are revealed."""

    return _triangle_edge_l2_gain(p, time, 1.0 / 3.0)


def triangle_palm_negative_syndrome_probability(p: float, time: float) -> float:
    """Syndrome-minus mass under the Palm law of a second merge at ``time``."""

    if time < 0.0:
        raise ValueError("time must be non-negative")
    q = 1.0 - p
    return q / (q + p * exp(-time * coupling(p)))


def triangle_palm_dendrogram_entropy_bits(p: float, time: float) -> float:
    """Average posterior entropy under the full critical-merge Palm oracle.

    The observation includes the face syndrome.  Given a full non-marked
    dendrogram with second merge at ``time``, the plus syndrome selects the
    unique all-satisfied state, while the minus syndrome leaves two equally
    likely states compatible with the first merged pair.  Their entropies are
    respectively zero and one bit.
    """

    return triangle_palm_negative_syndrome_probability(p, time)


def triangle_critical_palm_entropy_bits(p: float) -> float:
    """Closed form of the Palm entropy at the triangular percolation time."""

    beta_critical(p)
    return (1.0 - p) / (1.0 - Q_CRITICAL)


def main() -> None:
    p = conjectured_nishimori_root()
    critical_time = beta_critical(p)
    print(f"p_N^(0)                              = {p:.12f}")
    print(f"beta_c                               = {critical_time:.12f}")
    print(
        "I(Z; connection by beta_c | S)     = "
        f"{triangle_connection_information_bits(p, critical_time):.12f} bits"
    )
    print(
        "I(Z; connected first merge | S)    = "
        f"{triangle_connected_first_merge_information_bits(p, critical_time):.12f} bits"
    )
    print(
        "I(Z_1; connection by beta_c | S)   = "
        f"{triangle_edge_connection_information_bits(p, critical_time):.12f} bits"
    )
    print(
        "I(Z_1; connected first merge | S)  = "
        f"{triangle_edge_connected_first_merge_information_bits(p, critical_time):.12f} bits"
    )
    print(
        "L2 gain, connection only            = "
        f"{triangle_edge_connection_l2_gain(p, critical_time):.12f}"
    )
    print(
        "L2 gain, connected first merge      = "
        f"{triangle_edge_connected_first_merge_l2_gain(p, critical_time):.12f}"
    )
    print(
        "Palm entropy at second merge beta_c = "
        f"{triangle_critical_palm_entropy_bits(p):.12f} bits"
    )
    information_threshold = (1.0 + sqrt(Q_CRITICAL)) / 2.0
    information_critical_time = beta_critical(information_threshold)
    print("\ncomparison with the information-percolation threshold")
    print("p                 pair I(Y_beta|S)   pair L2 gain     Palm entropy")
    for parameter, time in (
        (information_threshold, information_critical_time),
        (p, critical_time),
    ):
        print(
            f"{parameter:.12f}    "
            f"{triangle_edge_connected_first_merge_information_bits(parameter, time):.12f}    "
            f"{triangle_edge_connected_first_merge_l2_gain(parameter, time):.12f}    "
            f"{triangle_critical_palm_entropy_bits(parameter):.12f}"
        )


if __name__ == "__main__":
    main()
