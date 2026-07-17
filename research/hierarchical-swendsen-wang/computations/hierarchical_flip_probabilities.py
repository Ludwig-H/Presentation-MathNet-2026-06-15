"""Exact flip probabilities from roots to internal nodes and path oracles.

The module keeps three objects separate:

* a two-state heat bath (root orientation or one-site/leaf orbit);
* the genuine four-state heat bath at an internal merger;
* an explicitly factorized path oracle, which is not the joint hierarchical
  sweep unless an independence theorem is supplied;
* an exact finite-state twisted transfer that retains path dependence.

All functions are dependency-free and use stable log-domain calculations.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, exp, fsum, inf, lgamma, log, log1p, prod, tanh
from typing import Hashable, Iterable, Mapping, Sequence


@dataclass(frozen=True)
class NodeClockParameters:
    """Residual mark parameters at a merger time."""

    p: float
    time: float
    coupling: float
    residual_log_odds: float
    satisfaction: float
    bias: float


@dataclass(frozen=True)
class WalshParameters:
    """Three nonconstant coefficients of four positive orbit weights."""

    field_left: float
    field_right: float
    coupling: float


@dataclass(frozen=True)
class FourStateProbabilities:
    """Probabilities of (0,0), (0,1), (1,0), and (1,1)."""

    p00: float
    p01: float
    p10: float
    p11: float

    @property
    def even(self) -> float:
        return self.p00 + self.p11

    @property
    def flip_left(self) -> float:
        return self.p10 + self.p11

    @property
    def flip_right(self) -> float:
        return self.p01 + self.p11


def _validate_probability_parameter(p: float) -> None:
    if not 0.5 < p < 1.0:
        raise ValueError("p must satisfy 1/2 < p < 1")


def node_clock_parameters(p: float, time: float) -> NodeClockParameters:
    """Return the residual Bernoulli channel at ``time`` in the GSBM."""

    _validate_probability_parameter(p)
    if not 0.0 <= time <= 1.0:
        raise ValueError("time must belong to [0, 1]")
    coupling = log(p / (1.0 - p))
    residual_log_odds = coupling * (1.0 - time)
    satisfaction = 1.0 / (1.0 + exp(-residual_log_odds))
    bias = 2.0 * satisfaction - 1.0
    return NodeClockParameters(
        p=p,
        time=time,
        coupling=coupling,
        residual_log_odds=residual_log_odds,
        satisfaction=satisfaction,
        bias=bias,
    )


def logistic(log_odds: float) -> float:
    """Stable logistic map on the extended real line."""

    if log_odds == inf:
        return 1.0
    if log_odds == -inf:
        return 0.0
    if log_odds >= 0.0:
        return 1.0 / (1.0 + exp(-log_odds))
    exponential = exp(log_odds)
    return exponential / (1.0 + exponential)


def root_flip_probability(log_prior_ratio: float = 0.0) -> float:
    """Two-state heat-bath flip probability of a final component."""

    return logistic(log_prior_ratio)


def leaf_heat_bath_flip_probability(log_weight_ratio: float) -> float:
    """Glauber/Barker flip probability on a one-site orbit."""

    return logistic(log_weight_ratio)


def leaf_metropolis_acceptance_probability(log_weight_ratio: float) -> float:
    """Metropolis acceptance for the deterministic one-site flip proposal."""

    return 1.0 if log_weight_ratio >= 0.0 else exp(log_weight_ratio)


def normalize_four_log_weights(
    log_weights: Sequence[float],
) -> FourStateProbabilities:
    """Normalize four log weights ordered as 00, 01, 10, 11."""

    if len(log_weights) != 4:
        raise ValueError("exactly four log weights are required")
    maximum = max(log_weights)
    if maximum == -inf:
        raise ValueError("at least one orbit weight must be positive")
    if maximum == inf:
        winners = [value == inf for value in log_weights]
        normalizer = sum(winners)
        probabilities = [float(winner) / normalizer for winner in winners]
        return FourStateProbabilities(*probabilities)
    masses = [
        0.0 if value == -inf else exp(value - maximum)
        for value in log_weights
    ]
    normalizer = fsum(masses)
    probabilities = [mass / normalizer for mass in masses]
    return FourStateProbabilities(*probabilities)


def walsh_parameters(log_weights: Sequence[float]) -> WalshParameters:
    """Return h_left, h_right, J for four strictly positive weights."""

    if len(log_weights) != 4:
        raise ValueError("exactly four log weights are required")
    if any(value in (-inf, inf) for value in log_weights):
        raise ValueError("Walsh parameters require finite log weights")
    g00, g01, g10, g11 = log_weights
    return WalshParameters(
        field_left=(g00 + g01 - g10 - g11) / 4.0,
        field_right=(g00 + g10 - g01 - g11) / 4.0,
        coupling=(g00 + g11 - g01 - g10) / 4.0,
    )


def four_state_probabilities_from_walsh(
    field_left: float, field_right: float, coupling: float
) -> FourStateProbabilities:
    """Return the internal-node heat bath from its three Walsh coefficients."""

    return normalize_four_log_weights(
        (
            field_left + field_right + coupling,
            field_left - field_right - coupling,
            -field_left + field_right - coupling,
            -field_left - field_right + coupling,
        )
    )


def _log_cosh(value: float) -> float:
    absolute = abs(value)
    return absolute + log1p(exp(-2.0 * absolute)) - log(2.0)


def even_log_odds_from_walsh(
    field_left: float, field_right: float, coupling: float
) -> float:
    """Log odds of (00 or 11) against (01 or 10)."""

    return (
        2.0 * coupling
        + _log_cosh(field_left + field_right)
        - _log_cosh(field_left - field_right)
    )


def internal_even_probability(
    field_left: float, field_right: float, coupling: float
) -> float:
    """Probability that the two child flips have even parity."""

    return logistic(
        even_log_odds_from_walsh(field_left, field_right, coupling)
    )


def node_plus_count_pmf(size: int, p: float, time: float) -> dict[int, float]:
    """Law of K=1+Bin(size-1,s_p(time)) under the true parity."""

    if size <= 0:
        raise ValueError("size must be positive")
    satisfaction = node_clock_parameters(p, time).satisfaction
    result: dict[int, float] = {}
    for count in range(1, size + 1):
        successes = count - 1
        failures = size - count
        log_mass = (
            lgamma(size)
            - lgamma(successes + 1)
            - lgamma(failures + 1)
            + successes * log(satisfaction)
            + failures * log1p(-satisfaction)
        )
        result[count] = exp(log_mass)
    return result


def node_local_log_odds(size: int, count: int, p: float, time: float) -> float:
    """Local parity LLR of an unmarked homogeneous merger bucket."""

    if size <= 0:
        raise ValueError("size must be positive")
    if not 0 <= count <= size:
        raise ValueError("count must belong to {0, ..., size}")
    if count == 0:
        return -inf
    if count == size:
        return inf
    residual = node_clock_parameters(p, time).residual_log_odds
    return log(count / (size - count)) + residual * (2 * count - size)


def node_oracle_reliability(size: int, p: float, time: float) -> float:
    """Mean squared posterior bias of the ancestor-neutral local bucket."""

    plus = node_plus_count_pmf(size, p, time)
    return fsum(
        mass * tanh(node_local_log_odds(size, count, p, time) / 2.0) ** 2
        for count, mass in plus.items()
    )


def node_mean_even_probability(size: int, p: float, time: float) -> float:
    """Mean heat-bath probability of preserving the true relative parity."""

    plus = node_plus_count_pmf(size, p, time)
    return fsum(
        mass * logistic(node_local_log_odds(size, count, p, time))
        for count, mass in plus.items()
    )


def node_strict_majority_probability(size: int, p: float, time: float) -> float:
    """Probability that the planted-winner bucket has a strict true majority."""

    if size <= 0:
        raise ValueError("size must be positive")
    satisfaction = node_clock_parameters(p, time).satisfaction
    first_success_count = size // 2
    return fsum(
        comb(size - 1, successes)
        * satisfaction**successes
        * (1.0 - satisfaction) ** (size - 1 - successes)
        for successes in range(first_success_count, size)
    )


def node_mean_error_exponent(p: float, time: float) -> float:
    """Return D(1/2 || s_p(time)) for the large-cut mean error."""

    bias = node_clock_parameters(p, time).bias
    return -0.5 * log(1.0 - bias * bias)


def factorized_path_same_parity_probability(
    signed_reliabilities: Iterable[float],
) -> float:
    """Same-parity probability for explicitly independent binary channels."""

    values = tuple(signed_reliabilities)
    if any(not -1.0 <= value <= 1.0 for value in values):
        raise ValueError("every signed reliability must belong to [-1, 1]")
    return 0.5 * (1.0 + prod(values))


def finite_state_path_correlation(
    initial_distribution: Mapping[Hashable, float],
    kernels: Sequence[
        Mapping[Hashable, Sequence[tuple[int, Hashable, float]]]
    ],
    separates_pair: Sequence[bool] | None = None,
) -> float:
    """Evaluate the exact twisted transfer recursion on a finite state.

    A transition is ``(decision, next_state, probability)``.  The decision
    belongs to ``{0, 1}``; it changes the pair parity when the corresponding
    value of ``separates_pair`` is true.  State dependence retains the joint
    correlations that the factorized path oracle deliberately discards.
    """

    if separates_pair is None:
        separates_pair = (True,) * len(kernels)
    if len(separates_pair) != len(kernels):
        raise ValueError("one path-incidence flag is required per kernel")
    if any(
        probability < 0.0
        for probability in initial_distribution.values()
    ):
        raise ValueError("initial probabilities must be nonnegative")
    if abs(fsum(initial_distribution.values()) - 1.0) > 1e-12:
        raise ValueError("initial probabilities must sum to one")

    signed_measure = dict(initial_distribution)
    for kernel, separates in zip(kernels, separates_pair):
        next_measure: dict[Hashable, float] = {}
        for state, signed_mass in signed_measure.items():
            if state not in kernel:
                raise ValueError(f"kernel is missing reachable state {state!r}")
            transitions = kernel[state]
            if any(
                decision not in (0, 1) or probability < 0.0
                for decision, _, probability in transitions
            ):
                raise ValueError("invalid decision or transition probability")
            total_probability = fsum(
                probability for _, _, probability in transitions
            )
            if abs(total_probability - 1.0) > 1e-12:
                raise ValueError("transition probabilities must sum to one")
            for decision, next_state, probability in transitions:
                sign = -1.0 if separates and decision == 1 else 1.0
                contribution = signed_mass * probability * sign
                next_measure[next_state] = (
                    next_measure.get(next_state, 0.0) + contribution
                )
        signed_measure = next_measure
    return fsum(signed_measure.values())


def symmetric_pair_status_probabilities(
    signed_correlation: float,
) -> FourStateProbabilities:
    """Joint conformity states after averaging a fair global orientation."""

    if not -1.0 <= signed_correlation <= 1.0:
        raise ValueError("signed_correlation must belong to [-1, 1]")
    same = (1.0 + signed_correlation) / 4.0
    different = (1.0 - signed_correlation) / 4.0
    return FourStateProbabilities(same, different, different, same)


def main() -> None:
    """Print a small reproducible table for the three hierarchy levels."""

    print("two-state endpoints")
    print(f"uniform root flip: {root_flip_probability():.6f}")
    log_ratio = log(3.0)
    print(
        "leaf with target ratio 3: "
        f"heat-bath={leaf_heat_bath_flip_probability(log_ratio):.6f}, "
        f"Metropolis={leaf_metropolis_acceptance_probability(log_ratio):.6f}"
    )
    print("\nlocal homogeneous node, p=0.8, m=8")
    print("time       majority       mean-even      reliability")
    for time in (0.0, 0.25, 0.5, 0.75, 1.0):
        print(
            f"{time:.2f}       "
            f"{node_strict_majority_probability(8, 0.8, time):.8f}     "
            f"{node_mean_even_probability(8, 0.8, time):.8f}     "
            f"{node_oracle_reliability(8, 0.8, time):.8f}"
        )


if __name__ == "__main__":
    main()
