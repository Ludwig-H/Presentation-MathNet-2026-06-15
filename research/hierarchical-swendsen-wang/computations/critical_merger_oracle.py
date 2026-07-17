"""Exact local calculus for a Kruskal merger at the percolation threshold.

The module is dependency-free.  It deliberately computes the reliability in
two ways:

* from the hierarchical local log-likelihood ratio;
* directly from the two mirrored binary experiments.

Their agreement is an independent audit of the latent-winner correction.
All results are conditional on an unlabelled critical bucket and set the
ancestral message B_u to zero.  They are not a non-oracle recovery threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atanh, exp, fsum, inf, lgamma, log, log1p, pi, sqrt, tanh

from critical_band_thresholds import (
    Q_CRITICAL,
    beta_critical,
    open_probability,
)


P_SW = (1.0 + Q_CRITICAL) / 2.0
P_INFO = (1.0 + sqrt(Q_CRITICAL)) / 2.0


@dataclass(frozen=True)
class CriticalParameters:
    """Closed-edge parameters evaluated at q_p(beta_c) = q_c."""

    p: float
    beta: float
    satisfaction: float
    bias: float
    residual_log_odds: float


def _validate_admissible_p(p: float) -> None:
    if not P_SW <= p < 1.0:
        raise ValueError("p must satisfy p_SW <= p < 1")


def critical_parameters(p: float) -> CriticalParameters:
    """Return beta_c, s_c, h_c, and a_c using their closed forms."""

    _validate_admissible_p(p)
    beta = beta_critical(p)
    bias = 2.0 * (p - P_SW) / (1.0 - Q_CRITICAL)
    satisfaction = (1.0 + bias) / 2.0
    residual_log_odds = log((p - Q_CRITICAL) / (1.0 - p))
    if abs(residual_log_odds) < 1e-15:
        residual_log_odds = 0.0
    return CriticalParameters(
        p=p,
        beta=beta,
        satisfaction=satisfaction,
        bias=bias,
        residual_log_odds=residual_log_odds,
    )


def critical_plus_count_pmf(size: int, p: float) -> dict[int, float]:
    """Law of K = 1 + Bin(size-1, s_c) under the true parity."""

    if size <= 0:
        raise ValueError("size must be positive")
    s_c = critical_parameters(p).satisfaction
    result: dict[int, float] = {}
    for count in range(1, size + 1):
        successes = count - 1
        failures = size - count
        log_mass = (
            lgamma(size)
            - lgamma(successes + 1)
            - lgamma(failures + 1)
            + successes * log(s_c)
            + failures * log1p(-s_c)
        )
        result[count] = exp(log_mass)
    return result


def critical_minus_count_pmf(size: int, p: float) -> dict[int, float]:
    """Mirrored count law under the opposite relative parity."""

    plus = critical_plus_count_pmf(size, p)
    return {count: plus.get(size - count, 0.0) for count in range(size + 1)}


def critical_local_log_odds(size: int, count: int, p: float) -> float:
    """Return log(P_+(K=count) / P_-(K=count))."""

    if size <= 0:
        raise ValueError("size must be positive")
    if not 0 <= count <= size:
        raise ValueError("count must belong to {0, ..., size}")
    if count == 0:
        return -inf
    if count == size:
        return inf
    a_c = critical_parameters(p).residual_log_odds
    return log(count / (size - count)) + a_c * (2 * count - size)


def _squared_posterior_bias(log_odds: float) -> float:
    if log_odds in (-inf, inf):
        return 1.0
    value = tanh(log_odds / 2.0)
    return value * value


def posterior_same_parity_probability(log_odds: float) -> float:
    """Return the heat-bath probability of the even pair of flip states."""

    if log_odds == inf:
        return 1.0
    if log_odds == -inf:
        return 0.0
    if log_odds >= 0.0:
        return 1.0 / (1.0 + exp(-log_odds))
    exponential = exp(log_odds)
    return exponential / (1.0 + exponential)


def critical_same_parity_probability_given_count(
    size: int, count: int, p: float
) -> float:
    """Return P((a,b) is (0,0) or (1,1) | K=count) for B_u=0."""

    return posterior_same_parity_probability(
        critical_local_log_odds(size, count, p)
    )


def critical_mean_opposite_parity_probability(size: int, p: float) -> float:
    """Return the mean local probability of choosing an odd flip parity.

    Computing the deficit directly avoids subtracting two floating-point
    numbers close to one when the critical bucket is large.
    """

    plus = critical_plus_count_pmf(size, p)
    return fsum(
        mass
        * (
            1.0
            - critical_same_parity_probability_given_count(
                size, count, p
            )
        )
        for count, mass in plus.items()
    )


def critical_mean_same_parity_probability(size: int, p: float) -> float:
    """Return the mean local probability of an even flip parity."""

    return 1.0 - critical_mean_opposite_parity_probability(size, p)


def critical_bayes_error(size: int, p: float) -> float:
    """Bayes error for the two mirrored critical count experiments."""

    plus = critical_plus_count_pmf(size, p)
    error = 0.0
    for count, mass in plus.items():
        log_odds = critical_local_log_odds(size, count, p)
        if log_odds < 0.0:
            error += mass
        elif log_odds == 0.0:
            error += 0.5 * mass
    return error


def critical_mean_error_exponent(p: float) -> float:
    """Return D(1/2 || s_c) = -log(1-h_c^2)/2."""

    bias = critical_parameters(p).bias
    return -0.5 * log(1.0 - bias * bias)


def _hyperbolic_secant(value: float) -> float:
    """Evaluate sech(value) without overflowing for large arguments."""

    exponential = exp(-abs(value))
    return 2.0 * exponential / (1.0 + exponential * exponential)


def critical_mean_error_prefactor(p: float, size_parity: int) -> float:
    """Return the sharp m^(-1/2) prefactor on an even/odd subsequence.

    For fixed ``p > P_SW`` and ``m`` with parity ``size_parity``, the mean
    opposite-parity probability is asymptotic to

        prefactor * exp(-m * critical_mean_error_exponent(p)) / sqrt(m).

    The defining shifted-sech series is evaluated directly for large
    residual log-odds and through Poisson summation for small log-odds.
    """

    parameters = critical_parameters(p)
    if parameters.residual_log_odds <= 0.0:
        raise ValueError("the sharp prefactor is defined only for p > p_SW")
    if size_parity not in (0, 1):
        raise ValueError("size_parity must be 0 (even) or 1 (odd)")

    a_c = parameters.residual_log_odds
    tolerance = 1e-16

    if a_c <= pi:
        # Poisson summation:
        # sum_j sech(a(j+eps)) = (pi/a) sum_k
        # exp(2 pi i k eps) sech(pi^2 k/a).
        transformed = 1.0
        frequency = 1
        while True:
            term = _hyperbolic_secant(pi * pi * frequency / a_c)
            sign = -1.0 if size_parity == 1 and frequency % 2 else 1.0
            transformed += 2.0 * sign * term
            if term < tolerance:
                break
            frequency += 1
        shifted_series = pi * transformed / a_c
    elif size_parity == 0:
        shifted_series = 1.0
        index = 1
        while True:
            term = _hyperbolic_secant(a_c * index)
            shifted_series += 2.0 * term
            if term < tolerance:
                break
            index += 1
    else:
        shifted_series = 0.0
        index = 0
        while True:
            term = _hyperbolic_secant(a_c * (index + 0.5))
            shifted_series += 2.0 * term
            if term < tolerance:
                break
            index += 1

    return shifted_series / (
        2.0 * parameters.satisfaction * sqrt(2.0 * pi)
    )


def critical_oracle_reliability(size: int, p: float) -> float:
    """Return Gamma_m^c from the local hierarchical likelihood ratio."""

    plus = critical_plus_count_pmf(size, p)
    return fsum(
        mass * _squared_posterior_bias(critical_local_log_odds(size, count, p))
        for count, mass in plus.items()
    )


def symmetric_experiment_reliability(size: int, p: float) -> float:
    """Compute the same L2 contraction directly from P_+ and P_-."""

    plus_sparse = critical_plus_count_pmf(size, p)
    minus = critical_minus_count_pmf(size, p)
    terms = []
    for count in range(size + 1):
        plus_mass = plus_sparse.get(count, 0.0)
        minus_mass = minus[count]
        total = plus_mass + minus_mass
        if total > 0.0:
            terms.append((plus_mass - minus_mass) ** 2 / total)
    return 0.5 * fsum(terms)


def reliability_deficit_bound(size: int, p: float) -> float:
    """Return the proved Hoeffding upper bound on 1 - Gamma_m^c."""

    if size < 2:
        raise ValueError("the exponential bound is stated for size >= 2")
    parameters = critical_parameters(p)
    h_c = parameters.bias
    a_c = parameters.residual_log_odds
    bound = exp(-(size - 1) * h_c * h_c / 8.0) + 4.0 * exp(
        -(size - 1) * a_c * h_c / 2.0
    )
    return min(1.0, bound)


def critical_scaling_p(size: int, alpha: float) -> float:
    """Return p_SW + (1-q_c) alpha / (2 sqrt(size))."""

    if size <= 0:
        raise ValueError("size must be positive")
    if alpha < 0.0:
        raise ValueError("alpha must be nonnegative")
    p = P_SW + (1.0 - Q_CRITICAL) * alpha / (2.0 * sqrt(size))
    if p >= 1.0:
        raise ValueError("size is too small for this alpha")
    return p


def gaussian_crossover(alpha: float, intervals: int = 20000) -> float:
    """Numerically integrate E[tanh^2(alpha Z + alpha^2)].

    Composite Simpson integration on [-9, 9] is used only as a numerical
    audit of the proved limit.  It is not part of the proof.
    """

    if alpha < 0.0:
        raise ValueError("alpha must be nonnegative")
    if intervals <= 0 or intervals % 2:
        raise ValueError("intervals must be a positive even integer")
    if alpha == 0.0:
        return 0.0

    radius = 9.0
    step = 2.0 * radius / intervals
    normalizer = sqrt(2.0 * pi)

    def integrand(z: float) -> float:
        value = tanh(alpha * z + alpha * alpha)
        return value * value * exp(-z * z / 2.0) / normalizer

    total = integrand(-radius) + integrand(radius)
    total += 4.0 * fsum(
        integrand(-radius + index * step)
        for index in range(1, intervals, 2)
    )
    total += 2.0 * fsum(
        integrand(-radius + index * step)
        for index in range(2, intervals, 2)
    )
    return total * step / 3.0


def _print_parameter_table() -> None:
    print("critical parameters")
    print("p              beta_c        s_c           h_c           a_c")
    for p in (P_SW, P_INFO, 0.8358058):
        values = critical_parameters(p)
        assert abs(open_probability(p, values.beta) - Q_CRITICAL) < 1e-12
        assert abs(values.residual_log_odds - 2.0 * atanh(values.bias)) < 1e-12
        print(
            f"{p:.10f}   {values.beta:.10f}   "
            f"{values.satisfaction:.10f}   {values.bias:.10f}   "
            f"{values.residual_log_odds:.10f}"
        )


def _print_reliability_table() -> None:
    print("\ncritical local oracle reliability Gamma_m^c")
    print("m        p_SW           p_info         0.8358058")
    for size in (1, 4, 16, 64, 256):
        values = (
            critical_oracle_reliability(size, P_SW),
            critical_oracle_reliability(size, P_INFO),
            critical_oracle_reliability(size, 0.8358058),
        )
        print(
            f"{size:<8d} "
            f"{values[0]:.10f}   {values[1]:.10f}   {values[2]:.10f}"
        )


def _print_scaling_table() -> None:
    print("\nfinite-size crossover")
    print("alpha     limit          m=64          m=256         m=1024")
    for alpha in (0.5, 1.0, 2.0):
        limit = gaussian_crossover(alpha)
        finite = []
        for size in (64, 256, 1024):
            p = critical_scaling_p(size, alpha)
            finite.append(critical_oracle_reliability(size, p))
        print(
            f"{alpha:<9.1f} {limit:.10f}   {finite[0]:.10f}   "
            f"{finite[1]:.10f}   {finite[2]:.10f}"
        )


def main() -> None:
    _print_parameter_table()
    _print_reliability_table()
    _print_scaling_table()


if __name__ == "__main__":
    main()
