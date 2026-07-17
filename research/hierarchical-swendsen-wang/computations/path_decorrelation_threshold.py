"""Decorrelation thresholds for the factorized hierarchical path channel.

The exact factorized-path identity is

    P(same relative orientation) = (1 + exp(-A)) / 2,
    A = sum_w -log(Gamma_w).

This module evaluates that attenuation, the logarithmic-cut transition on
the triangular lattice, and the high-p window for fixed bucket sizes.  It
does not identify the factorized oracle with the joint hierarchical sweep.
"""

from __future__ import annotations

from math import comb, exp, fsum, inf, log, log1p, sqrt
from typing import Iterable

from critical_band_thresholds import Q_CRITICAL
from critical_merger_oracle import (
    P_SW,
    critical_mean_error_exponent,
    critical_mean_error_prefactor,
    critical_minus_count_pmf,
    critical_parameters,
    critical_plus_count_pmf,
)
from hierarchical_flip_probabilities import (
    node_mean_error_exponent,
    node_oracle_reliability,
)
from nishimori_hierarchical_entropy import conjectured_nishimori_root


def attenuation_from_reliabilities(
    reliabilities: Iterable[float],
) -> float:
    """Return A = -sum log(Gamma), allowing the exact endpoints 0 and 1."""

    attenuation = 0.0
    for reliability in reliabilities:
        if not 0.0 <= reliability <= 1.0:
            raise ValueError("every reliability must belong to [0, 1]")
        if reliability == 0.0:
            return inf
        attenuation -= log(reliability)
    return attenuation


def same_relation_probability_from_attenuation(
    attenuation: float,
) -> float:
    """Return (1 + exp(-A))/2 from a nonnegative attenuation A."""

    if attenuation < 0.0:
        raise ValueError("attenuation must be nonnegative")
    if attenuation == inf:
        return 0.5
    return 0.5 * (1.0 + exp(-attenuation))


def regular_path_same_relation_probability(
    length: int,
    reliability: float,
) -> float:
    """Same-relation probability for ``length`` identical PATH-FAC nodes."""

    if length < 0:
        raise ValueError("length must be nonnegative")
    if not 0.0 <= reliability <= 1.0:
        raise ValueError("reliability must belong to [0, 1]")
    if length == 0:
        return 1.0
    if reliability == 0.0:
        return 0.5
    return same_relation_probability_from_attenuation(
        -length * log(reliability)
    )


def required_attenuation_for_accuracy(tolerance: float) -> float:
    """Attenuation needed to be within ``tolerance`` of random parity."""

    if not 0.0 < tolerance < 0.5:
        raise ValueError("tolerance must belong to (0, 1/2)")
    return log(1.0 / (2.0 * tolerance))


def critical_reliability_deficit(size: int, p: float) -> float:
    """Compute 1-Gamma_m^c stably through the two mirrored experiments."""

    if size <= 0:
        raise ValueError("size must be positive")
    plus = critical_plus_count_pmf(size, p)
    minus = critical_minus_count_pmf(size, p)
    harmonic_terms = []
    for count in range(size + 1):
        plus_mass = plus.get(count, 0.0)
        minus_mass = minus[count]
        total = plus_mass + minus_mass
        if total > 0.0:
            harmonic_terms.append(plus_mass * minus_mass / total)
    return 2.0 * fsum(harmonic_terms)


def regular_critical_path_same_relation_probability(
    length: int,
    size: int,
    p: float,
) -> float:
    """Exact PATH-FAC probability for identical critical buckets."""

    if length < 0:
        raise ValueError("length must be nonnegative")
    deficit = critical_reliability_deficit(size, p)
    if deficit == 1.0:
        return 0.5 if length > 0 else 1.0
    attenuation = -length * log1p(-deficit)
    return same_relation_probability_from_attenuation(attenuation)


def heterogeneous_critical_path_same_relation_probability(
    sizes: Iterable[int],
    p: float,
) -> float:
    """Exact PATH-FAC probability for a heterogeneous critical path."""

    attenuation = 0.0
    for size in sizes:
        deficit = critical_reliability_deficit(size, p)
        if deficit == 1.0:
            return 0.5
        attenuation -= log1p(-deficit)
    return same_relation_probability_from_attenuation(attenuation)


def heterogeneous_descendant_path_same_relation_probability(
    sizes: Iterable[int],
    times: Iterable[float],
    p: float,
) -> float:
    """Exact PATH-FAC probability for sizes and descendant clock levels."""

    size_values = tuple(sizes)
    time_values = tuple(times)
    if len(size_values) != len(time_values):
        raise ValueError("sizes and times must have the same length")
    beta = critical_parameters(p).beta
    attenuation = 0.0
    for size, time in zip(size_values, time_values, strict=True):
        if not 0.0 <= time <= beta:
            raise ValueError("every descendant time must belong to [0, beta_c]")
        reliability = node_oracle_reliability(size, p, time)
        if reliability == 0.0:
            return 0.5
        attenuation -= log(reliability)
    return same_relation_probability_from_attenuation(attenuation)


def critical_geometry_log_partition(
    sizes: Iterable[int],
    p: float,
) -> float:
    """Return log sum m^(-1/2) exp(-m I_c(p)) stably.

    The empty path has partition function zero and logarithm ``-inf``.
    """

    exponent = critical_mean_error_exponent(p)
    log_terms = []
    for size in sizes:
        if size <= 0:
            raise ValueError("every size must be positive")
        log_terms.append(-0.5 * log(size) - size * exponent)
    if not log_terms:
        return -inf
    maximum = max(log_terms)
    return maximum + log(fsum(exp(term - maximum) for term in log_terms))


def descendant_geometry_log_partition(
    sizes: Iterable[int],
    times: Iterable[float],
    p: float,
) -> float:
    """Return log Phi_desc for a path below its critical LCA."""

    size_values = tuple(sizes)
    time_values = tuple(times)
    if len(size_values) != len(time_values):
        raise ValueError("sizes and times must have the same length")
    beta = critical_parameters(p).beta
    log_terms = []
    for size, time in zip(size_values, time_values, strict=True):
        if size <= 0:
            raise ValueError("every size must be positive")
        if not 0.0 <= time <= beta:
            raise ValueError("every descendant time must belong to [0, beta_c]")
        exponent = node_mean_error_exponent(p, time)
        log_terms.append(-0.5 * log(size) - size * exponent)
    if not log_terms:
        return -inf
    maximum = max(log_terms)
    return maximum + log(fsum(exp(term - maximum) for term in log_terms))


def critical_p_from_error_exponent(exponent: float) -> float:
    """Invert I_c(p)=exponent on [p_SW, 1]."""

    if exponent < 0.0:
        raise ValueError("exponent must be nonnegative")
    if exponent == inf:
        return 1.0
    critical_bias = sqrt(1.0 - exp(-2.0 * exponent))
    return 0.5 * (
        1.0
        + Q_CRITICAL
        + (1.0 - Q_CRITICAL) * critical_bias
    )


def relative_level_error_exponent(p: float, theta: float) -> float:
    """Return I(theta beta_c(p); p) for a relative descendant level."""

    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must belong to [0, 1]")
    beta = critical_parameters(p).beta
    return node_mean_error_exponent(p, theta * beta)


def regular_relative_level_threshold(
    alpha: float,
    theta: float,
) -> float | None:
    """Solve alpha I_theta(p)=1 above p_SW, if an interior root exists."""

    if alpha <= 0.0:
        raise ValueError("alpha must be positive")
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must belong to [0, 1]")
    target = 1.0 / alpha
    if relative_level_error_exponent(P_SW, theta) >= target:
        return None

    lower = P_SW
    upper = 1.0 - 2.0**-52
    if relative_level_error_exponent(upper, theta) < target:
        return 1.0
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        if relative_level_error_exponent(midpoint, theta) < target:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def regular_log_cut_threshold(alpha: float) -> float:
    """Solve alpha * I_c(p) = 1 for the triangular critical channel."""

    if alpha <= 0.0:
        raise ValueError("alpha must be positive")
    return critical_p_from_error_exponent(1.0 / alpha)


def critical_log_cut_coordinate(
    length: int,
    size: int,
    p: float,
) -> float:
    """Return z=m I_c(p)-log(H)+log(m)/2 for the sharp path window."""

    if length <= 0 or size <= 0:
        raise ValueError("length and size must be positive")
    return (
        size * critical_mean_error_exponent(p)
        - log(length)
        + 0.5 * log(size)
    )


def critical_log_cut_limit_probability(
    p: float,
    size_parity: int,
    coordinate: float,
) -> float:
    """Sharp PATH-FAC limit at a finite logarithmic-cut coordinate."""

    prefactor = critical_mean_error_prefactor(p, size_parity)
    correlation = exp(-2.0 * prefactor * exp(-coordinate))
    return 0.5 * (1.0 + correlation)


def high_p_deficit_exponent(size: int) -> int:
    """Power d_m=ceil(m/2) in the fixed-bucket p-to-one window."""

    if size < 2:
        raise ValueError("the nontrivial fixed-size window requires size >= 2")
    return (size + 1) // 2


def high_p_deficit_constant(size: int) -> int:
    """Leading constant D_m in 1-Gamma_m^c ~ D_m epsilon^d_m."""

    if size < 2:
        raise ValueError("the nontrivial fixed-size window requires size >= 2")
    if size % 2 == 0:
        half = size // 2
        return comb(size - 1, half - 1)
    half = (size - 1) // 2
    return 4 * comb(size - 1, half - 1)


def high_p_critical_deficit_equivalent(size: int, p: float) -> float:
    """Leading fixed-size critical deficit as p tends to one."""

    if not P_SW <= p < 1.0:
        raise ValueError("p must satisfy p_SW <= p < 1")
    epsilon = (1.0 - p) / (1.0 - Q_CRITICAL)
    return high_p_deficit_constant(size) * epsilon ** high_p_deficit_exponent(
        size
    )


def high_p_correlation_length_equivalent(size: int, p: float) -> float:
    """Leading critical path correlation length for fixed bucket size."""

    deficit = high_p_critical_deficit_equivalent(size, p)
    return inf if deficit == 0.0 else 1.0 / deficit


def log_cut_coefficient_for_threshold(p: float) -> float:
    """Return alpha=1/I_c(p) for a requested first-order threshold."""

    if not P_SW < p < 1.0:
        raise ValueError("p must satisfy p_SW < p < 1")
    return 1.0 / critical_mean_error_exponent(p)


def nishimori_log_cut_coefficient() -> float:
    """Alpha for which the PATH-FAC first-order threshold equals p_N."""

    p_nishimori = conjectured_nishimori_root()
    return log_cut_coefficient_for_threshold(p_nishimori)


def main() -> None:
    """Print the two candidate decorrelation windows."""

    p_nishimori = conjectured_nishimori_root()
    p_information = 0.5 * (1.0 + sqrt(Q_CRITICAL))
    print(f"q_triangle = {Q_CRITICAL:.12f}")
    print(f"p_SW       = {P_SW:.12f}")
    print(f"p_N        = {p_nishimori:.12f}")
    print(
        "alpha matching p_N = "
        f"{nishimori_log_cut_coefficient():.12f}"
    )
    print(
        "alpha matching p_info = "
        f"{log_cut_coefficient_for_threshold(p_information):.12f}"
    )

    print("\nregular critical cuts m ~ alpha log H")
    print("alpha        p_path(alpha)")
    for alpha in (1.0, 2.0, 4.0, 7.0535961929, 10.0, 20.0):
        print(f"{alpha:10.6f}   {regular_log_cut_threshold(alpha):.12f}")

    print("\nrelative descendant level, alpha matching p_N at theta=1")
    print("theta        p_path(alpha_N, theta)")
    alpha_nishimori = nishimori_log_cut_coefficient()
    for theta in (0.0, 0.25, 0.5, 0.75, 1.0):
        threshold = regular_relative_level_threshold(alpha_nishimori, theta)
        rendered = "no loss phase" if threshold is None else f"{threshold:.12f}"
        print(f"{theta:5.2f}        {rendered}")

    print("\nfixed critical bucket, p -> 1")
    print("m       d_m       D_m")
    for size in range(2, 9):
        print(
            f"{size:1d}       {high_p_deficit_exponent(size):3d}       "
            f"{high_p_deficit_constant(size):3d}"
        )


if __name__ == "__main__":
    main()
