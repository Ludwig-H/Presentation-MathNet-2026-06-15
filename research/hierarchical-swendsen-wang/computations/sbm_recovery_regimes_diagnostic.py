"""Exact scalar benchmarks for recovery in the symmetric two-class SBM.

Write the within- and between-community edge probabilities as ``a_n / n``
and ``b_n / n``.  This module deliberately separates three questions.

* Weak recovery at bounded degree is governed by the Kesten--Stigum number

      lambda_KS = (a_n - b_n)^2 / (2(a_n + b_n)) = d theta^2.

* A finite-``n`` oracle degree-profile test has an exact Bhattacharyya
  affinity.  Conditional on fixed reference labels for the other ``n-1``
  vertices, the two hypotheses swap Bernoulli probabilities ``a_n/n`` and
  ``b_n/n``.  If the reference groups have sizes ``m_+`` and ``m_-``, their
  two binomial counts have product affinity

      [sqrt(p q) + sqrt((1-p)(1-q))] ** (m_+ + m_-).

  Thus it depends only on ``m_+ + m_- = n-1``.  This is an oracle benchmark,
  not an achievability result.  In an exactly balanced planted bisection,
  revealing every label except one would reveal the last label by the global
  count, so the experiment must instead be interpreted under an i.i.d. prior
  or a composition-blind genie (or replaced by a balance-compatible paired
  test).

* The Poisson/Chernoff constant

      C_n = (sqrt(a_n) - sqrt(b_n))^2 / 2

  is the asymptotic exponent of that oracle affinity.  The benchmark for
  almost exact recovery is ``C_n -> infinity``.  In the logarithmic regime
  ``a_n=A log n``, ``b_n=B log n``, the strict information-theoretic exact
  recovery condition is ``(sqrt(A)-sqrt(B))^2 > 2``.

None of the oracle calculations proves mixing, concentration, or
achievability for a hierarchical Gibbs dynamics.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass


DIAGNOSTIC_STATUS = "EXACT_SCALAR_FORMULAS_ORACLE_BENCHMARKS_ONLY"
THRESHOLD_TOLERANCE = 1e-12


@dataclass(frozen=True)
class WeakKSDiagnostic:
    """Kesten--Stigum parameters for the symmetric two-class SBM."""

    average_degree_d: float
    edge_correlation_theta: float
    signal_to_noise_lambda: float
    threshold: float
    regime: str
    weak_recovery_threshold_is_known_for_symmetric_two_class_sbm: bool
    hierarchical_dynamics_used_in_calculation: bool


@dataclass(frozen=True)
class OracleBinomialDiagnostic:
    """Exact finite-size affinity for the composition-blind oracle test."""

    n: int
    within_parameter_a_n: float
    between_parameter_b_n: float
    within_probability: float
    between_probability: float
    first_reference_group_size: int
    second_reference_group_size: int
    exposed_vertices: int
    single_edge_bhattacharyya_coefficient: float
    log_bhattacharyya_affinity: float
    bhattacharyya_affinity: float
    squared_hellinger_distance: float
    exact_bhattacharyya_exponent: float
    poisson_chernoff_constant_c_n: float
    exponent_minus_poisson_constant: float
    oracle_benchmark_only: bool
    exactly_balanced_single_vertex_oracle_valid: bool
    balance_warning: str


@dataclass(frozen=True)
class AlmostExactBenchmark:
    """Scope-correct statement of the almost-exact oracle benchmark."""

    poisson_chernoff_constant_c_n: float
    asymptotic_condition: str
    condition_can_be_decided_from_one_finite_n_instance: bool
    oracle_benchmark_only: bool
    hierarchical_achievability_claimed: bool
    interpretation: str


@dataclass(frozen=True)
class ExactLogarithmicBenchmark:
    """Strict logarithmic exact-recovery threshold and its scope."""

    within_coefficient_a: float
    between_coefficient_b: float
    squared_root_separation: float
    chernoff_coefficient: float
    threshold_squared_root_separation: float
    regime: str
    numerical_strict_inequality_met: bool
    fixed_coefficients_explicitly_supplied: bool
    strict_information_theoretic_condition_met: bool
    asymptotic_threshold_classification_claimed: bool
    equality_requires_separate_analysis: bool
    hierarchical_achievability_claimed: bool
    interpretation: str


@dataclass(frozen=True)
class SBMRecoveryRegimesDiagnostic:
    """Combined reproducible diagnostic for all three recovery regimes."""

    weak: WeakKSDiagnostic
    finite_n_oracle: OracleBinomialDiagnostic
    almost_exact: AlmostExactBenchmark
    exact_logarithmic: ExactLogarithmicBenchmark
    logarithmic_coefficients_source: str
    diagnostic_status: str
    graph_recovery_theorem_reproved: bool
    hierarchical_achievability_claimed: bool
    interpretation: str


def _validate_nonnegative_pair(first: float, second: float) -> None:
    if not math.isfinite(first) or not math.isfinite(second):
        raise ValueError("parameters must be finite")
    if first < 0.0 or second < 0.0:
        raise ValueError("parameters must be non-negative")


def _threshold_regime(
    value: float,
    threshold: float,
    *,
    below: str,
    at: str,
    above: str,
) -> str:
    if math.isclose(
        value,
        threshold,
        rel_tol=THRESHOLD_TOLERANCE,
        abs_tol=THRESHOLD_TOLERANCE,
    ):
        return at
    return below if value < threshold else above


def weak_ks_parameters(a_n: float, b_n: float) -> WeakKSDiagnostic:
    """Return ``d``, ``theta``, and ``lambda=d theta^2`` exactly."""

    _validate_nonnegative_pair(a_n, b_n)
    if a_n + b_n <= 0.0:
        raise ValueError("a_n + b_n must be positive")
    average_degree = 0.5 * (a_n + b_n)
    theta = (a_n - b_n) / (a_n + b_n)
    signal_to_noise = average_degree * theta * theta
    regime = _threshold_regime(
        signal_to_noise,
        1.0,
        below="BELOW_KS",
        at="AT_KS",
        above="ABOVE_KS",
    )
    return WeakKSDiagnostic(
        average_degree_d=average_degree,
        edge_correlation_theta=theta,
        signal_to_noise_lambda=signal_to_noise,
        threshold=1.0,
        regime=regime,
        weak_recovery_threshold_is_known_for_symmetric_two_class_sbm=True,
        hierarchical_dynamics_used_in_calculation=False,
    )


def bernoulli_bhattacharyya_coefficient(
    first_probability: float,
    second_probability: float,
) -> float:
    """Return the exact Bhattacharyya coefficient of two Bernoulli laws."""

    if not 0.0 <= first_probability <= 1.0:
        raise ValueError("first_probability must belong to [0, 1]")
    if not 0.0 <= second_probability <= 1.0:
        raise ValueError("second_probability must belong to [0, 1]")
    deficit = 0.5 * (
        (math.sqrt(first_probability) - math.sqrt(second_probability)) ** 2
        + (
            math.sqrt(1.0 - first_probability)
            - math.sqrt(1.0 - second_probability)
        )
        ** 2
    )
    return 1.0 - deficit


def binomial_bhattacharyya_affinity(
    trials: int,
    first_probability: float,
    second_probability: float,
) -> float:
    """Return the exact affinity of ``Bin(trials,p)`` and ``Bin(trials,q)``."""

    if trials < 0:
        raise ValueError("trials must be non-negative")
    coefficient = bernoulli_bhattacharyya_coefficient(
        first_probability,
        second_probability,
    )
    return coefficient**trials


def poisson_chernoff_constant(a_n: float, b_n: float) -> float:
    """Return ``C_n=(sqrt(a_n)-sqrt(b_n))^2/2``."""

    _validate_nonnegative_pair(a_n, b_n)
    return 0.5 * (math.sqrt(a_n) - math.sqrt(b_n)) ** 2


def oracle_binomial_diagnostic(
    n: int,
    a_n: float,
    b_n: float,
    *,
    first_reference_group_size: int | None = None,
) -> OracleBinomialDiagnostic:
    """Compute the exact finite-``n`` oracle affinity stably.

    The reference labels are held fixed under both hypotheses.  This is
    compatible with an i.i.d.-label model or a composition-blind genie.  It
    is not the literal single-vertex genie in an exactly balanced bisection.
    """

    if n < 2:
        raise ValueError("n must be at least 2")
    _validate_nonnegative_pair(a_n, b_n)
    if a_n >= n or b_n >= n:
        raise ValueError("a_n/n and b_n/n must be strictly smaller than 1")
    if first_reference_group_size is None:
        first_group_size = (n - 1) // 2
    else:
        first_group_size = first_reference_group_size
    if not 0 <= first_group_size <= n - 1:
        raise ValueError(
            "first_reference_group_size must belong to [0, n-1]"
        )
    second_group_size = n - 1 - first_group_size
    within_probability = a_n / n
    between_probability = b_n / n
    coefficient = bernoulli_bhattacharyya_coefficient(
        within_probability,
        between_probability,
    )
    coefficient_deficit = 0.5 * (
        (
            math.sqrt(within_probability)
            - math.sqrt(between_probability)
        )
        ** 2
        + (
            math.sqrt(1.0 - within_probability)
            - math.sqrt(1.0 - between_probability)
        )
        ** 2
    )
    if coefficient_deficit >= 1.0:
        raise ValueError("the Bhattacharyya coefficient must be positive")
    log_affinity = (n - 1) * math.log1p(-coefficient_deficit)
    affinity = math.exp(log_affinity)
    squared_hellinger = -math.expm1(log_affinity)
    exact_exponent = -log_affinity
    poisson_constant = poisson_chernoff_constant(a_n, b_n)
    return OracleBinomialDiagnostic(
        n=n,
        within_parameter_a_n=a_n,
        between_parameter_b_n=b_n,
        within_probability=within_probability,
        between_probability=between_probability,
        first_reference_group_size=first_group_size,
        second_reference_group_size=second_group_size,
        exposed_vertices=n - 1,
        single_edge_bhattacharyya_coefficient=coefficient,
        log_bhattacharyya_affinity=log_affinity,
        bhattacharyya_affinity=affinity,
        squared_hellinger_distance=squared_hellinger,
        exact_bhattacharyya_exponent=exact_exponent,
        poisson_chernoff_constant_c_n=poisson_constant,
        exponent_minus_poisson_constant=(
            exact_exponent - poisson_constant
        ),
        oracle_benchmark_only=True,
        exactly_balanced_single_vertex_oracle_valid=False,
        balance_warning=(
            "Revealing all other labels in an exactly balanced bisection "
            "determines the final label by the global count. Use an "
            "i.i.d.-label/composition-blind oracle or a balance-compatible "
            "paired test instead."
        ),
    )


def almost_exact_oracle_benchmark(
    a_n: float,
    b_n: float,
) -> AlmostExactBenchmark:
    """Return the oracle benchmark without inferring a limit from one point."""

    constant = poisson_chernoff_constant(a_n, b_n)
    return AlmostExactBenchmark(
        poisson_chernoff_constant_c_n=constant,
        asymptotic_condition="C_n -> infinity along the parameter sequence",
        condition_can_be_decided_from_one_finite_n_instance=False,
        oracle_benchmark_only=True,
        hierarchical_achievability_claimed=False,
        interpretation=(
            "Divergence of C_n is the oracle degree-profile benchmark for "
            "vanishing local error. Achievability requires a global "
            "initializer and a robust refinement argument."
        ),
    )


def logarithmic_exact_recovery_benchmark(
    within_coefficient_a: float,
    between_coefficient_b: float,
    *,
    fixed_coefficients_explicitly_supplied: bool = True,
) -> ExactLogarithmicBenchmark:
    """Classify the strict logarithmic exact-recovery benchmark."""

    _validate_nonnegative_pair(
        within_coefficient_a,
        between_coefficient_b,
    )
    root_separation = (
        math.sqrt(within_coefficient_a)
        - math.sqrt(between_coefficient_b)
    ) ** 2
    regime = _threshold_regime(
        root_separation,
        2.0,
        below="BELOW_EXACT_RECOVERY_THRESHOLD",
        at="AT_FIRST_ORDER_BOUNDARY",
        above="ABOVE_EXACT_RECOVERY_THRESHOLD",
    )
    if not fixed_coefficients_explicitly_supplied:
        regime = _threshold_regime(
            root_separation,
            2.0,
            below="FINITE_N_EFFECTIVE_COMPARISON_BELOW_TWO",
            at="FINITE_N_EFFECTIVE_COMPARISON_AT_TWO",
            above="FINITE_N_EFFECTIVE_COMPARISON_ABOVE_TWO",
        )
    numerical_strict_inequality = (
        root_separation > 2.0
        and not math.isclose(
            root_separation,
            2.0,
            rel_tol=THRESHOLD_TOLERANCE,
            abs_tol=THRESHOLD_TOLERANCE,
        )
    )
    return ExactLogarithmicBenchmark(
        within_coefficient_a=within_coefficient_a,
        between_coefficient_b=between_coefficient_b,
        squared_root_separation=root_separation,
        chernoff_coefficient=0.5 * root_separation,
        threshold_squared_root_separation=2.0,
        regime=regime,
        numerical_strict_inequality_met=numerical_strict_inequality,
        fixed_coefficients_explicitly_supplied=(
            fixed_coefficients_explicitly_supplied
        ),
        strict_information_theoretic_condition_met=(
            numerical_strict_inequality
            and fixed_coefficients_explicitly_supplied
        ),
        asymptotic_threshold_classification_claimed=(
            fixed_coefficients_explicitly_supplied
        ),
        equality_requires_separate_analysis=(
            regime == "AT_FIRST_ORDER_BOUNDARY"
        ),
        hierarchical_achievability_claimed=False,
        interpretation=(
            (
                "The strict condition (sqrt(A)-sqrt(B))^2>2 is the "
                "classical information-theoretic benchmark for fixed "
                "logarithmic coefficients. "
            )
            if fixed_coefficients_explicitly_supplied
            else (
                "The displayed A=a_n/log(n), B=b_n/log(n) are effective "
                "finite-n coefficients only; their numerical comparison "
                "with 2 is not an asymptotic recovery classification. "
            )
        )
        + (
            "This calculation does not prove hierarchical Gibbs "
            "achievability."
        ),
    )


def run_diagnostic(
    *,
    n: int,
    a_n: float,
    b_n: float,
    first_reference_group_size: int | None = None,
    logarithmic_within_coefficient: float | None = None,
    logarithmic_between_coefficient: float | None = None,
) -> SBMRecoveryRegimesDiagnostic:
    """Combine weak, finite-oracle, almost, and exact-log benchmarks."""

    if (
        logarithmic_within_coefficient is None
        and logarithmic_between_coefficient is None
    ):
        if n <= 1:
            raise ValueError("n must exceed 1 to infer logarithmic coefficients")
        log_n = math.log(n)
        log_within = a_n / log_n
        log_between = b_n / log_n
        source = "finite-n effective coefficients a_n/log(n), b_n/log(n)"
    elif (
        logarithmic_within_coefficient is not None
        and logarithmic_between_coefficient is not None
    ):
        log_within = logarithmic_within_coefficient
        log_between = logarithmic_between_coefficient
        source = "explicit fixed logarithmic coefficients"
    else:
        raise ValueError(
            "provide both logarithmic coefficients or neither of them"
        )

    return SBMRecoveryRegimesDiagnostic(
        weak=weak_ks_parameters(a_n, b_n),
        finite_n_oracle=oracle_binomial_diagnostic(
            n,
            a_n,
            b_n,
            first_reference_group_size=first_reference_group_size,
        ),
        almost_exact=almost_exact_oracle_benchmark(a_n, b_n),
        exact_logarithmic=logarithmic_exact_recovery_benchmark(
            log_within,
            log_between,
            fixed_coefficients_explicitly_supplied=(
                logarithmic_within_coefficient is not None
            ),
        ),
        logarithmic_coefficients_source=source,
        diagnostic_status=DIAGNOSTIC_STATUS,
        graph_recovery_theorem_reproved=False,
        hierarchical_achievability_claimed=False,
        interpretation=(
            "All displayed scalar formulas are exact at their stated level. "
            "The finite-n degree-profile calculation is oracle-only, and no "
            "mixing, graph-to-tree transfer, or hierarchical achievability "
            "is inferred."
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=100_000)
    parser.add_argument("--a", dest="a_n", type=float, default=30.0)
    parser.add_argument("--b", dest="b_n", type=float, default=10.0)
    parser.add_argument(
        "--first-reference-group-size",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--log-within-coefficient",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--log-between-coefficient",
        type=float,
        default=None,
    )
    arguments = parser.parse_args()
    diagnostic = run_diagnostic(
        n=arguments.n,
        a_n=arguments.a_n,
        b_n=arguments.b_n,
        first_reference_group_size=arguments.first_reference_group_size,
        logarithmic_within_coefficient=(
            arguments.log_within_coefficient
        ),
        logarithmic_between_coefficient=(
            arguments.log_between_coefficient
        ),
    )
    print(
        json.dumps(
            asdict(diagnostic),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
