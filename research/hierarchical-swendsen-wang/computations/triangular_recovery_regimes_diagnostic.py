"""Oracle recovery obstructions for the degree-six triangular GSBM.

Reveal every hidden label except the label at one vertex.  Its six incident
observations are then independent binary symmetric channels with reliability
``p``.  The optimal rule is majority vote, with a fair tie break.  Its exact
error is

    epsilon_6(p)
      = P(Bin(6,p) <= 2) + 0.5 P(Bin(6,p) = 3)
      = 10 delta^3 - 15 delta^4 + 6 delta^5,

where ``delta=1-p``.

Consequences are oracle lower bounds, hence necessary conditions.

* For every fixed ``p<1``, the error is positive, so almost exact recovery is
  impossible.
* A linear packing of vertices with disjoint incident-edge stars gives
  independent oracle tests.  Exact recovery therefore requires
  ``n epsilon_6(p_n) -> 0``.  Since ``epsilon_6(1-delta) ~ 10 delta^3``, this
  implies ``1-p_n=o(n^{-1/3})``.

No sufficiency or hierarchical-dynamics claim is made.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass


DIAGNOSTIC_STATUS = "EXACT_DEGREE_SIX_ORACLE_NECESSARY_CONDITIONS"


@dataclass(frozen=True)
class TriangularRecoveryRegimesDiagnostic:
    """Exact local error and finite-size versions of the oracle obstructions."""

    vertex_count: int
    p: float
    noise_delta: float
    degree: int
    exact_majority_oracle_error: float
    polynomial_oracle_error: float
    polynomial_identity_error: float
    leading_cubic_approximation: float
    error_over_leading_cubic: float | None
    scaled_exact_recovery_obstruction_n_epsilon: float
    scaled_cubic_proxy_10_n_delta_cubed: float
    fixed_p_less_than_one: bool
    almost_exact_recovery_impossible_for_this_fixed_p_sequence: bool
    exact_recovery_impossible_for_this_fixed_p_sequence: bool
    almost_exact_necessary_condition: str
    exact_recovery_necessary_condition: str
    equivalent_exact_scaling_necessary_condition: str
    oracle_benchmark_only: bool
    sufficiency_claimed: bool
    hierarchical_achievability_claimed: bool
    diagnostic_status: str
    interpretation: str


def majority_oracle_error_degree_six(p: float) -> float:
    """Return the Bayes error of six BSC observations with fair tie break."""

    if not 0.5 <= p <= 1.0 or not math.isfinite(p):
        raise ValueError("p must belong to [1/2, 1]")
    wrong = 1.0 - p
    strict_error = math.fsum(
        math.comb(6, correct)
        * p**correct
        * wrong ** (6 - correct)
        for correct in range(3)
    )
    tie_error = 0.5 * math.comb(6, 3) * p**3 * wrong**3
    return strict_error + tie_error


def majority_oracle_error_polynomial(p: float) -> float:
    """Return ``10 delta^3-15 delta^4+6 delta^5``."""

    if not 0.5 <= p <= 1.0 or not math.isfinite(p):
        raise ValueError("p must belong to [1/2, 1]")
    delta = 1.0 - p
    return 10.0 * delta**3 - 15.0 * delta**4 + 6.0 * delta**5


def run_diagnostic(
    *,
    vertex_count: int,
    p: float,
) -> TriangularRecoveryRegimesDiagnostic:
    """Evaluate the exact finite-size oracle obstruction."""

    if vertex_count <= 0:
        raise ValueError("vertex_count must be positive")
    exact_error = majority_oracle_error_degree_six(p)
    polynomial_error = majority_oracle_error_polynomial(p)
    identity_error = abs(exact_error - polynomial_error)
    delta = 1.0 - p
    leading = 10.0 * delta**3
    ratio = None if leading == 0.0 else exact_error / leading
    fixed_noisily_observed = p < 1.0
    return TriangularRecoveryRegimesDiagnostic(
        vertex_count=vertex_count,
        p=p,
        noise_delta=delta,
        degree=6,
        exact_majority_oracle_error=exact_error,
        polynomial_oracle_error=polynomial_error,
        polynomial_identity_error=identity_error,
        leading_cubic_approximation=leading,
        error_over_leading_cubic=ratio,
        scaled_exact_recovery_obstruction_n_epsilon=(
            vertex_count * exact_error
        ),
        scaled_cubic_proxy_10_n_delta_cubed=vertex_count * leading,
        fixed_p_less_than_one=fixed_noisily_observed,
        almost_exact_recovery_impossible_for_this_fixed_p_sequence=(
            fixed_noisily_observed
        ),
        exact_recovery_impossible_for_this_fixed_p_sequence=(
            fixed_noisily_observed
        ),
        almost_exact_necessary_condition="p_n -> 1",
        exact_recovery_necessary_condition="n epsilon_6(p_n) -> 0",
        equivalent_exact_scaling_necessary_condition=(
            "1-p_n = o(n^(-1/3))"
        ),
        oracle_benchmark_only=True,
        sufficiency_claimed=False,
        hierarchical_achievability_claimed=False,
        diagnostic_status=DIAGNOSTIC_STATUS,
        interpretation=(
            "Revealing all neighbouring labels can only help recovery. "
            "Positive oracle error rules out almost exact recovery at fixed "
            "p<1; a linear disjoint-star packing yields the exact-recovery "
            "necessary condition. Neither condition is asserted sufficient."
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vertices", type=int, default=1_000_000)
    parser.add_argument("--p", type=float, default=0.81)
    arguments = parser.parse_args()
    result = run_diagnostic(
        vertex_count=arguments.vertices,
        p=arguments.p,
    )
    print(
        json.dumps(
            asdict(result),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
