"""Exact one-edge audit of the critical-cut oracle on the symmetric SBM.

For a broadcast edge of correlation ``theta``, put

    p = (1 + theta) / 2,
    u = log(p / (1 - p)),
    q(beta) = p (1 - exp(-u beta)).

The geometric PGW critical cut is defined by ``q(beta_c) = 1 / d``.  If its
bit ``B`` is retained as oracle information, then

    E[ST | B=1] = 1,
    E[ST | B=0] = theta_res
                  = (theta - 1/d) / (1 - 1/d).

Consequently the oracle two-replica edge factor is

    E[E[ST | B]^2]
      = 1/d + (1 - 1/d) theta_res^2
      = theta^2 + (1/d) (1-theta)^2 / (1 - 1/d).

Marginalizing ``B`` first restores the original edge channel and its
two-replica Jacobian ``theta^2``.  This exact finite-edge distinction explains
why a permanently frozen critical cut can inflate the branching factor.  It
is neither a graph theorem nor a proof of a hierarchical threshold.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass


EXACT_STATUS = "EXACT_FINITE_EDGE_ORACLE_INFLATION_IDENTITY"


@dataclass(frozen=True)
class CriticalCutReplicaDiagnostic:
    """All scalar identities for one PGW broadcast edge."""

    degree: float
    theta: float
    p_equal: float
    clock_rate_u: float
    critical_open_probability: float
    beta_c: float
    critical_open_probability_from_beta: float
    critical_cut_within_final_time: bool
    open_given_equal_probability: float
    residual_correlation_theta_res: float
    marginalized_edge_correlation: float
    marginal_replica_jacobian: float
    oracle_replica_jacobian: float
    oracle_replica_jacobian_closed_form: float
    lambda_kesten_stigum: float
    oracle_branching_factor: float
    oracle_branching_factor_inflation: float
    oracle_branching_factor_inflation_closed_form: float
    critical_equation_error: float
    marginalization_identity_error: float
    oracle_identity_error: float
    inflation_identity_error: float
    exact_finite_edge_identity: bool
    graph_theorem_claimed: bool
    hierarchical_threshold_proof_claimed: bool
    diagnostic_status: str
    interpretation: str


def _validate_inputs(degree: float, theta: float) -> None:
    if degree <= 1.0:
        raise ValueError("degree must be strictly greater than 1")
    if not 0.0 < theta < 1.0:
        raise ValueError("theta must satisfy 0 < theta < 1")
    p_equal = 0.5 * (1.0 + theta)
    if degree * p_equal <= 1.0:
        raise ValueError(
            "q(beta)=1/d has no finite solution unless d p > 1"
        )


def critical_beta(degree: float, theta: float) -> float:
    """Solve ``p(1-exp(-u beta_c))=1/d`` in closed form."""

    _validate_inputs(degree, theta)
    p_equal = 0.5 * (1.0 + theta)
    clock_rate = math.log(p_equal / (1.0 - p_equal))
    return -math.log1p(-1.0 / (degree * p_equal)) / clock_rate


def residual_correlation(degree: float, theta: float) -> float:
    """Return the correlation conditional on the critical edge staying shut."""

    _validate_inputs(degree, theta)
    critical_probability = 1.0 / degree
    return (
        (theta - critical_probability)
        / (1.0 - critical_probability)
    )


def oracle_replica_jacobian(degree: float, theta: float) -> float:
    """Return ``E[E[ST|B]^2]`` when the critical-cut bit is retained."""

    theta_res = residual_correlation(degree, theta)
    critical_probability = 1.0 / degree
    return (
        critical_probability
        + (1.0 - critical_probability) * theta_res * theta_res
    )


def marginal_replica_jacobian(theta: float) -> float:
    """Return the two-replica Jacobian after marginalizing the cut bit."""

    if not 0.0 < theta < 1.0:
        raise ValueError("theta must satisfy 0 < theta < 1")
    return theta * theta


def diagnose_critical_cut_replica(
    degree: float,
    theta: float,
) -> CriticalCutReplicaDiagnostic:
    """Evaluate and cross-audit the exact finite-edge identities."""

    _validate_inputs(degree, theta)
    p_equal = 0.5 * (1.0 + theta)
    clock_rate = math.log(p_equal / (1.0 - p_equal))
    critical_probability = 1.0 / degree
    beta = critical_beta(degree, theta)
    open_given_equal = -math.expm1(-clock_rate * beta)
    critical_from_beta = p_equal * open_given_equal

    theta_res = residual_correlation(degree, theta)
    marginalized_correlation = (
        critical_probability
        + (1.0 - critical_probability) * theta_res
    )
    marginal_jacobian = marginal_replica_jacobian(theta)
    oracle_jacobian = oracle_replica_jacobian(degree, theta)
    oracle_closed_form = (
        marginal_jacobian
        + critical_probability
        * (1.0 - theta) ** 2
        / (1.0 - critical_probability)
    )
    lambda_kesten_stigum = degree * marginal_jacobian
    oracle_branching = degree * oracle_jacobian
    inflation = oracle_branching - lambda_kesten_stigum
    inflation_closed_form = (
        (1.0 - theta) ** 2
        / (1.0 - critical_probability)
    )

    critical_error = abs(critical_from_beta - critical_probability)
    marginalization_error = abs(marginalized_correlation - theta)
    oracle_error = abs(oracle_jacobian - oracle_closed_form)
    inflation_error = abs(inflation - inflation_closed_form)
    maximum_error = max(
        critical_error,
        marginalization_error,
        oracle_error,
        inflation_error,
    )
    exact_identity = maximum_error <= 1e-12

    return CriticalCutReplicaDiagnostic(
        degree=degree,
        theta=theta,
        p_equal=p_equal,
        clock_rate_u=clock_rate,
        critical_open_probability=critical_probability,
        beta_c=beta,
        critical_open_probability_from_beta=critical_from_beta,
        critical_cut_within_final_time=degree * theta >= 1.0,
        open_given_equal_probability=open_given_equal,
        residual_correlation_theta_res=theta_res,
        marginalized_edge_correlation=marginalized_correlation,
        marginal_replica_jacobian=marginal_jacobian,
        oracle_replica_jacobian=oracle_jacobian,
        oracle_replica_jacobian_closed_form=oracle_closed_form,
        lambda_kesten_stigum=lambda_kesten_stigum,
        oracle_branching_factor=oracle_branching,
        oracle_branching_factor_inflation=inflation,
        oracle_branching_factor_inflation_closed_form=inflation_closed_form,
        critical_equation_error=critical_error,
        marginalization_identity_error=marginalization_error,
        oracle_identity_error=oracle_error,
        inflation_identity_error=inflation_error,
        exact_finite_edge_identity=exact_identity,
        graph_theorem_claimed=False,
        hierarchical_threshold_proof_claimed=False,
        diagnostic_status=(
            EXACT_STATUS
            if exact_identity
            else "NUMERICAL_IDENTITY_AUDIT_FAILED"
        ),
        interpretation=(
            "Retaining the critical-cut bit B changes theta^2 into the larger "
            "oracle factor E[E[ST|B]^2]. Marginalizing B restores theta^2. "
            "This is an exact edge calculation, not a graph theorem or a "
            "hierarchical weak-recovery threshold proof."
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=float, default=3.0)
    parser.add_argument("--theta", type=float, default=0.5)
    arguments = parser.parse_args()
    diagnostic = diagnose_critical_cut_replica(
        degree=arguments.degree,
        theta=arguments.theta,
    )
    print(json.dumps(asdict(diagnostic), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
