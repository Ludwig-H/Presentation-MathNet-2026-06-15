"""Tests for the exact SBM critical-cut replica audit."""

from __future__ import annotations

import math
import unittest

from sbm_critical_cut_replica_diagnostic import (
    EXACT_STATUS,
    critical_beta,
    diagnose_critical_cut_replica,
    marginal_replica_jacobian,
    oracle_replica_jacobian,
    residual_correlation,
)


class CriticalCutReplicaDiagnosticTests(unittest.TestCase):
    def test_reference_example_exhibits_oracle_inflation(self) -> None:
        result = diagnose_critical_cut_replica(degree=3.0, theta=0.5)
        self.assertAlmostEqual(result.residual_correlation_theta_res, 0.25)
        self.assertAlmostEqual(result.marginal_replica_jacobian, 0.25)
        self.assertAlmostEqual(result.oracle_replica_jacobian, 0.375)
        self.assertAlmostEqual(result.lambda_kesten_stigum, 0.75)
        self.assertAlmostEqual(result.oracle_branching_factor, 1.125)
        self.assertAlmostEqual(
            result.oracle_branching_factor_inflation,
            0.375,
        )

    def test_beta_solves_the_critical_cut_equation(self) -> None:
        degree = 3.0
        theta = 0.5
        p_equal = 0.5 * (1.0 + theta)
        clock_rate = math.log(p_equal / (1.0 - p_equal))
        beta = critical_beta(degree, theta)
        q_beta = p_equal * (1.0 - math.exp(-clock_rate * beta))
        self.assertAlmostEqual(q_beta, 1.0 / degree)
        self.assertLess(beta, 1.0)

    def test_finite_edge_identities_hold_on_a_parameter_grid(self) -> None:
        for degree, theta in (
            (2.0, 0.4),
            (2.5, 0.6),
            (3.0, 0.5),
            (5.0, 0.3),
        ):
            with self.subTest(degree=degree, theta=theta):
                result = diagnose_critical_cut_replica(degree, theta)
                theta_res = (
                    (theta - 1.0 / degree)
                    / (1.0 - 1.0 / degree)
                )
                oracle_closed_form = (
                    theta * theta
                    + (1.0 / degree)
                    * (1.0 - theta) ** 2
                    / (1.0 - 1.0 / degree)
                )
                self.assertAlmostEqual(
                    residual_correlation(degree, theta),
                    theta_res,
                )
                self.assertAlmostEqual(
                    oracle_replica_jacobian(degree, theta),
                    oracle_closed_form,
                )
                self.assertAlmostEqual(
                    result.marginalized_edge_correlation,
                    theta,
                )
                self.assertAlmostEqual(
                    result.marginal_replica_jacobian,
                    theta * theta,
                )
                self.assertGreater(
                    result.oracle_replica_jacobian,
                    result.marginal_replica_jacobian,
                )
                self.assertTrue(result.exact_finite_edge_identity)

    def test_cut_time_flag_distinguishes_extended_beta(self) -> None:
        observed = diagnose_critical_cut_replica(3.0, 0.5)
        extended = diagnose_critical_cut_replica(2.0, 0.4)
        self.assertTrue(observed.critical_cut_within_final_time)
        self.assertLessEqual(observed.beta_c, 1.0)
        self.assertFalse(extended.critical_cut_within_final_time)
        self.assertGreater(extended.beta_c, 1.0)

    def test_status_flags_make_the_scope_explicit(self) -> None:
        result = diagnose_critical_cut_replica(3.0, 0.5)
        self.assertEqual(result.diagnostic_status, EXACT_STATUS)
        self.assertTrue(result.exact_finite_edge_identity)
        self.assertFalse(result.graph_theorem_claimed)
        self.assertFalse(result.hierarchical_threshold_proof_claimed)
        self.assertIn("not a graph theorem", result.interpretation)

    def test_marginal_jacobian_is_theta_squared(self) -> None:
        self.assertAlmostEqual(marginal_replica_jacobian(0.7), 0.49)

    def test_invalid_parameters_are_rejected(self) -> None:
        for degree, theta in (
            (1.0, 0.5),
            (3.0, 0.0),
            (3.0, 1.0),
            (1.1, 0.1),
        ):
            with self.subTest(degree=degree, theta=theta):
                with self.assertRaises(ValueError):
                    diagnose_critical_cut_replica(degree, theta)


if __name__ == "__main__":
    unittest.main()
