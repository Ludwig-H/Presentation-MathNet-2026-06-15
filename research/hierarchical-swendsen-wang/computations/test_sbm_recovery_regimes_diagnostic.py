"""Tests for the classical-SBM recovery-regime diagnostic."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

from sbm_recovery_regimes_diagnostic import (
    DIAGNOSTIC_STATUS,
    almost_exact_oracle_benchmark,
    bernoulli_bhattacharyya_coefficient,
    binomial_bhattacharyya_affinity,
    logarithmic_exact_recovery_benchmark,
    oracle_binomial_diagnostic,
    poisson_chernoff_constant,
    run_diagnostic,
    weak_ks_parameters,
)


def _brute_force_binomial_affinity(
    trials: int,
    first_probability: float,
    second_probability: float,
) -> float:
    affinity = 0.0
    for successes in range(trials + 1):
        multiplicity = math.comb(trials, successes)
        first_mass = (
            multiplicity
            * first_probability**successes
            * (1.0 - first_probability) ** (trials - successes)
        )
        second_mass = (
            multiplicity
            * second_probability**successes
            * (1.0 - second_probability) ** (trials - successes)
        )
        affinity += math.sqrt(first_mass * second_mass)
    return affinity


class WeakKSTests(unittest.TestCase):
    def test_ks_number_is_d_theta_squared(self) -> None:
        diagnostic = weak_ks_parameters(5.0, 1.0)
        self.assertEqual(diagnostic.average_degree_d, 3.0)
        self.assertAlmostEqual(diagnostic.edge_correlation_theta, 2.0 / 3.0)
        self.assertAlmostEqual(
            diagnostic.signal_to_noise_lambda,
            4.0 / 3.0,
        )
        self.assertEqual(diagnostic.regime, "ABOVE_KS")
        self.assertFalse(diagnostic.hierarchical_dynamics_used_in_calculation)

    def test_all_three_ks_regimes_are_distinguished(self) -> None:
        self.assertEqual(weak_ks_parameters(3.0, 1.0).regime, "BELOW_KS")
        self.assertEqual(weak_ks_parameters(2.0, 0.0).regime, "AT_KS")
        self.assertEqual(weak_ks_parameters(5.0, 1.0).regime, "ABOVE_KS")


class FiniteBinomialOracleTests(unittest.TestCase):
    def test_closed_form_matches_direct_binomial_sum(self) -> None:
        trials = 7
        first_probability = 0.2
        second_probability = 0.05
        expected = _brute_force_binomial_affinity(
            trials,
            first_probability,
            second_probability,
        )
        self.assertAlmostEqual(
            binomial_bhattacharyya_affinity(
                trials,
                first_probability,
                second_probability,
            ),
            expected,
        )

    def test_two_group_oracle_affinity_is_exact(self) -> None:
        n = 20
        a_n = 4.0
        b_n = 1.0
        first_group_size = 8
        result = oracle_binomial_diagnostic(
            n,
            a_n,
            b_n,
            first_reference_group_size=first_group_size,
        )
        coefficient = bernoulli_bhattacharyya_coefficient(
            a_n / n,
            b_n / n,
        )
        first_affinity = _brute_force_binomial_affinity(
            first_group_size,
            a_n / n,
            b_n / n,
        )
        second_affinity = _brute_force_binomial_affinity(
            n - 1 - first_group_size,
            b_n / n,
            a_n / n,
        )
        self.assertAlmostEqual(
            result.bhattacharyya_affinity,
            first_affinity * second_affinity,
        )
        self.assertAlmostEqual(
            result.bhattacharyya_affinity,
            coefficient ** (n - 1),
        )
        self.assertAlmostEqual(
            result.squared_hellinger_distance,
            1.0 - result.bhattacharyya_affinity,
        )
        self.assertAlmostEqual(
            result.exact_bhattacharyya_exponent,
            -math.log(result.bhattacharyya_affinity),
        )

    def test_poisson_constant_is_the_sparse_limit(self) -> None:
        a_n = 5.0
        b_n = 1.0
        expected = 0.5 * (math.sqrt(a_n) - math.sqrt(b_n)) ** 2
        self.assertAlmostEqual(
            poisson_chernoff_constant(a_n, b_n),
            expected,
        )
        finite = oracle_binomial_diagnostic(2_000_000, a_n, b_n)
        self.assertAlmostEqual(
            finite.exact_bhattacharyya_exponent,
            expected,
            delta=1e-5,
        )

    def test_exponent_is_stable_when_single_edge_affinity_rounds_to_one(
        self,
    ) -> None:
        a_n = 5.0
        b_n = 1.0
        result = oracle_binomial_diagnostic(10**16, a_n, b_n)
        self.assertGreaterEqual(
            result.single_edge_bhattacharyya_coefficient,
            1.0 - 2e-16,
        )
        self.assertAlmostEqual(
            result.exact_bhattacharyya_exponent,
            poisson_chernoff_constant(a_n, b_n),
            delta=1e-12,
        )

    def test_scope_flags_reject_the_balanced_single_vertex_oracle(self) -> None:
        result = oracle_binomial_diagnostic(1000, 10.0, 2.0)
        self.assertTrue(result.oracle_benchmark_only)
        self.assertFalse(result.exactly_balanced_single_vertex_oracle_valid)
        self.assertIn("global count", result.balance_warning)


class AlmostAndExactRecoveryTests(unittest.TestCase):
    def test_almost_exact_is_only_an_asymptotic_oracle_benchmark(self) -> None:
        result = almost_exact_oracle_benchmark(100.0, 25.0)
        self.assertEqual(result.poisson_chernoff_constant_c_n, 12.5)
        self.assertIn("C_n -> infinity", result.asymptotic_condition)
        self.assertFalse(
            result.condition_can_be_decided_from_one_finite_n_instance
        )
        self.assertTrue(result.oracle_benchmark_only)
        self.assertFalse(result.hierarchical_achievability_claimed)

    def test_exact_logarithmic_threshold_has_three_regimes(self) -> None:
        below = logarithmic_exact_recovery_benchmark(4.0, 1.0)
        boundary = logarithmic_exact_recovery_benchmark(
            (1.0 + math.sqrt(2.0)) ** 2,
            1.0,
        )
        above = logarithmic_exact_recovery_benchmark(9.0, 1.0)
        self.assertEqual(below.regime, "BELOW_EXACT_RECOVERY_THRESHOLD")
        self.assertEqual(boundary.regime, "AT_FIRST_ORDER_BOUNDARY")
        self.assertEqual(above.regime, "ABOVE_EXACT_RECOVERY_THRESHOLD")
        self.assertFalse(below.strict_information_theoretic_condition_met)
        self.assertTrue(boundary.equality_requires_separate_analysis)
        self.assertTrue(above.strict_information_theoretic_condition_met)
        self.assertTrue(above.numerical_strict_inequality_met)
        self.assertTrue(above.fixed_coefficients_explicitly_supplied)
        self.assertTrue(above.asymptotic_threshold_classification_claimed)
        self.assertFalse(above.hierarchical_achievability_claimed)

    def test_combined_summary_preserves_all_nonclaims(self) -> None:
        result = run_diagnostic(
            n=100_000,
            a_n=30.0,
            b_n=10.0,
            logarithmic_within_coefficient=9.0,
            logarithmic_between_coefficient=1.0,
        )
        self.assertEqual(result.diagnostic_status, DIAGNOSTIC_STATUS)
        self.assertEqual(
            result.logarithmic_coefficients_source,
            "explicit fixed logarithmic coefficients",
        )
        self.assertFalse(result.graph_recovery_theorem_reproved)
        self.assertFalse(result.hierarchical_achievability_claimed)
        self.assertTrue(result.finite_n_oracle.oracle_benchmark_only)
        self.assertFalse(
            result.almost_exact.hierarchical_achievability_claimed
        )
        self.assertFalse(
            result.exact_logarithmic.hierarchical_achievability_claimed
        )

    def test_inferred_effective_coefficients_do_not_claim_a_limit(self) -> None:
        result = run_diagnostic(
            n=100_000,
            a_n=90.0,
            b_n=10.0,
        )
        self.assertTrue(
            result.exact_logarithmic.numerical_strict_inequality_met
        )
        self.assertFalse(
            result
            .exact_logarithmic
            .fixed_coefficients_explicitly_supplied
        )
        self.assertFalse(
            result
            .exact_logarithmic
            .strict_information_theoretic_condition_met
        )
        self.assertFalse(
            result
            .exact_logarithmic
            .asymptotic_threshold_classification_claimed
        )
        self.assertEqual(
            result.exact_logarithmic.regime,
            "FINITE_N_EFFECTIVE_COMPARISON_ABOVE_TWO",
        )
        self.assertIn(
            "effective finite-n",
            result.exact_logarithmic.interpretation,
        )

    def test_cli_emits_machine_readable_json(self) -> None:
        module_path = Path(__file__).with_name(
            "sbm_recovery_regimes_diagnostic.py"
        )
        completed = subprocess.run(
            (
                sys.executable,
                str(module_path),
                "--n",
                "10000",
                "--a",
                "20",
                "--b",
                "4",
                "--log-within-coefficient",
                "9",
                "--log-between-coefficient",
                "1",
            ),
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["diagnostic_status"], DIAGNOSTIC_STATUS)
        self.assertIn("signal_to_noise_lambda", payload["weak"])
        self.assertIn(
            "bhattacharyya_affinity",
            payload["finite_n_oracle"],
        )
        self.assertFalse(payload["hierarchical_achievability_claimed"])


class ValidationTests(unittest.TestCase):
    def test_invalid_parameters_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            weak_ks_parameters(0.0, 0.0)
        with self.assertRaises(ValueError):
            oracle_binomial_diagnostic(1, 0.5, 0.1)
        with self.assertRaises(ValueError):
            oracle_binomial_diagnostic(10, 10.0, 1.0)
        with self.assertRaises(ValueError):
            binomial_bhattacharyya_affinity(-1, 0.2, 0.1)
        with self.assertRaises(ValueError):
            run_diagnostic(
                n=100,
                a_n=5.0,
                b_n=1.0,
                logarithmic_within_coefficient=4.0,
            )


if __name__ == "__main__":
    unittest.main()
