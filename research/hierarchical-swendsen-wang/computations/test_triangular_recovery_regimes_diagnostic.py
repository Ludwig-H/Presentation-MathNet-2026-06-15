"""Tests for the triangular degree-six oracle recovery diagnostic."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

from triangular_recovery_regimes_diagnostic import (
    DIAGNOSTIC_STATUS,
    majority_oracle_error_degree_six,
    majority_oracle_error_polynomial,
    run_diagnostic,
)


class MajorityOracleErrorTests(unittest.TestCase):
    def test_polynomial_identity_matches_binomial_sum(self) -> None:
        for p in (0.5, 0.7, 0.81, 0.835805792367, 0.99, 1.0):
            self.assertAlmostEqual(
                majority_oracle_error_degree_six(p),
                majority_oracle_error_polynomial(p),
                places=14,
            )

    def test_reference_values(self) -> None:
        self.assertAlmostEqual(
            majority_oracle_error_degree_six(0.81),
            0.0505275094,
            places=10,
        )
        self.assertAlmostEqual(
            majority_oracle_error_degree_six(0.835805792367),
            0.0340799611,
            places=10,
        )

    def test_cubic_asymptotic(self) -> None:
        diagnostic = run_diagnostic(vertex_count=10**9, p=0.9999)
        self.assertIsNotNone(diagnostic.error_over_leading_cubic)
        self.assertAlmostEqual(
            diagnostic.error_over_leading_cubic,
            1.0,
            delta=2e-4,
        )


class ScopeAndCliTests(unittest.TestCase):
    def test_fixed_noise_flags_only_impossibility(self) -> None:
        diagnostic = run_diagnostic(vertex_count=1000, p=0.81)
        self.assertEqual(diagnostic.diagnostic_status, DIAGNOSTIC_STATUS)
        self.assertTrue(
            diagnostic
            .almost_exact_recovery_impossible_for_this_fixed_p_sequence
        )
        self.assertTrue(
            diagnostic.exact_recovery_impossible_for_this_fixed_p_sequence
        )
        self.assertEqual(
            diagnostic.scaled_exact_recovery_obstruction_n_epsilon,
            1000 * diagnostic.exact_majority_oracle_error,
        )
        self.assertTrue(diagnostic.oracle_benchmark_only)
        self.assertFalse(diagnostic.sufficiency_claimed)
        self.assertFalse(diagnostic.hierarchical_achievability_claimed)

    def test_noiseless_endpoint_has_zero_oracle_error(self) -> None:
        diagnostic = run_diagnostic(vertex_count=1000, p=1.0)
        self.assertEqual(diagnostic.exact_majority_oracle_error, 0.0)
        self.assertIsNone(diagnostic.error_over_leading_cubic)
        self.assertFalse(
            diagnostic
            .almost_exact_recovery_impossible_for_this_fixed_p_sequence
        )
        self.assertFalse(
            diagnostic.exact_recovery_impossible_for_this_fixed_p_sequence
        )

    def test_cli_emits_json(self) -> None:
        module = Path(__file__).with_name(
            "triangular_recovery_regimes_diagnostic.py"
        )
        completed = subprocess.run(
            (
                sys.executable,
                str(module),
                "--vertices",
                "10000",
                "--p",
                "0.81",
            ),
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["diagnostic_status"], DIAGNOSTIC_STATUS)
        self.assertEqual(payload["degree"], 6)
        self.assertFalse(payload["sufficiency_claimed"])

    def test_invalid_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            majority_oracle_error_degree_six(0.49)
        with self.assertRaises(ValueError):
            majority_oracle_error_polynomial(1.01)
        with self.assertRaises(ValueError):
            run_diagnostic(vertex_count=0, p=0.81)


if __name__ == "__main__":
    unittest.main()
