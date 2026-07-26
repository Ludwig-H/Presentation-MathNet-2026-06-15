"""Tests for the symmetric-SBM PGW broadcast calibration."""

from __future__ import annotations

import math
import unittest

from sbm_broadcast_density_evolution import (
    PASS_STATUS,
    deterministic_regime,
    information_percolation_trajectory,
    linear_score_lower_bound,
    posterior_update,
    run_diagnostic,
)


class DeterministicBroadcastBoundsTests(unittest.TestCase):
    def test_linear_score_formula_in_all_three_regimes(self) -> None:
        self.assertEqual(linear_score_lower_bound(0.0, 0), 1.0)
        self.assertEqual(linear_score_lower_bound(0.0, 3), 0.0)
        self.assertAlmostEqual(
            linear_score_lower_bound(0.8, 5),
            1.0 / sum(0.8 ** (-step) for step in range(6)),
        )
        self.assertAlmostEqual(
            linear_score_lower_bound(1.0, 9),
            0.1,
        )
        self.assertAlmostEqual(
            linear_score_lower_bound(1.2, 5),
            1.0 / sum(1.2 ** (-step) for step in range(6)),
        )

    def test_information_percolation_recurrence(self) -> None:
        trajectory = information_percolation_trajectory(1.0, 3)
        self.assertEqual(trajectory[0], 1.0)
        self.assertAlmostEqual(trajectory[1], 1.0 - math.exp(-1.0))
        self.assertAlmostEqual(
            trajectory[2],
            1.0 - math.exp(-trajectory[1]),
        )
        self.assertTrue(
            all(
                following < preceding
                for preceding, following in zip(trajectory, trajectory[1:])
            )
        )

    def test_bounds_certify_the_exact_transition(self) -> None:
        self.assertLess(
            information_percolation_trajectory(0.95, 400)[-1],
            1e-8,
        )
        self.assertLess(
            information_percolation_trajectory(1.0, 300)[-1],
            0.01,
        )
        self.assertAlmostEqual(
            linear_score_lower_bound(1.2, 500),
            (1.2 - 1.0) / 1.2,
        )
        self.assertIn("NONRECONSTRUCTION_CERTIFIED", deterministic_regime(1.0))
        self.assertIn("RECONSTRUCTION_CERTIFIED", deterministic_regime(1.01))

    def test_exact_posterior_update(self) -> None:
        messages = (0.2, -0.4, 0.1)
        expected = math.tanh(
            sum(math.atanh(0.5 * message) for message in messages)
        )
        self.assertAlmostEqual(posterior_update(messages, 0.5), expected)
        self.assertEqual(posterior_update((), 0.5), 0.0)
        self.assertEqual(posterior_update((1.0,), 1.0), 1.0)


class PopulationBroadcastDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.arguments = dict(
            degree=3.0,
            signal_to_noise_values=(0.8, 1.2),
            max_depth=3,
            particles=800,
            batches=4,
            seed=20260726,
        )
        cls.summary = run_diagnostic(**cls.arguments)

    def test_seeded_population_dynamics_is_exactly_reproducible(self) -> None:
        self.assertEqual(self.summary, run_diagnostic(**self.arguments))

    def test_output_contains_requested_audits_and_nonclaim(self) -> None:
        self.assertEqual(self.summary.diagnostic_status, PASS_STATUS)
        self.assertFalse(self.summary.graph_weak_recovery_claimed)
        self.assertEqual(len(self.summary.lambda_diagnostics), 2)
        for diagnostic in self.summary.lambda_diagnostics:
            self.assertEqual(len(diagnostic.estimates), 4)
            self.assertEqual(
                diagnostic.empirical_consistency_status,
                PASS_STATUS,
            )
            for estimate in diagnostic.estimates:
                self.assertGreaterEqual(estimate.q_hat, 0.0)
                self.assertLessEqual(estimate.q_hat, 1.0)
                self.assertLessEqual(estimate.ell_t, estimate.r_t)
                self.assertTrue(
                    estimate.nishimori_consistent_up_to_tolerance
                )
                self.assertTrue(estimate.sandwich_consistent_up_to_tolerance)

    def test_parameter_validation(self) -> None:
        with self.assertRaises(ValueError):
            run_diagnostic(
                degree=3.0,
                signal_to_noise_values=(3.1,),
                max_depth=1,
                particles=10,
                batches=2,
                seed=1,
            )
        with self.assertRaises(ValueError):
            run_diagnostic(
                degree=3.0,
                signal_to_noise_values=(1.0,),
                max_depth=1,
                particles=2,
                batches=3,
                seed=1,
            )


if __name__ == "__main__":
    unittest.main()
