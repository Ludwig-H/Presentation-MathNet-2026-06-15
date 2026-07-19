"""Tests for the unselected exact two-step population diagnostic."""

from __future__ import annotations

import unittest

from two_step_l2_population_diagnostic import run_population_diagnostic


class TwoStepL2PopulationDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.summary, cls.cells = run_population_diagnostic(
            repetitions=1,
            p=0.805,
            distance_fraction=0.5,
            seed=20260729,
        )

    def test_pair_sampling_and_cells_are_unselected(self) -> None:
        self.assertEqual(self.summary.connected_pair_count, 1)
        self.assertGreater(self.summary.eligible_cell_count, 0)
        self.assertEqual(len(self.cells), self.summary.eligible_cell_count)
        self.assertTrue(self.summary.pair_sample_is_independent_of_dendrogram)
        self.assertFalse(self.summary.witness_selected_after_scan)
        self.assertTrue(all(cell.pair_distance == 2 for cell in self.cells))
        self.assertEqual(self.summary.critical_rank_window, 0.02)

    def test_energy_statistics_are_probability_like(self) -> None:
        for value in (
            self.summary.energy_weighted_second_relative_loss,
            self.summary.positive_second_loss_fraction,
            self.summary.second_loss_above_one_percent_fraction,
            self.summary.top_ten_percent_second_absolute_loss_share,
            self.summary.top_twenty_percent_second_absolute_loss_share,
            self.summary.energy_weighted_strictly_positive_exterior_share,
            self.summary.boundary_zero_margin_cell_fraction,
            self.summary.near_critical_cell_fraction,
            self.summary.near_critical_incoming_energy_share,
            self.summary.near_critical_second_absolute_loss_share,
            self.summary.near_critical_boundary_zero_margin_cell_fraction,
            self.summary.near_critical_energy_weighted_strictly_positive_exterior_share,
        ):
            self.assertIsNotNone(value)
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)

    def test_exact_audits_hold_for_every_cell(self) -> None:
        self.assertTrue(self.summary.exact_full_posterior_enumeration)
        self.assertLess(self.summary.maximum_pythagorean_error, 1e-11)
        self.assertLess(self.summary.maximum_projective_factorization_error, 1e-11)
        self.assertTrue(
            all(cell.maximum_pythagorean_error < 1e-11 for cell in self.cells)
        )
        self.assertTrue(
            all(
                cell.maximum_projective_factorization_error < 1e-11
                for cell in self.cells
            )
        )
        self.assertFalse(self.summary.weak_recovery_claimed)


if __name__ == "__main__":
    unittest.main()
