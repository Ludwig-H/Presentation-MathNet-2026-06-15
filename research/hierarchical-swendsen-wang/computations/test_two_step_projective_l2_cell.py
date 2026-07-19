"""Tests for the exact two-step projective L2 cell."""

from __future__ import annotations

import math
import unittest
from math import fsum

from two_step_projective_l2_cell import build_preselected_witness


class TwoStepProjectiveL2CellTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_preselected_witness()

    def test_declared_three_bit_geometry_is_real_and_reproducible(self) -> None:
        geometry = self.result.geometry
        self.assertEqual((geometry.pair_first, geometry.pair_second), (0, 13))
        self.assertEqual((geometry.lower_node, geometry.upper_node), (21, 23))
        self.assertEqual(
            (geometry.spine_size, geometry.attachment_size),
            (1, 4),
        )
        self.assertEqual(geometry.ancestral_sibling_size, 1)
        self.assertEqual(geometry.first_flip_rank, 2)
        self.assertEqual(geometry.cumulative_flip_rank, 3)
        self.assertLess(geometry.lower_rank, geometry.upper_rank)
        self.assertTrue(geometry.upper_is_strictly_before_lca)
        self.assertTrue(geometry.ancestral_sibling_contains_no_endpoint)
        self.assertTrue(self.result.witness_selected_after_exploratory_scan)
        self.assertFalse(self.result.unbiased_palm_sample)

    def test_second_loss_is_on_the_propagated_function(self) -> None:
        result = self.result
        self.assertAlmostEqual(result.initial_norm_square, 1.0)
        self.assertAlmostEqual(result.first_absolute_loss, 0.030873433894622826)
        self.assertAlmostEqual(result.second_absolute_loss, 0.13225611899017553)
        self.assertAlmostEqual(result.second_relative_loss, 0.13646939792567278)
        self.assertAlmostEqual(
            result.first_absolute_loss,
            result.first_difference_norm_square,
        )
        self.assertAlmostEqual(
            result.second_absolute_loss,
            result.second_difference_norm_square,
        )
        self.assertLess(result.maximum_pythagorean_error, 1e-12)
        self.assertLess(
            result.direct_ancestral_sibling_loss_on_pair_character,
            1e-12,
        )
        self.assertTrue(result.second_projection_applied_to_first_output)

    def test_exterior_potentials_are_exact_attained_coset_potentials(self) -> None:
        result = self.result
        orbits = result.potential_orbits
        self.assertEqual(result.attainable_potential_count, 128)
        self.assertEqual(result.strictly_positive_exterior_potential_count, 64)
        self.assertEqual(result.boundary_exterior_potential_count, 64)
        self.assertAlmostEqual(fsum(orbit.posterior_mass for orbit in orbits), 1.0)
        self.assertTrue(all(orbit.posterior_mass > 0.0 for orbit in orbits))
        self.assertTrue(
            all(max(orbit.exterior_projective_log_weights) == 0.0 for orbit in orbits)
        )
        self.assertTrue(
            all(
                all(
                    math.isfinite(value)
                    for value in orbit.exterior_projective_log_weights
                )
                for orbit in orbits
                if orbit.strictly_positive_exterior
            )
        )
        self.assertTrue(
            all(
                any(
                    not math.isfinite(value)
                    for value in orbit.exterior_projective_log_weights
                )
                for orbit in orbits
                if not orbit.strictly_positive_exterior
            )
        )
        self.assertLess(result.maximum_projective_factorization_error, 1e-12)
        self.assertLess(result.orbit_second_loss_audit_error, 1e-12)
        self.assertTrue(result.exterior_potentials_are_attained)
        self.assertTrue(result.full_posterior_enumerated)

    def test_interior_margin_is_finite_cell_only(self) -> None:
        result = self.result
        self.assertGreater(result.strictly_positive_exterior_posterior_mass, 0.9)
        self.assertGreater(result.strictly_positive_exterior_energy_ratio, 0.1)
        self.assertGreater(
            result.minimum_strictly_positive_exterior_second_relative_loss,
            0.003,
        )
        self.assertEqual(result.minimum_all_attainable_second_relative_loss, 0.0)
        self.assertGreater(result.fraction_incoming_energy_on_positive_second_loss, 0.9)
        self.assertFalse(result.weak_recovery_claimed)


if __name__ == "__main__":
    unittest.main()
