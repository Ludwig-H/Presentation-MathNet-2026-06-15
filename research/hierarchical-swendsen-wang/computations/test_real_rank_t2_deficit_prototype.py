"""Tests for the faithful/projected real-rank T2 deficit audit."""

from __future__ import annotations

import unittest
from math import exp, isinf

from critical_band_thresholds import Q_CRITICAL
from real_rank_t2_deficit_prototype import (
    REPLICATED_RELATIVE_STATES,
    find_sample_prototype,
)


class RealRankT2DeficitPrototypeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = find_sample_prototype(
            side_length=8,
            p=0.805,
            seed=20260724,
        )

    def test_real_ranks_and_full_exterior_vectors_are_retained(self) -> None:
        result = self.result
        geometry = result.geometry
        self.assertFalse(geometry.criticalization_used)
        self.assertFalse(geometry.winning_edge_identity_used)
        self.assertTrue(geometry.all_factor_ranks_are_realized)
        self.assertLess(geometry.actual_rank, 2.0 * geometry.p - 1.0)
        self.assertTrue(geometry.is_postcritical)
        self.assertGreaterEqual(geometry.actual_rank, Q_CRITICAL)
        self.assertEqual(geometry.bucket_size, 5)
        self.assertEqual(geometry.attachment_child_size, 1)

        for boundary in result.replica_exterior_states:
            self.assertEqual(len(boundary.exterior_projective_log_weights), 4)
            self.assertEqual(len(boundary.full_projective_log_weights), 4)
            self.assertEqual(
                tuple(factor.node for factor in boundary.strict_ancestor_factors),
                geometry.strict_ancestor_nodes,
            )
            self.assertLess(boundary.heat_bath_audit_error, 3e-12)
            for factor in (
                boundary.local_factor,
                *boundary.strict_ancestor_factors,
            ):
                self.assertEqual(len(factor.satisfied_counts), 4)
                self.assertEqual(len(factor.log_weights), 4)
                self.assertGreaterEqual(factor.rank, 0.0)
                self.assertLess(factor.rank, geometry.p)

        # Equal reference replicas share the complete observed environment;
        # only their subsequent heat-bath draws are independent.
        self.assertEqual(
            result.replica_exterior_states[0],
            result.replica_exterior_states[1],
        )
        self.assertTrue(result.shared_nishimori_environment)
        self.assertTrue(result.replica_draws_conditionally_independent)

    def test_faithful_transfer_has_exactly_zero_one_step_deficit(self) -> None:
        result = self.result
        mass = result.faithful.transfer.mass[0]
        twisted = result.faithful.transfer.twisted[0]
        ratio = result.faithful.transfer.ratio[0]
        self.assertAlmostEqual(sum(mass), 1.0, places=14)
        self.assertAlmostEqual(result.faithful.feynman_kac_envelope, 1.0, places=14)
        self.assertAlmostEqual(result.faithful.logarithmic_attenuation, 0.0, places=14)
        self.assertLess(result.faithful.finite_horizon_doob_envelope_error, 1e-15)
        self.assertTrue(result.faithful.finite_horizon_doob_inequality_holds)
        self.assertEqual(len(result.faithful_targets), 16)

        for target, k_value, u_value, local_ratio, deficit in zip(
            result.faithful_targets,
            mass,
            twisted,
            ratio,
            result.faithful.transition_deficits,
            strict=True,
        ):
            self.assertAlmostEqual(
                u_value, target.replicated_twist * k_value, places=15
            )
            if k_value > 0.0:
                self.assertEqual(local_ratio, 1.0)
                self.assertAlmostEqual(deficit, 0.0, places=15)
            else:
                self.assertEqual(local_ratio, 0.0)
                self.assertTrue(isinf(deficit))

    def test_two_replicas_use_product_draws_in_one_environment(self) -> None:
        result = self.result
        first_probabilities = result.replica_exterior_states[0].heat_bath_probabilities
        second_probabilities = result.replica_exterior_states[1].heat_bath_probabilities
        for target, mass in zip(
            result.faithful_targets,
            result.faithful.transfer.mass[0],
            strict=True,
        ):
            expected = (
                first_probabilities[target.first_orbit_index]
                * second_probabilities[target.second_orbit_index]
            )
            self.assertAlmostEqual(mass, expected, places=15)

    def test_projected_numbers_are_exact_aggregates_but_not_closed(self) -> None:
        result = self.result
        expected_mass = {state: 0.0 for state in REPLICATED_RELATIVE_STATES}
        expected_twisted = {state: 0.0 for state in REPLICATED_RELATIVE_STATES}
        for target, mass, twisted in zip(
            result.faithful_targets,
            result.faithful.transfer.mass[0],
            result.faithful.transfer.twisted[0],
            strict=True,
        ):
            state = (
                target.first_relative_orientation,
                target.second_relative_orientation,
            )
            expected_mass[state] += mass
            expected_twisted[state] += twisted

        for index, state in enumerate(REPLICATED_RELATIVE_STATES):
            self.assertAlmostEqual(
                result.projected.transfer.mass[0][index],
                expected_mass[state],
                places=14,
            )
            self.assertAlmostEqual(
                result.projected.transfer.twisted[0][index],
                expected_twisted[state],
                places=14,
            )
            self.assertLessEqual(
                abs(result.projected.transfer.twisted[0][index]),
                result.projected.transfer.mass[0][index] + 1e-15,
            )

        # Reproducible, falsifiable one-cell diagnostic.  Its nonzero value
        # cannot be composed because the projection discarded the global
        # orientations that transform the exterior potential.
        self.assertAlmostEqual(
            result.projected.feynman_kac_envelope,
            0.9917743479761344,
            places=13,
        )
        self.assertAlmostEqual(
            result.projected.logarithmic_attenuation,
            0.008259669371149435,
            places=13,
        )
        self.assertAlmostEqual(
            exp(-result.projected.logarithmic_attenuation),
            result.projected.feynman_kac_envelope,
            places=14,
        )
        self.assertLess(result.projected.finite_horizon_doob_envelope_error, 1e-15)
        self.assertTrue(result.projected.finite_horizon_doob_inequality_holds)
        self.assertFalse(result.projected_boundary_is_markov_closed)
        self.assertFalse(result.composable_t2_deficit_certified)

    def test_scope_flags_rule_out_palm_and_weak_recovery_overclaims(self) -> None:
        result = self.result
        self.assertTrue(result.full_boundary_retained_by_faithful_transfer)
        self.assertFalse(result.palm_abundance_estimated)
        self.assertFalse(result.weak_recovery_claimed)
        self.assertIn("zero deficit", result.interpretation)
        self.assertIn("not composable", result.interpretation)

    def test_search_input_validation(self) -> None:
        with self.assertRaises(ValueError):
            find_sample_prototype(side_length=3)
        with self.assertRaises(ValueError):
            find_sample_prototype(p=0.5)
        with self.assertRaises(ValueError):
            find_sample_prototype(maximum_environments=0)
        with self.assertRaises(ValueError):
            find_sample_prototype(minimum_rank=1.0)


if __name__ == "__main__":
    unittest.main()
