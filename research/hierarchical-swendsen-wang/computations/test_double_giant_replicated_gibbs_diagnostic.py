"""Tests for the exact double-giant replicated Gibbs diagnostic."""

from __future__ import annotations

import random
import unittest
from math import fsum

from critical_pair_path_geometry import triangular_torus_edges
from double_giant_replicated_gibbs_diagnostic import (
    BASELINE_P,
    ReplicaHierarchyDiagnostic,
    _ReplicaWork,
    analyze_replica_pair,
    conditional_rank_audit,
    edge_relation,
    enumerate_observation_posterior,
    pair_persistence,
    physical_pair_correlation_matrix,
    run_diagnostic,
    sample_conditional_ranked_edges,
    walsh_transform,
)


class DoubleGiantReplicatedGibbsDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.summary, cls.observations = run_diagnostic(
            side_length=4,
            p=BASELINE_P,
            observation_count=1,
            replica_pairs_per_observation=1,
            seed=3801,
        )
        cls.observation = cls.observations[0]
        cls.pair = cls.observation.replica_pairs[0]

    def test_exact_observation_posterior_and_direct_target(self) -> None:
        signs = tuple(
            1 if index % 5 else -1
            for index, _ in enumerate(triangular_torus_edges(4))
        )
        posterior = enumerate_observation_posterior(4, signs, BASELINE_P)
        self.assertEqual(len(posterior), 1 << 16)
        self.assertAlmostEqual(fsum(posterior), 1.0)
        complement = (1 << 16) - 1
        self.assertLess(
            max(
                abs(posterior[state] - posterior[state ^ complement])
                for state in range(1 << 16)
            ),
            1e-15,
        )
        matrix = physical_pair_correlation_matrix(posterior, 0)
        self.assertTrue(
            all(matrix[index][index] == 1.0 for index in range(16))
        )
        self.assertGreaterEqual(pair_persistence(matrix), 1.0 / 16.0)
        self.assertLessEqual(pair_persistence(matrix), 1.0)

    def test_conditional_clock_ranks_encode_the_reference_gauge(self) -> None:
        edges = triangular_torus_edges(4)
        signs = tuple(1 if index % 3 else -1 for index, _ in enumerate(edges))
        reference_state = 0b1010010110010110
        first = sample_conditional_ranked_edges(
            4,
            signs,
            reference_state,
            BASELINE_P,
            random.Random(91),
        )
        second = sample_conditional_ranked_edges(
            4,
            signs,
            reference_state,
            BASELINE_P,
            random.Random(91),
        )
        self.assertEqual(first, second)
        self.assertTrue(
            conditional_rank_audit(
                signs,
                reference_state,
                first,
                BASELINE_P,
            )
        )
        for sign, (rank, vertex, neighbour) in zip(signs, first, strict=True):
            satisfied = (
                sign * edge_relation(reference_state, vertex, neighbour) == 1
            )
            self.assertEqual(rank <= BASELINE_P, satisfied)

    def test_walsh_pair_matrix_matches_direct_character_sum(self) -> None:
        probabilities = (0.1, 0.2, 0.3, 0.4)
        transform = walsh_transform(probabilities)
        self.assertAlmostEqual(transform[0], 1.0)
        matrix = physical_pair_correlation_matrix(probabilities, 0)
        expected = probabilities[0] - probabilities[1] - probabilities[2] + probabilities[3]
        self.assertAlmostEqual(matrix[0][1], expected)

    def test_independent_hierarchies_close_the_root_decomposition(self) -> None:
        self.assertFalse(self.pair.shared_hierarchy_used)
        self.assertTrue(self.pair.independent_hierarchy_seeds)
        self.assertNotEqual(
            self.pair.first.hierarchy_seed,
            self.pair.second.hierarchy_seed,
        )
        self.assertTrue(
            self.pair.first.conditional_ranks_match_observation_in_reference_gauge
        )
        self.assertTrue(
            self.pair.second.conditional_ranks_match_observation_in_reference_gauge
        )
        self.assertTrue(self.pair.exact_root_intersection_decomposition_passed)
        self.assertLess(self.pair.root_intersection_decomposition_error, 5e-11)
        self.assertTrue(
            self.pair.separated_pair_products_zero_within_tolerance
        )
        self.assertLess(
            self.pair.maximum_separated_pair_product_absolute_value,
            5e-11,
        )
        self.assertAlmostEqual(
            self.pair.normalized_pair_correlation_product,
            self.pair.root_intersection_contribution_sum,
            places=10,
        )

    def test_double_giant_and_critical_refinement_statistics_are_exact(self) -> None:
        largest = next(
            item
            for item in self.pair.root_intersection_contributions
            if (
                item.first_root_index
                == self.pair.first.largest_final_component_index
                and item.second_root_index
                == self.pair.second.largest_final_component_index
            )
        )
        self.assertEqual(
            largest.vertices,
            self.pair.largest_root_intersection_vertices,
        )
        self.assertAlmostEqual(
            largest.normalized_pair_product_contribution,
            self.pair.largest_root_intersection_contribution,
        )
        self.assertEqual(
            sum(self.pair.common_critical_refinement_cell_sizes),
            16,
        )
        self.assertEqual(
            sum(
                self.pair.largest_root_intersection_refinement_cell_sizes
            ),
            self.pair.largest_root_intersection_size,
        )
        self.assertGreaterEqual(
            self.pair.common_critical_refinement_diagonal_mass,
            1.0 / 16.0,
        )
        self.assertLessEqual(
            self.pair.common_critical_refinement_diagonal_mass,
            1.0,
        )

    def test_disjoint_largest_roots_give_an_empty_double_giant(self) -> None:
        def work(
            replica_index: int,
            final_components: tuple[tuple[int, ...], ...],
        ) -> _ReplicaWork:
            diagnostic = ReplicaHierarchyDiagnostic(
                replica_index=replica_index,
                posterior_spin_seed=10 + replica_index,
                hierarchy_seed=20 + replica_index,
                reference_state=0,
                reference_spins=(1, 1),
                exact_gibbs_state_count=4,
                positive_gibbs_state_count=4,
                final_component_sizes=(1, 1),
                largest_final_component_index=0,
                largest_final_component_vertices=final_components[0],
                critical_component_sizes=(1, 1),
                pair_correlation_matrix=((1.0, 0.0), (0.0, 1.0)),
                maximum_matrix_diagonal_error=0.0,
                maximum_matrix_symmetry_error=0.0,
                maximum_cross_final_root_correlation=0.0,
                conditional_ranks_match_observation_in_reference_gauge=True,
            )
            return _ReplicaWork(
                diagnostic=diagnostic,
                final_components=final_components,
                critical_components=final_components,
            )

        pair = analyze_replica_pair(
            observation_index=0,
            replica_pair_index=0,
            first=work(0, ((0,), (1,))),
            second=work(1, ((1,), (0,))),
        )
        self.assertEqual(pair.largest_root_intersection_vertices, ())
        self.assertEqual(pair.largest_root_intersection_size, 0)
        self.assertEqual(pair.largest_root_intersection_fraction, 0.0)
        self.assertEqual(pair.largest_root_intersection_contribution, 0.0)
        self.assertEqual(
            pair.largest_root_intersection_refinement_cell_sizes,
            (),
        )

    def test_direct_q_is_exposed_as_the_unbiased_monte_carlo_target(self) -> None:
        direct_matrix = self.observation.direct_posterior_pair_correlation_matrix
        self.assertAlmostEqual(
            self.observation.direct_posterior_pair_persistence,
            pair_persistence(direct_matrix),
        )
        self.assertTrue(
            self.observation.independent_hierarchy_estimator_targets_direct_persistence
        )
        self.assertFalse(
            self.observation.equality_with_direct_target_expected_per_sample
        )
        self.assertTrue(
            self.summary.independent_hierarchy_estimator_targets_q_direct
        )
        self.assertFalse(self.summary.shared_hierarchy_used)
        self.assertFalse(self.summary.weak_recovery_claimed)
        self.assertFalse(
            self.summary.asymptotic_double_giant_reduction_claimed
        )

    def test_caps_and_non_l4_requests_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=5)
        with self.assertRaises(ValueError):
            run_diagnostic(maximum_state_count=(1 << 16) - 1)
        with self.assertRaises(ValueError):
            run_diagnostic(observation_count=2, maximum_observation_count=1)
        with self.assertRaises(ValueError):
            run_diagnostic(
                replica_pairs_per_observation=2,
                maximum_replica_pairs_per_observation=1,
            )


if __name__ == "__main__":
    unittest.main()
