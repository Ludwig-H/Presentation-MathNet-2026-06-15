from __future__ import annotations

import random
import unittest
from math import sqrt

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    triangular_torus_distance,
    triangular_torus_edges,
)
from giant_component_quotient_diagnostic import (
    BASELINE_P,
    P_SW,
    critical_partition_of_largest_tree,
    diagnose_environment,
    mean_with_environment_standard_error,
    run_diagnostic,
    sample_uniform_far_cross_block_pair,
    upper_quotient_path_geometry,
)


def _three_block_forest() -> CriticalKruskalForest:
    """Build {0,1}, {2,3}, {4,5}, then merge above q_c."""

    selected_ranks = {
        (0, 1): 0.10,
        (2, 3): 0.20,
        (4, 5): 0.30,
        (1, 2): 0.40,
        (0, 4): 0.50,
    }
    ranked_edges = tuple(
        (selected_ranks.get(edge, 0.99), *edge)
        for edge in triangular_torus_edges(4)
    )
    return CriticalKruskalForest(
        4,
        ranked_edges,
        critical_rank=2.0 * BASELINE_P - 1.0,
    )


class GiantComponentQuotientDiagnosticTests(unittest.TestCase):
    def test_largest_final_tree_is_cut_into_expected_critical_blocks(self) -> None:
        partition = critical_partition_of_largest_tree(
            _three_block_forest()
        )
        self.assertEqual(
            partition.blocks,
            ((0, 1), (2, 3), (4, 5)),
        )
        self.assertEqual(len(partition.component), 6)
        self.assertEqual(partition.block_size(0), 2)
        self.assertEqual(partition.block_size(5), 2)

    def test_quotient_path_keeps_only_strictly_postcritical_nodes(self) -> None:
        forest = _three_block_forest()
        partition = critical_partition_of_largest_tree(forest)
        path = upper_quotient_path_geometry(
            forest, partition, first=0, second=4
        )
        self.assertEqual(path.bucket_ranks, (0.40, 0.50))
        self.assertEqual(path.path_length, 2)
        self.assertAlmostEqual(path.lca_rank, 0.50)
        self.assertTrue(all(rank > Q_CRITICAL for rank in path.bucket_ranks))
        self.assertEqual(path.first_critical_block_size, 2)
        self.assertEqual(path.second_critical_block_size, 2)

    def test_pair_sampler_enforces_distance_and_distinct_blocks(self) -> None:
        partition = critical_partition_of_largest_tree(
            _three_block_forest()
        )
        rng = random.Random(123)
        for _ in range(100):
            first, second = sample_uniform_far_cross_block_pair(
                partition,
                side_length=4,
                distance_fraction=0.25,
                rng=rng,
            )
            self.assertNotEqual(
                partition.block_index(first),
                partition.block_index(second),
            )
            self.assertGreaterEqual(
                triangular_torus_distance(first, second, 4),
                1.0,
            )

    def test_environment_diagnostic_labels_a_finite_local_oracle(self) -> None:
        forest = _three_block_forest()
        result = diagnose_environment(
            forest=forest,
            repetition=0,
            p=BASELINE_P,
            pairs_per_environment=12,
            distance_fraction=0.25,
            rng=random.Random(7),
        )
        self.assertEqual(result.pair_count, 12)
        self.assertGreaterEqual(
            result.second_largest_final_component_fraction, 0.0
        )
        self.assertGreater(
            result.same_critical_block_pair_probability_in_giant, 0.0
        )
        self.assertAlmostEqual(
            result.same_critical_block_pair_probability_in_giant,
            0.2,
        )
        self.assertGreater(result.mean_upper_quotient_node_count, 0.0)
        self.assertGreaterEqual(
            result.mean_path_fac_local_oracle_attenuation, 0.0
        )
        self.assertGreaterEqual(
            result.mean_path_fac_local_oracle_correlation, 0.0
        )
        self.assertLessEqual(
            result.mean_path_fac_local_oracle_correlation, 1.0
        )

    def test_standard_error_clusters_environment_means(self) -> None:
        estimate = mean_with_environment_standard_error((1.0, 3.0))
        self.assertEqual(estimate.mean, 2.0)
        self.assertAlmostEqual(estimate.standard_error, 1.0)
        singleton = mean_with_environment_standard_error((sqrt(2.0),))
        self.assertIsNone(singleton.standard_error)

    def test_seeded_summary_is_deterministic_and_explicitly_nonproof(self) -> None:
        arguments = dict(
            side_length=8,
            repetitions=3,
            pairs_per_environment=5,
            p=BASELINE_P,
            distance_fraction=0.25,
            seed=20260726,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(first, second)
        self.assertEqual(first.pair_count, 15)
        self.assertEqual(first.eligible_environment_count, 3)
        self.assertEqual(first.ineligible_environment_count, 0)
        self.assertTrue(first.estimates_conditioned_on_eligible_environments)
        self.assertFalse(first.weak_recovery_claimed)
        self.assertIn("non-proof", first.path_fac_status)
        self.assertIsNotNone(
            first.geometry_estimates[
                "mean_upper_quotient_node_count"
            ].standard_error
        )

    def test_nonsupercritical_p_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=8,
                repetitions=1,
                pairs_per_environment=1,
                p=P_SW,
                distance_fraction=0.25,
                seed=1,
            )

    def test_near_critical_ineligible_environments_are_counted(self) -> None:
        result = run_diagnostic(
            side_length=8,
            repetitions=1,
            pairs_per_environment=1,
            p=P_SW + 1e-8,
            distance_fraction=0.25,
            seed=1,
        )
        self.assertEqual(result.eligible_environment_count, 0)
        self.assertEqual(result.ineligible_environment_count, 1)
        self.assertEqual(result.pair_count, 0)
        self.assertEqual(result.geometry_estimates, {})
        self.assertEqual(
            result.path_fac_local_oracle_nonproof_estimates,
            {},
        )


if __name__ == "__main__":
    unittest.main()
