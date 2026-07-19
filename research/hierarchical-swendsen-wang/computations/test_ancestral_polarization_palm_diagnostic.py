from __future__ import annotations

import math
import unittest
from dataclasses import asdict

from ancestral_polarization_palm_diagnostic import (
    analyze_node_polarization,
    parity_log_odds,
    run_diagnostic,
)
from critical_pair_path_geometry import triangular_torus_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


def _ranked_hierarchy_with_ancestor(
    side_length: int = 4,
) -> tuple[tuple[float, int, int], ...]:
    """Create {0,1}, then add 4, then add 5 as a strict ancestor."""

    selected = {
        tuple(sorted((0, 1))): 0.01,
        tuple(sorted((0, side_length))): 0.50,
        tuple(sorted((side_length, side_length + 1))): 0.55,
    }
    return tuple(
        (selected.get(edge, 0.99), *edge)
        for edge in triangular_torus_edges(side_length)
    )


class AncestralPolarizationPalmDiagnosticTests(unittest.TestCase):
    def test_parity_log_odds_matches_direct_weight_ratio(self) -> None:
        weights = (2.0, 3.0, 5.0, 7.0)
        logs = tuple(math.log(value) for value in weights)
        self.assertAlmostEqual(
            parity_log_odds(logs),
            math.log((weights[0] + weights[3]) / (weights[1] + weights[2])),
        )
        self.assertEqual(
            parity_log_odds((0.0, float("-inf"), float("-inf"), 0.0)),
            float("inf"),
        )

    def test_node_factorization_recovers_total_log_odds(self) -> None:
        environment = HierarchicalSweepEnvironment(
            4, _ranked_hierarchy_with_ancestor(), 0.805
        )
        first_merge = environment.forest.tree_lca(0, 1)
        second_merge = environment.forest.tree_lca(0, 4)
        self.assertNotEqual(first_merge, second_merge)
        self.assertGreater(len(environment.factor_chains[second_merge]), 1)
        result = analyze_node_polarization(environment, second_merge)
        self.assertGreaterEqual(result.strict_ancestor_count, 1)
        self.assertIsNotNone(result.actual_identity_error)
        self.assertIsNotNone(result.favourable_identity_error)
        self.assertLess(result.actual_identity_error, 2e-14)
        self.assertLess(result.favourable_identity_error, 2e-14)
        self.assertAlmostEqual(
            result.actual_total_log_odds,
            result.actual_external_message + result.actual_local_log_odds,
        )
        self.assertAlmostEqual(
            result.favourable_total_log_odds,
            result.favourable_external_message + result.favourable_local_log_odds,
        )

    def test_root_node_has_zero_external_message(self) -> None:
        environment = HierarchicalSweepEnvironment(
            4, _ranked_hierarchy_with_ancestor(), 0.805
        )
        root = environment.forest.tree_lca(0, 5)
        result = analyze_node_polarization(environment, root)
        self.assertEqual(result.strict_ancestor_count, 0)
        self.assertAlmostEqual(result.actual_external_message, 0.0)
        self.assertAlmostEqual(result.favourable_external_message, 0.0)

    def test_small_palm_run_is_reproducible_and_thresholds_are_nested(
        self,
    ) -> None:
        arguments = dict(
            side_length=6,
            repetitions=3,
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=8,
            maximum_charge=1.0,
            message_thresholds=(1.0, 2.0, 4.0),
            seed=7291,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first), asdict(second))
        self.assertGreater(first.total_palm_weight, 0)
        self.assertGreater(first.weighted_mean_candidate_count, 0.0)
        self.assertLess(first.maximum_finite_log_odds_identity_error, 5e-13)
        for counts in (
            first.weighted_mean_actual_bounded_message_counts,
            first.weighted_mean_favourable_bounded_message_counts,
        ):
            self.assertLessEqual(counts[0], counts[1])
            self.assertLessEqual(counts[1], counts[2])
            self.assertLessEqual(counts[2], first.weighted_mean_candidate_count)
        self.assertIn("not lateral screening", first.interpretation)

    def test_invalid_inputs_are_rejected(self) -> None:
        environment = HierarchicalSweepEnvironment(
            4, _ranked_hierarchy_with_ancestor(), 0.805
        )
        with self.assertRaises(ValueError):
            analyze_node_polarization(environment, 0)
        with self.assertRaises(ValueError):
            parity_log_odds((0.0, 0.0))
        with self.assertRaises(ValueError):
            parity_log_odds(
                (
                    float("-inf"),
                    float("-inf"),
                    float("-inf"),
                    float("-inf"),
                )
            )
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=6,
                repetitions=1,
                p=0.805,
                distance_fraction=0.2,
                maximum_bucket_size=8,
                maximum_charge=1.0,
                message_thresholds=(2.0, 1.0),
                seed=1,
            )


if __name__ == "__main__":
    unittest.main()
