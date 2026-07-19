"""Tests for exact nested odd-sector projection dissipation."""

from __future__ import annotations

import random
import unittest
from dataclasses import asdict

from critical_pair_path_geometry import sample_ranked_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from nested_projection_l2_diagnostic import (
    analyze_pair,
    collapsed_projection,
    pair_character,
    posterior_weights,
    run_diagnostic,
    weighted_norm_square,
)


class NestedProjectionL2DiagnosticTests(unittest.TestCase):
    def test_one_heat_bath_projection_is_idempotent_and_contracting(self) -> None:
        environment = HierarchicalSweepEnvironment(
            4, sample_ranked_edges(4, random.Random(81)), 0.805
        )
        weights = posterior_weights(environment)
        values = pair_character(len(weights), 0, 1)
        node = next(iter(environment.forest.internal_nodes))
        generators = tuple(mask for mask in environment.proposal_masks(node) if mask)
        projected = collapsed_projection(values, weights, generators)
        projected_twice = collapsed_projection(projected, weights, generators)
        self.assertLessEqual(
            weighted_norm_square(projected, weights),
            weighted_norm_square(values, weights) + 1e-12,
        )
        self.assertLess(
            max(
                abs(first - second) for first, second in zip(projected, projected_twice)
            ),
            1e-12,
        )

    def test_connected_pair_satisfies_pythagoras_and_posterior_bound(self) -> None:
        rng = random.Random(991)
        for _ in range(200):
            environment = HierarchicalSweepEnvironment(
                4, sample_ranked_edges(4, rng), 0.805
            )
            connected_pair = next(
                (
                    (first, second)
                    for first in range(environment.vertex_count)
                    for second in range(first + 1, environment.vertex_count)
                    if environment.forest.connected(first, second)
                ),
                None,
            )
            if connected_pair is not None:
                break
        else:
            self.fail("unable to sample a connected pair")
        result = analyze_pair(environment, *connected_pair, repetition=0)
        self.assertTrue(result.connected_in_final_forest)
        self.assertLess(result.pythagorean_error, 1e-10)
        self.assertTrue(all(step.pythagorean_error < 1e-10 for step in result.steps))
        self.assertTrue(result.posterior_below_collapsed)
        self.assertLessEqual(
            result.posterior_mean_square,
            result.final_collapsed_persistence + 1e-12,
        )
        self.assertGreater(result.path_node_count, 0)
        self.assertGreater(len(result.steps), 0)

    def test_run_is_reproducible_and_rejects_large_tori(self) -> None:
        arguments = dict(
            side_length=4,
            repetitions=1,
            p=0.805,
            distance_fraction=0.25,
            seed=2126,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first[0]), asdict(second[0]))
        self.assertEqual(
            tuple(asdict(item) for item in first[1]),
            tuple(asdict(item) for item in second[1]),
        )
        self.assertTrue(first[0].every_posterior_below_collapsed)
        self.assertLess(first[0].maximum_pythagorean_error, 1e-10)
        self.assertLessEqual(
            first[0].positive_dissipation_environment_count,
            first[0].connected_environment_count,
        )
        if first[0].second_pre_lca_energy_weighted_relative_loss is not None:
            self.assertGreaterEqual(
                first[0].second_pre_lca_energy_weighted_relative_loss, 0.0
            )
            self.assertLessEqual(
                first[0].second_pre_lca_energy_weighted_relative_loss, 1.0
            )
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=5, repetitions=1)


if __name__ == "__main__":
    unittest.main()
