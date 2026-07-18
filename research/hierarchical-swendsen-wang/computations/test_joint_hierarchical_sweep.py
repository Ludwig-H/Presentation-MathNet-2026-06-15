from __future__ import annotations

import random
import unittest
from math import exp, fsum

from joint_hierarchical_sweep import (
    HierarchicalSweepEnvironment,
    unbiased_rademacher_mean_square,
)
from critical_pair_path_geometry import sample_ranked_edges


class JointHierarchicalSweepTests(unittest.TestCase):
    def test_every_kruskal_bucket_is_assigned_exactly(self) -> None:
        rng = random.Random(4)
        environment = HierarchicalSweepEnvironment(
            4, sample_ranked_edges(4, rng), 0.81
        )
        for node in environment.forest.internal_nodes:
            self.assertEqual(
                len(environment.bucket_edges[node]),
                environment.forest.bucket_size[node],
            )

    def test_planted_state_has_positive_mass_in_every_orbit(self) -> None:
        rng = random.Random(7)
        environment = HierarchicalSweepEnvironment(
            4, sample_ranked_edges(4, rng), 0.81
        )
        for node in environment.bottom_up_order:
            self.assertGreater(
                environment.conditional_log_weight(node, 0),
                float("-inf"),
            )

    def test_local_heat_bath_agrees_with_full_conditional_density(self) -> None:
        rng = random.Random(19)
        environment = HierarchicalSweepEnvironment(
            4, sample_ranked_edges(4, rng), 0.81
        )
        states = (0,) + tuple(
            environment.sweep("bottom-up", rng) for _ in range(4)
        )
        for node in environment.bottom_up_order:
            for state in states:
                orbit = environment.proposal_probabilities(node, state)
                full_logs = tuple(
                    environment.full_log_weight(state ^ mask)
                    for mask, _ in orbit
                )
                maximum = max(full_logs)
                weights = tuple(
                    0.0 if value == float("-inf") else exp(value - maximum)
                    for value in full_logs
                )
                normalizer = fsum(weights)
                expected = tuple(weight / normalizer for weight in weights)
                for (_, actual), target in zip(orbit, expected, strict=True):
                    self.assertAlmostEqual(actual, target)

    def test_final_root_update_contains_fair_global_recolouring(self) -> None:
        rng = random.Random(23)
        environment = HierarchicalSweepEnvironment(
            4, sample_ranked_edges(4, rng), 0.81
        )
        roots = [
            node
            for node in environment.forest.internal_nodes
            if environment.forest.tree_parent[node] == -1
        ]
        self.assertTrue(roots)
        for root in roots:
            probabilities = dict(environment.proposal_probabilities(root, 0))
            left = environment.cluster_mask[
                environment.forest.left_child[root]
            ]
            right = environment.cluster_mask[
                environment.forest.right_child[root]
            ]
            self.assertAlmostEqual(probabilities[0], probabilities[left | right])
            self.assertAlmostEqual(probabilities[left], probabilities[right])

    def test_rademacher_square_estimator_is_exact_on_constant_samples(self) -> None:
        self.assertEqual(unbiased_rademacher_mean_square([1] * 10), 1.0)
        self.assertEqual(unbiased_rademacher_mean_square([-1] * 10), 1.0)

    def test_rademacher_square_estimator_removes_diagonal_terms(self) -> None:
        samples = [1, -1, 1, -1]
        self.assertAlmostEqual(
            unbiased_rademacher_mean_square(samples), -1.0 / 3.0
        )


if __name__ == "__main__":
    unittest.main()
