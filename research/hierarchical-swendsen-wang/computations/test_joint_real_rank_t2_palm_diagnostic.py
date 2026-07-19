from __future__ import annotations

import random
import unittest
from dataclasses import asdict

from ancestral_polarization_palm_diagnostic import analyze_node_polarization
from corridor_t2_signature_diagnostic import event_palm_corridor_observations
from critical_pair_path_geometry import triangular_torus_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from joint_real_rank_t2_palm_diagnostic import (
    analyze_joint_pair,
    run_diagnostic,
    run_size_series,
)


def _deterministic_comb_edges(
    side_length: int = 4,
) -> tuple[tuple[float, int, int], ...]:
    """Create a rank-0.4 strict-arm merger under a rank-0.55 LCA."""

    special_ranks = {
        tuple(sorted((0, 1))): 0.01,
        tuple(sorted((0, side_length))): 0.4,
        tuple(sorted((1, side_length + 1))): 0.55,
    }
    result = []
    for index, edge in enumerate(triangular_torus_edges(side_length)):
        result.append((special_ranks.get(edge, 0.7 + 0.003 * index), *edge))
    if len({rank for rank, _, _ in result}) != len(result):
        raise AssertionError("deterministic ranks must be distinct")
    return tuple(result)


class JointRealRankT2PalmDiagnosticTests(unittest.TestCase):
    def test_same_oriented_node_uses_real_charge_and_actual_message(self) -> None:
        ranked_edges = _deterministic_comb_edges()
        environment = HierarchicalSweepEnvironment(4, ranked_edges, 0.805)
        pairs, audit = event_palm_corridor_observations(
            environment.forest,
            ranked_edges,
            repetition=0,
            p=0.805,
            distance_fraction=0.1,
            rng=random.Random(5),
        )
        self.assertTrue(audit.passed)
        final_pair = next(
            pair
            for pair in pairs
            if environment.forest.merge_rank[pair.lca_node] == 0.55
        )
        result = analyze_joint_pair(
            environment,
            final_pair,
            maximum_bucket_size=8,
            maximum_attachment_size=2,
            maximum_charge=0.28,
            message_thresholds=(1.0, 2.0, 4.0),
        )
        target = next(node for node in result.nodes if node.actual_rank == 0.4)
        self.assertEqual(target.node, environment.forest.tree_lca(0, 4))
        self.assertEqual(target.bucket_size, 2)
        self.assertEqual(target.attachment_child_size, 2)
        self.assertAlmostEqual(target.actual_charge, 0.245)
        self.assertGreater(target.criticalized_charge_proxy, 0.28)
        self.assertTrue(target.passes_actual_charge)
        self.assertFalse(target.passes_criticalized_charge_proxy)
        direct_message = analyze_node_polarization(environment, target.node)
        self.assertEqual(
            target.actual_external_message,
            direct_message.actual_external_message,
        )
        self.assertEqual(
            target.actual_identity_error,
            direct_message.actual_identity_error,
        )
        self.assertEqual(result.palm_weight, final_pair.n_rho)
        self.assertTrue(result.port_count_is_geometric_proxy)
        self.assertFalse(result.screening_computed)
        self.assertFalse(result.blackwell_domination_assumed)

    def test_reproducible_palm_run_has_nested_counts_and_all_audits(self) -> None:
        arguments = dict(
            side_length=6,
            repetitions=3,
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=8,
            maximum_attachment_size=4,
            maximum_charge=1.0,
            message_thresholds=(1.0, 2.0, 4.0),
            seed=5821,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first), asdict(second))
        self.assertGreater(first.total_palm_weight, 0)
        self.assertTrue(first.reconstruction_audit.passed)
        self.assertTrue(first.palm_partition_audit.passed)
        self.assertTrue(first.joint_consistency_audit.passed)
        self.assertIsNotNone(
            first.jackknife_standard_error_actual_charge_candidate_count
        )
        self.assertIsNotNone(
            first.jackknife_standard_error_global_port_count_per_corridor_node
        )
        actual = first.weighted_mean_actual_charge_and_bounded_message_counts
        proxy = first.weighted_mean_criticalized_charge_proxy_and_bounded_message_counts
        self.assertLessEqual(actual[0], actual[1])
        self.assertLessEqual(actual[1], actual[2])
        self.assertLessEqual(
            actual[-1], first.weighted_mean_actual_charge_candidate_count
        )
        self.assertLessEqual(proxy[0], proxy[1])
        self.assertLessEqual(proxy[1], proxy[2])
        self.assertLessEqual(
            proxy[-1],
            first.weighted_mean_criticalized_charge_proxy_candidate_count,
        )
        for proxy_value, actual_value in zip(proxy, actual, strict=True):
            self.assertLessEqual(proxy_value, actual_value)
        self.assertFalse(first.blackwell_domination_assumed)
        self.assertFalse(first.transfer_deficit_computed)
        self.assertFalse(first.weak_recovery_claimed)
        self.assertIn("no domination", first.interpretation)

    def test_size_series_uses_explicit_efforts_and_distinct_streams(self) -> None:
        result = run_size_series(
            (6, 6),
            (2, 2),
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=6,
            maximum_attachment_size=3,
            maximum_charge=1.0,
            seed=19,
        )
        self.assertEqual(tuple(item.repetitions for item in result), (2, 2))
        self.assertEqual(result[0].seed, 19)
        self.assertEqual(result[1].seed, 1_000_022)
        self.assertNotEqual(asdict(result[0]), asdict(result[1]))

    def test_invalid_filters_and_effort_vectors_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=6, repetitions=0)
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=6,
                repetitions=1,
                maximum_attachment_size=0,
            )
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=6,
                repetitions=1,
                message_thresholds=(2.0, 1.0),
            )
        with self.assertRaises(ValueError):
            run_size_series((6, 8), (2,))


if __name__ == "__main__":
    unittest.main()
