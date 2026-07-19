"""Tests for incidence last-use bounds on event-Palm attachments."""

from __future__ import annotations

import random
import unittest
from dataclasses import asdict

from corridor_t2_signature_diagnostic import event_palm_corridor_observations
from critical_pair_path_geometry import triangular_torus_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from last_use_attachment_palm_diagnostic import (
    analyze_attachment_last_use,
    run_diagnostic,
)


def _deterministic_comb_edges(
    side_length: int = 4,
) -> tuple[tuple[float, int, int], ...]:
    special_ranks = {
        tuple(sorted((0, 1))): 0.01,
        tuple(sorted((0, side_length))): 0.4,
        tuple(sorted((1, side_length + 1))): 0.55,
    }
    result = []
    for index, edge in enumerate(triangular_torus_edges(side_length)):
        result.append((special_ranks.get(edge, 0.7 + 0.003 * index), *edge))
    if len({rank for rank, _, _ in result}) != len(result):
        raise AssertionError("test ranks must be distinct")
    return tuple(result)


class LastUseAttachmentPalmDiagnosticTests(unittest.TestCase):
    def test_comb_incidence_bound_is_reconstructed_from_edge_lcas(self) -> None:
        ranked_edges = _deterministic_comb_edges()
        environment = HierarchicalSweepEnvironment(4, ranked_edges, 0.805)
        observations, geometry_audit = event_palm_corridor_observations(
            environment.forest,
            ranked_edges,
            repetition=0,
            p=0.805,
            distance_fraction=0.1,
            rng=random.Random(5),
        )
        self.assertTrue(geometry_audit.passed)
        pair = next(
            item
            for item in observations
            if environment.forest.merge_rank[item.lca_node] == 0.55
        )
        signature = next(
            item
            for item in pair.signatures
            if environment.forest.merge_rank[item.node] == 0.4
        )
        result, audit = analyze_attachment_last_use(environment, pair, signature)
        self.assertTrue(audit.passed)
        self.assertEqual(audit.checked_candidate_count, 1)
        self.assertEqual(result.attachment_child_size, 2)
        self.assertEqual(result.pair_lca_depth, 1)
        self.assertEqual(result.last_attachment_use_depth, 1)
        self.assertEqual(result.last_attachment_use_rank, 0.55)
        self.assertTrue(result.last_attachment_use_no_later_than_pair_lca)
        self.assertTrue(result.incidence_upper_bound_only)
        self.assertFalse(result.attachment_forgettable_within(0))
        self.assertTrue(result.attachment_forgettable_within(1))
        self.assertEqual(
            result.last_attachment_use_depth,
            max(
                item.depth_from_update
                for item in result.incidences
                if item.attachment_orientation_sensitive
            ),
        )
        self.assertTrue(
            all(
                not item.attachment_orientation_sensitive
                for item in result.incidences[result.last_attachment_use_depth + 1 :]
            )
        )
        with self.assertRaises(ValueError):
            result.attachment_forgettable_within(-1)

    def test_summary_is_reproducible_palm_weighted_and_jackknifed(self) -> None:
        arguments = dict(
            side_length=6,
            repetitions=3,
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=8,
            maximum_attachment_size=4,
            minimum_rank=0.0,
            ancestor_windows=(0, 1, 2, 4),
            seed=981,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first), asdict(second))
        self.assertTrue(first.reconstruction_audit.passed)
        self.assertTrue(first.palm_partition_audit.passed)
        self.assertGreater(first.weighted_candidate_mass, 0)
        self.assertIsNotNone(first.jackknife_standard_error_last_attachment_use_depth)
        self.assertIsNotNone(
            first.jackknife_standard_errors_attachment_forgettable_fractions[0]
        )
        self.assertEqual(
            tuple(sorted(first.weighted_attachment_forgettable_fractions)),
            first.weighted_attachment_forgettable_fractions,
        )
        for attachment, union in zip(
            first.weighted_attachment_forgettable_fractions,
            first.weighted_union_forgettable_fractions,
            strict=True,
        ):
            self.assertGreaterEqual(attachment + 1e-15, union)
        self.assertIn("N_rho only", first.palm_weight_convention)
        self.assertFalse(first.transfer_deficit_computed)
        self.assertFalse(first.weak_recovery_claimed)
        self.assertIn("upper bound", first.exact_criterion)
        self.assertIn("no projected Markov closure", first.interpretation)

    def test_invalid_windows_and_cutoffs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=6, repetitions=1, ancestor_windows=(1, 1))
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=6, repetitions=1, maximum_attachment_size=0)
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=3, repetitions=1)


if __name__ == "__main__":
    unittest.main()
