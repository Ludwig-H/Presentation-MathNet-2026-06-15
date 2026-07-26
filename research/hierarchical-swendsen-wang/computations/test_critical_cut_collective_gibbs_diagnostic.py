"""Tests for the exact collective critical-cut Gibbs diagnostic."""

from __future__ import annotations

import unittest
from dataclasses import asdict

from critical_cut_collective_gibbs_diagnostic import (
    ROOT_CAP_SKIP,
    TOTAL_CAP_SKIP,
    critical_cut_partition,
    diagnose_environment,
    run_diagnostic,
    summarize,
)
from critical_pair_path_geometry import triangular_torus_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


SIDE_LENGTH = 4
P = 0.809439


def _one_postcritical_edge() -> tuple[tuple[float, int, int], ...]:
    """Leave all critical blocks isolated, then merge vertices 0 and 1."""

    return tuple(
        (0.50 if edge == (0, 1) else 0.99, *edge)
        for edge in triangular_torus_edges(SIDE_LENGTH)
    )


def _one_critical_component() -> tuple[tuple[float, int, int], ...]:
    """Put the entire torus in one critical and final component."""

    return tuple(
        (0.01, *edge)
        for edge in triangular_torus_edges(SIDE_LENGTH)
    )


class CriticalCutCollectiveGibbsDiagnosticTests(unittest.TestCase):
    def test_partition_groups_critical_blocks_by_final_root(self) -> None:
        ranked_edges = _one_postcritical_edge()
        environment = HierarchicalSweepEnvironment(
            SIDE_LENGTH,
            ranked_edges,
            P,
        )
        partition = critical_cut_partition(environment, ranked_edges)
        self.assertEqual(len(partition.blocks), 16)
        self.assertTrue(all(size == 1 for size in partition.block_sizes))
        block_counts = sorted(
            len(indices)
            for indices in partition.final_component_block_indices
        )
        self.assertEqual(block_counts, [1] * 14 + [2])

    def test_exact_root_matrices_exclude_cross_root_correlations(self) -> None:
        result = diagnose_environment(
            side_length=SIDE_LENGTH,
            ranked_edges=_one_postcritical_edge(),
            p=P,
            maximum_block_count=16,
        )
        self.assertTrue(result.used)
        self.assertIsNone(result.skip_reason)
        coupled_root = next(
            root for root in result.roots if root.final_component_size == 2
        )
        self.assertEqual(
            coupled_root.correlation_matrix,
            ((1.0, 1.0), (1.0, 1.0)),
        )
        self.assertEqual(coupled_root.orientation_state_count, 4)
        self.assertEqual(coupled_root.positive_orientation_state_count, 2)

        # The size-two root contributes four ordered pairs.  Each of the
        # fourteen other roots contributes only its diagonal singleton.
        self.assertAlmostEqual(result.collective_persistence, 18.0 / 256.0)
        self.assertAlmostEqual(
            result.critical_diagonal_persistence,
            16.0 / 256.0,
        )
        self.assertAlmostEqual(result.off_diagonal_persistence, 2.0 / 256.0)
        self.assertLess(result.maximum_decomposition_error, 1e-14)
        self.assertLess(result.maximum_diagonal_error, 1e-14)
        self.assertLess(result.maximum_symmetry_error, 1e-14)
        self.assertLess(result.maximum_global_flip_log_weight_error, 1e-14)

    def test_cap_exclusions_are_classified_and_warn_about_bias(self) -> None:
        root_excluded = diagnose_environment(
            side_length=SIDE_LENGTH,
            ranked_edges=_one_postcritical_edge(),
            p=P,
            maximum_block_count=1,
            repetition=0,
        )
        total_excluded = diagnose_environment(
            side_length=SIDE_LENGTH,
            ranked_edges=_one_postcritical_edge(),
            p=P,
            maximum_block_count=2,
            repetition=1,
        )
        used = diagnose_environment(
            side_length=SIDE_LENGTH,
            ranked_edges=_one_critical_component(),
            p=P,
            maximum_block_count=2,
            repetition=2,
        )
        self.assertEqual(root_excluded.skip_reason, ROOT_CAP_SKIP)
        self.assertEqual(total_excluded.skip_reason, TOTAL_CAP_SKIP)
        self.assertTrue(used.used)

        summary = summarize(
            (root_excluded, total_excluded, used),
            side_length=SIDE_LENGTH,
            p=P,
            maximum_block_count=2,
            seed=7,
        )
        self.assertEqual(summary.used_environment_count, 1)
        self.assertEqual(summary.skipped_environment_count, 2)
        self.assertEqual(summary.skipped_for_root_block_cap_count, 1)
        self.assertEqual(summary.skipped_for_total_block_cap_count, 1)
        self.assertTrue(summary.maximum_block_count_selection_bias)
        self.assertIn("WARNING", summary.selection_bias_warning)
        self.assertIn("selection-biased", summary.selection_bias_warning)
        self.assertIsNone(
            summary.unconditional_normalized_off_diagonal_persistence
        )

    def test_uncapped_normalized_off_diagonal_formula_is_exact(self) -> None:
        result = diagnose_environment(
            side_length=SIDE_LENGTH,
            ranked_edges=_one_postcritical_edge(),
            p=P,
            maximum_block_count=16,
        )
        summary = summarize(
            (result,),
            side_length=SIDE_LENGTH,
            p=P,
            maximum_block_count=16,
            seed=11,
        )
        self.assertFalse(summary.maximum_block_count_selection_bias)
        self.assertEqual(summary.used_environment_count, 1)
        expected = (2.0 / 256.0) / (1.0 - 16.0 / 256.0)
        self.assertAlmostEqual(
            summary.used_sample_normalized_off_diagonal_persistence,
            expected,
        )
        self.assertAlmostEqual(
            summary.unconditional_normalized_off_diagonal_persistence,
            expected,
        )
        self.assertFalse(summary.weak_recovery_claimed)

    def test_small_seeded_run_is_reproducible(self) -> None:
        arguments = dict(
            side_length=4,
            repetitions=2,
            p=P,
            maximum_block_count=16,
            seed=20260726,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first[0]), asdict(second[0]))
        self.assertEqual(
            tuple(asdict(item) for item in first[1]),
            tuple(asdict(item) for item in second[1]),
        )
        self.assertEqual(first[0].used_environment_count, 2)
        self.assertEqual(first[0].skipped_environment_count, 0)
        self.assertTrue(
            first[0].exact_orientation_enumeration_on_used_environments
        )

    def test_invalid_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=3, repetitions=1)
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=4, repetitions=0)
        with self.assertRaises(ValueError):
            run_diagnostic(side_length=4, repetitions=1, maximum_block_count=0)


if __name__ == "__main__":
    unittest.main()
