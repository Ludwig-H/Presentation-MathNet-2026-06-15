"""Tests for exact critical-cut quotient factor elimination."""

from __future__ import annotations

import random
import unittest
from math import exp, fsum, log

from critical_cut_collective_gibbs_diagnostic import (
    analyze_final_root,
    critical_cut_partition,
)
from critical_cut_quotient_elimination import (
    BinaryFactor,
    analyze_root_by_elimination,
    eliminate_log_partition,
    run_elimination_diagnostic,
)
from critical_pair_path_geometry import (
    sample_ranked_edges,
    triangular_torus_edges,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


P = 0.81


class CriticalCutQuotientEliminationTests(unittest.TestCase):
    def test_generic_elimination_matches_direct_constrained_sum(self) -> None:
        factors = (
            BinaryFactor((0, 1), (1.0, 2.0, 3.0, 5.0)),
            BinaryFactor((1, 2), (7.0, 11.0, 13.0, 17.0)),
        )
        for fixed in ({}, {0: 1}, {0: 0, 2: 1}):
            audit = eliminate_log_partition(factors, 3, fixed)
            direct = fsum(
                factors[0].values[
                    ((state >> 0) & 1)
                    | (((state >> 1) & 1) << 1)
                ]
                * factors[1].values[
                    ((state >> 1) & 1)
                    | (((state >> 2) & 1) << 1)
                ]
                for state in range(1 << 3)
                if all(
                    ((state >> variable) & 1) == bit
                    for variable, bit in fixed.items()
                )
            )
            self.assertAlmostEqual(exp(audit.log_partition), direct)
        self.assertEqual(
            eliminate_log_partition(factors, 3).induced_width,
            1,
        )

    def test_free_variable_and_impossible_constraint_are_handled(self) -> None:
        free = eliminate_log_partition((), 3, {0: 1})
        self.assertAlmostEqual(free.log_partition, log(4.0))
        impossible = eliminate_log_partition(
            (BinaryFactor((0,), (1.0, 0.0)),),
            1,
            {0: 1},
        )
        self.assertEqual(impossible.log_partition, float("-inf"))

    def test_one_critical_block_has_no_nonconstant_factor(self) -> None:
        side_length = 4
        ranked_edges = tuple(
            (0.01, *edge)
            for edge in triangular_torus_edges(side_length)
        )
        environment = HierarchicalSweepEnvironment(
            side_length,
            ranked_edges,
            P,
        )
        partition = critical_cut_partition(environment, ranked_edges)
        self.assertEqual(len(partition.final_components), 1)
        result = analyze_root_by_elimination(environment, partition, 0)
        self.assertEqual(result.block_count, 1)
        self.assertEqual(result.nonconstant_factor_count, 0)
        self.assertEqual(result.correlation_matrix, ((1.0,),))
        self.assertAlmostEqual(
            result.scaled_unconstrained_log_partition,
            log(2.0),
        )

    def test_factor_elimination_matches_exhaustive_orientation_law(self) -> None:
        for side_length, seed in ((4, 19), (5, 23), (6, 29)):
            ranked_edges = sample_ranked_edges(
                side_length,
                random.Random(seed),
            )
            environment = HierarchicalSweepEnvironment(
                side_length,
                ranked_edges,
                P,
            )
            partition = critical_cut_partition(environment, ranked_edges)
            for component_index in range(len(partition.final_components)):
                exact = analyze_final_root(
                    environment,
                    partition,
                    component_index,
                )
                eliminated = analyze_root_by_elimination(
                    environment,
                    partition,
                    component_index,
                )
                self.assertEqual(
                    eliminated.block_indices,
                    exact.block_indices,
                )
                self.assertEqual(
                    eliminated.orientation_state_count,
                    exact.orientation_state_count,
                )
                for first_row, second_row in zip(
                    eliminated.correlation_matrix,
                    exact.correlation_matrix,
                    strict=True,
                ):
                    for first, second in zip(
                        first_row,
                        second_row,
                        strict=True,
                    ):
                        self.assertAlmostEqual(first, second, places=11)
                self.assertLess(
                    eliminated.maximum_global_flip_log_partition_error,
                    2e-12,
                )
                self.assertTrue(eliminated.exact_factor_elimination)
                self.assertTrue(
                    eliminated.conditional_on_fixed_internal_representative
                )
                self.assertTrue(
                    eliminated.factor_log_normalization_constants_omitted
                )
                self.assertIn(
                    "divided by its maximum",
                    eliminated.scaled_log_partition_definition,
                )

    def test_width_is_strictly_below_enumeration_dimension_on_seeded_l6(self) -> None:
        side_length = 6
        ranked_edges = sample_ranked_edges(
            side_length,
            random.Random(20260726),
        )
        environment = HierarchicalSweepEnvironment(
            side_length,
            ranked_edges,
            P,
        )
        partition = critical_cut_partition(environment, ranked_edges)
        largest = max(
            range(len(partition.final_components)),
            key=lambda index: len(
                partition.final_component_block_indices[index]
            ),
        )
        result = analyze_root_by_elimination(
            environment,
            partition,
            largest,
        )
        self.assertLessEqual(
            result.maximum_initial_factor_scope_size,
            result.block_count,
        )
        self.assertLess(
            result.unconstrained_induced_width,
            result.block_count,
        )
        self.assertEqual(
            result.constrained_partition_evaluation_count,
            4 * result.block_count * (result.block_count - 1) // 2,
        )

    def test_uncapped_summary_is_reproducible_and_reports_no_selection(self) -> None:
        arguments = dict(
            side_length=5,
            repetitions=3,
            p=P,
            seed=41,
        )
        first = run_elimination_diagnostic(**arguments)
        second = run_elimination_diagnostic(**arguments)
        self.assertEqual(first, second)
        summary, environments = first
        self.assertEqual(len(environments), 3)
        self.assertEqual(summary.cap_exclusion_count, 0)
        self.assertFalse(summary.selection_bias_from_complexity_cap)
        self.assertTrue(summary.every_environment_eliminated_exactly)
        self.assertTrue(summary.factor_log_normalization_constants_omitted)
        self.assertFalse(summary.weak_recovery_claimed)
        self.assertGreaterEqual(
            summary.maximum_root_block_count,
            summary.maximum_unconstrained_induced_width,
        )

    def test_invalid_factor_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            BinaryFactor((1, 0), (1.0, 1.0, 1.0, 1.0))
        with self.assertRaises(ValueError):
            BinaryFactor((0,), (1.0,))
        with self.assertRaises(ValueError):
            BinaryFactor((0,), (1.0, -1.0))
        with self.assertRaises(ValueError):
            eliminate_log_partition((), 1, {1: 0})
        with self.assertRaises(ValueError):
            eliminate_log_partition((), 1, {0: 2})
        with self.assertRaises(ValueError):
            run_elimination_diagnostic(side_length=3, repetitions=1)
        with self.assertRaises(ValueError):
            run_elimination_diagnostic(side_length=4, repetitions=0)
        with self.assertRaises(ValueError):
            run_elimination_diagnostic(side_length=4, repetitions=1, p=1.0)


if __name__ == "__main__":
    unittest.main()
