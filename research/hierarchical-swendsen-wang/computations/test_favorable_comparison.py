from __future__ import annotations

import unittest
from math import isclose, log

from critical_pair_path_geometry import (
    CriticalKruskalForest,
    sample_ranked_edges,
    triangular_torus_edges,
)
from favorable_time_comparison import (
    anti_alignment_violation,
    bucket_binary_experiment,
    bucket_blackwell_minimum_call_gap,
    critical_beta,
    local_comparison_error_bound,
    local_parity_log_odds,
    local_persistence,
    p_eight_counterexample,
    screened_two_edge_contraction,
    two_edge_blackwell_erasure_probability,
)
from pair_favorability_diagnostic import (
    PAIR_CLASSES,
    _jackknife_ratio_difference_standard_error,
    _weighted_ratio,
    classify_pair,
    far_pair_classes,
)


class FavorableTimeComparisonTests(unittest.TestCase):
    def test_heat_bath_sweep_has_the_palindromic_l2_identity(self) -> None:
        stationary = (0.1, 0.2, 0.3, 0.4)

        def projection(bit: int) -> tuple[tuple[float, ...], ...]:
            rows = []
            for state in range(4):
                first = state & ~(1 << bit)
                second = first | (1 << bit)
                normalizer = stationary[first] + stationary[second]
                row = [0.0] * 4
                row[first] = stationary[first] / normalizer
                row[second] = stationary[second] / normalizer
                rows.append(tuple(row))
            return tuple(rows)

        def multiply(
            first: tuple[tuple[float, ...], ...],
            second: tuple[tuple[float, ...], ...],
        ) -> tuple[tuple[float, ...], ...]:
            return tuple(
                tuple(
                    sum(first[i][k] * second[k][j] for k in range(4))
                    for j in range(4)
                )
                for i in range(4)
            )

        def apply(
            matrix: tuple[tuple[float, ...], ...],
            values: tuple[float, ...],
        ) -> tuple[float, ...]:
            return tuple(
                sum(matrix[i][j] * values[j] for j in range(4))
                for i in range(4)
            )

        first = projection(0)
        second = projection(1)
        for projection_matrix in (first, second):
            square = multiply(projection_matrix, projection_matrix)
            for i in range(4):
                for j in range(4):
                    self.assertAlmostEqual(square[i][j], projection_matrix[i][j])
                    self.assertAlmostEqual(
                        stationary[i] * projection_matrix[i][j],
                        stationary[j] * projection_matrix[j][i],
                    )

        sweep = multiply(first, second)
        reverse_sweep = multiply(second, first)
        palindrome = multiply(reverse_sweep, sweep)
        pair_function = (1.0, -1.0, -1.0, 1.0)
        swept = apply(sweep, pair_function)
        conditional_transfer = tuple(
            pair_function[state] * swept[state] for state in range(4)
        )
        replicated_second_moment = sum(
            stationary[state] * conditional_transfer[state] ** 2
            for state in range(4)
        )
        l2_norm = sum(
            stationary[state] * swept[state] ** 2 for state in range(4)
        )
        palindromic_correlation = sum(
            stationary[state]
            * pair_function[state]
            * apply(palindrome, pair_function)[state]
            for state in range(4)
        )
        self.assertAlmostEqual(replicated_second_moment, l2_norm)
        self.assertAlmostEqual(l2_norm, palindromic_correlation)

    def test_anti_alignment_condition_matches_direct_comparison(self) -> None:
        p = 0.8
        early = critical_beta(p)
        late = 0.8
        for size in range(2, 8):
            for satisfied in range(1, size):
                for outside in (-4.0, -2.0, -0.5, 0.0, 0.5, 2.0, 4.0):
                    violation = anti_alignment_violation(
                        outside, p, size, satisfied, early, late
                    )
                    early_value = local_persistence(
                        outside, p, size, satisfied, early
                    )
                    late_value = local_persistence(
                        outside, p, size, satisfied, late
                    )
                    self.assertEqual(violation, late_value > early_value)
                    self.assertLessEqual(
                        late_value - early_value,
                        local_comparison_error_bound(
                            outside, p, size, satisfied, early, late
                        )
                        + 1e-15,
                    )

    def test_aligned_messages_make_the_critical_node_most_persistent(self) -> None:
        p = 0.8
        early = critical_beta(p)
        late = 0.8
        for size in range(2, 8):
            for satisfied in range(1, size):
                difference = 2 * satisfied - size
                outside = 2.0 if difference >= 0 else -2.0
                self.assertGreaterEqual(
                    local_persistence(outside, p, size, satisfied, early),
                    local_persistence(outside, p, size, satisfied, late),
                )

    def test_p_eight_counterexample_is_admissible_and_strict(self) -> None:
        example = p_eight_counterexample()
        self.assertLess(example.critical_level, example.late_level)
        self.assertLess(example.late_level, 0.81)
        self.assertLess(example.critical_persistence, 1e-6)
        self.assertGreater(example.late_persistence, 0.06)
        self.assertAlmostEqual(
            example.late_persistence, 0.06928167028453398
        )

    def test_screened_two_edge_coefficient_has_the_stated_p_eight_values(
        self,
    ) -> None:
        p = 0.8
        level = critical_beta(p)
        expected = (0.6935822227524088, 0.7590184337431733)
        self.assertAlmostEqual(
            screened_two_edge_contraction(p, level, 0.0), expected[0]
        )
        self.assertAlmostEqual(
            screened_two_edge_contraction(p, level, 1.0), expected[1]
        )
        self.assertLess(expected[0], expected[1])

    def test_bucket_likelihood_matches_the_local_message(self) -> None:
        p = 0.8
        level = critical_beta(p)
        for size in range(2, 10):
            plus, minus = bucket_binary_experiment(p, size, level)
            self.assertAlmostEqual(sum(plus), 1.0)
            self.assertAlmostEqual(sum(minus), 1.0)
            for count in range(1, size):
                likelihood = plus[count] / minus[count]
                self.assertAlmostEqual(
                    log(likelihood),
                    local_parity_log_odds(p, size, count, level),
                )

    def test_critical_bucket_blackwell_dominates_later_buckets(self) -> None:
        p = 0.8
        critical = critical_beta(p)
        for size in range(1, 25):
            for late in (critical, 0.5, 0.7, 0.9, 1.0):
                if late < critical:
                    continue
                self.assertGreaterEqual(
                    bucket_blackwell_minimum_call_gap(
                        p, size, critical, late
                    ),
                    -2e-15,
                )

    def test_two_edge_blackwell_kernel_is_an_explicit_erasure(self) -> None:
        p = 0.8
        early = critical_beta(p)
        late = 1.0
        erasure = two_edge_blackwell_erasure_probability(p, early, late)
        early_plus, early_minus = bucket_binary_experiment(p, 2, early)
        late_plus, late_minus = bucket_binary_experiment(p, 2, late)
        for source, target in (
            (early_plus, late_plus),
            (early_minus, late_minus),
        ):
            degraded = (
                (1.0 - erasure) * source[0],
                source[1] + erasure * (source[0] + source[2]),
                (1.0 - erasure) * source[2],
            )
            for actual, expected in zip(degraded, target, strict=True):
                self.assertAlmostEqual(actual, expected)
        self.assertAlmostEqual(erasure, 0.2791049372404585)


class PairFavorabilityDiagnosticTests(unittest.TestCase):
    def test_pair_classes_cover_early_critical_late_and_separate(self) -> None:
        side = 4
        special = {
            (0, 1): 0.1,
            (4, 5): 0.34,
            (8, 9): 0.36,
        }
        ranked = tuple(
            (special.get(edge, 0.99), *edge)
            for edge in triangular_torus_edges(side)
        )
        forest = CriticalKruskalForest(side, ranked, critical_rank=0.6)
        self.assertEqual(classify_pair(forest, 0, 1, 0.02), "early")
        self.assertEqual(classify_pair(forest, 4, 5, 0.02), "critical")
        self.assertEqual(classify_pair(forest, 8, 9, 0.02), "late")
        self.assertEqual(classify_pair(forest, 12, 13, 0.02), "separate")

    def test_far_pair_classes_form_an_exact_partition(self) -> None:
        import random

        rng = random.Random(31)
        forest = CriticalKruskalForest(
            4, sample_ranked_edges(4, rng), critical_rank=0.6
        )
        classes = far_pair_classes(forest, 0.25, 0.05)
        self.assertEqual(tuple(classes), PAIR_CLASSES)
        pairs = [pair for values in classes.values() for pair in values]
        self.assertEqual(len(pairs), len(set(pairs)))
        self.assertTrue(pairs)

    def test_weighted_ratio_uses_pair_counts(self) -> None:
        value = _weighted_ratio(((1.0, 0.0), (3.0, 1.0)))
        self.assertTrue(isclose(value, 0.75))

    def test_paired_jackknife_keeps_zero_weight_environments(self) -> None:
        first = ((1.0, 0.0), (0.0, 0.0), (3.0, 1.0))
        second = ((2.0, 0.5), (1.0, 0.0), (0.0, 0.0))
        paired = _jackknife_ratio_difference_standard_error(first, second)
        self.assertGreater(paired, 0.0)
        with self.assertRaises(ValueError):
            _jackknife_ratio_difference_standard_error(first, second[:-1])


if __name__ == "__main__":
    unittest.main()
