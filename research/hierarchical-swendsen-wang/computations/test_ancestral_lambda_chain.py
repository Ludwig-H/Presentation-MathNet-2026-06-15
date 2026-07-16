from __future__ import annotations

import unittest
from math import log

from ancestral_lambda_chain import (
    AncestorBucket,
    AncestorGeometry,
    ancestral_log_weights,
    ancestral_message,
    closed_satisfaction_probability,
    conditional_total_count_pmf,
    exact_reliability_given_skeleton,
    grouped_count_pmf,
    local_result,
    message_from_coefficients,
    walsh_coefficients,
)


class AncestorLambdaTests(unittest.TestCase):
    def test_four_rates_match_direct_complementation(self) -> None:
        bucket = AncestorBucket(
            beta=0.6,
            total=(11.0, 7.0, 5.0),
            satisfied=(4.0, 2.0, 3.0),
        )
        expected = {
            (0, 0): 9.0,
            (0, 1): 8.0,
            (1, 0): 12.0,
            (1, 1): 11.0,
        }
        self.assertEqual(bucket.four_lambdas(), expected)

    def test_message_matches_direct_parity_sums(self) -> None:
        ancestors = (
            AncestorBucket(0.55, (2.0, 3.0, 4.0), (1.0, 2.0, 1.0)),
            AncestorBucket(0.80, (5.0, 2.0, 1.0), (2.0, 1.0, 1.0)),
        )
        direct = local_result(3.0, 7.0, 0.4, ancestors).log_odds
        separated = ancestral_message(ancestors)
        self.assertAlmostEqual(
            direct,
            separated + log(3.0 / 4.0) + 0.6 * (3.0 - 4.0),
        )

    def test_grouped_law_is_normalized_and_has_total_binomial_law(self) -> None:
        group_sizes = (3, 2, 1)
        p = 0.79
        beta = 0.43
        grouped = grouped_count_pmf(group_sizes, p, beta)
        total = conditional_total_count_pmf(sum(group_sizes), p, beta)
        self.assertAlmostEqual(sum(grouped.values()), 1.0)
        collapsed: dict[int, float] = {}
        for counts, mass in grouped.items():
            collapsed[sum(counts)] = collapsed.get(sum(counts), 0.0) + mass
        self.assertEqual(set(collapsed), set(total))
        for count in total:
            self.assertAlmostEqual(collapsed[count], total[count])

    def test_walsh_collapse_matches_four_state_message(self) -> None:
        ancestors = (
            AncestorBucket(0.55, (2.0, 3.0, 4.0), (1.0, 2.0, 1.0)),
            AncestorBucket(0.80, (5.0, 2.0, 1.0), (2.0, 1.0, 1.0)),
        )
        coefficients = walsh_coefficients(ancestral_log_weights(ancestors))
        self.assertAlmostEqual(
            message_from_coefficients(coefficients),
            ancestral_message(ancestors),
        )

    def test_group_means_include_uniform_winner_correction(self) -> None:
        group_sizes = (4, 3, 2)
        p = 0.83
        beta = 0.37
        grouped = grouped_count_pmf(group_sizes, p, beta)
        s = closed_satisfaction_probability(p, beta)
        total_size = sum(group_sizes)
        for group, size in enumerate(group_sizes):
            empirical = sum(counts[group] * mass for counts, mass in grouped.items())
            expected = size * s + size / total_size * (1.0 - s)
            self.assertAlmostEqual(empirical, expected)

    def test_terminal_identity(self) -> None:
        for size in range(1, 8):
            value = exact_reliability_given_skeleton(size, 1.0, (), 0.8)
            self.assertAlmostEqual(value, 1.0 / size)

    def test_small_ancestor_enumeration_stays_in_unit_interval(self) -> None:
        value = exact_reliability_given_skeleton(
            local_size=3,
            local_beta=0.45,
            ancestor_geometry=(
                AncestorGeometry(0.65, (2, 1, 1)),
                AncestorGeometry(0.85, (1, 1, 2)),
            ),
            p=0.8,
        )
        self.assertGreaterEqual(value, 0.0)
        self.assertLessEqual(value, 1.0)

    def test_local_zero_rate_is_handled_without_nan(self) -> None:
        result = local_result(1.0, 1.0, 0.5)
        self.assertEqual(result.persistence, 1.0)
        self.assertEqual(result.reliability, 1.0)


if __name__ == "__main__":
    unittest.main()
