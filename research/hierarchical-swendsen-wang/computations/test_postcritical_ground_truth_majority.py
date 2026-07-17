from __future__ import annotations

import itertools
import unittest
from math import comb, exp, log, tanh

from critical_band_thresholds import Q_CRITICAL, open_probability
from postcritical_ground_truth_majority import (
    P_INFO,
    P_LATE,
    P_SW,
    bucket_vote_moments,
    clock_factor_weights,
    closed_categories,
    critical_closed_categories,
    critical_masses,
    critical_time_untruncated,
    grouped_rates,
    late_pool_true_fraction,
    margin_failure_bound,
    margin_success_probability,
    multiply_four_weights,
    parity_keep_probability,
    parity_log_odds,
    strict_majority_probability,
    walsh_coefficients,
)


class PostcriticalGroundTruthMajorityTests(unittest.TestCase):
    def test_untruncated_critical_time_and_censoring_domain(self) -> None:
        for p in (0.51, 0.6, P_SW, P_INFO, 0.9):
            beta = critical_time_untruncated(p)
            self.assertAlmostEqual(open_probability(p, beta), Q_CRITICAL)
            self.assertEqual(beta <= 1.0 + 1e-14, p >= P_SW - 1e-14)

    def test_closed_category_decomposition_matches_direct_masses(self) -> None:
        for p, time in ((0.7, 0.2), (0.8, 0.6), (0.9, 1.0)):
            result = closed_categories(p, time)
            q = 1.0 - p
            rate = log(p / q)
            denominator = q + p * exp(-rate * time)
            direct = (
                p * (exp(-rate * time) - exp(-rate)) / denominator,
                p * exp(-rate) / denominator,
                q / denominator,
            )
            self.assertAlmostEqual(result.late_true, direct[0])
            self.assertAlmostEqual(result.censored_true, direct[1])
            self.assertAlmostEqual(result.false, direct[2])
            self.assertAlmostEqual(
                result.late_true + result.censored_true + result.false, 1.0
            )
            self.assertAlmostEqual(result.censored_true, result.false)
            self.assertAlmostEqual(result.signed_margin, result.late_true)
            self.assertAlmostEqual(
                result.signed_margin, tanh(rate * (1.0 - time) / 2.0)
            )

    def test_critical_unconditional_masses_and_closed_forms(self) -> None:
        for p in (P_SW, P_LATE, P_INFO, 0.835805792367, 0.95):
            masses = critical_masses(p)
            categories = critical_closed_categories(p)
            self.assertAlmostEqual(
                masses.early_true
                + masses.late_true
                + masses.censored_true
                + masses.false,
                1.0,
            )
            self.assertAlmostEqual(masses.censored_true, masses.false)
            self.assertAlmostEqual(
                categories.true_probability, (p - Q_CRITICAL) / (1.0 - Q_CRITICAL)
            )
            self.assertAlmostEqual(
                categories.signed_margin,
                (2.0 * p - 1.0 - Q_CRITICAL) / (1.0 - Q_CRITICAL),
            )

    def test_two_majority_thresholds_are_distinct(self) -> None:
        self.assertAlmostEqual(late_pool_true_fraction(P_LATE), 0.5)
        self.assertAlmostEqual(critical_closed_categories(P_SW).true_probability, 0.5)
        self.assertLess(P_SW, P_LATE)
        self.assertLess(P_LATE, P_INFO)
        self.assertGreater(late_pool_true_fraction(P_INFO), 0.5)

    def test_critical_time_maximizes_postcritical_quality(self) -> None:
        for p in (0.7, P_INFO, 0.835805792367, 0.95):
            beta = critical_time_untruncated(p)
            values = [
                closed_categories(p, beta + index * (1.0 - beta) / 20.0)
                for index in range(21)
            ]
            for left, right in zip(values, values[1:]):
                self.assertGreaterEqual(
                    left.true_probability + 1e-15, right.true_probability
                )
                self.assertGreaterEqual(
                    left.signed_margin + 1e-15, right.signed_margin
                )

    def test_majority_probability_matches_full_enumeration(self) -> None:
        size = 6
        p = 0.82
        time = critical_time_untruncated(p)
        s = closed_categories(p, time).true_probability
        direct = 0.0
        for word in itertools.product((0, 1), repeat=size - 1):
            probability = 1.0
            for bit in word:
                probability *= s if bit else 1.0 - s
            true_count = 1 + sum(word)
            if 2 * true_count > size:
                direct += probability
        self.assertAlmostEqual(
            strict_majority_probability(size, p, time), direct
        )

    def test_bucket_vote_moments_match_full_enumeration(self) -> None:
        size = 7
        p = 0.81
        time = 0.55
        s = closed_categories(p, time).true_probability
        probabilities: list[tuple[float, float]] = []
        for word in itertools.product((0, 1), repeat=size - 1):
            probability = 1.0
            for bit in word:
                probability *= s if bit else 1.0 - s
            vote = 2.0 * (1 + sum(word)) - size
            probabilities.append((probability, vote))
        direct_mean = sum(probability * vote for probability, vote in probabilities)
        direct_variance = sum(
            probability * (vote - direct_mean) ** 2
            for probability, vote in probabilities
        )
        mean, variance = bucket_vote_moments(size, p, time)
        self.assertAlmostEqual(mean, direct_mean)
        self.assertAlmostEqual(variance, direct_variance)

    def test_hoeffding_margin_bound_covers_exact_binomial_tail(self) -> None:
        for residual_count in (1, 2, 5, 12):
            for winner_indicator in (0, 1):
                for h in (0.0, 0.15, 0.5, 0.9):
                    for require_strict in (False, True):
                        s = (1.0 + h) / 2.0
                        exact = 0.0
                        for successes in range(residual_count + 1):
                            margin = (
                                winner_indicator + 2 * successes - residual_count
                            )
                            fails = margin <= 0 if require_strict else margin < 0
                            if fails:
                                exact += (
                                    comb(residual_count, successes)
                                    * s**successes
                                    * (1.0 - s)
                                    ** (residual_count - successes)
                                )
                        self.assertLessEqual(
                            exact,
                            margin_failure_bound(
                                residual_count,
                                winner_indicator,
                                h,
                                require_strict,
                            )
                            + 1e-14,
                        )

    def test_exact_margin_tail_matches_full_enumeration(self) -> None:
        for residual_count in range(7):
            for winner_indicator in (0, 1):
                for s in (0.2, 0.5, 0.83):
                    for require_strict in (False, True):
                        direct = 0.0
                        for word in itertools.product(
                            (0, 1), repeat=residual_count
                        ):
                            probability = 1.0
                            for bit in word:
                                probability *= s if bit else 1.0 - s
                            margin = (
                                winner_indicator
                                + 2 * sum(word)
                                - residual_count
                            )
                            succeeds = (
                                margin > 0 if require_strict else margin >= 0
                            )
                            if succeeds:
                                direct += probability
                        self.assertAlmostEqual(
                            margin_success_probability(
                                residual_count,
                                winner_indicator,
                                s,
                                require_strict,
                            ),
                            direct,
                        )

    def test_grouped_rates_match_direct_flip_enumeration(self) -> None:
        groups = ((1, 0), (1, 1), (2, 0), (0, 1))
        sizes = (
            sum(group == 0 for group, _ in groups),
            sum(group == 1 for group, _ in groups),
            sum(group == 2 for group, _ in groups),
        )
        counts = (
            sum(group == 0 and truth == 1 for group, truth in groups),
            sum(group == 1 and truth == 1 for group, truth in groups),
            sum(group == 2 and truth == 1 for group, truth in groups),
        )
        weight = 1.7
        rates = grouped_rates(sizes, counts, weight)
        for a, b in itertools.product((0, 1), repeat=2):
            direct = 0.0
            for group, truth in groups:
                flipped = (group == 1 and a == 1) or (group == 2 and b == 1)
                satisfied = truth if not flipped else 1 - truth
                direct += weight * satisfied
            self.assertAlmostEqual(rates[(a, b)], direct)

        m1 = 2 * counts[1] - sizes[1]
        m2 = 2 * counts[2] - sizes[2]
        self.assertAlmostEqual(rates[(0, 0)] - rates[(1, 0)], weight * m1)
        self.assertAlmostEqual(rates[(0, 0)] - rates[(0, 1)], weight * m2)

    def test_exact_four_weight_parity_criterion(self) -> None:
        weights = {(0, 0): 8.0, (1, 0): 2.0, (0, 1): 3.0, (1, 1): 1.0}
        probability = parity_keep_probability(weights)
        log_odds = parity_log_odds(weights)
        self.assertAlmostEqual(probability, 9.0 / 14.0)
        self.assertAlmostEqual(log_odds, log(9.0 / 5.0))
        self.assertEqual(probability > 0.5, log_odds > 0.0)

        counterexample = {
            (0, 0): 1.01,
            (1, 0): 1.0,
            (0, 1): 1.0,
            (1, 1): 0.01,
        }
        self.assertGreater(counterexample[(0, 0)], counterexample[(1, 0)])
        self.assertGreater(counterexample[(0, 0)], counterexample[(0, 1)])
        self.assertLess(parity_keep_probability(counterexample), 0.5)

    def test_group_majorities_generate_the_positive_walsh_cone(self) -> None:
        factors: list[dict[tuple[int, int], float]] = []
        for beta in (0.0, 0.4, 1.0):
            for m0 in range(3):
                for k0 in range(m0 + 1):
                    for m1 in range(4):
                        for k1 in range(m1 + 1):
                            if 2 * k1 < m1:
                                continue
                            for m2 in range(4):
                                for k2 in range(m2 + 1):
                                    if 2 * k2 < m2:
                                        continue
                                    rates = grouped_rates(
                                        (m0, m1, m2), (k0, k1, k2), 1.3
                                    )
                                    if max(rates.values()) == 0.0:
                                        continue
                                    factor = clock_factor_weights(rates, beta)
                                    coefficients = walsh_coefficients(factor)
                                    for coefficient in coefficients.values():
                                        self.assertGreaterEqual(
                                            coefficient + 1e-12, 0.0
                                        )
                                    factors.append(factor)

        for left, right in zip(factors[::19], factors[7::19]):
            product = multiply_four_weights(left, right)
            coefficients = walsh_coefficients(product)
            for coefficient in coefficients.values():
                self.assertGreaterEqual(coefficient + 1e-10, 0.0)

        for total, margin_one, margin_two, beta in (
            (7.3, 1.1, 2.4, 0.2),
            (4.8, 0.0, 1.7, 0.65),
            (9.0, 3.2, 0.4, 1.0),
        ):
            weighted_rectangle = {
                (0, 0): total,
                (1, 0): total - margin_one,
                (0, 1): total - margin_two,
                (1, 1): total - margin_one - margin_two,
            }
            factor = clock_factor_weights(weighted_rectangle, beta)
            for coefficient in walsh_coefficients(factor).values():
                self.assertGreaterEqual(coefficient + 1e-12, 0.0)

    def test_local_and_grouped_majorities_certify_correct_parity(self) -> None:
        local_rates = {
            (0, 0): 5.0,
            (1, 1): 5.0,
            (1, 0): 3.0,
            (0, 1): 3.0,
        }
        local = clock_factor_weights(local_rates, 0.35)
        first = clock_factor_weights(
            grouped_rates((2, 5, 4), (1, 4, 3), 0.7), 0.55
        )
        second = clock_factor_weights(
            grouped_rates((3, 4, 7), (2, 2, 4), 0.7), 0.85
        )
        weights = multiply_four_weights(local, first, second)
        self.assertGreater(parity_keep_probability(weights), 0.5)
        self.assertGreater(walsh_coefficients(weights)[(1, 2)], 0.0)


if __name__ == "__main__":
    unittest.main()
