from __future__ import annotations

import itertools
import unittest
from math import comb, exp, log, tanh

from critical_band_thresholds import Q_CRITICAL, open_probability
from critical_component_boundary import (
    P_BOUNDARY_LATE,
    P_INFO,
    P_SW,
    boundary_late_true_fraction,
    bucket_vote_moments,
    clock_factor_weights,
    closed_categories,
    connected_components,
    critical_closed_categories,
    critical_masses,
    critical_time_untruncated,
    direct_cut_merger_hazard,
    future_activation_crossover_time,
    future_true_minus_false_mass,
    grouped_rates,
    internal_proportions,
    margin_failure_bound,
    margin_success_probability,
    multiply_four_weights,
    pair_palm_weights,
    parity_keep_probability,
    parity_log_odds,
    partition_boundary,
    strict_majority_probability,
    snapshot_cut_information_load,
    snapshot_cut_reliability,
    snapshot_cut_signal_to_noise,
    snapshot_cut_vote_moments,
    unactivated_true_minus_false_mass,
    walsh_coefficients,
)


class CriticalComponentBoundaryTests(unittest.TestCase):
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

    def test_two_residual_balance_crossovers_are_distinct(self) -> None:
        p = 0.8
        beta_active = future_activation_crossover_time(p)
        self.assertAlmostEqual(beta_active, 0.5)
        self.assertAlmostEqual(
            future_true_minus_false_mass(p, beta_active), 0.0
        )
        for time in (0.0, 0.2, 0.5, 0.9):
            self.assertGreater(unactivated_true_minus_false_mass(p, time), 0.0)
        self.assertAlmostEqual(unactivated_true_minus_false_mass(p, 1.0), 0.0)

    def test_boundary_late_threshold_aligns_the_two_times(self) -> None:
        self.assertAlmostEqual(
            critical_time_untruncated(P_BOUNDARY_LATE),
            future_activation_crossover_time(P_BOUNDARY_LATE),
        )

    def test_internal_proportions_separate_geometry_from_residual_marks(
        self,
    ) -> None:
        for p, time in ((0.7, 0.2), (0.8, 0.6), (0.9, 1.0)):
            rate = log(p / (1.0 - p))
            open_density = p * (1.0 - exp(-rate * time))
            baseline = internal_proportions(p, time, open_density)
            self.assertAlmostEqual(baseline.open_true, open_density)
            self.assertAlmostEqual(
                baseline.late_true,
                p * (exp(-rate * time) - exp(-rate)),
            )
            self.assertAlmostEqual(baseline.censored_true, 1.0 - p)
            self.assertAlmostEqual(baseline.false, 1.0 - p)
            self.assertAlmostEqual(
                baseline.open_true
                + baseline.late_true
                + baseline.censored_true
                + baseline.false,
                1.0,
            )

            more_connected = internal_proportions(
                p, time, min(1.0, open_density + 0.1)
            )
            self.assertLessEqual(more_connected.false, baseline.false)
            self.assertLessEqual(
                more_connected.unactivated_true,
                baseline.unactivated_true,
            )

    def test_critical_unconditional_masses_and_closed_forms(self) -> None:
        for p in (P_SW, P_BOUNDARY_LATE, P_INFO, 0.835805792367, 0.95):
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
        self.assertAlmostEqual(
            boundary_late_true_fraction(P_BOUNDARY_LATE), 0.5
        )
        self.assertAlmostEqual(critical_closed_categories(P_SW).true_probability, 0.5)
        self.assertLess(P_SW, P_BOUNDARY_LATE)
        self.assertLess(P_BOUNDARY_LATE, P_INFO)
        self.assertGreater(boundary_late_true_fraction(P_INFO), 0.5)

    def test_internal_edges_cancel_and_only_the_current_cut_is_retained(self) -> None:
        left_child = frozenset({0, 1})
        right_child = frozenset({2, 3})
        partition = (left_child, right_child)
        marked_edges = (
            (0, 1, -1),  # false and internal to the left child
            (2, 3, 1),   # true and internal to the right child
            (0, 2, -1),  # false and crossing the current cut
            (1, 3, 1),   # true and crossing the current cut
        )
        boundary = partition_boundary(
            ((left, right) for left, right, _ in marked_edges), partition
        )
        self.assertEqual(set(boundary), {(0, 2), (1, 3)})

        satisfaction: dict[tuple[int, int], tuple[bool, ...]] = {}
        for a, b in itertools.product((0, 1), repeat=2):
            values: list[bool] = []
            for left, right, truth_sign in marked_edges:
                left_spin = -1 if (left in left_child and a) or (
                    left in right_child and b
                ) else 1
                right_spin = -1 if (right in left_child and a) or (
                    right in right_child and b
                ) else 1
                values.append(truth_sign * left_spin * right_spin == 1)
            satisfaction[(a, b)] = tuple(values)

        for edge_index in (0, 1):
            self.assertEqual(
                {values[edge_index] for values in satisfaction.values()},
                {satisfaction[(0, 0)][edge_index]},
            )
        for edge_index in (2, 3):
            self.assertEqual(
                {values[edge_index] for values in satisfaction.values()},
                {False, True},
            )

    def test_boundary_marks_factor_after_conditioning_on_triangle_partition(
        self,
    ) -> None:
        vertices = (0, 1, 2)
        edges = ((0, 1), (0, 2), (1, 2))
        target = (frozenset({0, 1}), frozenset({2}))
        self.assertEqual(partition_boundary(edges, target), edges[1:])

        p = 0.82
        time = 0.47
        q = 1.0 - p
        rate = log(p / q)
        category_probabilities = (
            p * (1.0 - exp(-rate * time)),  # early true
            p * (exp(-rate * time) - exp(-rate)),  # late true
            p * exp(-rate),  # censored true
            q,  # false
        )
        total = 0.0
        marginal = [[0.0] * 4 for _ in range(2)]
        joint = [[0.0] * 4 for _ in range(4)]
        internal = [0.0] * 4
        for categories in itertools.product(range(4), repeat=len(edges)):
            probability = 1.0
            early_edges = []
            for edge, category in zip(edges, categories):
                probability *= category_probabilities[category]
                if category == 0:
                    early_edges.append(edge)
            if connected_components(vertices, early_edges) != target:
                continue
            total += probability
            internal[categories[0]] += probability
            marginal[0][categories[1]] += probability
            marginal[1][categories[2]] += probability
            joint[categories[1]][categories[2]] += probability

        self.assertGreater(total, 0.0)
        residual = closed_categories(p, time)
        expected = (0.0, residual.late_true, residual.censored_true, residual.false)
        for edge_marginal in marginal:
            for observed, target_probability in zip(edge_marginal, expected):
                self.assertAlmostEqual(observed / total, target_probability)
        for first_category in range(4):
            for second_category in range(4):
                self.assertAlmostEqual(
                    joint[first_category][second_category] / total,
                    expected[first_category] * expected[second_category],
                )
        self.assertAlmostEqual(internal[0] / total, 1.0)
        self.assertAlmostEqual(sum(internal[1:]) / total, 0.0)

    def test_pair_palm_selection_is_weighted_by_ordered_pairs(self) -> None:
        large = frozenset({0, 1, 2})
        medium = frozenset({4, 5})
        singleton = frozenset({9})
        weights = pair_palm_weights(
            (large, medium, singleton),
            distance=lambda left, right: abs(left - right),
            minimum_distance=0.0,
        )
        self.assertAlmostEqual(weights[large], 6.0 / 8.0)
        self.assertAlmostEqual(weights[medium], 2.0 / 8.0)
        self.assertEqual(weights[singleton], 0.0)

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

    def test_snapshot_cut_moments_match_full_enumeration(self) -> None:
        size = 7
        p = 0.81
        time = 0.55
        s = closed_categories(p, time).true_probability
        probabilities: list[tuple[float, float]] = []
        for word in itertools.product((0, 1), repeat=size):
            probability = 1.0
            for bit in word:
                probability *= s if bit else 1.0 - s
            vote = 2.0 * sum(word) - size
            probabilities.append((probability, vote))
        direct_mean = sum(probability * vote for probability, vote in probabilities)
        direct_variance = sum(
            probability * (vote - direct_mean) ** 2
            for probability, vote in probabilities
        )
        mean, variance = snapshot_cut_vote_moments(size, p, time)
        self.assertAlmostEqual(mean, direct_mean)
        self.assertAlmostEqual(variance, direct_variance)
        self.assertAlmostEqual(
            snapshot_cut_signal_to_noise(size, p, time),
            direct_mean * direct_mean / direct_variance,
        )

    def test_snapshot_reliability_matches_two_independent_calculations(self) -> None:
        size = 6
        p = 0.82
        time = 0.63
        categories = closed_categories(p, time)
        s = categories.true_probability
        residual_log_odds = log(s / (1.0 - s))
        plus: dict[int, float] = {}
        minus: dict[int, float] = {}
        direct_from_plus = 0.0
        for count in range(size + 1):
            plus[count] = (
                comb(size, count) * s**count * (1.0 - s) ** (size - count)
            )
            minus[count] = (
                comb(size, count)
                * (1.0 - s) ** count
                * s ** (size - count)
            )
            bias = tanh(residual_log_odds * (2 * count - size) / 2.0)
            direct_from_plus += plus[count] * bias * bias
        symmetric = 0.5 * sum(
            (plus[count] - minus[count]) ** 2
            / (plus[count] + minus[count])
            for count in range(size + 1)
        )
        result = snapshot_cut_reliability(size, p, time)
        self.assertAlmostEqual(result, direct_from_plus)
        self.assertAlmostEqual(result, symmetric)

    def test_snapshot_terminal_cut_has_zero_information(self) -> None:
        for size in (1, 2, 7, 20):
            self.assertAlmostEqual(snapshot_cut_reliability(size, 0.8, 1.0), 0.0)
            self.assertAlmostEqual(
                snapshot_cut_signal_to_noise(size, 0.8, 1.0), 0.0
            )
            self.assertAlmostEqual(snapshot_cut_information_load(size, 0.8, 1.0), 0.0)

    def test_information_load_has_the_expected_weak_bias_equivalent(self) -> None:
        p = 0.8
        size = 37
        for time in (0.9, 0.99, 0.999):
            h = closed_categories(p, time).signed_margin
            load = snapshot_cut_information_load(size, p, time)
            ratio = load / (0.5 * size * h * h)
            self.assertLess(abs(ratio - 1.0), 0.01)

    def test_merger_hazard_is_exactly_size_biased(self) -> None:
        p = 0.8
        time = 0.61
        unit = direct_cut_merger_hazard(1, p, time)
        for size in (2, 3, 11):
            self.assertAlmostEqual(
                direct_cut_merger_hazard(size, p, time), size * unit
            )

    def test_merger_compensator_matches_direct_cut_activation(self) -> None:
        p = 0.8
        q = 1.0 - p
        rate = log(p / q)
        for size in (1, 2, 5, 9):
            # Independent integration of
            # P(cut still closed at t) * conditional cut hazard at t.
            compensator = 0.0
            for power in range(size):
                compensator += (
                    size
                    * comb(size - 1, power)
                    * q ** (size - 1 - power)
                    * p ** (power + 1)
                    * (1.0 - exp(-(power + 1) * rate))
                    / (power + 1)
                )
            direct = 1.0 - (2.0 * q) ** size
            self.assertAlmostEqual(compensator, direct)

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
