from __future__ import annotations

import itertools
import unittest
from math import exp, log, prod, tanh

from critical_merger_oracle import (
    critical_mean_same_parity_probability,
    critical_oracle_reliability,
    critical_parameters,
)
from hierarchical_flip_probabilities import (
    even_log_odds_from_walsh,
    factorized_path_same_parity_probability,
    finite_state_path_correlation,
    four_state_probabilities_from_walsh,
    internal_even_probability,
    leaf_heat_bath_flip_probability,
    leaf_metropolis_acceptance_probability,
    node_clock_parameters,
    node_local_log_odds,
    node_mean_error_exponent,
    node_mean_even_probability,
    node_oracle_reliability,
    node_plus_count_pmf,
    node_strict_majority_probability,
    normalize_four_log_weights,
    root_flip_probability,
    symmetric_pair_status_probabilities,
    walsh_parameters,
)


class HierarchicalFlipProbabilityTests(unittest.TestCase):
    def test_uniform_root_recolouring_is_fair(self) -> None:
        self.assertEqual(root_flip_probability(), 0.5)
        self.assertAlmostEqual(root_flip_probability(log(3.0)), 0.75)

    def test_leaf_heat_bath_and_metropolis_have_the_same_target_ratio(
        self,
    ) -> None:
        for log_ratio in (-3.0, -0.4, 0.0, 0.7, 2.0):
            target_ratio = exp(log_ratio)
            heat_forward = leaf_heat_bath_flip_probability(log_ratio)
            heat_reverse = leaf_heat_bath_flip_probability(-log_ratio)
            metro_forward = leaf_metropolis_acceptance_probability(log_ratio)
            metro_reverse = leaf_metropolis_acceptance_probability(-log_ratio)
            self.assertAlmostEqual(heat_forward / heat_reverse, target_ratio)
            self.assertAlmostEqual(
                metro_forward / metro_reverse, target_ratio
            )
            if log_ratio != 0.0:
                self.assertNotAlmostEqual(heat_forward, metro_forward)

    def test_walsh_coordinates_reproduce_arbitrary_four_weights(self) -> None:
        log_weights = (0.3, -1.2, 0.8, 2.1)
        direct = normalize_four_log_weights(log_weights)
        parameters = walsh_parameters(log_weights)
        reconstructed = four_state_probabilities_from_walsh(
            parameters.field_left,
            parameters.field_right,
            parameters.coupling,
        )
        for left, right in zip(
            direct.__dict__.values(), reconstructed.__dict__.values()
        ):
            self.assertAlmostEqual(left, right)
        self.assertAlmostEqual(
            direct.even,
            internal_even_probability(
                parameters.field_left,
                parameters.field_right,
                parameters.coupling,
            ),
        )
        self.assertAlmostEqual(
            log(direct.even / (1.0 - direct.even)),
            even_log_odds_from_walsh(
                parameters.field_left,
                parameters.field_right,
                parameters.coupling,
            ),
        )
        even_split = 1.0 / (
            1.0
            + exp(
                -2.0
                * (parameters.field_left + parameters.field_right)
            )
        )
        odd_split = 1.0 / (
            1.0
            + exp(
                -2.0
                * (parameters.field_left - parameters.field_right)
            )
        )
        self.assertAlmostEqual(direct.p00 / direct.even, even_split)
        self.assertAlmostEqual(
            direct.p01 / (direct.p01 + direct.p10), odd_split
        )

    def test_ancestor_neutral_internal_node_splits_each_parity_equally(
        self,
    ) -> None:
        for local_log_odds in (-2.0, 0.0, 1.4):
            probabilities = four_state_probabilities_from_walsh(
                0.0, 0.0, local_log_odds / 2.0
            )
            expected_even = 1.0 / (1.0 + exp(-local_log_odds))
            self.assertAlmostEqual(probabilities.even, expected_even)
            self.assertAlmostEqual(probabilities.p00, probabilities.p11)
            self.assertAlmostEqual(probabilities.p01, probabilities.p10)

    def test_walsh_reconstruction_on_a_four_weight_grid(self) -> None:
        for log_weights in itertools.product((-1.7, 0.2, 2.1), repeat=4):
            direct = normalize_four_log_weights(log_weights)
            parameters = walsh_parameters(log_weights)
            reconstructed = four_state_probabilities_from_walsh(
                parameters.field_left,
                parameters.field_right,
                parameters.coupling,
            )
            with self.subTest(log_weights=log_weights):
                for left, right in zip(
                    direct.__dict__.values(),
                    reconstructed.__dict__.values(),
                ):
                    self.assertAlmostEqual(left, right)

    def test_clock_parameters_match_the_exponential_conditioning(self) -> None:
        for p, time in ((0.7, 0.0), (0.8, 0.4), (0.9, 1.0)):
            parameters = node_clock_parameters(p, time)
            direct = p * exp(-parameters.coupling * time)
            direct /= 1.0 - p + p * exp(-parameters.coupling * time)
            self.assertAlmostEqual(parameters.satisfaction, direct)
            self.assertAlmostEqual(
                parameters.bias,
                tanh(parameters.residual_log_odds / 2.0),
            )

    def test_node_log_odds_matches_the_mirrored_count_experiment(self) -> None:
        size = 8
        p = 0.79
        time = 0.37
        plus = node_plus_count_pmf(size, p, time)
        for count in range(1, size):
            minus_mass = plus[size - count]
            self.assertAlmostEqual(
                node_local_log_odds(size, count, p, time),
                log(plus[count] / minus_mass),
            )

    def test_mean_even_probability_is_one_plus_reliability_over_two(
        self,
    ) -> None:
        for size in (1, 2, 5, 20):
            for p, time in ((0.7, 0.0), (0.8, 0.4), (0.9, 1.0)):
                reliability = node_oracle_reliability(size, p, time)
                statuses = symmetric_pair_status_probabilities(reliability)
                self.assertAlmostEqual(
                    node_mean_even_probability(size, p, time),
                    (1.0 + reliability) / 2.0,
                )
                self.assertAlmostEqual(statuses.p00, statuses.p11)
                self.assertAlmostEqual(
                    statuses.p00,
                    node_mean_even_probability(size, p, time) / 2.0,
                )

    def test_terminal_level_has_exact_one_over_size_reliability(self) -> None:
        for size in range(1, 20):
            self.assertAlmostEqual(
                node_oracle_reliability(size, 0.81, 1.0), 1.0 / size
            )

    def test_earlier_level_improves_majority_and_large_cut_exponent(self) -> None:
        p = 0.8
        times = (0.1, 0.4, 0.8)
        majorities = [
            node_strict_majority_probability(12, p, time) for time in times
        ]
        exponents = [node_mean_error_exponent(p, time) for time in times]
        self.assertGreater(majorities[0], majorities[1])
        self.assertGreater(majorities[1], majorities[2])
        self.assertGreater(exponents[0], exponents[1])
        self.assertGreater(exponents[1], exponents[2])

    def test_critical_specialization_matches_the_critical_oracle(self) -> None:
        for p in (0.7, 0.75, 0.8358058):
            beta = critical_parameters(p).beta
            for size in (1, 4, 16):
                self.assertAlmostEqual(
                    node_oracle_reliability(size, p, beta),
                    critical_oracle_reliability(size, p),
                )
                self.assertAlmostEqual(
                    node_mean_even_probability(size, p, beta),
                    critical_mean_same_parity_probability(size, p),
                )

    def test_factorized_path_formula_matches_full_sign_enumeration(self) -> None:
        reliabilities = (0.2, 0.6, 0.85)
        direct = 0.0
        for signs in itertools.product((-1, 1), repeat=len(reliabilities)):
            probability = 1.0
            for sign, reliability in zip(signs, reliabilities):
                probability *= (1.0 + sign * reliability) / 2.0
            if prod(signs) == 1:
                direct += probability
        self.assertAlmostEqual(
            factorized_path_same_parity_probability(reliabilities), direct
        )

    def test_twisted_transfer_keeps_dependence_between_path_decisions(
        self,
    ) -> None:
        initial = {0: 1.0}
        kernels = (
            {0: ((0, 0, 0.5), (1, 1, 0.5))},
            {
                0: ((0, 0, 1.0),),
                1: ((1, 1, 1.0),),
            },
        )
        self.assertEqual(
            finite_state_path_correlation(initial, kernels), 1.0
        )
        self.assertEqual(
            factorized_path_same_parity_probability((0.0, 0.0)), 0.5
        )

    def test_fair_global_orientation_splits_same_and_opposite_states(
        self,
    ) -> None:
        correlation = 0.37
        probabilities = symmetric_pair_status_probabilities(correlation)
        self.assertAlmostEqual(probabilities.p00, (1.0 + correlation) / 4.0)
        self.assertAlmostEqual(probabilities.p11, probabilities.p00)
        self.assertAlmostEqual(probabilities.p01, (1.0 - correlation) / 4.0)
        self.assertAlmostEqual(probabilities.p10, probabilities.p01)
        self.assertAlmostEqual(probabilities.even, (1.0 + correlation) / 2.0)


if __name__ == "__main__":
    unittest.main()
