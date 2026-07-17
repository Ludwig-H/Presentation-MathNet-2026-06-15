from __future__ import annotations

import unittest
from math import atanh, cosh, exp, log, pi, sqrt

from critical_band_thresholds import Q_CRITICAL, open_probability
from critical_merger_oracle import (
    P_INFO,
    P_SW,
    critical_bayes_error,
    critical_local_log_odds,
    critical_mean_error_exponent,
    critical_mean_error_prefactor,
    critical_mean_opposite_parity_probability,
    critical_mean_same_parity_probability,
    critical_minus_count_pmf,
    critical_oracle_reliability,
    critical_parameters,
    critical_plus_count_pmf,
    critical_scaling_p,
    gaussian_crossover,
    reliability_deficit_bound,
    symmetric_experiment_reliability,
)


class CriticalMergerOracleTests(unittest.TestCase):
    def test_closed_forms_match_clock_definitions(self) -> None:
        for p in (P_SW, P_INFO, 0.8358058, 0.9):
            parameters = critical_parameters(p)
            coupling = log(p / (1.0 - p))
            exponential = exp(-coupling * parameters.beta)
            satisfaction = p * exponential / (1.0 - p + p * exponential)

            self.assertAlmostEqual(
                open_probability(p, parameters.beta), Q_CRITICAL
            )
            self.assertAlmostEqual(parameters.satisfaction, satisfaction)
            self.assertAlmostEqual(parameters.bias, 2.0 * satisfaction - 1.0)
            self.assertAlmostEqual(
                parameters.residual_log_odds,
                coupling * (1.0 - parameters.beta),
            )
            self.assertAlmostEqual(
                parameters.residual_log_odds, 2.0 * atanh(parameters.bias)
            )

    def test_sw_endpoint_is_exactly_terminal(self) -> None:
        parameters = critical_parameters(P_SW)
        self.assertAlmostEqual(parameters.beta, 1.0)
        self.assertEqual(parameters.satisfaction, 0.5)
        self.assertEqual(parameters.bias, 0.0)
        self.assertEqual(parameters.residual_log_odds, 0.0)

    def test_two_count_experiments_are_normalized_and_mirrored(self) -> None:
        for size in (1, 2, 7, 20):
            plus = critical_plus_count_pmf(size, 0.81)
            minus = critical_minus_count_pmf(size, 0.81)
            self.assertAlmostEqual(sum(plus.values()), 1.0)
            self.assertAlmostEqual(sum(minus.values()), 1.0)
            for count in range(size + 1):
                self.assertAlmostEqual(
                    minus[count], plus.get(size - count, 0.0)
                )

    def test_heat_bath_log_odds_matches_direct_likelihood_ratio(self) -> None:
        for size in (2, 5, 11):
            for p in (P_SW, P_INFO, 0.8358058):
                plus = critical_plus_count_pmf(size, p)
                minus = critical_minus_count_pmf(size, p)
                for count in range(1, size):
                    direct = log(plus[count] / minus[count])
                    self.assertAlmostEqual(
                        critical_local_log_odds(size, count, p), direct
                    )

    def test_two_reliability_representations_agree(self) -> None:
        for size in (1, 2, 5, 16, 64):
            for p in (P_SW, 0.72, P_INFO, 0.8358058):
                self.assertAlmostEqual(
                    critical_oracle_reliability(size, p),
                    symmetric_experiment_reliability(size, p),
                )

    def test_exact_sw_reliability(self) -> None:
        for size in range(1, 21):
            self.assertAlmostEqual(
                critical_oracle_reliability(size, P_SW), 1.0 / size
            )

    def test_mean_even_probability_is_one_plus_reliability_over_two(
        self,
    ) -> None:
        for size in (1, 2, 5, 16, 64):
            for p in (P_SW, 0.72, P_INFO, 0.8358058):
                probability = critical_mean_same_parity_probability(size, p)
                reliability = critical_oracle_reliability(size, p)
                self.assertAlmostEqual(
                    probability, (1.0 + reliability) / 2.0
                )

    def test_exact_sw_even_probability(self) -> None:
        for size in range(1, 21):
            self.assertAlmostEqual(
                critical_mean_same_parity_probability(size, P_SW),
                0.5 + 0.5 / size,
            )

    def test_even_deficit_is_comparable_to_bayes_error(self) -> None:
        for size in (1, 2, 7, 32):
            for p in (P_SW, 0.72, P_INFO, 0.8358058):
                deficit = critical_mean_opposite_parity_probability(size, p)
                bayes = critical_bayes_error(size, p)
                self.assertGreaterEqual(deficit + 1e-15, bayes)
                self.assertLessEqual(deficit, 2.0 * bayes + 1e-15)

    def test_mean_even_deficit_has_the_claimed_logarithmic_rate(self) -> None:
        p = 0.75
        exponent = critical_mean_error_exponent(p)
        errors = []
        for size in (128, 256, 512):
            deficit = critical_mean_opposite_parity_probability(size, p)
            errors.append(abs(-log(deficit) / size - exponent))
        self.assertLess(errors[-1], errors[0])
        self.assertLess(errors[-1], 0.02)

    def test_sharp_error_prefactor_matches_its_defining_series(self) -> None:
        p = 0.75
        parameters = critical_parameters(p)
        a_c = parameters.residual_log_odds
        for parity, shift in ((0, 0.0), (1, 0.5)):
            direct_series = sum(
                1.0 / cosh(a_c * (index + shift))
                for index in range(-100, 101)
            )
            direct = direct_series / (
                2.0 * parameters.satisfaction * sqrt(2.0 * pi)
            )
            self.assertAlmostEqual(
                critical_mean_error_prefactor(p, parity), direct
            )

    def test_mean_even_deficit_has_the_sharp_prefactor(self) -> None:
        p = 0.75
        exponent = critical_mean_error_exponent(p)
        for parity in (0, 1):
            ratios = []
            for base_size in (128, 256, 512):
                size = base_size + parity
                asymptotic = (
                    critical_mean_error_prefactor(p, parity)
                    * exp(-size * exponent)
                    / sqrt(size)
                )
                exact = critical_mean_opposite_parity_probability(size, p)
                ratios.append(exact / asymptotic)
            self.assertLess(abs(1.0 - ratios[-1]), abs(1.0 - ratios[0]))
            self.assertGreater(ratios[-1], 0.94)

    def test_sharp_prefactor_rejects_the_critical_endpoint(self) -> None:
        with self.assertRaises(ValueError):
            critical_mean_error_prefactor(P_SW, 0)

    def test_even_and_odd_prefactors_are_genuinely_distinct(self) -> None:
        self.assertNotAlmostEqual(
            critical_mean_error_prefactor(0.9, 0),
            critical_mean_error_prefactor(0.9, 1),
            places=4,
        )

    def test_single_edge_oracle_is_perfect(self) -> None:
        for p in (P_SW, 0.72, P_INFO, 0.8358058, 0.95):
            self.assertEqual(critical_oracle_reliability(1, p), 1.0)

    def test_proved_exponential_deficit_bound(self) -> None:
        for size in (2, 8, 64, 256):
            for p in (P_SW, 0.72, P_INFO, 0.8358058):
                deficit = 1.0 - critical_oracle_reliability(size, p)
                self.assertLessEqual(
                    deficit, reliability_deficit_bound(size, p) + 1e-14
                )

    def test_fixed_supercritical_parameter_tends_to_one(self) -> None:
        self.assertGreater(critical_oracle_reliability(256, P_INFO), 0.999)
        self.assertGreater(critical_oracle_reliability(64, 0.8358058), 0.999)

    def test_finite_size_scaling_matches_gaussian_limit(self) -> None:
        size = 1024
        for alpha in (0.5, 1.0, 2.0):
            p = critical_scaling_p(size, alpha)
            finite = critical_oracle_reliability(size, p)
            limit = gaussian_crossover(alpha, intervals=12000)
            self.assertLess(abs(finite - limit), 0.03)


if __name__ == "__main__":
    unittest.main()
