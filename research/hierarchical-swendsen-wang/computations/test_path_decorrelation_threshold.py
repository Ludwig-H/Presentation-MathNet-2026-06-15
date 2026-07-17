from __future__ import annotations

import unittest
from math import exp, log, sqrt

from critical_band_thresholds import Q_CRITICAL
from critical_merger_oracle import (
    critical_mean_error_exponent,
    critical_parameters,
)
from hierarchical_flip_probabilities import (
    factorized_path_same_parity_probability,
    node_oracle_reliability,
)
from nishimori_hierarchical_entropy import conjectured_nishimori_root
from path_decorrelation_threshold import (
    attenuation_from_reliabilities,
    critical_geometry_log_partition,
    critical_log_cut_coordinate,
    critical_log_cut_limit_probability,
    critical_p_from_error_exponent,
    critical_reliability_deficit,
    descendant_geometry_log_partition,
    heterogeneous_critical_path_same_relation_probability,
    heterogeneous_descendant_path_same_relation_probability,
    high_p_critical_deficit_equivalent,
    high_p_deficit_constant,
    high_p_deficit_exponent,
    log_cut_coefficient_for_threshold,
    nishimori_log_cut_coefficient,
    regular_critical_path_same_relation_probability,
    regular_log_cut_threshold,
    regular_path_same_relation_probability,
    regular_relative_level_threshold,
    relative_level_error_exponent,
    required_attenuation_for_accuracy,
    same_relation_probability_from_attenuation,
)


class PathDecorrelationThresholdTests(unittest.TestCase):
    def test_attenuation_is_exactly_the_logarithm_of_the_path_product(
        self,
    ) -> None:
        reliabilities = (0.23, 0.61, 0.88, 1.0)
        attenuation = attenuation_from_reliabilities(reliabilities)
        self.assertAlmostEqual(
            same_relation_probability_from_attenuation(attenuation),
            factorized_path_same_parity_probability(reliabilities),
        )

    def test_accuracy_threshold_is_exact(self) -> None:
        for tolerance in (0.2, 0.05, 1e-4):
            attenuation = required_attenuation_for_accuracy(tolerance)
            self.assertAlmostEqual(
                same_relation_probability_from_attenuation(attenuation),
                0.5 + tolerance,
            )

    def test_fixed_nontrivial_channel_loses_correlation(self) -> None:
        self.assertLess(
            regular_path_same_relation_probability(10_000, 0.99),
            0.5000000000000001,
        )
        self.assertEqual(
            regular_path_same_relation_probability(10_000, 1.0), 1.0
        )

    def test_log_cut_threshold_inverts_the_critical_rate(self) -> None:
        for alpha in (0.5, 1.0, 2.0, 7.0, 20.0):
            threshold = regular_log_cut_threshold(alpha)
            self.assertAlmostEqual(
                alpha * critical_mean_error_exponent(threshold), 1.0
            )

    def test_critical_exponent_inverse_includes_both_endpoints(self) -> None:
        self.assertAlmostEqual(
            critical_p_from_error_exponent(0.0),
            0.5 * (1.0 + Q_CRITICAL),
        )
        self.assertEqual(critical_p_from_error_exponent(float("inf")), 1.0)
        for exponent in (0.01, 0.2, 1.0, 4.0):
            self.assertAlmostEqual(
                critical_mean_error_exponent(
                    critical_p_from_error_exponent(exponent)
                ),
                exponent,
            )

    def test_heterogeneous_path_and_geometric_partition(self) -> None:
        p = 0.82
        sizes = (7, 11, 7, 20)
        reliabilities = []
        beta = critical_parameters(p).beta
        for size in sizes:
            reliabilities.append(node_oracle_reliability(size, p, beta))
        self.assertAlmostEqual(
            heterogeneous_critical_path_same_relation_probability(sizes, p),
            factorized_path_same_parity_probability(reliabilities),
        )

        alpha = 4.0
        exponent = critical_mean_error_exponent(p)
        below_length = 1_000
        below_size = round(alpha * log(below_length))
        direct = critical_geometry_log_partition(
            [below_size] * below_length, p
        )
        expected = (
            log(below_length)
            - 0.5 * log(below_size)
            - exponent * below_size
        )
        self.assertAlmostEqual(direct, expected)

    def test_time_decorated_partition_and_exact_path_at_critical_level(
        self,
    ) -> None:
        p = 0.82
        beta = critical_parameters(p).beta
        sizes = (8, 13, 21)
        critical_times = (beta,) * len(sizes)
        self.assertAlmostEqual(
            descendant_geometry_log_partition(sizes, critical_times, p),
            critical_geometry_log_partition(sizes, p),
        )
        self.assertAlmostEqual(
            heterogeneous_descendant_path_same_relation_probability(
                sizes, critical_times, p
            ),
            heterogeneous_critical_path_same_relation_probability(sizes, p),
        )

        early_times = tuple(0.5 * beta for _ in sizes)
        self.assertLess(
            descendant_geometry_log_partition(sizes, early_times, p),
            descendant_geometry_log_partition(sizes, critical_times, p),
        )
        self.assertGreater(
            heterogeneous_descendant_path_same_relation_probability(
                sizes, early_times, p
            ),
            heterogeneous_descendant_path_same_relation_probability(
                sizes, critical_times, p
            ),
        )

    def test_relative_level_threshold_is_monotone_and_has_exact_endpoints(
        self,
    ) -> None:
        alpha = nishimori_log_cut_coefficient()
        thresholds = [
            regular_relative_level_threshold(alpha, theta)
            for theta in (0.0, 0.25, 0.5, 0.75, 1.0)
        ]
        self.assertTrue(all(value is not None for value in thresholds))
        numeric = [float(value) for value in thresholds if value is not None]
        self.assertEqual(numeric, sorted(numeric))
        self.assertAlmostEqual(numeric[-1], regular_log_cut_threshold(alpha))

        bias = sqrt(1.0 - exp(-2.0 / alpha))
        self.assertAlmostEqual(numeric[0], 0.5 * (1.0 + bias))
        for theta, threshold in zip(
            (0.0, 0.25, 0.5, 0.75, 1.0), numeric, strict=True
        ):
            self.assertAlmostEqual(
                alpha * relative_level_error_exponent(threshold, theta),
                1.0,
            )

        self.assertIsNone(regular_relative_level_threshold(20.0, 0.0))

    def test_near_critical_descendant_band_has_inverse_size_width(
        self,
    ) -> None:
        p = 0.82
        parameters = critical_parameters(p)
        epsilon = 1e-6
        finite_difference = (
            relative_level_error_exponent(p, 1.0 - epsilon)
            - relative_level_error_exponent(p, 1.0)
        ) / epsilon
        expected = (
            0.5
            * parameters.bias
            * parameters.beta
            * log(p / (1.0 - p))
        )
        self.assertAlmostEqual(finite_difference, expected, delta=1e-6)

    def test_nishimori_value_only_calibrates_one_geometric_coefficient(
        self,
    ) -> None:
        p_nishimori = conjectured_nishimori_root()
        alpha = nishimori_log_cut_coefficient()
        self.assertAlmostEqual(regular_log_cut_threshold(alpha), p_nishimori)
        self.assertAlmostEqual(alpha, 7.0535961929, places=9)

        p_information = 0.5 * (1.0 + sqrt(Q_CRITICAL))
        alpha_information = log_cut_coefficient_for_threshold(p_information)
        self.assertAlmostEqual(
            regular_log_cut_threshold(alpha_information), p_information
        )
        self.assertAlmostEqual(alpha_information, 13.5216281646, places=9)

    def test_fixed_bucket_high_p_equivalent_and_constants(self) -> None:
        expected = {
            2: (1, 1),
            3: (2, 4),
            4: (2, 3),
            5: (3, 16),
            6: (3, 10),
            7: (4, 60),
            8: (4, 35),
        }
        p = 0.999
        beta = critical_parameters(p).beta
        for size, (exponent, constant) in expected.items():
            with self.subTest(size=size):
                self.assertEqual(high_p_deficit_exponent(size), exponent)
                self.assertEqual(high_p_deficit_constant(size), constant)
                exact = 1.0 - node_oracle_reliability(size, p, beta)
                equivalent = high_p_critical_deficit_equivalent(size, p)
                self.assertAlmostEqual(exact / equivalent, 1.0, delta=0.004)

    def test_stable_deficit_matches_one_minus_direct_reliability(self) -> None:
        for size in (2, 5, 20, 80):
            p = 0.8
            beta = critical_parameters(p).beta
            direct = 1.0 - node_oracle_reliability(size, p, beta)
            stable = critical_reliability_deficit(size, p)
            self.assertAlmostEqual(direct, stable, places=13)

    def test_sharp_log_cut_coordinate_predicts_the_path_limit(self) -> None:
        p = 0.8
        size = 320
        target_coordinate = 0.0
        exponent = critical_mean_error_exponent(p)
        length = round(exp(size * exponent) * sqrt(size))
        coordinate = critical_log_cut_coordinate(length, size, p)
        exact = regular_critical_path_same_relation_probability(
            length, size, p
        )
        predicted = critical_log_cut_limit_probability(
            p, size % 2, target_coordinate
        )
        self.assertAlmostEqual(coordinate, target_coordinate, places=10)
        self.assertAlmostEqual(exact, predicted, delta=0.005)

    def test_descendant_high_p_residual_error_scale(self) -> None:
        p = 0.9999
        beta = critical_parameters(p).beta
        for theta in (0.0, 0.4, 1.0):
            time = theta * beta
            direct_error = 1.0 / (
                1.0
                + exp(
                    log(p / (1.0 - p))
                    * (1.0 - time)
                )
            )
            equivalent = (1.0 - p) / (1.0 - Q_CRITICAL) ** theta
            with self.subTest(theta=theta):
                self.assertAlmostEqual(
                    direct_error / equivalent, 1.0, delta=0.002
                )


if __name__ == "__main__":
    unittest.main()
