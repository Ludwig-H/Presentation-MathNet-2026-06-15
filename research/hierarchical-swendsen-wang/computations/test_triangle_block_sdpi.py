from __future__ import annotations

import random
import unittest
from math import sqrt

from triangle_block_sdpi import (
    P_FIXED_PRIOR_TRIANGLE,
    P_INFORMATION,
    P_SCALAR_TRIANGLE,
    conditional_candidate_polynomial,
    conditional_candidate_q,
    direct_chi_square_ratio,
    erasure_envelope,
    less_noisy_gap_matrix,
    prior_profile,
    profile_chi_square_ratio,
    projection_diagonal_lower_bound,
    projection_variance_sum,
    sampled_line_defect,
    triangle_channel,
    triangle_uniform_contraction,
    two_observation_contraction,
)


class TriangleBlockSdpiTests(unittest.TestCase):
    def test_channel_rows_are_probabilities(self) -> None:
        for q in (0.0, 0.2, 0.6, 0.9):
            for row in triangle_channel(q):
                self.assertTrue(all(value >= 0.0 for value in row))
                self.assertAlmostEqual(sum(row), 1.0)

    def test_fixed_prior_profile_matches_direct_enumeration(self) -> None:
        generator = random.Random(20260716)
        for q in (0.15, 0.4, 0.62, 0.9):
            for _ in range(80):
                raw = [generator.expovariate(1.0) for _ in range(4)]
                total = sum(raw)
                prior = tuple(value / total for value in raw)
                function = tuple(generator.uniform(-2.0, 2.0) for _ in range(4))
                direct = direct_chi_square_ratio(q, prior, function)
                closed = profile_chi_square_ratio(q, prior, function)
                self.assertAlmostEqual(direct, closed, places=12)

    def test_uniform_and_binary_calibrations(self) -> None:
        characters = (1.0, 1.0, -1.0, -1.0)
        for q in (0.1, 0.4, 0.62, 0.9):
            uniform = direct_chi_square_ratio(
                q, (0.25, 0.25, 0.25, 0.25), characters
            )
            binary = direct_chi_square_ratio(
                q, (0.5, 0.5, 0.0, 0.0), (1.0, -1.0, 0.0, 0.0)
            )
            self.assertAlmostEqual(uniform, triangle_uniform_contraction(q))
            self.assertAlmostEqual(binary, two_observation_contraction(q))

    def test_side_information_strictly_increases_contraction(self) -> None:
        for q in (0.1, 0.4, 0.62, 0.9):
            self.assertGreater(
                two_observation_contraction(q), triangle_uniform_contraction(q)
            )

    def test_scalar_triangle_route_is_worse_than_edge_information(self) -> None:
        self.assertLess(P_SCALAR_TRIANGLE, P_INFORMATION)
        self.assertGreater(P_FIXED_PRIOR_TRIANGLE, P_INFORMATION)
        q_naive = 2.0 * P_FIXED_PRIOR_TRIANGLE - 1.0
        self.assertAlmostEqual(triangle_uniform_contraction(q_naive), 0.5)

    def test_projection_lower_bound(self) -> None:
        generator = random.Random(17072026)
        for _ in range(300):
            raw = [generator.expovariate(1.0) for _ in range(4)]
            total = sum(raw)
            prior = tuple(value / total for value in raw)
            function = tuple(generator.uniform(-3.0, 3.0) for _ in range(4))
            self.assertGreaterEqual(
                projection_variance_sum(prior, function) + 1e-13,
                projection_diagonal_lower_bound(prior, function),
            )

    def test_gap_matrix_matches_both_quadratic_forms(self) -> None:
        generator = random.Random(18072026)
        for q in (0.2, 0.45, 0.62):
            candidate = erasure_envelope(q)
            for _ in range(100):
                raw = [generator.expovariate(1.0) for _ in range(4)]
                total = sum(raw)
                prior = tuple(value / total for value in raw)
                function = tuple(
                    generator.uniform(-2.0, 2.0) for _ in range(4)
                )
                mean = sum(
                    mass * value for mass, value in zip(prior, function)
                )
                centered = tuple(value - mean for value in function)
                variance = sum(
                    mass * value * value
                    for mass, value in zip(prior, centered)
                )
                physical = sum(
                    mass * prior_profile(q, mass) * value * value
                    for mass, value in zip(prior, centered)
                )
                erasure = (
                    candidate.full * variance
                    + candidate.single
                    * projection_variance_sum(prior, function)
                )
                coordinates = tuple(
                    prior[index] * centered[index]
                    for index in range(1, 4)
                )
                matrix = less_noisy_gap_matrix(
                    q,
                    prior,
                    candidate.full,
                    candidate.single,
                )
                matrix_gap = sum(
                    coordinates[row]
                    * matrix[row][column]
                    * coordinates[column]
                    for row in range(3)
                    for column in range(3)
                )
                self.assertAlmostEqual(erasure - physical, matrix_gap, places=11)

    def test_endpoint_only_erasure_parameters_fail_the_full_profile(self) -> None:
        q = 2.0 * P_INFORMATION - 1.0
        uniform = triangle_uniform_contraction(q)
        global_coefficient = two_observation_contraction(q)
        naive_single = global_coefficient - uniform
        naive_full = 2.0 * uniform - global_coefficient
        defect = sampled_line_defect(q, naive_full, naive_single)
        self.assertLess(defect, -1e-4)

    def test_affine_envelope_and_conditional_candidate(self) -> None:
        q = conditional_candidate_q()
        candidate = erasure_envelope(q)
        square = q * q
        self.assertAlmostEqual(q, 0.6198185785038384, places=11)
        self.assertAlmostEqual((1.0 + q) / 2.0, 0.8099092892519192, places=11)
        self.assertAlmostEqual(
            candidate.full,
            2.0 * (1.0 - 2.0 * square) / (1.0 + square),
            places=10,
        )
        self.assertAlmostEqual(
            candidate.single,
            (3.0 * square - 1.0) / (1.0 + square),
            places=10,
        )
        self.assertAlmostEqual(
            candidate.tangent_mass,
            0.5
            - square * (1.0 - square) / (8.0 * (3.0 * square - 1.0)),
            places=10,
        )
        self.assertAlmostEqual(candidate.full, candidate.empty, places=10)
        self.assertAlmostEqual(candidate.self_dual_score, 1.0, places=10)
        self.assertGreaterEqual(
            sampled_line_defect(q, candidate.full, candidate.single),
            -1e-11,
        )
        self.assertLess(abs(conditional_candidate_polynomial(q)), 1e-10)

    def test_candidate_survives_randomized_polarized_psd_audit(self) -> None:
        """Search for PSD counterexamples to P; this finite audit is not a proof."""

        generator = random.Random(16072026)
        q_star = conditional_candidate_q()
        for q in (0.2, 0.45, 0.6, q_star):
            candidate = erasure_envelope(q)
            for _ in range(1000):
                dominant_mass = generator.uniform(0.500001, 0.999999)
                tail = [generator.expovariate(1.0) for _ in range(3)]
                tail_total = sum(tail)
                prior = [dominant_mass]
                prior.extend(
                    (1.0 - dominant_mass) * value / tail_total
                    for value in tail
                )
                generator.shuffle(prior)
                matrix = less_noisy_gap_matrix(
                    q,
                    prior,
                    candidate.full,
                    candidate.single,
                )
                diagonal = tuple(matrix[index][index] for index in range(3))
                two_by_two = tuple(
                    diagonal[first] * diagonal[second]
                    - matrix[first][second] ** 2
                    for first, second in ((0, 1), (0, 2), (1, 2))
                )
                determinant = (
                    matrix[0][0]
                    * (matrix[1][1] * matrix[2][2] - matrix[1][2] ** 2)
                    - matrix[0][1]
                    * (
                        matrix[0][1] * matrix[2][2]
                        - matrix[0][2] * matrix[1][2]
                    )
                    + matrix[0][2]
                    * (
                        matrix[0][1] * matrix[1][2]
                        - matrix[0][2] * matrix[1][1]
                    )
                )
                scale = diagonal[0] * diagonal[1] * diagonal[2]
                self.assertTrue(all(value >= -1e-11 for value in diagonal))
                self.assertTrue(all(value >= -1e-11 for value in two_by_two))
                self.assertGreaterEqual(determinant / scale, -1e-11)

    def test_chayes_lei_side_conditions_hold_at_candidate(self) -> None:
        candidate = erasure_envelope(conditional_candidate_q())
        threshold = 2.0 * sqrt(2.0) / (3.0 + 2.0 * sqrt(2.0))
        self.assertGreaterEqual(
            candidate.full * candidate.empty,
            2.0 * candidate.single * candidate.single,
        )
        self.assertGreater(candidate.full + candidate.empty, threshold)


if __name__ == "__main__":
    unittest.main()
