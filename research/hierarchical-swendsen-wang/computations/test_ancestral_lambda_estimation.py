from __future__ import annotations

import random
import unittest
from itertools import product
from math import atanh, exp, log, prod, sqrt

from ancestral_lambda_chain import (
    AncestorBucket,
    STATES,
    ancestral_log_weights,
    ancestral_message,
    closed_satisfaction_probability as homogeneous_closed_probability,
    coupling,
    grouped_count_pmf,
    walsh_coefficients,
)
from ancestral_lambda_estimation import (
    RELIABILITY_LIPSCHITZ_CONSTANT,
    WeightedEdge,
    ancestral_tail_bound,
    closed_satisfaction_probability,
    conditional_group_moments,
    contrast_envelope,
    four_rate_moments,
    homogeneous_four_rate_moments,
    lca_reliability,
    reliability_tail_bound,
    winner_probabilities,
)


def brute_force_group_moments(
    edges: tuple[WeightedEdge, ...], beta: float
) -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """Independent enumeration of the winner and every residual Bernoulli."""

    winner = winner_probabilities(edges, beta)
    probabilities = tuple(
        closed_satisfaction_probability(edge, beta) for edge in edges
    )
    atoms: list[tuple[tuple[float, float, float], float]] = []
    for winner_index, winner_mass in enumerate(winner):
        others = tuple(index for index in range(len(edges)) if index != winner_index)
        for marks in product((0, 1), repeat=len(others)):
            values = [0.0, 0.0, 0.0]
            winning_edge = edges[winner_index]
            values[winning_edge.group] += winning_edge.weight
            mass = winner_mass
            for index, mark in zip(others, marks):
                edge = edges[index]
                probability = probabilities[index]
                mass *= probability if mark else 1.0 - probability
                if mark:
                    values[edge.group] += edge.weight
            atoms.append(((values[0], values[1], values[2]), mass))

    mean = tuple(sum(value[r] * mass for value, mass in atoms) for r in range(3))
    covariance = tuple(
        tuple(
            sum(
                (value[r] - mean[r]) * (value[s] - mean[s]) * mass
                for value, mass in atoms
            )
            for s in range(3)
        )
        for r in range(3)
    )
    return mean, covariance


class WeightedAncestorEstimationTests(unittest.TestCase):
    def test_nishimori_closed_probability_simplifies_to_logistic(self) -> None:
        for weight in (0.2, 1.3, 4.0):
            for beta in (0.0, 0.4, 1.0):
                edge = WeightedEdge.nishimori(weight, 0)
                expected = 1.0 / (1.0 + exp(-weight * (1.0 - beta)))
                self.assertAlmostEqual(
                    closed_satisfaction_probability(edge, beta), expected
                )

    def test_weighted_winner_is_conditional_hazard_not_uniform(self) -> None:
        edges = (
            WeightedEdge(0.4, 0, 0.7),
            WeightedEdge(1.1, 1, 0.6),
            WeightedEdge(2.0, 2, 0.8),
        )
        beta = 0.55
        probabilities = tuple(
            closed_satisfaction_probability(edge, beta) for edge in edges
        )
        hazards = tuple(
            edge.weight * probability
            for edge, probability in zip(edges, probabilities)
        )
        expected = tuple(value / sum(hazards) for value in hazards)
        actual = winner_probabilities(edges, beta)
        for left, right in zip(actual, expected):
            self.assertAlmostEqual(left, right)
        self.assertGreater(max(actual) - min(actual), 0.1)

        # Counter-audit the hazard factorization against the original race
        # density, before any algebraic cancellation.
        survival = tuple(
            1.0
            - edge.prior_satisfied
            + edge.prior_satisfied * exp(-edge.weight * beta)
            for edge in edges
        )
        densities = tuple(
            edge.prior_satisfied
            * edge.weight
            * exp(-edge.weight * beta)
            * prod(
                survival[other]
                for other in range(len(edges))
                if other != index
            )
            for index, edge in enumerate(edges)
        )
        direct = tuple(value / sum(densities) for value in densities)
        for left, right in zip(actual, direct):
            self.assertAlmostEqual(left, right)

    def test_equal_weights_reduce_to_uniform_winner(self) -> None:
        weight = coupling(0.81)
        edges = tuple(
            WeightedEdge.nishimori(weight, index % 3) for index in range(7)
        )
        for probability in winner_probabilities(edges, 0.63):
            self.assertAlmostEqual(probability, 1.0 / len(edges))

    def test_homogeneous_generating_polynomial_matches_grouped_law(self) -> None:
        group_sizes = (3, 2, 1)
        p = 0.79
        beta = 0.43
        z = (0.7, 1.2, 0.9)
        distribution = grouped_count_pmf(group_sizes, p, beta)
        enumerated = sum(
            mass * prod(z[group] ** counts[group] for group in range(3))
            for counts, mass in distribution.items()
        )
        s = homogeneous_closed_probability(p, beta)
        size = sum(group_sizes)
        polynomial = 0.0
        for winner_group, winner_size in enumerate(group_sizes):
            if winner_size == 0:
                continue
            term = winner_size / size * z[winner_group]
            for group, group_size in enumerate(group_sizes):
                exponent = group_size - int(group == winner_group)
                term *= (1.0 - s + s * z[group]) ** exponent
            polynomial += term
        self.assertAlmostEqual(enumerated, polynomial)

    def test_weighted_moments_match_independent_full_enumeration(self) -> None:
        edges = (
            WeightedEdge(0.4, 0, 0.68),
            WeightedEdge(0.9, 0, 0.73),
            WeightedEdge(1.2, 1, 0.61),
            WeightedEdge(1.8, 2, 0.82),
        )
        expected_mean, expected_covariance = brute_force_group_moments(edges, 0.47)
        actual = conditional_group_moments(edges, 0.47)
        for left, right in zip(actual.mean, expected_mean):
            self.assertAlmostEqual(left, right)
        for r in range(3):
            for s in range(3):
                self.assertAlmostEqual(
                    actual.covariance[r][s], expected_covariance[r][s]
                )

    def test_random_weighted_moments_match_full_enumeration(self) -> None:
        generator = random.Random(104729)
        for _ in range(100):
            edges = tuple(
                WeightedEdge(
                    weight=generator.uniform(0.1, 2.5),
                    group=generator.randrange(3),
                    prior_satisfied=generator.uniform(0.52, 0.93),
                )
                for _ in range(generator.randint(1, 6))
            )
            beta = generator.random()
            expected_mean, expected_covariance = brute_force_group_moments(
                edges, beta
            )
            actual = conditional_group_moments(edges, beta)
            for left, right in zip(actual.mean, expected_mean):
                self.assertAlmostEqual(left, right, places=11)
            for r in range(3):
                for s in range(3):
                    self.assertAlmostEqual(
                        actual.covariance[r][s],
                        expected_covariance[r][s],
                        places=11,
                    )

    def test_stated_hoeffding_and_bernstein_radii_cover_exact_tails(self) -> None:
        generator = random.Random(65537)
        for _ in range(100):
            edges = tuple(
                WeightedEdge(
                    weight=generator.uniform(0.1, 1.8),
                    group=generator.randrange(3),
                    prior_satisfied=generator.uniform(0.55, 0.9),
                )
                for _ in range(generator.randint(2, 7))
            )
            beta = generator.random()
            probabilities = tuple(
                closed_satisfaction_probability(edge, beta) for edge in edges
            )
            winner_index = generator.randrange(len(edges))
            group = generator.randrange(3)
            residual = tuple(
                index
                for index, edge in enumerate(edges)
                if index != winner_index and edge.group == group
            )
            atoms: list[tuple[float, float]] = []
            for marks in product((0, 1), repeat=len(residual)):
                centered = 0.0
                mass = 1.0
                for index, mark in zip(residual, marks):
                    edge = edges[index]
                    probability = probabilities[index]
                    centered += edge.weight * (mark - probability)
                    mass *= probability if mark else 1.0 - probability
                atoms.append((centered, mass))
            variance = sum(
                edges[index].weight**2
                * probabilities[index]
                * (1.0 - probabilities[index])
                for index in residual
            )
            square_range = sum(edges[index].weight**2 for index in residual)
            maximum = max((edges[index].weight for index in residual), default=0.0)
            if not residual:
                self.assertEqual(atoms, [(0.0, 1.0)])
                continue
            for x in (1.5, 3.0, 5.0):
                hoeffding_radius = (square_range * x / 2.0) ** 0.5
                bernstein_radius = (
                    (2.0 * variance * x) ** 0.5 + 2.0 * maximum * x / 3.0
                )
                bound = 2.0 * exp(-x)
                hoeffding_tail = sum(
                    mass
                    for centered, mass in atoms
                    if abs(centered) + 1e-13 >= hoeffding_radius
                )
                bernstein_tail = sum(
                    mass
                    for centered, mass in atoms
                    if abs(centered) + 1e-13 >= bernstein_radius
                )
                self.assertLessEqual(hoeffding_tail, bound + 1e-12)
                self.assertLessEqual(bernstein_tail, bound + 1e-12)

    def test_affine_four_rate_moments_match_enumeration(self) -> None:
        edges = (
            WeightedEdge(0.5, 0, 0.7),
            WeightedEdge(0.8, 1, 0.65),
            WeightedEdge(1.3, 1, 0.75),
            WeightedEdge(1.6, 2, 0.8),
        )
        group_mean, group_covariance = brute_force_group_moments(edges, 0.51)
        totals = tuple(
            sum(edge.weight for edge in edges if edge.group == group)
            for group in range(3)
        )
        expected_mean: list[float] = []
        expected_covariance: list[list[float]] = []
        coefficients = tuple(
            (1.0, 1.0 - 2.0 * a, 1.0 - 2.0 * b) for a, b in STATES
        )
        for (a, b), coefficient in zip(STATES, coefficients):
            expected_mean.append(
                a * totals[1]
                + b * totals[2]
                + sum(coefficient[r] * group_mean[r] for r in range(3))
            )
        for first in coefficients:
            expected_covariance.append(
                [
                    sum(
                        first[r] * group_covariance[r][s] * second[s]
                        for r in range(3)
                        for s in range(3)
                    )
                    for second in coefficients
                ]
            )
        actual = four_rate_moments(totals, conditional_group_moments(edges, 0.51))
        for left, right in zip(actual.mean, expected_mean):
            self.assertAlmostEqual(left, right)
        for r in range(4):
            for s in range(4):
                self.assertAlmostEqual(
                    actual.covariance[r][s], expected_covariance[r][s]
                )

    def test_homogeneous_closed_forms_for_mean_and_variance(self) -> None:
        group_sizes = (5, 3, 2)
        p = 0.82
        beta = 0.64
        weight = coupling(p)
        s = 1.0 / (1.0 + exp(-weight * (1.0 - beta)))
        h = 2.0 * s - 1.0
        size = sum(group_sizes)
        actual = homogeneous_four_rate_moments(group_sizes, p, beta)
        for index, (a, b) in enumerate(STATES):
            signed_size = (
                group_sizes[0]
                + (1 - 2 * a) * group_sizes[1]
                + (1 - 2 * b) * group_sizes[2]
            )
            alpha = h + (1.0 - h) / size
            expected_mean = 0.5 * weight * (size + alpha * signed_size)
            expected_variance = weight**2 * (
                s * (1.0 - s) * (size - 1)
                + (1.0 - s) ** 2
                * (1.0 - (signed_size / size) ** 2)
            )
            self.assertAlmostEqual(actual.mean[index], expected_mean)
            self.assertAlmostEqual(
                actual.covariance[index][index], expected_variance
            )

    def test_homogeneous_closed_form_matches_generic_weighted_kernel(self) -> None:
        for group_sizes in ((1, 0, 0), (2, 1, 3), (4, 4, 2)):
            p = 0.78
            beta = 0.57
            weight = coupling(p)
            edges = tuple(
                WeightedEdge.nishimori(weight, group)
                for group, size in enumerate(group_sizes)
                for _ in range(size)
            )
            generic = four_rate_moments(
                tuple(weight * size for size in group_sizes),
                conditional_group_moments(edges, beta),
            )
            closed = homogeneous_four_rate_moments(group_sizes, p, beta)
            for left, right in zip(closed.mean, generic.mean):
                self.assertAlmostEqual(left, right)
            for r in range(4):
                for s in range(4):
                    self.assertAlmostEqual(
                        closed.covariance[r][s], generic.covariance[r][s]
                    )

    def test_contrast_envelope_bounds_direct_walsh_coefficients(self) -> None:
        generator = random.Random(20260716)
        for _ in range(250):
            totals = tuple(generator.uniform(2.0, 12.0) for _ in range(3))
            satisfied = tuple(
                generator.uniform(0.2 * total, 0.8 * total) for total in totals
            )
            bucket = AncestorBucket(
                beta=generator.uniform(0.05, 0.98),
                total=totals,
                satisfied=satisfied,
            )
            rates = bucket.four_lambdas()
            if min(rates.values()) <= 0.0:
                continue
            coefficients = walsh_coefficients(ancestral_log_weights((bucket,)))
            envelope = contrast_envelope(bucket)
            self.assertLessEqual(abs(coefficients.field_1), envelope.field_1 + 1e-12)
            self.assertLessEqual(abs(coefficients.field_2), envelope.field_2 + 1e-12)
            self.assertLessEqual(abs(coefficients.coupling), envelope.coupling + 1e-12)

    def test_tail_bound_controls_the_actual_message_change(self) -> None:
        buckets = (
            AncestorBucket(0.60, (8.0, 6.0, 5.0), (4.5, 3.4, 2.8)),
            AncestorBucket(0.83, (7.0, 4.0, 6.0), (3.8, 2.1, 3.2)),
        )
        coefficients = walsh_coefficients(ancestral_log_weights(buckets))
        actual_message = (
            2.0 * coefficients.coupling
            + log(
                (
                    exp(coefficients.field_1 + coefficients.field_2)
                    + exp(-coefficients.field_1 - coefficients.field_2)
                )
                / (
                    exp(coefficients.field_1 - coefficients.field_2)
                    + exp(-coefficients.field_1 + coefficients.field_2)
                )
            )
        )
        self.assertLessEqual(abs(actual_message), ancestral_tail_bound(buckets))

    def test_tail_bound_is_valid_around_an_arbitrary_retained_chain(self) -> None:
        generator = random.Random(99991)
        for _ in range(250):
            retained: list[AncestorBucket] = []
            omitted: list[AncestorBucket] = []
            for destination in (retained, retained, omitted, omitted):
                totals = tuple(generator.uniform(3.0, 15.0) for _ in range(3))
                satisfied = tuple(
                    generator.uniform(0.15 * total, 0.85 * total)
                    for total in totals
                )
                destination.append(
                    AncestorBucket(
                        beta=generator.uniform(0.1, 0.95),
                        total=totals,
                        satisfied=satisfied,
                    )
                )
            log_prior = {state: generator.uniform(-2.0, 2.0) for state in STATES}
            actual_change = abs(
                ancestral_message((*retained, *omitted), log_prior)
                - ancestral_message(retained, log_prior)
            )
            self.assertLessEqual(
                actual_change, ancestral_tail_bound(omitted) + 1e-11
            )

    def test_reliability_transport_bounds_are_sharp_and_quadratic(self) -> None:
        generator = random.Random(314159)
        for _ in range(1000):
            first = generator.uniform(-20.0, 20.0)
            second = generator.uniform(-20.0, 20.0)
            actual_error = abs(lca_reliability(first) - lca_reliability(second))
            self.assertLessEqual(
                actual_error,
                reliability_tail_bound(abs(first - second)) + 1e-14,
            )
            self.assertLessEqual(
                lca_reliability(first),
                min(1.0, first**2 / 4.0) + 1e-14,
            )

        # The derivative reaches its global maximum when
        # tanh(x / 2) = 1 / sqrt(3).
        maximizer = 2.0 * atanh(1.0 / sqrt(3.0))
        increment = 1e-6
        numerical_derivative = (
            lca_reliability(maximizer + increment)
            - lca_reliability(maximizer - increment)
        ) / (2.0 * increment)
        self.assertAlmostEqual(
            numerical_derivative,
            RELIABILITY_LIPSCHITZ_CONSTANT,
            places=9,
        )


if __name__ == "__main__":
    unittest.main()
