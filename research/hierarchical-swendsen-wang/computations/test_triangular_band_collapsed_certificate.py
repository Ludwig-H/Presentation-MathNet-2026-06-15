from __future__ import annotations

import unittest
from fractions import Fraction
from itertools import product
from math import exp, fsum, log

from critical_band_thresholds import Q_CRITICAL
from triangular_band_collapsed_certificate import (
    BOUNDARY_SPINS,
    ENVIRONMENTS,
    PARITY_PAIRS,
    REPLICATED_BOUNDARIES,
    boundary_polarization_diagnostic,
    cell_edge_relations,
    collapsed_cell_kernel,
    independent_environment_replicated_kernel,
    p_805_uniform_contraction_certificate,
    polarized_cell_second_moment,
    replicated_fourier_sector,
    residual_environment_probability,
    residual_satisfaction_probability,
    shared_environment_replicated_kernel,
    strict_decimal_upper_bound,
    triangular_cell_transfer,
    twisted_chain_second_moment,
    weighted_absolute_contraction,
)


class TriangularBandCollapsedCertificateTests(unittest.TestCase):
    def test_residual_environment_and_each_heat_bath_are_normalized(self) -> None:
        p = 0.805
        self.assertAlmostEqual(
            fsum(
                residual_environment_probability(environment, p, Q_CRITICAL)
                for environment in ENVIRONMENTS
            ),
            1.0,
        )
        for environment in ENVIRONMENTS:
            kernel = collapsed_cell_kernel(environment, p, Q_CRITICAL)
            for row in kernel:
                self.assertAlmostEqual(fsum(row), 1.0)
                self.assertGreater(min(row), 0.0)

    def test_documented_edges_contain_one_triangle_and_two_ports(self) -> None:
        left = (1, -1)
        right = (-1, -1)
        self.assertEqual(
            cell_edge_relations(left, right),
            (-1, 1, 1, 1),
        )

    def test_shared_replicated_kernel_is_stochastic_and_equivariant(self) -> None:
        kernel = shared_environment_replicated_kernel(0.805, Q_CRITICAL)
        for row in kernel:
            self.assertAlmostEqual(fsum(row), 1.0)
            self.assertGreaterEqual(min(row), 0.0)

        index = {state: position for position, state in enumerate(REPLICATED_BOUNDARIES)}

        def flip_first(state):
            first, second = state
            return ((-first[0], -first[1]), second)

        def flip_second(state):
            first, second = state
            return (first, (-second[0], -second[1]))

        for left, right in product(REPLICATED_BOUNDARIES, repeat=2):
            original = kernel[index[left]][index[right]]
            self.assertAlmostEqual(
                original,
                kernel[index[flip_first(left)]][index[flip_first(right)]],
            )
            self.assertAlmostEqual(
                original,
                kernel[index[flip_second(left)]][index[flip_second(right)]],
            )

    def test_trivial_sector_is_the_normalized_mass_transfer(self) -> None:
        transfer = triangular_cell_transfer()
        direct = replicated_fourier_sector(
            shared_environment_replicated_kernel(0.805, Q_CRITICAL), 0, 0
        )
        self.assertEqual(transfer.mass, direct)
        self.assertLess(transfer.mass_row_sum_error, 2e-15)
        for row in transfer.mass:
            self.assertAlmostEqual(fsum(row), 1.0)
            self.assertGreaterEqual(min(row), 0.0)

    def test_shared_chi_tensor_chi_has_a_neutral_uniform_margin(self) -> None:
        transfer = triangular_cell_transfer()
        self.assertEqual(transfer.worst_boundary_state, (1, 1))
        self.assertAlmostEqual(
            transfer.uniform_weight_contraction,
            0.2939937883401138,
        )
        self.assertTrue(transfer.has_strict_uniform_contraction)
        self.assertLess(transfer.uniform_weight_contraction, 0.3)
        self.assertIn("E1+", transfer.scope_label)
        self.assertIn("critical Palm law", transfer.missing_features)

    def test_shared_environment_is_not_replaced_by_two_independent_ones(
        self,
    ) -> None:
        transfer = triangular_cell_transfer()
        independent = independent_environment_replicated_kernel(
            0.805, Q_CRITICAL
        )
        self.assertEqual(
            transfer.independent_chi_tensor_chi,
            replicated_fourier_sector(independent, 1, 1),
        )
        independent_coefficient = weighted_absolute_contraction(
            transfer.independent_chi_tensor_chi, (1.0,) * 4
        )
        self.assertAlmostEqual(independent_coefficient, 0.08643234758257165)
        self.assertGreater(
            transfer.uniform_weight_contraction, independent_coefficient
        )

    def test_transfer_matches_an_independent_global_environment_enumeration(
        self,
    ) -> None:
        p = 0.805
        for initial_parity in (-1, 1):
            for depth in (0, 1, 2):
                with self.subTest(parity=initial_parity, depth=depth):
                    self.assertAlmostEqual(
                        twisted_chain_second_moment(
                            p, depth, initial_parity, Q_CRITICAL
                        ),
                        _brute_chain_second_moment(
                            p, Q_CRITICAL, depth, initial_parity
                        ),
                    )

    def test_reference_chain_values_decay_in_the_neutral_cell_model(self) -> None:
        values = tuple(
            twisted_chain_second_moment(0.805, depth)
            for depth in (1, 2, 3)
        )
        expected = (
            0.2939937883401139,
            0.06497530380616247,
            0.01357548484723256,
        )
        for actual, reference in zip(values, expected, strict=True):
            self.assertAlmostEqual(actual, reference)
        self.assertGreater(values[0], values[1])
        self.assertGreater(values[1], values[2])

    def test_rational_intervals_certify_the_neutral_uniform_weight(self) -> None:
        transfer = triangular_cell_transfer()
        certificate = p_805_uniform_contraction_certificate()
        self.assertTrue(certificate.is_strict)
        self.assertLess(certificate.contraction_upper_bound, 3 / 10)
        self.assertLessEqual(float(certificate.rank[0]), Q_CRITICAL)
        self.assertLessEqual(Q_CRITICAL, float(certificate.rank[1]))
        self.assertLessEqual(
            float(certificate.residual_satisfaction[0]),
            transfer.residual_satisfaction,
        )
        self.assertLessEqual(
            transfer.residual_satisfaction,
            float(certificate.residual_satisfaction[1]),
        )
        for row, interval_row in zip(
            transfer.chi_tensor_chi,
            certificate.twisted_entries,
            strict=True,
        ):
            for value, interval in zip(row, interval_row, strict=True):
                self.assertGreaterEqual(value + 2e-15, float(interval[0]))
                self.assertLessEqual(value - 2e-15, float(interval[1]))

    def test_displayed_decimal_certificate_rounds_strictly_upward(self) -> None:
        certificate = p_805_uniform_contraction_certificate()
        displayed = strict_decimal_upper_bound(
            certificate.contraction_upper_bound, 12
        )
        self.assertEqual(displayed, "0.293993788341")
        whole, fractional = displayed.split(".")
        displayed_fraction = Fraction(
            int(whole) * 10**len(fractional) + int(fractional),
            10**len(fractional),
        )
        self.assertGreater(
            displayed_fraction, certificate.contraction_upper_bound
        )

    def test_unbounded_boundary_polarization_is_an_exact_no_go(self) -> None:
        neutral = twisted_chain_second_moment(0.805, 1)
        fields = (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0)
        values = tuple(
            polarized_cell_second_moment(0.805, field) for field in fields
        )
        self.assertAlmostEqual(values[0], neutral)
        for first, second in zip(values[:-1], values[1:], strict=True):
            self.assertLess(first, second)
        self.assertGreater(values[-1], 0.9999999)
        diagnostic = boundary_polarization_diagnostic(0.805, fields)
        self.assertEqual(diagnostic.second_moments, values)
        self.assertEqual(diagnostic.limiting_second_moment, 1.0)
        self.assertFalse(diagnostic.has_uniform_absolute_contraction)

    def test_input_validation(self) -> None:
        with self.assertRaises(ValueError):
            residual_satisfaction_probability(0.5, Q_CRITICAL)
        with self.assertRaises(ValueError):
            residual_satisfaction_probability(0.805, 0.7)
        with self.assertRaises(ValueError):
            replicated_fourier_sector(((1.0,),), 0, 0)
        with self.assertRaises(ValueError):
            weighted_absolute_contraction(((1.0,),), (0.0,))
        with self.assertRaises(ValueError):
            twisted_chain_second_moment(0.805, -1)
        with self.assertRaises(ValueError):
            polarized_cell_second_moment(0.805, float("inf"))
        with self.assertRaises(ValueError):
            boundary_polarization_diagnostic(0.805, ())


def _brute_chain_second_moment(
    p: float, rank: float, depth: int, initial_parity: int
) -> float:
    """Independent enumeration without using the replicated transfer code."""

    spins = tuple(product((-1, 1), repeat=2))
    environments = tuple(product((-1, 1), repeat=4))
    satisfaction = (p - rank) / (1.0 - rank)

    def sign_probability(signs) -> float:
        answer = 1.0
        for sign in signs:
            answer *= satisfaction if sign == 1 else 1.0 - satisfaction
        return answer

    def local_kernel(environment, left):
        log_weights = []
        for right in spins:
            relations = (
                left[0] * right[0],
                left[1] * right[0],
                left[1] * right[1],
                right[0] * right[1],
            )
            log_weights.append(
                fsum(
                    log(satisfaction if observed == relation else 1.0 - satisfaction)
                    for observed, relation in zip(
                        environment, relations, strict=True
                    )
                )
            )
        maximum = max(log_weights)
        weights = tuple(exp(value - maximum) for value in log_weights)
        normalizer = fsum(weights)
        return {
            right: weight / normalizer
            for right, weight in zip(spins, weights, strict=True)
        }

    answer = 0.0
    for environment_sequence in product(environments, repeat=depth):
        environment_mass = 1.0
        for environment in environment_sequence:
            environment_mass *= sign_probability(environment)
        law = {(1, initial_parity): 1.0}
        for environment in environment_sequence:
            next_law = {right: 0.0 for right in spins}
            for left, left_mass in law.items():
                kernel = local_kernel(environment, left)
                for right in spins:
                    next_law[right] += left_mass * kernel[right]
            law = next_law
        correlation = fsum(mass * right[0] for right, mass in law.items())
        answer += environment_mass * correlation * correlation
    return answer


if __name__ == "__main__":
    unittest.main()
