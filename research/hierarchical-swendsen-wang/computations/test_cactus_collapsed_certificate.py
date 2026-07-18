from __future__ import annotations

import unittest
from fractions import Fraction

from cactus_collapsed_certificate import (
    bsc_kernel,
    connected_cactus_conformity_probability,
    connected_cactus_second_moment,
    connected_triangle_decomposition,
    connected_triangle_reliability,
    connected_triangle_reliability_by_quadrature,
    connected_triangle_replicated_kernel,
    critical_to_late_merger_flux_garbling,
    critical_to_late_replicated_garbling,
    direct_first_connection_density,
    direct_first_connection_mass,
    fixed_path_first_cactus_direct_second_moment,
    fixed_path_first_cactus_transfer_second_moment,
    lca_rank_cactus_conformity_probability,
    lca_rank_cactus_second_moment,
    merger_flux_triangle_reliability,
    multiply_two_by_two,
    p_eight_critical_cactus_interval,
    p_eight_critical_lca_interval,
    path_first_connection_mass,
    path_first_connection_density,
    residual_satisfaction_probability,
    screened_connected_triangle_reliability,
    triangle_connection_probability,
    triangle_connection_density,
    unconditional_cactus_second_moment,
)
from critical_band_thresholds import Q_CRITICAL


class CactusCollapsedCertificateTests(unittest.TestCase):
    def test_history_masses_partition_local_connection(self) -> None:
        for rank in (0.05, 0.2, Q_CRITICAL, 0.55):
            self.assertAlmostEqual(
                direct_first_connection_mass(rank)
                + path_first_connection_mass(rank),
                triangle_connection_probability(rank),
            )
            self.assertAlmostEqual(
                direct_first_connection_density(rank)
                + path_first_connection_density(rank),
                triangle_connection_density(rank),
            )

    def test_closed_form_matches_history_integral(self) -> None:
        for p, rank in ((0.7, 0.2), (0.8, Q_CRITICAL), (0.9, 0.6)):
            closed = connected_triangle_reliability(p, rank)
            decomposed = connected_triangle_decomposition(
                p, rank
            ).replicated_reliability
            integrated = connected_triangle_reliability_by_quadrature(
                p, rank, subdivisions=100
            )
            self.assertAlmostEqual(closed, decomposed)
            self.assertAlmostEqual(closed, integrated)

    def test_fixed_global_enumeration_matches_local_transfer(self) -> None:
        p = 0.8
        for ranks in ((0.2,), (0.18, 0.27), (0.12, 0.24, 0.34)):
            direct = fixed_path_first_cactus_direct_second_moment(p, ranks)
            transfer = fixed_path_first_cactus_transfer_second_moment(p, ranks)
            self.assertAlmostEqual(direct, transfer)

    def test_satellite_levels_cancel_from_endpoint_transfer(self) -> None:
        p = 0.8
        ranks = (0.2, 0.34)
        expected = fixed_path_first_cactus_transfer_second_moment(p, ranks)
        for satellites in ((0.01, 0.02), (0.19, 0.33), (0.05, 0.2)):
            actual = fixed_path_first_cactus_direct_second_moment(
                p, ranks, satellites
            )
            self.assertAlmostEqual(actual, expected)

    def test_replicated_kernel_has_the_reliability_eigenvalue(self) -> None:
        p = 0.8
        kernel = connected_triangle_replicated_kernel(p, Q_CRITICAL)
        reliability = connected_triangle_reliability(p, Q_CRITICAL)
        for row in kernel:
            self.assertAlmostEqual(sum(row), 1.0)
        self.assertAlmostEqual(kernel[0][0] - kernel[0][1], reliability)

    def test_late_connected_block_is_an_explicit_bsc_degradation(self) -> None:
        p = 0.8
        early_rank = Q_CRITICAL
        late_rank = 2.0 * p - 1.0
        early = connected_triangle_replicated_kernel(p, early_rank)
        late = connected_triangle_replicated_kernel(p, late_rank)
        garbling = critical_to_late_replicated_garbling(
            p, early_rank, late_rank
        )
        composed = multiply_two_by_two(early, garbling)
        for actual_row, expected_row in zip(composed, late, strict=True):
            for actual, expected in zip(
                actual_row, expected_row, strict=True
            ):
                self.assertAlmostEqual(actual, expected)

    def test_late_pivotal_block_is_an_explicit_bsc_degradation(self) -> None:
        p = 0.8
        early_rank = Q_CRITICAL
        late_rank = 2.0 * p - 1.0
        early_reliability = merger_flux_triangle_reliability(p, early_rank)
        late_reliability = merger_flux_triangle_reliability(p, late_rank)
        early = bsc_kernel(early_reliability)
        late = bsc_kernel(late_reliability)
        garbling = critical_to_late_merger_flux_garbling(
            p, early_rank, late_rank
        )
        composed = multiply_two_by_two(early, garbling)
        for actual_row, expected_row in zip(composed, late, strict=True):
            for actual, expected in zip(
                actual_row, expected_row, strict=True
            ):
                self.assertAlmostEqual(actual, expected)

    def test_critical_connection_is_the_favorable_postcritical_case(self) -> None:
        p = 0.8
        critical = connected_triangle_reliability(p, Q_CRITICAL)
        critical_flux = merger_flux_triangle_reliability(p, Q_CRITICAL)
        for rank in (0.4, 0.5, 2.0 * p - 1.0):
            self.assertGreater(critical, connected_triangle_reliability(p, rank))
            self.assertGreater(
                critical_flux, merger_flux_triangle_reliability(p, rank)
            )

    def test_screened_boundary_coefficient_is_strict_below_one(self) -> None:
        p = 0.8
        neutral = connected_triangle_reliability(p, Q_CRITICAL)
        self.assertAlmostEqual(
            screened_connected_triangle_reliability(p, Q_CRITICAL, 0.0),
            neutral,
        )
        self.assertLess(
            screened_connected_triangle_reliability(p, Q_CRITICAL, 2.0),
            1.0,
        )

    def test_p_eight_reference_values(self) -> None:
        p = 0.8
        coefficient = connected_triangle_reliability(p, Q_CRITICAL)
        self.assertAlmostEqual(coefficient, 0.8867525668568206)
        self.assertAlmostEqual(
            connected_cactus_second_moment(p, Q_CRITICAL, 20),
            0.09037516135894078,
        )
        self.assertAlmostEqual(
            connected_cactus_conformity_probability(p, Q_CRITICAL, 40),
            0.5040838348953273,
        )
        self.assertAlmostEqual(
            merger_flux_triangle_reliability(p, Q_CRITICAL),
            0.7915307368661747,
        )
        self.assertAlmostEqual(
            lca_rank_cactus_second_moment(p, Q_CRITICAL, 20),
            0.08067043811150573,
        )
        self.assertAlmostEqual(
            lca_rank_cactus_conformity_probability(p, Q_CRITICAL, 40),
            0.5036453019306119,
        )

    def test_lca_palm_formula_matches_cumulative_flux_ratio(self) -> None:
        p = 0.8
        rank = Q_CRITICAL
        step = 1e-6
        for blocks in (1, 2, 3):
            def informative_cumulative(value: float) -> float:
                decomposition = connected_triangle_decomposition(p, value)
                informative = (
                    decomposition.connection_mass
                    * decomposition.replicated_reliability
                )
                return informative**blocks

            def connection_cumulative(value: float) -> float:
                return triangle_connection_probability(value) ** blocks

            informative_density = (
                informative_cumulative(rank + step)
                - informative_cumulative(rank - step)
            ) / (2.0 * step)
            connection_density = (
                connection_cumulative(rank + step)
                - connection_cumulative(rank - step)
            ) / (2.0 * step)
            self.assertAlmostEqual(
                informative_density / connection_density,
                lca_rank_cactus_second_moment(p, rank, blocks),
                places=8,
            )

    def test_rational_interval_certifies_p_eight_values(self) -> None:
        lower, upper = p_eight_critical_cactus_interval(3)
        actual = connected_cactus_second_moment(0.8, Q_CRITICAL, 3)
        self.assertLess(float(lower), actual)
        self.assertLess(actual, float(upper))
        self.assertIsInstance(lower, Fraction)
        lca_lower, lca_upper = p_eight_critical_lca_interval(20)
        lca_actual = lca_rank_cactus_second_moment(0.8, Q_CRITICAL, 20)
        self.assertLess(float(lca_lower), lca_actual)
        self.assertLess(lca_actual, float(lca_upper))

    def test_unconditional_value_includes_exact_connection_mass(self) -> None:
        p = 0.8
        rank = Q_CRITICAL
        blocks = 3
        conditional = connected_cactus_second_moment(p, rank, blocks)
        connection = triangle_connection_probability(rank) ** blocks
        self.assertAlmostEqual(
            unconditional_cactus_second_moment(p, rank, blocks),
            connection * conditional,
        )

    def test_input_validation(self) -> None:
        with self.assertRaises(ValueError):
            residual_satisfaction_probability(0.5, 0.1)
        with self.assertRaises(ValueError):
            connected_triangle_reliability(0.8, 0.7)
        with self.assertRaises(ValueError):
            connected_triangle_reliability_by_quadrature(0.8, 0.2, 3)
        with self.assertRaises(ValueError):
            fixed_path_first_cactus_direct_second_moment(0.8, (0.3, 0.2))
        with self.assertRaises(ValueError):
            bsc_kernel(1.1)
        with self.assertRaises(ValueError):
            lca_rank_cactus_second_moment(0.8, 0.2, 0)


if __name__ == "__main__":
    unittest.main()
