"""Tests of the exact certificates for note SBM/08."""

from __future__ import annotations

import random
import unittest
from math import log

from fractions import Fraction

from sbm_tree_threshold_proofs import (
    frozen_probability_cross_class,
    frozen_probability_same_class,
    pair_affinity,
    pair_affinity_from_channel,
    psi,
    psi_series,
    residual_channel_identity,
    tree_reconstruction_quantities,
)


class TestPsiLemmas(unittest.TestCase):
    def test_series_matches_closed_form(self) -> None:
        rng = random.Random(3)
        for _ in range(500):
            m = rng.uniform(-0.95, 0.95)
            self.assertAlmostEqual(psi_series(m), psi(m), places=11)

    def test_sdpi_and_quadratic_minoration(self) -> None:
        rng = random.Random(5)
        for _ in range(5000):
            m = rng.uniform(-0.999, 0.999)
            theta = rng.uniform(0.0, 1.0)
            self.assertLessEqual(psi(theta * m), theta * theta * psi(m) + 1e-12)
            self.assertGreaterEqual(psi(m), m * m / 2.0 - 1e-12)


class TestTreeTheorems(unittest.TestCase):
    def test_information_bound_theorem_I18(self) -> None:
        for arity, depth, theta in [(2, 2, 0.4), (2, 3, 0.6), (3, 2, 0.7), (2, 2, 0.9)]:
            q_t, info, bound, _ = tree_reconstruction_quantities(depth, arity, theta)
            self.assertLessEqual(info, bound + 1e-12)
            self.assertLessEqual(q_t, 2.0 * info + 1e-12)

    def test_second_moment_bound_theorem_I17(self) -> None:
        # weakened corollary on the regular tree (v_regular <= v_Poisson)
        for arity, depth, theta in [(2, 2, 0.4), (2, 3, 0.6), (3, 2, 0.7), (2, 2, 0.9)]:
            q_t, _, _, ell_t = tree_reconstruction_quantities(depth, arity, theta)
            self.assertGreaterEqual(q_t, ell_t - 1e-12)

    def test_poisson_recursion_exact_rational(self) -> None:
        """The Poisson recursion of I.17 equals its closed form exactly."""

        from sbm_tree_threshold_proofs import poisson_second_moment_defect

        zero = Fraction(0)
        for d, theta in [
            (Fraction(3), Fraction(1, 2)),
            (Fraction(5, 2), Fraction(3, 5)),
            (Fraction(4), Fraction(7, 10)),
            (Fraction(2), Fraction(9, 10)),
        ]:
            for depth in range(1, 13):
                self.assertEqual(
                    poisson_second_moment_defect(d, theta, depth), zero
                )

    def test_subcritical_information_decays(self) -> None:
        # arity 2, theta = 0.4: lambda = 0.32 < 1 — the bound decays with depth.
        bounds = [
            tree_reconstruction_quantities(t, 2, 0.4)[2] for t in (1, 2, 3)
        ]
        self.assertLess(bounds[1], bounds[0])
        self.assertLess(bounds[2], bounds[1])


class TestScalarIdentities(unittest.TestCase):
    def test_frozen_probability_is_a_minus_b_over_n(self) -> None:
        for n, a, b in [(100, 8.0, 2.0), (1000, 30.0, 10.0), (50, 5.0, 1.0)]:
            expected = (a - b) / n
            self.assertAlmostEqual(
                frozen_probability_same_class(n, a, b), expected, places=12
            )
            self.assertAlmostEqual(
                frozen_probability_cross_class(n, a, b), expected, places=12
            )

    def test_residual_channel_identities(self) -> None:
        rng = random.Random(7)
        for _ in range(2000):
            theta = rng.uniform(0.01, 0.99)
            q = rng.uniform(0.0, theta)  # cut level below the SW endpoint
            marg, shared = residual_channel_identity(theta, q)
            self.assertAlmostEqual(marg, 0.0, places=12)
            self.assertAlmostEqual(shared, 0.0, places=12)

    def test_pair_affinity_value(self) -> None:
        for n, a, b in [(200, 10.0, 3.0), (5000, 40.0, 12.0)]:
            self.assertAlmostEqual(
                pair_affinity(n, a, b), pair_affinity_from_channel(n, a, b), places=14
            )


if __name__ == "__main__":
    unittest.main()
