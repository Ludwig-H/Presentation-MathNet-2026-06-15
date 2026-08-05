"""Tests of the beta = 0 stability benchmark (note SBM/07)."""

from __future__ import annotations

import unittest
from math import exp, log, sqrt

from sbm_glauber_stability_benchmark import (
    affinity,
    exp_half_llr,
    expected_sweep_flips,
    flip_probability,
    log_regime_crossing,
)


class TestAffinityIdentity(unittest.TestCase):
    def test_exp_half_llr_equals_affinity(self) -> None:
        """Proposition 3.1: E[exp(-Delta/2)] == rho ** (n-1) exactly."""

        for n, a, b in [(60, 8.0, 2.0), (120, 12.0, 3.0), (200, 20.0, 5.0)]:
            self.assertAlmostEqual(
                exp_half_llr(n, a, b) / affinity(n, a, b), 1.0, places=10
            )

    def test_flip_probability_below_affinity(self) -> None:
        """1 / (1 + e^x) <= e^{-x/2} transfers to the expectations."""

        for n, a, b in [(60, 8.0, 2.0), (150, 15.0, 4.0)]:
            self.assertLessEqual(flip_probability(n, a, b), affinity(n, a, b))

    def test_flip_probability_positive(self) -> None:
        self.assertGreater(flip_probability(80, 10.0, 3.0), 0.0)


class TestExactRecoveryThreshold(unittest.TestCase):
    def test_stable_side_has_vanishing_sweep_flips(self) -> None:
        """(sqrt A - sqrt B)^2 > 2: log(n * flip) / log n stays negative."""

        big_a, big_b = 9.0, 1.0  # (3 - 1)^2 = 4 > 2
        values = [log_regime_crossing(big_a, big_b, n) for n in (200, 400, 800)]
        for value in values:
            self.assertLess(value, 0.0)

    def test_unstable_side_has_growing_sweep_flips(self) -> None:
        """(sqrt A - sqrt B)^2 < 2: expected sweep flips grow with n."""

        big_a, big_b = 4.0, 1.0  # (2 - 1)^2 = 1 < 2
        values = [log_regime_crossing(big_a, big_b, n) for n in (200, 400, 800)]
        self.assertGreater(values[-1], 0.0)

    def test_crossing_approaches_first_order_theory(self) -> None:
        """The empirical exponent moves monotonically toward the predicted
        1 - (sqrt A - sqrt B)^2 / 2 as n grows (slow log corrections)."""

        big_a, big_b = 9.0, 1.0
        predicted = 1.0 - (sqrt(big_a) - sqrt(big_b)) ** 2 / 2.0
        gaps = [
            abs(log_regime_crossing(big_a, big_b, n) - predicted)
            for n in (200, 400, 800, 1600)
        ]
        for earlier, later in zip(gaps, gaps[1:]):
            self.assertLess(later, earlier)
        self.assertLess(gaps[-1], 0.5)


class TestAlmostExactRegime(unittest.TestCase):
    def test_per_site_error_decreases_with_lambda(self) -> None:
        """Larger lambda_n = (a-b)^2 / (2(a+b)) gives a smaller flip rate."""

        n = 300
        errors = []
        for a, b in [(6.0, 3.0), (10.0, 3.0), (16.0, 3.0)]:
            errors.append(flip_probability(n, a, b))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])

    def test_affinity_matches_sparse_exponent(self) -> None:
        """-log affinity ~ (sqrt a - sqrt b)^2 / 2 in the sparse regime."""

        n, a, b = 4000, 12.0, 3.0
        exponent = -log(affinity(n, a, b))
        predicted = (sqrt(a) - sqrt(b)) ** 2 / 2.0
        self.assertLess(abs(exponent / predicted - 1.0), 0.05)


if __name__ == "__main__":
    unittest.main()
