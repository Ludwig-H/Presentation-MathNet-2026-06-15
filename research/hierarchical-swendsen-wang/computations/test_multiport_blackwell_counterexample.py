from __future__ import annotations

import unittest
from fractions import Fraction

from multiport_blackwell_counterexample import (
    LATE_SATISFACTION,
    PRIOR,
    STATES,
    TARGET,
    candidate_covariant_degradation_noise,
    marginalized_winner_channel,
    marginalized_winner_noise,
    p_805_multiport_blackwell_counterexample,
    posterior_variance,
    symmetrize_postprocessing,
    winner_noise_fourier,
)


class MultiportBlackwellCounterexampleTests(unittest.TestCase):
    def test_marginalized_winner_channel_is_exactly_additive(self) -> None:
        satisfaction = Fraction(17, 30)
        noise = marginalized_winner_noise(satisfaction)
        self.assertEqual(
            noise,
            (
                Fraction(17, 30),
                Fraction(13, 60),
                Fraction(13, 60),
                Fraction(0),
            ),
        )
        channel = marginalized_winner_channel(satisfaction)
        self.assertTrue(all(sum(row) == 1 for row in channel))
        self.assertEqual(channel[0], noise)
        for latent_index, latent in enumerate(STATES):
            for output_index, output in enumerate(STATES):
                relation = output[0] * latent[0], output[1] * latent[1]
                self.assertEqual(
                    channel[latent_index][output_index],
                    noise[STATES.index(relation)],
                )

    def test_noise_fourier_coefficients(self) -> None:
        satisfaction = Fraction(17, 30)
        coefficients = winner_noise_fourier(satisfaction)
        self.assertEqual(
            coefficients,
            (
                Fraction(1),
                satisfaction,
                satisfaction,
                Fraction(2, 15),
            ),
        )

    def test_symmetrization_produces_a_covariant_stochastic_kernel(self) -> None:
        kernel = (
            (Fraction(1, 2), Fraction(1, 2), Fraction(0), Fraction(0)),
            (Fraction(1, 3), Fraction(0), Fraction(2, 3), Fraction(0)),
            (Fraction(0), Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)),
            (Fraction(1, 5), Fraction(1, 5), Fraction(1, 5), Fraction(2, 5)),
        )
        symmetric = symmetrize_postprocessing(kernel)
        self.assertTrue(all(sum(row) == 1 for row in symmetric))
        self.assertTrue(all(entry >= 0 for row in symmetric for entry in row))
        index = {state: position for position, state in enumerate(STATES)}
        for shift in STATES:
            for source in STATES:
                for target in STATES:
                    shifted_source = (
                        shift[0] * source[0],
                        shift[1] * source[1],
                    )
                    shifted_target = (
                        shift[0] * target[0],
                        shift[1] * target[1],
                    )
                    self.assertEqual(
                        symmetric[index[source]][index[target]],
                        symmetric[index[shifted_source]][index[shifted_target]],
                    )

    def test_unique_covariant_candidate_has_negative_mass(self) -> None:
        candidate = candidate_covariant_degradation_noise(
            Fraction(7, 10), Fraction(17, 30)
        )
        self.assertEqual(sum(candidate), 1)
        self.assertLess(candidate[STATES.index((-1, -1))], 0)
        dominant_margin = Fraction(2, 5)
        degraded_margin = Fraction(2, 15)
        expected = (
            (1 - dominant_margin)
            * (degraded_margin - dominant_margin)
            / (4 * dominant_margin * (1 + dominant_margin))
        )
        self.assertEqual(candidate[STATES.index((-1, -1))], expected)

    def test_exact_prior_and_target_moments(self) -> None:
        self.assertEqual(sum(PRIOR), 1)
        mean = sum(mass * value for mass, value in zip(PRIOR, TARGET, strict=True))
        variance = sum(
            mass * (value - mean) ** 2
            for mass, value in zip(PRIOR, TARGET, strict=True)
        )
        self.assertEqual(mean, Fraction(91, 100))
        self.assertEqual(variance, Fraction(1719, 10000))

    def test_late_variance_is_exact_and_reproducible(self) -> None:
        value = posterior_variance(
            marginalized_winner_channel(LATE_SATISFACTION), PRIOR, TARGET
        )
        self.assertEqual(
            value,
            Fraction(5961475663560973, 207013559198670000),
        )
        self.assertAlmostEqual(float(value), 0.028797513006574463)

    def test_actual_p805_interval_certificate(self) -> None:
        certificate = p_805_multiport_blackwell_counterexample()
        self.assertTrue(certificate.certifies_no_blackwell_degradation)
        self.assertLess(certificate.candidate_noise_minus_minus[1], 0)
        self.assertTrue(certificate.certifies_posterior_variance_reversal)
        self.assertGreater(certificate.variance_reversal_gap[0], 0)
        self.assertLess(
            certificate.critical_posterior_second_moment[1],
            certificate.late_posterior_second_moment,
        )
        self.assertAlmostEqual(
            float(certificate.critical_posterior_variance[0]),
            0.022535603548554394,
        )
        self.assertAlmostEqual(
            float(certificate.late_posterior_second_moment),
            0.8568975130065746,
        )


if __name__ == "__main__":
    unittest.main()
