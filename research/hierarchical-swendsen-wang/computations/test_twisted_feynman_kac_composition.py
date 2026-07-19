from __future__ import annotations

import unittest
from fractions import Fraction

from triangular_band_collapsed_certificate import twisted_chain_second_moment
from twisted_feynman_kac_composition import (
    OPEN_STEPS,
    brute_force_path_expansion,
    canonical_fraction_lift_audit,
    canonical_lift_from_mass_and_twist,
    direct_fraction_lift_audit,
    finite_composition_certificate,
    normalize_lifted_transfer,
    triangular_neutral_power_diagnostic,
)


class TwistedFeynmanKacCompositionTests(unittest.TestCase):
    def test_direct_fraction_lift_matches_exhaustive_paths_exactly(self) -> None:
        audit = direct_fraction_lift_audit()
        self.assertTrue(audit.matches_exactly)
        self.assertTrue(audit.certificate.is_exact)
        self.assertTrue(audit.certificate.inequality_holds)
        self.assertTrue(
            all(isinstance(value, Fraction) for value in audit.brute_signed_values)
        )
        self.assertTrue(all(slack >= 0 for slack in audit.certificate.slacks))

    def test_canonical_fraction_lift_is_an_independent_exact_audit(self) -> None:
        audit = canonical_fraction_lift_audit()
        self.assertTrue(audit.matches_exactly)
        self.assertTrue(audit.certificate.is_exact)
        self.assertTrue(audit.certificate.inequality_holds)
        self.assertNotEqual(
            audit.certificate.signed_values,
            direct_fraction_lift_audit().certificate.signed_values,
        )

    def test_zero_mass_has_zero_twist_and_zero_ratio(self) -> None:
        step = normalize_lifted_transfer(
            (
                (
                    (Fraction(0), Fraction(0)),
                    (Fraction(1, 4), Fraction(3, 4)),
                ),
            )
        )
        self.assertEqual(step.mass[0][0], 0)
        self.assertEqual(step.twisted[0][0], 0)
        self.assertEqual(step.ratio[0][0], 0)
        certificate = finite_composition_certificate(
            (step,), (Fraction(7, 9), Fraction(-2, 5))
        )
        brute = brute_force_path_expansion(
            (step,), (Fraction(7, 9), Fraction(-2, 5))
        )
        self.assertEqual(certificate.signed_values, brute[0])
        self.assertEqual(certificate.feynman_kac_envelope, brute[1])

    def test_empty_composition_is_identity_and_absolute_terminal(self) -> None:
        terminal = (Fraction(-2, 3), Fraction(5, 7))
        certificate = finite_composition_certificate((), terminal)
        self.assertEqual(certificate.signed_values, terminal)
        self.assertEqual(
            certificate.feynman_kac_envelope,
            (Fraction(2, 3), Fraction(5, 7)),
        )
        self.assertEqual(
            brute_force_path_expansion((), terminal),
            (terminal, (Fraction(2, 3), Fraction(5, 7))),
        )

    def test_entrywise_domination_is_equivalent_to_a_positive_binary_lift(
        self,
    ) -> None:
        mass = (
            (Fraction(1, 3), Fraction(2, 3)),
            (Fraction(1, 2), Fraction(1, 2)),
        )
        twist = (
            (Fraction(-1, 6), Fraction(1, 3)),
            (Fraction(0), Fraction(-1, 2)),
        )
        step = canonical_lift_from_mass_and_twist(mass, twist)
        self.assertEqual(step.mass, mass)
        self.assertEqual(step.twisted, twist)
        for mass_row, twist_row, ratio_row in zip(
            step.mass, step.twisted, step.ratio, strict=True
        ):
            for mass_value, twist_value, ratio in zip(
                mass_row, twist_row, ratio_row, strict=True
            ):
                self.assertEqual(mass_value * ratio, abs(twist_value))

    def test_integer_inputs_are_promoted_to_exact_rationals(self) -> None:
        direct = normalize_lifted_transfer(((((0, 1)),),))
        canonical = canonical_lift_from_mass_and_twist(((1,),), ((0,),))
        self.assertTrue(direct.is_exact)
        self.assertTrue(canonical.is_exact)
        self.assertEqual(direct.ratio, ((Fraction(1),),))
        self.assertEqual(
            canonical.lifted,
            (((Fraction(1, 2), Fraction(1, 2)),),),
        )

    def test_inhomogeneous_dimensions_are_checked(self) -> None:
        two_to_three = canonical_lift_from_mass_and_twist(
            ((Fraction(1, 3),) * 3, (Fraction(1, 3),) * 3),
            ((Fraction(0),) * 3, (Fraction(0),) * 3),
        )
        two_to_two = canonical_lift_from_mass_and_twist(
            ((Fraction(1, 2),) * 2, (Fraction(1, 2),) * 2),
            ((Fraction(0),) * 2, (Fraction(0),) * 2),
        )
        with self.assertRaises(ValueError):
            finite_composition_certificate(
                (two_to_three, two_to_two), (Fraction(1), Fraction(1))
            )
        with self.assertRaises(ValueError):
            brute_force_path_expansion(
                (two_to_three, two_to_two), (Fraction(1), Fraction(1))
            )

    def test_triangular_neutral_powers_obey_both_composition_bounds(self) -> None:
        depths = (0, 1, 2, 3, 5, 10)
        diagnostic = triangular_neutral_power_diagnostic(depths)
        self.assertTrue(diagnostic.all_composition_bounds_hold)
        self.assertLess(diagnostic.uniform_coefficient, 0.3)
        self.assertIn("not a real Palm-corridor certificate", diagnostic.scope_label)
        self.assertEqual(diagnostic.open_steps, OPEN_STEPS)
        for depth, value, envelope, power_bound in zip(
            depths,
            diagnostic.signed_values,
            diagnostic.feynman_kac_envelopes,
            diagnostic.uniform_power_bounds,
            strict=True,
        ):
            self.assertAlmostEqual(value, twisted_chain_second_moment(0.805, depth))
            self.assertLessEqual(abs(value), envelope + 1e-12)
            self.assertLessEqual(envelope, power_bound + 1e-12)

    def test_strict_validation_rejects_bad_inputs(self) -> None:
        with self.assertRaises(ValueError):
            normalize_lifted_transfer(())
        with self.assertRaises(ValueError):
            normalize_lifted_transfer((((Fraction(-1, 4), Fraction(5, 4)),),))
        with self.assertRaises(ValueError):
            normalize_lifted_transfer((((Fraction(1, 4), Fraction(1, 4)),),))
        with self.assertRaises(ValueError):
            normalize_lifted_transfer((((Fraction(1),),),))
        with self.assertRaises(ValueError):
            normalize_lifted_transfer(((1,),))
        with self.assertRaises(ValueError):
            canonical_lift_from_mass_and_twist(
                ((Fraction(1),),), ((Fraction(5, 4),),)
            )
        with self.assertRaises(ValueError):
            canonical_lift_from_mass_and_twist(
                ((Fraction(1),),), ((Fraction(0), Fraction(0)),)
            )
        with self.assertRaises(ValueError):
            finite_composition_certificate((), ())
        with self.assertRaises(ValueError):
            triangular_neutral_power_diagnostic(())
        with self.assertRaises(ValueError):
            triangular_neutral_power_diagnostic((True,))


if __name__ == "__main__":
    unittest.main()
