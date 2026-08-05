"""Tests of the ancestral Walsh identities (note 42, section 5.2)."""

from __future__ import annotations

import random
import unittest

from fractions import Fraction

from ancestral_walsh_identities import (
    WalshAncestor,
    ancestral_message,
    ancestral_message_walsh,
    exact_rational_defects,
    linear_coupling_defect,
    paired_product_defect,
    paired_sum_gap,
    paired_sum_value_defect,
)


def random_ancestor(rng: random.Random, beta: float | None = None) -> WalshAncestor:
    """Draw an ancestor whose four rates are strictly positive."""

    while True:
        common = rng.uniform(-4.0, 4.0)
        imbalance_1 = rng.uniform(-4.0, 4.0)
        imbalance_2 = rng.uniform(-4.0, 4.0)
        slack = rng.uniform(0.1, 6.0)
        total = abs(common) + abs(imbalance_1) + abs(imbalance_2) + slack
        candidate = WalshAncestor(
            beta=rng.uniform(0.0, 1.0) if beta is None else beta,
            common=common,
            imbalance_1=imbalance_1,
            imbalance_2=imbalance_2,
            total=total,
        )
        if min(candidate.four_rates().values()) > 1e-9:
            return candidate


class TestPairedIdentities(unittest.TestCase):
    def test_paired_sums_and_products(self) -> None:
        rng = random.Random(7)
        for _ in range(2000):
            ancestor = random_ancestor(rng)
            self.assertAlmostEqual(paired_sum_gap(ancestor), 0.0, places=10)
            self.assertAlmostEqual(paired_sum_value_defect(ancestor), 0.0, places=10)
            self.assertAlmostEqual(paired_product_defect(ancestor), 0.0, places=8)

    def test_exact_rational_certificates(self) -> None:
        """Identities (i)-(ii) hold in exact rational arithmetic."""

        rng = random.Random(29)
        zero = Fraction(0)
        for _ in range(500):
            imb_1 = Fraction(rng.randint(-40, 40), rng.randint(1, 9))
            imb_2 = Fraction(rng.randint(-40, 40), rng.randint(1, 9))
            common = Fraction(rng.randint(-40, 40), rng.randint(1, 9))
            total = abs(common) + abs(imb_1) + abs(imb_2) + Fraction(
                rng.randint(1, 50), rng.randint(1, 9)
            )
            gaps = exact_rational_defects(total, common, imb_1, imb_2)
            self.assertEqual(gaps, (zero, zero, zero))

    def test_coupling_sign_is_opposite_to_imbalance_product(self) -> None:
        rng = random.Random(11)
        for _ in range(500):
            ancestor = random_ancestor(rng)
            rates = ancestor.four_rates()
            product_gap = (
                rates[0, 0] * rates[1, 1] - rates[0, 1] * rates[1, 0]
            )
            imbalance_product = ancestor.imbalance_1 * ancestor.imbalance_2
            if abs(imbalance_product) > 1e-9:
                self.assertLess(product_gap * imbalance_product, 0.0)

    def test_linear_part_never_couples(self) -> None:
        rng = random.Random(13)
        for _ in range(500):
            ancestor = random_ancestor(rng)
            self.assertAlmostEqual(linear_coupling_defect(ancestor), 0.0, places=10)


class TestClosedWalshForm(unittest.TestCase):
    def test_direct_and_walsh_evaluations_agree(self) -> None:
        rng = random.Random(17)
        for _ in range(300):
            chain = [random_ancestor(rng) for _ in range(rng.randint(1, 5))]
            self.assertAlmostEqual(
                ancestral_message(chain),
                ancestral_message_walsh(chain),
                places=8,
            )

    def test_single_terminal_ancestor_is_mute(self) -> None:
        rng = random.Random(19)
        for _ in range(500):
            ancestor = random_ancestor(rng, beta=1.0)
            self.assertAlmostEqual(ancestral_message([ancestor]), 0.0, places=10)

    def test_majority_certificate_yields_nonnegative_message(self) -> None:
        rng = random.Random(23)
        for _ in range(1000):
            chain = []
            for _ in range(rng.randint(1, 4)):
                base = random_ancestor(rng)
                chain.append(
                    WalshAncestor(
                        beta=base.beta,
                        common=base.common,
                        imbalance_1=abs(base.imbalance_1),
                        imbalance_2=abs(base.imbalance_2),
                        total=base.total,
                    )
                )
            self.assertGreaterEqual(ancestral_message(chain), -1e-9)


if __name__ == "__main__":
    unittest.main()
