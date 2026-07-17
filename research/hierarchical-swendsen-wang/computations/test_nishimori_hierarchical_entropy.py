from __future__ import annotations

import random
import unittest
from itertools import product
from math import log, log2

from nishimori_hierarchical_entropy import (
    binary_entropy_bits,
    conditional_face_noise_laws,
    conjectured_nishimori_root,
    entropy_balance_derivative_nats,
    exponential_race_entropy_bits,
    face_noise_probability,
    face_residual_entropy_bits,
    face_syndrome_probability,
    nishimori_entropy_balance_bits,
    nishimori_ohzeki_equation_nats,
    race_winner_probabilities,
)


class NishimoriHierarchicalEntropyTests(unittest.TestCase):
    def test_conditional_laws_have_four_states_and_correct_parity(self) -> None:
        for p in (0.51, 0.7, 0.8358058, 0.95):
            laws = conditional_face_noise_laws(p)
            for syndrome in (-1, 1):
                self.assertEqual(len(laws[syndrome]), 4)
                self.assertAlmostEqual(sum(laws[syndrome].values()), 1.0)
                for word, probability in laws[syndrome].items():
                    self.assertGreater(probability, 0.0)
                    self.assertEqual(word[0] * word[1] * word[2], syndrome)

    def test_entropy_chain_rule_by_independent_enumeration(self) -> None:
        generator = random.Random(20260717)
        for _ in range(100):
            p = generator.uniform(0.02, 0.98)
            joint_entropy = 0.0
            syndrome_masses = {-1: 0.0, 1: 0.0}
            conditional_entropy = 0.0
            words = tuple(product((-1, 1), repeat=3))
            for word in words:
                probability = face_noise_probability(word, p)
                joint_entropy -= probability * log2(probability)
                syndrome_masses[word[0] * word[1] * word[2]] += probability
            for word in words:
                probability = face_noise_probability(word, p)
                syndrome = word[0] * word[1] * word[2]
                conditional = probability / syndrome_masses[syndrome]
                conditional_entropy -= probability * log2(conditional)
            syndrome_entropy = -sum(
                mass * log2(mass) for mass in syndrome_masses.values()
            )
            self.assertAlmostEqual(joint_entropy, 3.0 * binary_entropy_bits(p))
            self.assertAlmostEqual(
                syndrome_masses[1], face_syndrome_probability(p)
            )
            self.assertAlmostEqual(joint_entropy - syndrome_entropy, conditional_entropy)
            self.assertAlmostEqual(
                conditional_entropy, face_residual_entropy_bits(p)
            )

    def test_published_equation_is_exactly_the_entropy_balance(self) -> None:
        generator = random.Random(28032006)
        for _ in range(100):
            p = generator.uniform(0.02, 0.98)
            published = nishimori_ohzeki_equation_nats(p)
            entropy = -log(2.0) * nishimori_entropy_balance_bits(p)
            self.assertAlmostEqual(published, entropy, places=12)

    def test_root_and_strict_uniqueness_on_upper_half(self) -> None:
        root = conjectured_nishimori_root()
        self.assertAlmostEqual(root, 0.8358058, places=7)
        self.assertAlmostEqual(nishimori_entropy_balance_bits(root), 0.0, places=14)
        for index in range(1, 1000):
            p = 0.5 + 0.5 * index / 1000.0
            self.assertLess(entropy_balance_derivative_nats(p), 0.0)

    def test_four_state_race_reproduces_conditional_entropy(self) -> None:
        for p in (0.55, 0.7, 0.8358058, 0.93):
            syndrome_plus = face_syndrome_probability(p)
            laws = conditional_face_noise_laws(p)
            race_entropy = (
                syndrome_plus * exponential_race_entropy_bits(laws[1].values())
                + (1.0 - syndrome_plus)
                * exponential_race_entropy_bits(laws[-1].values())
            )
            self.assertAlmostEqual(race_entropy, face_residual_entropy_bits(p))

    def test_race_winner_law_is_normalized_rates(self) -> None:
        rates = (0.5, 1.5, 3.0, 5.0)
        probabilities = race_winner_probabilities(rates)
        self.assertEqual(probabilities, (0.05, 0.15, 0.3, 0.5))
        direct_entropy = -sum(value * log2(value) for value in probabilities)
        self.assertAlmostEqual(exponential_race_entropy_bits(rates), direct_entropy)


if __name__ == "__main__":
    unittest.main()
