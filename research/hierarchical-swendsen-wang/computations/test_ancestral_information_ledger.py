from __future__ import annotations

import random
import unittest
from itertools import combinations, product
from math import exp, isinf, log, log2, sqrt

from ancestral_information_ledger import (
    binary_censored_exponential_information_bound_nats,
    binary_filtration_ledger,
    censored_exponential_kl_nats,
    four_rate_censored_information_bound_nats,
    logit_increment_audit,
    satisfied_edge_open_probability,
    triangle_connected_first_merge_information_bits,
    triangle_connection_information_bits,
    triangle_critical_palm_entropy_bits,
    triangle_edge_connected_first_merge_information_bits,
    triangle_edge_connected_first_merge_l2_gain,
    triangle_edge_connection_information_bits,
    triangle_edge_connection_l2_gain,
    triangle_full_connection_probability,
    triangle_palm_dendrogram_entropy_bits,
    triangle_palm_negative_syndrome_probability,
)
from critical_band_thresholds import Q_CRITICAL, beta_critical, coupling
from nishimori_hierarchical_entropy import (
    binary_entropy_bits,
    conjectured_nishimori_root,
    face_noise_probability,
)


def _conditional_information_from_joint(
    joint: dict[tuple[object, int, object], float]
) -> float:
    syndrome_mass: dict[int, float] = {}
    syndrome_output_mass: dict[tuple[int, object], float] = {}
    state_syndrome_mass: dict[tuple[object, int], float] = {}
    for (word, syndrome, output), probability in joint.items():
        syndrome_mass[syndrome] = syndrome_mass.get(syndrome, 0.0) + probability
        key_output = (syndrome, output)
        syndrome_output_mass[key_output] = (
            syndrome_output_mass.get(key_output, 0.0) + probability
        )
        key_state = (word, syndrome)
        state_syndrome_mass[key_state] = (
            state_syndrome_mass.get(key_state, 0.0) + probability
        )
    information = 0.0
    for (word, syndrome, output), probability in joint.items():
        if probability == 0.0:
            continue
        numerator = probability * syndrome_mass[syndrome]
        denominator = (
            state_syndrome_mass[(word, syndrome)]
            * syndrome_output_mass[(syndrome, output)]
        )
        information += probability * log2(numerator / denominator)
    return information


def _direct_triangle_joint(
    p: float, time: float, first_edge: bool, edge_target: bool = False
) -> dict[tuple[object, int, object], float]:
    joint: dict[tuple[object, int, object], float] = {}
    for word in product((-1, 1), repeat=3):
        target: object = word[0] if edge_target else word
        probability = face_noise_probability(word, p)
        syndrome = word[0] * word[1] * word[2]
        positive_edges = tuple(index for index, sign in enumerate(word) if sign == 1)
        count = len(positive_edges)
        connection = triangle_full_connection_probability(count, p, time)
        null_output: object = None if first_edge else 0
        null_key = (target, syndrome, null_output)
        joint[null_key] = joint.get(null_key, 0.0) + probability * (
            1.0 - connection
        )
        if first_edge:
            if count >= 2:
                for edge in positive_edges:
                    key = (target, syndrome, edge)
                    joint[key] = joint.get(key, 0.0) + (
                        probability * connection / count
                    )
        else:
            key = (target, syndrome, 1)
            joint[key] = joint.get(key, 0.0) + probability * connection
    return joint


def _direct_triangle_information(
    p: float, time: float, first_edge: bool, edge_target: bool = False
) -> float:
    return _conditional_information_from_joint(
        _direct_triangle_joint(p, time, first_edge, edge_target)
    )


def _direct_triangle_edge_l2_gain(p: float, time: float, first_edge: bool) -> float:
    joint = _direct_triangle_joint(p, time, first_edge, edge_target=True)
    syndrome_mass: dict[int, float] = {}
    syndrome_signed_mass: dict[int, float] = {}
    output_mass: dict[tuple[int, object], float] = {}
    output_signed_mass: dict[tuple[int, object], float] = {}
    for (target, syndrome, output), probability in joint.items():
        signed_target = int(target)
        syndrome_mass[syndrome] = syndrome_mass.get(syndrome, 0.0) + probability
        syndrome_signed_mass[syndrome] = (
            syndrome_signed_mass.get(syndrome, 0.0)
            + probability * signed_target
        )
        key = (syndrome, output)
        output_mass[key] = output_mass.get(key, 0.0) + probability
        output_signed_mass[key] = (
            output_signed_mass.get(key, 0.0) + probability * signed_target
        )
    baseline = sum(
        syndrome_signed_mass[syndrome] ** 2 / mass
        for syndrome, mass in syndrome_mass.items()
    )
    terminal = sum(
        output_signed_mass[key] ** 2 / mass
        for key, mass in output_mass.items()
        if mass > 0.0
    )
    return terminal - baseline


def _censored_mixture_information_nats(
    prior_plus: float,
    plus_rates: tuple[float, ...],
    plus_weights: tuple[float, ...],
    minus_rates: tuple[float, ...],
    minus_weights: tuple[float, ...],
    horizon: float,
    steps: int = 40000,
) -> float:
    """Independent midpoint integration of a censored clock mixture."""

    spacing = horizon / steps
    information = 0.0
    for index in range(steps):
        time = (index + 0.5) * spacing
        plus_density = sum(
            weight * rate * exp(-rate * time)
            for rate, weight in zip(plus_rates, plus_weights)
        )
        minus_density = sum(
            weight * rate * exp(-rate * time)
            for rate, weight in zip(minus_rates, minus_weights)
        )
        mixture_density = (
            prior_plus * plus_density
            + (1.0 - prior_plus) * minus_density
        )
        if plus_density > 0.0:
            information += (
                spacing
                * prior_plus
                * plus_density
                * log(plus_density / mixture_density)
            )
        if minus_density > 0.0:
            information += (
                spacing
                * (1.0 - prior_plus)
                * minus_density
                * log(minus_density / mixture_density)
            )

    plus_atom = sum(
        weight * exp(-rate * horizon)
        for rate, weight in zip(plus_rates, plus_weights)
    )
    minus_atom = sum(
        weight * exp(-rate * horizon)
        for rate, weight in zip(minus_rates, minus_weights)
    )
    mixture_atom = prior_plus * plus_atom + (1.0 - prior_plus) * minus_atom
    if plus_atom > 0.0:
        information += prior_plus * plus_atom * log(plus_atom / mixture_atom)
    if minus_atom > 0.0:
        information += (
            (1.0 - prior_plus)
            * minus_atom
            * log(minus_atom / mixture_atom)
        )
    return information


class AncestralInformationLedgerTests(unittest.TestCase):
    def test_discrete_ledger_telescopes_and_is_monotone(self) -> None:
        joint: dict[tuple[int, int, int], float] = {}
        for hidden in (-1, 1):
            for first in (-1, 1):
                for second in (-1, 1):
                    probability = 0.5
                    probability *= 0.8 if first == hidden else 0.2
                    probability *= 0.7 if second == hidden else 0.3
                    joint[(hidden, first, second)] = probability
        ledger = binary_filtration_ledger(joint)
        self.assertAlmostEqual(ledger.second_moments[0], 0.0)
        self.assertAlmostEqual(ledger.conditional_entropies_bits[0], 1.0)
        self.assertAlmostEqual(ledger.second_moments[1], 0.6**2)
        self.assertAlmostEqual(
            ledger.conditional_entropies_bits[1], binary_entropy_bits(0.8)
        )
        self.assertTrue(all(gain >= -1e-14 for gain in ledger.l2_gains))
        self.assertTrue(
            all(gain >= -1e-14 for gain in ledger.information_gains_bits)
        )
        self.assertAlmostEqual(
            sum(ledger.l2_gains),
            ledger.second_moments[-1] - ledger.second_moments[0],
        )
        self.assertAlmostEqual(
            sum(ledger.information_gains_bits),
            ledger.conditional_entropies_bits[0]
            - ledger.conditional_entropies_bits[-1],
        )

    def test_logit_increment_bounds_hold(self) -> None:
        generator = random.Random(20260717)
        for _ in range(1000):
            prior = generator.uniform(-8.0, 8.0)
            posterior = generator.uniform(-8.0, 8.0)
            audit = logit_increment_audit(prior, posterior)
            self.assertLessEqual(audit.kl_nats, audit.kl_quadratic_bound + 1e-13)
            self.assertLessEqual(
                audit.squared_magnetization_change,
                audit.magnetization_quadratic_bound + 1e-13,
            )
            self.assertLessEqual(
                audit.squared_magnetization_change,
                audit.magnetization_pinsker_bound + 1e-13,
            )

    def test_censored_exponential_kl_and_binary_bound(self) -> None:
        rate = 1.7
        reference = 0.6
        horizon = 0.8
        steps = 100000
        spacing = horizon / steps
        direct_kl = 0.0
        for index in range(steps):
            time = (index + 0.5) * spacing
            density = rate * exp(-rate * time)
            direct_kl += spacing * density * (
                log(rate / reference) - (rate - reference) * time
            )
        direct_kl += exp(-rate * horizon) * (-(rate - reference) * horizon)
        self.assertAlmostEqual(
            direct_kl,
            censored_exponential_kl_nats(rate, reference, horizon),
            places=10,
        )
        self.assertAlmostEqual(
            censored_exponential_kl_nats(0.0, reference, horizon),
            reference * horizon,
        )
        self.assertTrue(
            isinf(censored_exponential_kl_nats(reference, 0.0, horizon))
        )

        prior = 0.37
        exact_information = _censored_mixture_information_nats(
            prior,
            (rate,),
            (1.0,),
            (reference,),
            (1.0,),
            horizon,
        )
        bound = binary_censored_exponential_information_bound_nats(
            prior, rate, reference, horizon
        )
        self.assertLessEqual(exact_information, bound + 1e-10)

    def test_four_rate_transport_bound_controls_exact_mixture(self) -> None:
        prior = 0.43
        even_rates = (0.7, 2.3)
        even_weights = (0.31, 0.69)
        odd_rates = (1.1, 3.0)
        odd_weights = (0.58, 0.42)
        horizon = 0.65
        exact_information = _censored_mixture_information_nats(
            prior,
            even_rates,
            even_weights,
            odd_rates,
            odd_weights,
            horizon,
        )
        bound = four_rate_censored_information_bound_nats(
            prior,
            even_rates,
            even_weights,
            odd_rates,
            odd_weights,
            horizon,
        )
        self.assertGreaterEqual(bound, 0.0)
        self.assertLessEqual(exact_information, bound + 1e-10)

    def test_triangle_connection_probabilities_by_subset_enumeration(self) -> None:
        for p in (0.6, 0.8, 0.95):
            for time in (0.0, 0.2, 0.7, 1.0):
                opened = satisfied_edge_open_probability(p, time)
                for count in range(4):
                    direct = 0.0
                    edges = tuple(range(count))
                    for size in range(count + 1):
                        for _subset in combinations(edges, size):
                            probability = opened**size * (1.0 - opened) ** (count - size)
                            if size >= 2:
                                direct += probability
                    self.assertAlmostEqual(
                        direct,
                        triangle_full_connection_probability(count, p, time),
                    )

    def test_information_formulas_match_eight_state_enumeration(self) -> None:
        for p in (0.55, 0.7, 0.835805792367, 0.93):
            for time in (0.0, 0.15, 0.5, 1.0):
                self.assertAlmostEqual(
                    triangle_connection_information_bits(p, time),
                    _direct_triangle_information(p, time, first_edge=False),
                    places=12,
                )
                self.assertAlmostEqual(
                    triangle_connected_first_merge_information_bits(p, time),
                    _direct_triangle_information(p, time, first_edge=True),
                    places=12,
                )
                self.assertAlmostEqual(
                    triangle_edge_connection_information_bits(p, time),
                    _direct_triangle_information(
                        p, time, first_edge=False, edge_target=True
                    ),
                    places=12,
                )
                self.assertAlmostEqual(
                    triangle_edge_connected_first_merge_information_bits(
                        p, time
                    ),
                    _direct_triangle_information(
                        p, time, first_edge=True, edge_target=True
                    ),
                    places=12,
                )
                self.assertAlmostEqual(
                    triangle_edge_connection_l2_gain(p, time),
                    _direct_triangle_edge_l2_gain(
                        p, time, first_edge=False
                    ),
                    places=12,
                )
                self.assertAlmostEqual(
                    triangle_edge_connected_first_merge_l2_gain(p, time),
                    _direct_triangle_edge_l2_gain(
                        p, time, first_edge=True
                    ),
                    places=12,
                )

    def test_palm_formula_matches_direct_clock_densities(self) -> None:
        generator = random.Random(314159)
        for _ in range(100):
            p = generator.uniform(0.51, 0.97)
            time = generator.uniform(0.02, 1.0)
            first_time = generator.uniform(0.0, time)
            rate = coupling(p)
            plus_density = (
                p**3
                * 2.0
                * rate**2
                * exp(-rate * first_time - 2.0 * rate * time)
            )
            minus_density = (
                2.0
                * p**2
                * (1.0 - p)
                * rate**2
                * exp(-rate * (first_time + time))
            )
            direct = minus_density / (plus_density + minus_density)
            self.assertAlmostEqual(
                direct, triangle_palm_negative_syndrome_probability(p, time)
            )
            self.assertAlmostEqual(
                direct, triangle_palm_dendrogram_entropy_bits(p, time)
            )

    def test_critical_palm_simplification_and_nishimori_audit(self) -> None:
        p = conjectured_nishimori_root()
        critical_time = beta_critical(p)
        direct = triangle_palm_dendrogram_entropy_bits(p, critical_time)
        closed = (1.0 - p) / (1.0 - Q_CRITICAL)
        self.assertAlmostEqual(direct, closed)
        self.assertAlmostEqual(triangle_critical_palm_entropy_bits(p), closed)
        self.assertAlmostEqual(closed, 0.251560120699, places=11)
        self.assertAlmostEqual(
            triangle_connection_information_bits(p, critical_time),
            0.043883918779,
            places=11,
        )
        self.assertAlmostEqual(
            triangle_connected_first_merge_information_bits(p, critical_time),
            0.078638140273,
            places=11,
        )
        self.assertAlmostEqual(
            triangle_edge_connection_information_bits(p, critical_time),
            0.027809400607,
            places=11,
        )
        self.assertAlmostEqual(
            triangle_edge_connected_first_merge_information_bits(
                p, critical_time
            ),
            0.042759377412,
            places=11,
        )
        self.assertAlmostEqual(
            triangle_edge_connection_l2_gain(p, critical_time),
            0.006320258880,
            places=11,
        )
        self.assertAlmostEqual(
            triangle_edge_connected_first_merge_l2_gain(p, critical_time),
            0.019523088673,
            places=11,
        )

        information_threshold = (1.0 + sqrt(Q_CRITICAL)) / 2.0
        information_critical_time = beta_critical(information_threshold)
        self.assertAlmostEqual(
            triangle_edge_connected_first_merge_information_bits(
                information_threshold, information_critical_time
            ),
            0.061716309354,
            places=11,
        )
        self.assertAlmostEqual(
            triangle_edge_connected_first_merge_l2_gain(
                information_threshold, information_critical_time
            ),
            0.031834915394,
            places=11,
        )


if __name__ == "__main__":
    unittest.main()
