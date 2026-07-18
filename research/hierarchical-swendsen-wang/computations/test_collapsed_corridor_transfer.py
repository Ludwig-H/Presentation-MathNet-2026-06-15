from __future__ import annotations

import unittest
from fractions import Fraction
from math import fsum

from collapsed_corridor_transfer import (
    binary_experiment_reliability,
    collapsed_corridor_second_moment,
    factorized_corridor_reliability,
    ising_chain_prior,
    screened_two_edge_corridor_bound,
    uniform_spin_prior,
)
from favorable_time_comparison import (
    bucket_binary_experiment,
    critical_beta,
)


class CollapsedCorridorTransferTests(unittest.TestCase):
    def test_two_edge_experiment_is_an_erasure_channel(self) -> None:
        p = 0.8
        level = critical_beta(p)
        plus, minus = bucket_binary_experiment(p, 2, level)
        reliability = binary_experiment_reliability(plus, minus)
        self.assertAlmostEqual(reliability, plus[2])
        self.assertAlmostEqual(reliability, 0.6935822227524088)

    def test_factorized_formula_matches_full_corridor_enumeration(self) -> None:
        p = 0.8
        levels = (critical_beta(p), 0.6, 0.9)
        sizes = (2, 3, 4)
        enumerated = collapsed_corridor_second_moment(p, sizes, levels)
        factorized = factorized_corridor_reliability(p, sizes, levels)
        self.assertAlmostEqual(enumerated, factorized)

    def test_blackwell_tensorization_survives_a_correlated_prior(self) -> None:
        p = 0.8
        critical = critical_beta(p)
        sizes = (2, 3, 2, 4)
        early = (critical,) * len(sizes)
        late = (0.55, 0.7, 0.85, 1.0)
        for interaction, field in ((-0.8, 0.0), (0.6, 0.0), (0.4, 0.2)):
            prior = ising_chain_prior(len(sizes), interaction, field)
            self.assertGreaterEqual(
                collapsed_corridor_second_moment(p, sizes, early, prior),
                collapsed_corridor_second_moment(p, sizes, late, prior)
                - 2e-15,
            )

    def test_blackwell_tensorization_applies_to_an_arbitrary_target(self) -> None:
        p = 0.8
        critical = critical_beta(p)
        sizes = (2, 3, 2)
        prior = ising_chain_prior(3, interaction=-0.5, field=0.15)
        target = {
            state: 0.3 * state[0] - 0.2 * state[1] * state[2]
            for state in prior
        }
        early = collapsed_corridor_second_moment(
            p, sizes, (critical,) * 3, prior, target
        )
        late = collapsed_corridor_second_moment(
            p, sizes, (0.55, 0.75, 1.0), prior, target
        )
        self.assertGreaterEqual(early, late - 2e-15)

    def test_collapsed_projection_is_below_a_systematic_sweep(self) -> None:
        stationary = (0.04, 0.06, 0.1, 0.2, 0.08, 0.12, 0.16, 0.24)

        def projection(bits: tuple[int, ...]) -> tuple[tuple[float, ...], ...]:
            rows = []
            masks = tuple(sum(1 << bit for bit in subset) for subset in _subsets(bits))
            for state in range(8):
                orbit = tuple(sorted({state ^ mask for mask in masks}))
                normalizer = fsum(stationary[item] for item in orbit)
                row = tuple(
                    stationary[item] / normalizer if item in orbit else 0.0
                    for item in range(8)
                )
                rows.append(row)
            return tuple(rows)

        def apply(matrix, values):
            return tuple(
                fsum(matrix[i][j] * values[j] for j in range(8))
                for i in range(8)
            )

        first = projection((0,))
        second = projection((1,))
        collapsed = projection((0, 1))
        values = (1.0, -0.4, 0.2, -1.0, -0.7, 0.3, -0.1, 0.8)
        lca_only = apply(first, values)
        sweep = apply(first, apply(second, values))
        block = apply(collapsed, values)
        sweep_norm = fsum(
            mass * value * value for mass, value in zip(stationary, sweep)
        )
        block_norm = fsum(
            mass * value * value for mass, value in zip(stationary, block)
        )
        lca_norm = fsum(
            mass * value * value
            for mass, value in zip(stationary, lca_only)
        )
        depth_gap_norm = fsum(
            mass * (lca - collapsed_value) ** 2
            for mass, lca, collapsed_value in zip(
                stationary, lca_only, block, strict=True
            )
        )
        original_mean = fsum(
            mass * value for mass, value in zip(stationary, values)
        )
        block_mean = fsum(
            mass * value for mass, value in zip(stationary, block)
        )
        self.assertAlmostEqual(original_mean, block_mean)
        self.assertAlmostEqual(lca_norm, block_norm + depth_gap_norm)
        self.assertLessEqual(original_mean * original_mean, block_norm)
        self.assertLessEqual(block_norm, sweep_norm + 1e-15)

    def test_top_down_insertion_has_no_abstract_projection_order(self) -> None:
        values = (
            Fraction(1),
            Fraction(-1),
            Fraction(0),
            Fraction(0),
        )

        def conditional_expectation(partition, vector):
            answer = [Fraction(0)] * len(vector)
            for cell in partition:
                average = sum((vector[index] for index in cell), Fraction(0))
                average /= len(cell)
                for index in cell:
                    answer[index] = average
            return tuple(answer)

        lca_partition = ((0, 1), (2,), (3,))
        descendant_partition = ((0,), (1, 2), (3,))
        lca_only = conditional_expectation(lca_partition, values)
        descendant_first = conditional_expectation(
            descendant_partition, values
        )
        top_down_operator = conditional_expectation(
            lca_partition, descendant_first
        )
        lca_norm = sum((value * value for value in lca_only), Fraction(0)) / 4
        top_down_norm = (
            sum((value * value for value in top_down_operator), Fraction(0))
            / 4
        )
        self.assertEqual(lca_norm, 0)
        self.assertEqual(
            top_down_operator,
            (
                Fraction(1, 4),
                Fraction(1, 4),
                Fraction(-1, 2),
                Fraction(0),
            ),
        )
        self.assertEqual(top_down_norm, Fraction(3, 32))

    def test_neutral_two_edge_blocks_contract_at_p_eight(self) -> None:
        p = 0.8
        critical = critical_beta(p)
        self.assertAlmostEqual(
            screened_two_edge_corridor_bound(p, critical, 0.0, 10),
            0.025761997386046745,
        )
        self.assertLess(
            screened_two_edge_corridor_bound(p, critical, 0.0, 40),
            5e-7,
        )

    def test_invalid_prior_and_target_are_rejected(self) -> None:
        p = 0.8
        level = critical_beta(p)
        prior = uniform_spin_prior(2)
        prior.pop((-1, -1))
        with self.assertRaises(ValueError):
            collapsed_corridor_second_moment(p, (2, 2), (level, level), prior)


def _subsets(bits: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    result = [()]
    for bit in bits:
        result += [subset + (bit,) for subset in result]
    return tuple(result)


if __name__ == "__main__":
    unittest.main()
