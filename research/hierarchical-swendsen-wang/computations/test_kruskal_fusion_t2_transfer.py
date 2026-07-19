from __future__ import annotations

import unittest
from math import fsum

from kruskal_fusion_t2_transfer import (
    REPLICATED_STATES,
    STATES,
    T2Geometry,
    closed_satisfaction_probability,
    compare_critical_and_late,
    edgewise_cell_channel_probability,
    edgewise_heat_bath,
    edgewise_merger_channel_probability,
    edgewise_t2_transfer,
    environment_counts,
    grouped_cell_channel_probability,
    grouped_heat_bath,
    grouped_merger_count_distribution,
    grouped_t2_transfer,
    p805_geometry,
    p805_scan,
    prior_distribution,
)


class KruskalFusionT2TransferTests(unittest.TestCase):
    def setUp(self) -> None:
        self.p = 0.805
        self.geometry = p805_geometry(attachment_beta=0.8)

    def _edgewise_environments(self, geometry: T2Geometry):
        from itertools import product

        first_size, second_size = geometry.attachment_groups
        return tuple(
            (target, (first, second))
            for target, first, second in product(
                product((-1, 1), repeat=geometry.target_size),
                product((-1, 1), repeat=first_size),
                product((-1, 1), repeat=second_size),
            )
        )

    def test_geometry_contains_target_fusion_and_two_relation_attachment(
        self,
    ) -> None:
        self.assertEqual(self.geometry.target_size, 2)
        self.assertEqual(self.geometry.attachment_groups, (1, 1))
        self.assertLess(self.geometry.beta_target, self.geometry.beta_attachment)
        with self.assertRaises(ValueError):
            T2Geometry(1, (1, 1), 0.3, 0.5)
        with self.assertRaises(ValueError):
            T2Geometry(2, (0, 1), 0.3, 0.5)
        with self.assertRaises(ValueError):
            T2Geometry(2, (1, 1), 0.7, 0.5)

    def test_latent_winner_is_summed_and_never_returned(self) -> None:
        beta = self.geometry.beta_attachment
        satisfaction = closed_satisfaction_probability(self.p, beta)
        plus_plus = edgewise_merger_channel_probability(
            ((1,), (1,)), (1, 1), self.p, beta
        )
        plus_minus = edgewise_merger_channel_probability(
            ((1,), (-1,)), (1, 1), self.p, beta
        )
        impossible = edgewise_merger_channel_probability(
            ((-1,), (-1,)), (1, 1), self.p, beta
        )
        self.assertAlmostEqual(plus_plus, satisfaction)
        self.assertAlmostEqual(plus_minus, (1.0 - satisfaction) / 2.0)
        self.assertEqual(impossible, 0.0)
        self.assertAlmostEqual(plus_plus + 2.0 * plus_minus + impossible, 1.0)

    def test_each_edgewise_channel_row_is_normalized(self) -> None:
        environments = self._edgewise_environments(self.geometry)
        for state in STATES:
            self.assertAlmostEqual(
                fsum(
                    edgewise_cell_channel_probability(
                        environment, state, self.geometry, self.p
                    )
                    for environment in environments
                ),
                1.0,
            )

    def test_grouped_elimination_matches_edgewise_channel_counts(self) -> None:
        environments = self._edgewise_environments(self.geometry)
        for state in STATES:
            collapsed: dict[tuple[tuple[int, ...], tuple[int, ...]], float] = {}
            for environment in environments:
                counts = environment_counts(environment)
                collapsed[counts] = collapsed.get(counts, 0.0) + (
                    edgewise_cell_channel_probability(
                        environment, state, self.geometry, self.p
                    )
                )
            for counts, mass in collapsed.items():
                self.assertAlmostEqual(
                    mass,
                    grouped_cell_channel_probability(
                        counts, state, self.geometry, self.p
                    ),
                )
            self.assertAlmostEqual(fsum(collapsed.values()), 1.0)

    def test_grouped_winner_law_has_exact_uniform_winner_correction(self) -> None:
        groups = (2, 1)
        beta = 0.61
        law = grouped_merger_count_distribution(groups, self.p, beta)
        satisfaction = closed_satisfaction_probability(self.p, beta)
        self.assertAlmostEqual(fsum(law.values()), 1.0)
        means = tuple(
            fsum(counts[group] * mass for counts, mass in law.items())
            for group in range(2)
        )
        total = sum(groups)
        for size, mean in zip(groups, means, strict=True):
            self.assertAlmostEqual(
                mean,
                size * satisfaction + size / total * (1.0 - satisfaction),
            )

    def test_lambda_factor_heat_bath_matches_bayes_posterior_edgewise(
        self,
    ) -> None:
        prior = prior_distribution(2.0, 1.0)
        for environment in self._edgewise_environments(self.geometry):
            likelihoods = tuple(
                edgewise_cell_channel_probability(
                    environment, state, self.geometry, self.p
                )
                for state in STATES
            )
            evidence = fsum(
                mass * likelihood
                for mass, likelihood in zip(prior, likelihoods, strict=True)
            )
            if evidence == 0.0:
                continue
            bayes = tuple(
                mass * likelihood / evidence
                for mass, likelihood in zip(prior, likelihoods, strict=True)
            )
            factor_heat_bath = edgewise_heat_bath(
                environment, self.geometry, self.p, 2.0, 1.0
            )
            for first, second in zip(bayes, factor_heat_bath, strict=True):
                self.assertAlmostEqual(first, second)

    def test_two_heat_bath_implementations_match_on_the_same_states(self) -> None:
        for environment in self._edgewise_environments(self.geometry):
            direct = edgewise_heat_bath(environment, self.geometry, self.p, -1.7, 0.8)
            eliminated = grouped_heat_bath(
                environment_counts(environment),
                self.geometry,
                self.p,
                -1.7,
                0.8,
            )
            for first, second in zip(direct, eliminated, strict=True):
                self.assertAlmostEqual(first, second)

    def test_ancestral_lambda_distinguishes_the_two_child_flips(self) -> None:
        environment = ((1, -1), ((1,), (-1,)))
        posterior = edgewise_heat_bath(environment, self.geometry, self.p, 0.0, 0.0)
        self.assertNotAlmostEqual(
            posterior[STATES.index((1, -1))],
            posterior[STATES.index((-1, 1))],
        )

    def test_edgewise_and_grouped_full_transfers_agree(self) -> None:
        direct = edgewise_t2_transfer(self.geometry, self.p, 2.0, 1.0)
        eliminated = grouped_t2_transfer(self.geometry, self.p, 2.0, 1.0)
        for first_matrix, second_matrix in (
            (direct.mass, eliminated.mass),
            (direct.shared_replicated_law, eliminated.shared_replicated_law),
            (
                direct.full_replicated_heat_bath,
                eliminated.full_replicated_heat_bath,
            ),
        ):
            for first_row, second_row in zip(first_matrix, second_matrix, strict=True):
                for first, second in zip(first_row, second_row, strict=True):
                    self.assertAlmostEqual(first, second, places=13)
        self.assertAlmostEqual(
            direct.chi_square_reliability,
            eliminated.chi_square_reliability,
            places=13,
        )

    def test_mass_is_reversible_and_full_shared_replica_is_stochastic(
        self,
    ) -> None:
        transfer = grouped_t2_transfer(self.geometry, self.p, 2.0, 1.0)
        self.assertLess(transfer.mass_row_sum_error, 2e-14)
        self.assertLess(transfer.detailed_balance_error, 2e-14)
        self.assertLess(transfer.shared_law_mass_error, 2e-14)
        for row in transfer.mass:
            self.assertAlmostEqual(fsum(row), 1.0)
            self.assertGreaterEqual(min(row), 0.0)
        self.assertEqual(len(transfer.full_replicated_heat_bath), 16)
        for row in transfer.full_replicated_heat_bath:
            self.assertEqual(len(row), 16)
            self.assertAlmostEqual(fsum(row), 1.0)
            self.assertGreaterEqual(min(row), 0.0)
        self.assertEqual(len(REPLICATED_STATES), 16)

    def test_shared_environment_is_not_two_independent_posteriors(self) -> None:
        transfer = grouped_t2_transfer(self.geometry, self.p, 0.0, 0.0)
        independent = tuple(
            tuple(first * second for second in transfer.prior)
            for first in transfer.prior
        )
        difference = max(
            abs(first - second)
            for first_row, second_row in zip(
                transfer.shared_replicated_law, independent, strict=True
            )
            for first, second in zip(first_row, second_row, strict=True)
        )
        self.assertGreater(difference, 1e-3)
        self.assertGreater(transfer.chi_square_reliability, 0.0)

    def test_neutral_prior_favors_the_critical_ancestor(self) -> None:
        comparison = compare_critical_and_late(
            outside_field=0.0,
            attachment_coupling=0.0,
            late_beta=0.8,
        )
        self.assertAlmostEqual(comparison.critical_reliability, 0.7496396953178438)
        self.assertAlmostEqual(comparison.late_reliability, 0.7071689378811097)
        self.assertLess(comparison.late_minus_critical, 0.0)
        self.assertFalse(comparison.refutes_uniform_critical_domination)

    def test_polarized_prior_is_an_explicit_criticalization_counterexample(
        self,
    ) -> None:
        grouped = compare_critical_and_late(
            outside_field=4.0,
            attachment_coupling=3.0,
            late_beta=0.8,
            route="grouped",
        )
        edgewise = compare_critical_and_late(
            outside_field=4.0,
            attachment_coupling=3.0,
            late_beta=0.8,
            route="edgewise",
        )
        self.assertAlmostEqual(grouped.critical_reliability, 0.73511220324355)
        self.assertAlmostEqual(grouped.late_reliability, 0.7556375347686635)
        self.assertAlmostEqual(grouped.late_minus_critical, 0.020525331525113488)
        self.assertGreater(grouped.late_minus_critical, 0.02)
        self.assertTrue(grouped.refutes_uniform_critical_domination)
        self.assertAlmostEqual(
            grouped.critical_reliability,
            edgewise.critical_reliability,
        )
        self.assertAlmostEqual(grouped.late_reliability, edgewise.late_reliability)

    def test_scan_is_reproducible_and_contains_both_orderings(self) -> None:
        first = p805_scan(
            fields=(0.0, 2.0, 4.0),
            attachment_couplings=(0.0, 3.0),
            attachment_groups=((1, 1),),
        )
        second = p805_scan(
            fields=(0.0, 2.0, 4.0),
            attachment_couplings=(0.0, 3.0),
            attachment_groups=((1, 1),),
        )
        self.assertEqual(first, second)
        self.assertTrue(any(row.late_minus_critical < 0.0 for row in first))
        self.assertTrue(any(row.late_minus_critical > 0.0 for row in first))
        self.assertIn(
            "no Palm",
            grouped_t2_transfer(self.geometry).scope_label,
        )

    def test_input_validation(self) -> None:
        with self.assertRaises(ValueError):
            closed_satisfaction_probability(0.5, 0.5)
        with self.assertRaises(ValueError):
            closed_satisfaction_probability(self.p, 1.1)
        with self.assertRaises(ValueError):
            prior_distribution(float("inf"), 0.0)
        with self.assertRaises(ValueError):
            grouped_merger_count_distribution((1, 0), self.p, 0.5)
        with self.assertRaises(ValueError):
            compare_critical_and_late(late_beta=0.2)
        with self.assertRaises(ValueError):
            compare_critical_and_late(route="missing")
        with self.assertRaises(ValueError):
            p805_scan(fields=())


if __name__ == "__main__":
    unittest.main()
