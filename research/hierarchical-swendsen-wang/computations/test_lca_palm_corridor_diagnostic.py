from __future__ import annotations

import random
import unittest
from dataclasses import asdict

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    triangular_torus_edges,
)
from lca_palm_corridor_diagnostic import (
    campbell_weight_counter_audit,
    charge_identity_counter_audit,
    critical_cut_intensity_observations,
    cut_intensity_palm_weight,
    final_realized_event_observations,
    geometric_charge,
    lca_pair_partition_counts,
    realized_event_palm_weight,
    residual_margin_from_rank,
    run_diagnostic,
)


def _two_merger_ranked_edges(
    *, side_length: int, second_rank: float
) -> tuple[tuple[float, int, int], ...]:
    """A deterministic hierarchy with bucket sizes one and two."""

    first = 0
    second = 1
    common_neighbour = side_length
    ranked_edges = []
    for edge in triangular_torus_edges(side_length):
        if edge == tuple(sorted((first, second))):
            rank = 0.01
        elif edge == tuple(sorted((first, common_neighbour))):
            rank = second_rank
        else:
            rank = 0.99
        ranked_edges.append((rank, *edge))
    return tuple(ranked_edges)


class LcaPalmCorridorDiagnosticTests(unittest.TestCase):
    def test_campbell_counter_audit_matches_m_n_and_detects_double_m(
        self,
    ) -> None:
        audit = campbell_weight_counter_audit(
            cut_sizes=(2, 5), distant_pair_counts=(3, 7)
        )
        self.assertEqual(audit.cut_intensity_contributions, (6, 35))
        self.assertEqual(audit.event_rate_factors, (2, 5))
        self.assertEqual(audit.realized_event_weights, (3, 7))
        self.assertEqual(audit.event_expected_contributions, (6, 35))
        self.assertEqual(
            audit.incorrectly_double_biased_contributions, (12, 175)
        )
        self.assertTrue(audit.representations_agree)
        self.assertTrue(audit.double_bias_detected)

        self.assertEqual(cut_intensity_palm_weight(5, 7), 35)
        self.assertEqual(realized_event_palm_weight(7), 7)

    def test_charge_rank_formula_has_an_independent_beta_counter_audit(
        self,
    ) -> None:
        p = 0.805
        final_rank = 2.0 * p - 1.0
        ranks = (0.0, 0.17, Q_CRITICAL, 0.5, final_rank)
        audit = charge_identity_counter_audit(p, ranks)
        self.assertTrue(audit.passed)
        self.assertLessEqual(audit.maximum_margin_error, 2e-12)
        for rank in ranks:
            direct = (final_rank - rank) / (1.0 - rank)
            self.assertAlmostEqual(residual_margin_from_rank(p, rank), direct)
            self.assertAlmostEqual(
                geometric_charge(6, p, rank), 6.0 * direct * direct
            )

    def test_fixed_level_oracle_uses_explicit_m_n_weight_and_disjoint_buckets(
        self,
    ) -> None:
        side_length = 4
        forest = CriticalKruskalForest(
            side_length,
            _two_merger_ranked_edges(
                side_length=side_length, second_rank=0.5
            ),
            critical_rank=Q_CRITICAL,
        )
        observations = critical_cut_intensity_observations(
            forest,
            repetition=0,
            p=0.805,
            distance_fraction=0.1,
            rng=random.Random(11),
        )
        self.assertTrue(observations)
        for observation in observations:
            self.assertEqual(
                observation.palm_weight,
                observation.lca_bucket_size * observation.n_rho,
            )
            self.assertEqual(
                observation.palm_weight_definition,
                "m(A,B) * N_rho(A,B)",
            )
            self.assertTrue(observation.physical_buckets_pairwise_disjoint)
            self.assertFalse(observation.screening_computed)
            self.assertIsNone(observation.screened_cut_count)
            self.assertGreaterEqual(observation.distance, 0.1 * side_length)
            lca_cuts = [cut for cut in observation.cuts if cut.is_lca]
            self.assertEqual(len(lca_cuts), 1)
            self.assertAlmostEqual(lca_cuts[0].actual_rank, Q_CRITICAL)

        # The component {0,1} and the singleton {4} share two physical edges
        # and have four eligible ordered endpoint pairs.
        self.assertTrue(
            any(
                item.lca_bucket_size == 2
                and item.n_rho == 4
                and item.palm_weight == 8
                for item in observations
            )
        )

    def test_realized_event_uses_n_only_and_criticalizes_without_reskeletonizing(
        self,
    ) -> None:
        side_length = 4
        p = 0.805
        second_rank = 0.5
        forest = CriticalKruskalForest(
            side_length,
            _two_merger_ranked_edges(
                side_length=side_length, second_rank=second_rank
            ),
            critical_rank=2.0 * p - 1.0,
        )
        observations = final_realized_event_observations(
            forest,
            repetition=0,
            p=p,
            distance_fraction=0.1,
            rng=random.Random(19),
        )
        target = next(
            item
            for item in observations
            if item.lca_bucket_size == 2 and item.n_rho == 4
        )
        self.assertEqual(target.palm_weight, 4)
        self.assertNotEqual(
            target.palm_weight, target.lca_bucket_size * target.n_rho
        )
        self.assertIn("already carries m", target.palm_weight_definition)
        self.assertTrue(target.physical_buckets_pairwise_disjoint)

        lca = next(cut for cut in target.cuts if cut.is_lca)
        self.assertAlmostEqual(lca.actual_rank, second_rank)
        self.assertAlmostEqual(lca.favourable_rank, Q_CRITICAL)
        self.assertEqual(lca.bucket_size, 2)
        self.assertGreater(lca.favourable_charge, lca.actual_charge)

        path_nodes = forest.pair_path_nodes(target.first, target.second)
        self.assertEqual(tuple(cut.node for cut in target.cuts), path_nodes)
        for cut in target.cuts:
            self.assertEqual(cut.bucket_size, forest.bucket_size[cut.node])
            self.assertAlmostEqual(
                cut.favourable_rank, min(cut.actual_rank, Q_CRITICAL)
            )

        event_pairs, connected_pairs = lca_pair_partition_counts(forest, 0.1)
        self.assertEqual(event_pairs, 6)
        self.assertEqual(connected_pairs, 6)

    def test_geometric_counts_never_claim_missing_screening(self) -> None:
        result = run_diagnostic(
            side_length=6,
            repetitions=2,
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=8,
            maximum_charge=2.0,
            seed=123,
        )
        for summary in (
            result.critical_snapshot_benchmark,
            result.final_corridor,
        ):
            self.assertGreater(summary.observation_count, 0)
            self.assertGreater(summary.total_palm_weight, 0)
            self.assertEqual(summary.repetition_count, 2)
            self.assertGreater(summary.effective_repetition_count, 0.0)
            self.assertIsNotNone(
                summary.jackknife_standard_error_corridor_cut_count
            )
            self.assertGreaterEqual(
                summary.weighted_mean_disjoint_bucket_size_two_count, 0.0
            )
            self.assertIsNotNone(
                summary.jackknife_standard_error_disjoint_bucket_size_two_count
            )
            self.assertLessEqual(
                summary.weighted_mean_disjoint_bucket_size_two_count,
                summary.weighted_mean_disjoint_bounded_size_count,
            )
            self.assertEqual(summary.weighted_pairwise_disjoint_fraction, 1.0)
            self.assertIn("geometric proxy only", summary.interpretation)
            self.assertIn("rank-environment", summary.uncertainty_note)
            self.assertFalse(summary.screening_computed)
            self.assertIsNone(summary.weighted_mean_screened_cut_count)
        self.assertIn(
            "not a Blackwell domination",
            result.critical_snapshot_benchmark.interpretation,
        )
        self.assertIn(
            "geometry-fixed rank criticalization",
            result.final_corridor.interpretation,
        )
        self.assertTrue(result.lca_pair_partition_audit.passed)
        self.assertEqual(
            result.lca_pair_partition_audit.realized_event_pair_total,
            result.lca_pair_partition_audit.connected_pair_total,
        )
        self.assertEqual(
            result.final_corridor.total_palm_weight,
            result.lca_pair_partition_audit.connected_pair_total,
        )
        self.assertLessEqual(
            result.final_corridor.weighted_mean_disjoint_favourable_bounded_charge_count,
            result.final_corridor.weighted_mean_disjoint_actual_bounded_charge_count,
        )

    def test_generation_and_pair_selection_are_reproducible(self) -> None:
        arguments = dict(
            side_length=6,
            repetitions=2,
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=8,
            maximum_charge=2.0,
            seed=9182,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first), asdict(second))

    def test_invalid_weight_and_proxy_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            cut_intensity_palm_weight(0, 3)
        with self.assertRaises(ValueError):
            realized_event_palm_weight(-1)
        with self.assertRaises(ValueError):
            residual_margin_from_rank(0.805, 0.7)
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=6,
                repetitions=1,
                p=1.1,
                distance_fraction=0.2,
                maximum_bucket_size=8,
                maximum_charge=1.0,
                seed=1,
            )
        postcritical_forest = CriticalKruskalForest(
            4,
            _two_merger_ranked_edges(side_length=4, second_rank=0.5),
            critical_rank=0.5,
        )
        with self.assertRaises(ValueError):
            critical_cut_intensity_observations(
                postcritical_forest,
                repetition=0,
                p=0.805,
                distance_fraction=0.1,
                rng=random.Random(1),
            )


if __name__ == "__main__":
    unittest.main()
