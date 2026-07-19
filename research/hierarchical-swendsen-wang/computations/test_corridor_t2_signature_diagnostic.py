from __future__ import annotations

import random
import unittest
from dataclasses import asdict

from critical_pair_path_geometry import (
    CriticalKruskalForest,
    triangular_torus_edges,
)
from corridor_t2_signature_diagnostic import (
    event_palm_corridor_observations,
    reconstruct_node_geometries,
    run_diagnostic,
)
from lca_palm_corridor_diagnostic import lca_pair_partition_counts


def _deterministic_comb_edges(
    side_length: int = 4,
) -> tuple[tuple[float, int, int], ...]:
    """Three mergers with cut sizes 1, 2, 2 and otherwise late ranks."""

    special_ranks = {
        tuple(sorted((0, 1))): 0.01,
        tuple(sorted((0, side_length))): 0.4,
        tuple(sorted((1, side_length + 1))): 0.55,
    }
    ranked_edges = []
    for index, edge in enumerate(triangular_torus_edges(side_length)):
        rank = special_ranks.get(edge, 0.7 + 0.003 * index)
        ranked_edges.append((rank, *edge))
    ranks = [rank for rank, _, _ in ranked_edges]
    if len(set(ranks)) != len(ranks):
        raise AssertionError("the deterministic test ranks must be distinct")
    return tuple(ranked_edges)


class CorridorT2SignatureDiagnosticTests(unittest.TestCase):
    def test_reconstruction_recovers_cut_boundary_and_port_partition(
        self,
    ) -> None:
        ranked_edges = _deterministic_comb_edges()
        forest = CriticalKruskalForest(4, ranked_edges, critical_rank=0.61)
        geometries, audit = reconstruct_node_geometries(
            forest, ranked_edges, censor_rank=0.61
        )
        self.assertTrue(audit.passed)
        self.assertEqual(audit.checked_internal_node_count, 3)
        self.assertEqual(audit.maximum_winning_rank_error, 0.0)
        self.assertEqual(audit.bucket_size_mismatch_count, 0)
        self.assertEqual(audit.external_boundary_identity_mismatch_count, 0)
        self.assertEqual(audit.port_partition_mismatch_count, 0)

        target_node = next(
            node for node in forest.internal_nodes if forest.merge_rank[node] == 0.4
        )
        target = geometries[target_node]
        self.assertEqual(target.bucket_size, 2)
        self.assertEqual(
            sorted((len(target.left_vertices), len(target.right_vertices))),
            [1, 2],
        )
        self.assertEqual(target.external_physical_edge_count, 12)
        self.assertEqual(target.port_count, 9)

    def test_counter_audit_detects_a_corrupted_stored_bucket_size(self) -> None:
        ranked_edges = _deterministic_comb_edges()
        forest = CriticalKruskalForest(4, ranked_edges, critical_rank=0.61)
        target_node = next(
            node for node in forest.internal_nodes if forest.merge_rank[node] == 0.4
        )
        forest.bucket_size[target_node] += 1
        geometries, audit = reconstruct_node_geometries(
            forest, ranked_edges, censor_rank=0.61
        )
        self.assertFalse(audit.passed)
        self.assertEqual(audit.bucket_size_mismatch_count, 1)
        self.assertEqual(geometries[target_node].bucket_size, 2)

    def test_event_palm_uses_n_rho_only_and_orients_comb_attachment(
        self,
    ) -> None:
        ranked_edges = _deterministic_comb_edges()
        forest = CriticalKruskalForest(4, ranked_edges, critical_rank=0.61)
        observations, audit = event_palm_corridor_observations(
            forest,
            ranked_edges,
            repetition=0,
            p=0.805,
            distance_fraction=0.1,
            rng=random.Random(5),
        )
        self.assertTrue(audit.passed)
        self.assertEqual(len(observations), 3)
        for observation in observations:
            self.assertEqual(observation.palm_weight, observation.n_rho)
            self.assertIn("N_rho only", observation.palm_weight_definition)
            self.assertFalse(observation.screening_computed)
            for signature in observation.signatures:
                self.assertTrue(signature.ports_are_geometric_proxy)
                self.assertFalse(signature.screening_computed)

        event_pairs, connected_pairs = lca_pair_partition_counts(forest, 0.1)
        self.assertEqual(sum(item.palm_weight for item in observations), 12)
        self.assertEqual(event_pairs, 12)
        self.assertEqual(connected_pairs, 12)

        final_observation = next(
            item for item in observations if forest.merge_rank[item.lca_node] == 0.55
        )
        strict_arm = next(
            signature
            for signature in final_observation.signatures
            if forest.merge_rank[signature.node] == 0.4
        )
        self.assertFalse(strict_arm.is_lca)
        self.assertEqual(strict_arm.bucket_size, 2)
        self.assertEqual(strict_arm.small_child_size, 1)
        self.assertEqual(strict_arm.spine_child_size, 1)
        self.assertEqual(strict_arm.attachment_child_size, 2)
        self.assertFalse(strict_arm.attachment_is_smaller_child)
        self.assertEqual(strict_arm.port_count, 9)
        self.assertLess(strict_arm.favourable_charge, 1.0)
        self.assertTrue(
            strict_arm.is_t2_proxy_candidate(
                maximum_bucket_size=6,
                maximum_charge=1.0,
                maximum_ports=9,
                maximum_attachment_size=2,
            )
        )
        self.assertFalse(
            strict_arm.is_t2_proxy_candidate(
                maximum_bucket_size=6,
                maximum_charge=1.0,
                maximum_ports=9,
                maximum_attachment_size=1,
            )
        )
        lca_signature = final_observation.signatures[-1]
        self.assertTrue(lca_signature.is_lca)
        self.assertIsNone(lca_signature.spine_child_size)
        self.assertIsNone(lca_signature.attachment_child_size)

    def test_summary_has_environment_jackknife_and_never_claims_screening(
        self,
    ) -> None:
        arguments = dict(
            side_length=6,
            repetitions=3,
            p=0.805,
            distance_fraction=0.2,
            maximum_bucket_size=6,
            maximum_charge=1.0,
            maximum_ports=6,
            maximum_attachment_size=4,
            seed=123,
        )
        first = run_diagnostic(**arguments)
        second = run_diagnostic(**arguments)
        self.assertEqual(asdict(first), asdict(second))
        self.assertTrue(first.reconstruction_audit.passed)
        self.assertTrue(first.palm_partition_audit.passed)
        self.assertEqual(
            first.total_palm_weight,
            first.palm_partition_audit.connected_distant_ordered_pair_total,
        )
        self.assertIsNotNone(first.jackknife_standard_error_t2_proxy_candidate_count)
        self.assertIsNotNone(
            first.jackknife_standard_error_fraction_with_t2_proxy_candidate
        )
        self.assertLessEqual(
            first.weighted_mean_t2_proxy_candidate_count,
            first.weighted_mean_small_attachment_candidate_count,
        )
        self.assertLessEqual(
            first.weighted_mean_t2_proxy_candidate_count,
            first.weighted_mean_low_port_candidate_count,
        )
        self.assertIn("N_rho only", first.palm_weight_convention)
        self.assertIn("proxies only", first.interpretation)
        self.assertIn("rank-environment", first.uncertainty_note)
        self.assertTrue(first.ports_are_geometric_proxy)
        self.assertFalse(first.screening_computed)

    def test_ambiguous_event_ranks_and_invalid_proxy_cutoffs_are_rejected(
        self,
    ) -> None:
        ranked_edges = list(_deterministic_comb_edges())
        ranked_edges[-1] = (
            ranked_edges[-2][0],
            ranked_edges[-1][1],
            ranked_edges[-1][2],
        )
        forest = CriticalKruskalForest(4, ranked_edges, critical_rank=0.61)
        with self.assertRaisesRegex(ValueError, "distinct ranks"):
            reconstruct_node_geometries(forest, ranked_edges, censor_rank=0.61)
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=6,
                repetitions=1,
                maximum_ports=-1,
            )
        with self.assertRaises(ValueError):
            run_diagnostic(
                side_length=6,
                repetitions=1,
                maximum_attachment_size=0,
            )


if __name__ == "__main__":
    unittest.main()
