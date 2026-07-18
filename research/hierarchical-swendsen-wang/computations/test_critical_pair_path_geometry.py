from __future__ import annotations

import random
import unittest
from math import exp, log

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    clock_time_from_rank,
    evaluate_ranked_geometry,
    ranked_pair_path_geometry,
    sample_ranked_edges,
    triangular_torus_distance,
    triangular_torus_edges,
)


class CriticalPairPathGeometryTests(unittest.TestCase):
    def test_triangular_torus_has_three_edges_per_vertex(self) -> None:
        for side_length in (4, 5, 8):
            self.assertEqual(
                len(triangular_torus_edges(side_length)),
                3 * side_length * side_length,
            )

    def test_torus_distance_matches_the_three_lattice_directions(self) -> None:
        side_length = 8
        origin = 0
        neighbours = (1, side_length, 1 + side_length * (side_length - 1))
        for neighbour in neighbours:
            self.assertEqual(
                triangular_torus_distance(origin, neighbour, side_length), 1
            )
        self.assertEqual(
            triangular_torus_distance(origin, 1 + side_length, side_length), 2
        )

    def test_kruskal_bucket_uses_all_crossing_physical_edges(self) -> None:
        side_length = 4
        first = 0
        second = 1
        common_neighbour = side_length
        ranked_edges = []
        for edge in triangular_torus_edges(side_length):
            if edge == tuple(sorted((first, second))):
                rank = 0.01
            elif edge == tuple(sorted((first, common_neighbour))):
                rank = 0.02
            else:
                rank = 0.99
            ranked_edges.append((rank, *edge))
        forest = CriticalKruskalForest(side_length, ranked_edges)
        geometry = ranked_pair_path_geometry(
            forest, second, common_neighbour
        )
        self.assertEqual(geometry.bucket_sizes, (1, 2))
        self.assertEqual(geometry.bucket_ranks, (0.01, 0.02))
        self.assertEqual(geometry.lca_rank, 0.02)

    def test_clock_coordinate_is_the_exact_inverse(self) -> None:
        for p in (0.75, 0.81, 0.9):
            coupling = log(p / (1.0 - p))
            for rank in (0.01, 0.2, Q_CRITICAL):
                if rank >= 2.0 * p - 1.0:
                    continue
                time = clock_time_from_rank(rank, p)
                reconstructed = p * (1.0 - exp(-coupling * time))
                self.assertAlmostEqual(reconstructed, rank)

    def test_same_ranked_hierarchy_can_be_evaluated_at_several_p(self) -> None:
        rng = random.Random(12)
        forest = CriticalKruskalForest(
            8, sample_ranked_edges(8, rng)
        )
        first, second = forest.sample_far_connected_pair(rng, 0.2)
        geometry = ranked_pair_path_geometry(forest, first, second)
        lower = evaluate_ranked_geometry(geometry, 0.80)
        higher = evaluate_ranked_geometry(geometry, 0.86)
        self.assertGreaterEqual(lower.attenuation, 0.0)
        self.assertGreaterEqual(higher.attenuation, 0.0)
        self.assertLessEqual(lower.correlation, 1.0)
        self.assertLessEqual(higher.correlation, 1.0)


if __name__ == "__main__":
    unittest.main()
