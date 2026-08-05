"""Tests of the replicated-target measurement (GSBM, note GSBM/03)."""

from __future__ import annotations

import random
import unittest

from gsbm_replicated_target_exact import (
    GsbmExact,
    beta_c,
    coupling,
    measure,
    triangular_torus,
)


class TestGeometry(unittest.TestCase):
    def test_torus_edges(self) -> None:
        n, edges = triangular_torus(3)
        self.assertEqual(n, 9)
        self.assertEqual(len(edges), 27)
        self.assertEqual(len(set(map(lambda e: tuple(sorted(e)), edges))), 27)
        for (i, j) in edges:
            self.assertNotEqual(i, j)

    def test_beta_c_values(self) -> None:
        self.assertAlmostEqual(beta_c(0.81), 0.386168, places=5)
        self.assertAlmostEqual(beta_c(0.8), 0.410717, places=5)


class TestExactInner(unittest.TestCase):
    def setUp(self) -> None:
        self.model = GsbmExact(3, 0.8, random.Random(11))

    def test_posterior_concentrates_when_all_edges_agree(self) -> None:
        obs = [1] * len(self.model.edges)
        sat = self.model.satisfied_tables(obs)
        weights = self.model.posterior_weights(sat)
        corr = self.model.pair_correlation(weights, 0, 1)
        self.assertGreater(corr, 0.99)

    def test_correlations_bounded(self) -> None:
        obs = self.model.sample_observation()
        sat = self.model.satisfied_tables(obs)
        weights = self.model.posterior_weights(sat)
        for (i, j), c in self.model.all_pair_correlations(weights).items():
            self.assertLessEqual(abs(c), 1.0 + 1e-12)

    def test_cut_edges_partition_same_root_edges(self) -> None:
        """Each same-final-root edge appears in exactly one merge cut
        (partition-by-LCA lemma, note 42 section 3b)."""

        obs = self.model.sample_observation()
        sat = self.model.satisfied_tables(obs)
        weights = self.model.posterior_weights(sat)
        sigma = self.model.sample_configuration(weights)
        merges, blocks, roots = self.model.sample_dendrogram(sigma, sat)
        counts = {}
        for (_, cut) in merges:
            for e in cut:
                counts[e] = counts.get(e, 0) + 1
        for e, (i, j) in enumerate(self.model.edges):
            expected = 1 if roots.find(i) == roots.find(j) else 0
            self.assertEqual(counts.get(e, 0), expected)

    def test_merge_times_increasing_and_winner_satisfied(self) -> None:
        obs = self.model.sample_observation()
        sat = self.model.satisfied_tables(obs)
        weights = self.model.posterior_weights(sat)
        sigma = self.model.sample_configuration(weights)
        merges, _, _ = self.model.sample_dendrogram(sigma, sat)
        times = [t for (t, _) in merges]
        self.assertEqual(times, sorted(times))
        for (_, cut) in merges:
            self.assertTrue(any(sat[e][sigma] for e in cut))

    def test_gibbs_given_dendrogram_supported_and_bounded(self) -> None:
        obs = self.model.sample_observation()
        sat = self.model.satisfied_tables(obs)
        weights = self.model.posterior_weights(sat)
        sigma = self.model.sample_configuration(weights)
        merges, _, _ = self.model.sample_dendrogram(sigma, sat)
        nu = self.model.gibbs_given_dendrogram(merges, sat)
        self.assertGreater(nu[sigma], 0.0)
        self.assertGreater(sum(nu), 0.0)
        corr = self.model.gibbs_pair_correlation(nu, 0, 5)
        self.assertLessEqual(abs(corr), 1.0 + 1e-12)


class TestEstimator(unittest.TestCase):
    def test_measure_runs_and_bounds(self) -> None:
        result = measure(size=3, p=0.8, n_obs=2, k_dendro=3, seed=7)
        self.assertGreaterEqual(result["q_full"], 1.0 / 9.0 - 1e-12)
        self.assertLessEqual(result["q_full"], 1.0 + 1e-12)
        self.assertLessEqual(result["d_cross"], result["j_cross"] + 1e-9)
        self.assertGreaterEqual(result["j_cross"], 0.0)
        self.assertGreaterEqual(result["s_crit"], 0.0)
        self.assertLessEqual(result["s_crit"], 1.0 + 1e-12)

    def test_unbiased_square_can_be_negative_but_dominated(self) -> None:
        # the U-statistic may fluctuate below zero for single pairs,
        # but the Jensen envelope always dominates the debiased square.
        result = measure(size=3, p=0.76, n_obs=3, k_dendro=4, seed=19)
        self.assertLessEqual(result["d_cross"], result["j_cross"] + 1e-9)


if __name__ == "__main__":
    unittest.main()
