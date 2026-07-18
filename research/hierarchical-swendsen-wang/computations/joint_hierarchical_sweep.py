"""Small-volume audit of the genuine joint hierarchical heat-bath sweep.

Unlike PATH-FAC, this module recomputes every factor affected by a proposed
cluster flip.  For a node ``u`` it therefore uses all buckets at ``u`` and
above ``u`` in the exact conditional law

    product_v Lambda_v * exp((1-beta_v) Lambda_v).

The implementation is deliberately exhaustive and intended for small tori.
It provides a counter-audit for factorization and a numerical prototype for
the finite-state twisted transfer; it is not a proof of an infinite-volume
threshold.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import exp, fsum, log, sqrt
from statistics import fmean
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    RankedEdge,
    clock_time_from_rank,
    sample_ranked_edges,
)


@dataclass(frozen=True)
class ConditionalCorrelationEstimate:
    """Monte-Carlo estimates conditional on one marked hierarchy."""

    correlation: float
    correlation_square_unbiased: float


@dataclass(frozen=True)
class JointSweepSummary:
    """Second-moment diagnostics for top-down and bottom-up sweeps."""

    side_length: int
    repetitions: int
    sweeps_per_environment: int
    passes_per_sweep: int
    distance_fraction: float
    results: dict[str, dict[str, dict[str, float]]]


class HierarchicalSweepEnvironment:
    """Exact conditional factors for one GSBM hierarchy at fixed ``p``."""

    def __init__(
        self,
        side_length: int,
        ranked_edges: Sequence[RankedEdge],
        p: float,
    ) -> None:
        if Q_CRITICAL >= 2.0 * p - 1.0:
            raise ValueError("p must put the critical cut strictly before time one")
        self.side_length = side_length
        self.p = p
        self.coupling = log(p / (1.0 - p))
        self.forest = CriticalKruskalForest(
            side_length, ranked_edges, critical_rank=2.0 * p - 1.0
        )
        self.vertex_count = self.forest.vertex_count
        self.edge_rank = {
            tuple(sorted((first, second))): rank
            for rank, first, second in ranked_edges
        }

        capacity = len(self.forest.tree_parent)
        self.cluster_mask = [0] * capacity
        for vertex in range(self.vertex_count):
            self.cluster_mask[vertex] = 1 << vertex
        for node in self.forest.internal_nodes:
            self.cluster_mask[node] = (
                self.cluster_mask[self.forest.left_child[node]]
                | self.cluster_mask[self.forest.right_child[node]]
            )

        bucket_edges: list[list[tuple[int, int]]] = [
            [] for _ in range(capacity)
        ]
        for edge in self.forest.edges:
            first, second = edge
            if self.forest.connected(first, second):
                bucket_edges[self.forest.tree_lca(first, second)].append(edge)
        self.bucket_edges = tuple(tuple(edges) for edges in bucket_edges)
        for node in self.forest.internal_nodes:
            if len(self.bucket_edges[node]) != self.forest.bucket_size[node]:
                raise AssertionError("bucket assignment disagrees with the Kruskal cut")

        factor_chains: list[tuple[int, ...]] = [tuple() for _ in range(capacity)]
        for update_node in range(self.vertex_count):
            factor_chains[update_node] = self._ancestors_from(
                self.forest.tree_parent[update_node]
            )
        for update_node in self.forest.internal_nodes:
            factor_chains[update_node] = self._ancestors_from(update_node)
        self.factor_chains = tuple(factor_chains)

        self.bottom_up_order = tuple(range(self.vertex_count)) + tuple(
            self.forest.internal_nodes
        )
        descending_nodes = tuple(reversed(tuple(self.forest.internal_nodes)))
        self.top_down_order = descending_nodes + tuple(range(self.vertex_count))

    def _ancestors_from(self, first_node: int) -> tuple[int, ...]:
        result = []
        node = first_node
        while node != -1:
            result.append(node)
            node = self.forest.tree_parent[node]
        return tuple(result)

    def _satisfied_count(self, factor_node: int, state: int) -> int:
        count = 0
        for first, second in self.bucket_edges[factor_node]:
            truth_satisfied = self.edge_rank[(first, second)] <= self.p
            relative_parity = ((state >> first) ^ (state >> second)) & 1
            count += truth_satisfied != bool(relative_parity)
        return count

    def conditional_log_weight(self, update_node: int, state: int) -> float:
        """Log of all factors changed by ``update_node``, up to constants."""

        result = 0.0
        for factor_node in self.factor_chains[update_node]:
            count = self._satisfied_count(factor_node, state)
            if count == 0:
                return float("-inf")
            time = clock_time_from_rank(
                self.forest.merge_rank[factor_node], self.p
            )
            result += log(count) + (1.0 - time) * self.coupling * count
        return result

    def full_log_weight(self, state: int) -> float:
        """Log of the complete conditional density, up to one constant.

        This slower expression is kept as an independent counter-audit of
        ``conditional_log_weight``.  It is useful on the small volumes for
        which this module is intended.
        """

        result = 0.0
        for factor_node in self.forest.internal_nodes:
            count = self._satisfied_count(factor_node, state)
            if count == 0:
                return float("-inf")
            time = clock_time_from_rank(
                self.forest.merge_rank[factor_node], self.p
            )
            result += log(count) + (1.0 - time) * self.coupling * count
        return result

    def proposal_masks(self, update_node: int) -> tuple[int, ...]:
        """Two leaf moves or the genuine four-state internal-node orbit."""

        if update_node < self.vertex_count:
            return (0, 1 << update_node)
        left = self.cluster_mask[self.forest.left_child[update_node]]
        right = self.cluster_mask[self.forest.right_child[update_node]]
        return (0, left, right, left | right)

    def proposal_probabilities(
        self, update_node: int, state: int
    ) -> tuple[tuple[int, float], ...]:
        """Return the exact heat-bath law on one update orbit."""

        masks = self.proposal_masks(update_node)
        log_weights = tuple(
            self.conditional_log_weight(update_node, state ^ mask)
            for mask in masks
        )
        maximum = max(log_weights)
        if maximum == float("-inf"):
            raise AssertionError("the current conditional orbit has zero mass")
        weights = tuple(
            0.0 if value == float("-inf") else exp(value - maximum)
            for value in log_weights
        )
        normalizer = fsum(weights)
        return tuple(
            (mask, weight / normalizer)
            for mask, weight in zip(masks, weights, strict=True)
        )

    def update(self, update_node: int, state: int, rng: random.Random) -> int:
        probabilities = self.proposal_probabilities(update_node, state)
        draw = rng.random()
        cumulative = 0.0
        for mask, probability in probabilities:
            cumulative += probability
            if draw <= cumulative:
                return state ^ mask
        return state ^ probabilities[-1][0]

    def sweep(
        self, order: str, rng: random.Random, pass_count: int = 1
    ) -> int:
        """Start at the truth and perform one complete systematic sweep.

        At an internal final root, the masks ``0`` and ``left | right``
        have equal weight, as do ``left`` and ``right``.  The root update
        therefore includes the fair global recolouring.  For an isolated
        root, the leaf update has two equal weights.  No extra recolouring
        step is needed; it would not change the relation of the selected
        pair, which already lies in one critical component.
        """

        if order == "bottom-up":
            update_order = self.bottom_up_order
        elif order == "top-down":
            update_order = self.top_down_order
        else:
            raise ValueError("order must be 'bottom-up' or 'top-down'")
        if pass_count <= 0:
            raise ValueError("pass_count must be positive")
        state = 0
        for _ in range(pass_count):
            for node in update_order:
                state = self.update(node, state, rng)
        return state


def pair_relation(state: int, first: int, second: int) -> int:
    """Return the final relation relative to the planted relation."""

    return 1 if (((state >> first) ^ (state >> second)) & 1) == 0 else -1


def unbiased_rademacher_mean_square(samples: Sequence[int]) -> float:
    """Unbiased estimate of the square of a Rademacher mean."""

    sample_count = len(samples)
    if sample_count < 2:
        raise ValueError("at least two samples are required")
    total = sum(samples)
    return (total * total - sample_count) / (
        sample_count * (sample_count - 1)
    )


def estimate_conditional_pair_correlation(
    environment: HierarchicalSweepEnvironment,
    first: int,
    second: int,
    order: str,
    sweep_count: int,
    pass_count: int,
    rng: random.Random,
) -> ConditionalCorrelationEstimate:
    """Estimate ``H_S(i,j)`` and its square without plug-in bias."""

    samples = [
        pair_relation(
            environment.sweep(order, rng, pass_count), first, second
        )
        for _ in range(sweep_count)
    ]
    return ConditionalCorrelationEstimate(
        correlation=fmean(samples),
        correlation_square_unbiased=unbiased_rademacher_mean_square(samples),
    )


def _standard_error(values: Sequence[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = fmean(values)
    variance = fsum((value - mean) ** 2 for value in values) / (
        len(values) - 1
    )
    return sqrt(variance / len(values))


def run_joint_sweep_diagnostic(
    side_length: int,
    repetitions: int,
    sweeps_per_environment: int,
    passes_per_sweep: int,
    p_values: Sequence[float],
    distance_fraction: float,
    seed: int,
) -> JointSweepSummary:
    """Compare first and second transfer moments under critical-pair Palm bias."""

    if repetitions <= 0 or sweeps_per_environment < 2 or passes_per_sweep <= 0:
        raise ValueError("positive repetitions and at least two sweeps are required")
    geometry_rng = random.Random(seed)
    stored: dict[tuple[float, str], list[ConditionalCorrelationEstimate]] = {
        (p, order): []
        for p in p_values
        for order in ("bottom-up", "top-down")
    }
    for repetition in range(repetitions):
        ranked_edges = sample_ranked_edges(side_length, geometry_rng)
        critical_forest = CriticalKruskalForest(side_length, ranked_edges)
        first, second = critical_forest.sample_far_connected_pair(
            geometry_rng, distance_fraction
        )
        for p_index, p in enumerate(p_values):
            environment = HierarchicalSweepEnvironment(
                side_length, ranked_edges, p
            )
            for order_index, order in enumerate(("bottom-up", "top-down")):
                sweep_seed = (
                    seed
                    ^ ((repetition + 1) * 0x9E3779B97F4A7C15)
                    ^ ((p_index + 1) * 0xBF58476D1CE4E5B9)
                    ^ (order_index * 0x94D049BB133111EB)
                ) & ((1 << 64) - 1)
                stored[(p, order)].append(
                    estimate_conditional_pair_correlation(
                        environment,
                        first,
                        second,
                        order,
                        sweeps_per_environment,
                        passes_per_sweep,
                        random.Random(sweep_seed),
                    )
                )

    results: dict[str, dict[str, dict[str, float]]] = {}
    for p in p_values:
        p_result: dict[str, dict[str, float]] = {}
        for order in ("bottom-up", "top-down"):
            estimates = stored[(p, order)]
            first_moments = [value.correlation for value in estimates]
            second_moments = [
                value.correlation_square_unbiased for value in estimates
            ]
            mean_correlation = fmean(first_moments)
            p_result[order] = {
                "annealed_correlation": mean_correlation,
                "annealed_correlation_standard_error": _standard_error(
                    first_moments
                ),
                "same_relation_probability": 0.5 * (1.0 + mean_correlation),
                "conditional_second_moment": fmean(second_moments),
                "conditional_second_moment_standard_error": _standard_error(
                    second_moments
                ),
            }
        results[f"{p:.12g}"] = p_result

    return JointSweepSummary(
        side_length=side_length,
        repetitions=repetitions,
        sweeps_per_environment=sweeps_per_environment,
        passes_per_sweep=passes_per_sweep,
        distance_fraction=distance_fraction,
        results=results,
    )


def _parse_p_values(raw: str) -> tuple[float, ...]:
    return tuple(float(value) for value in raw.split(","))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=6)
    parser.add_argument("--repetitions", type=int, default=40)
    parser.add_argument("--sweeps", type=int, default=200)
    parser.add_argument("--passes", type=int, default=1)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument(
        "--p-values", default="0.8,0.81,0.835805792367"
    )
    parser.add_argument("--seed", type=int, default=20260718)
    arguments = parser.parse_args()
    summary = run_joint_sweep_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        sweeps_per_environment=arguments.sweeps,
        passes_per_sweep=arguments.passes,
        p_values=_parse_p_values(arguments.p_values),
        distance_fraction=arguments.distance_fraction,
        seed=arguments.seed,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
