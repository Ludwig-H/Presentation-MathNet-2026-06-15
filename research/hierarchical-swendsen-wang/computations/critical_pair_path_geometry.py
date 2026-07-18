"""Critical-pair Palm diagnostics for the hierarchical Kruskal path.

The geometry below the critical cut is universal in the percolation
coordinate ``q``.  Each edge receives one uniform rank and is present at
level ``q`` when its rank is at most ``q``.  For a GSBM parameter ``p``, the
same rank is mapped to the exponential-clock time

    q = p * (1 - exp(-u_p * t)),  u_p = log(p / (1-p)).

This lets one reuse exactly the same critical hierarchy for several values
of ``p``.  Pairs are sampled with the correct finite-volume Palm bias:
conditional on being connected at ``q_c``, a component of size ``s`` has
weight ``s(s-1)``.  Rejection on the torus distance then conditions further
on the pair being macroscopically separated.

The reported attenuation is the explicitly factorized PATH-FAC oracle.  It
is a geometry diagnostic, not a replacement for the joint twisted transfer
of a full hierarchical sweep.
"""

from __future__ import annotations

import argparse
import json
import random
from bisect import bisect_left
from dataclasses import asdict, dataclass
from math import exp, log
from statistics import fmean, median
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from hierarchical_flip_probabilities import node_oracle_reliability


Edge = tuple[int, int]
RankedEdge = tuple[float, int, int]


@dataclass(frozen=True)
class PairPathGeometry:
    """Geometry of the two dendrogram arms joining a critical pair."""

    distance: int
    component_size: int
    lca_rank: float
    path_length: int
    bucket_sizes: tuple[int, ...]

    def bounded_bucket_count(self, maximum: int) -> int:
        """Count nontrivial buckets with size at most ``maximum``."""

        return sum(2 <= size <= maximum for size in self.bucket_sizes)


@dataclass(frozen=True)
class PairPathEvaluation:
    """One geometry evaluated at one signal parameter ``p``."""

    p: float
    attenuation: float
    correlation: float


@dataclass(frozen=True)
class DiagnosticSummary:
    """Compact reproducible summary of a diagnostic run."""

    side_length: int
    repetitions: int
    pairs_per_repetition: int
    distance_fraction: float
    pair_count: int
    mean_component_size: float
    mean_path_length: float
    median_path_length: float
    mean_lca_rank_fraction: float
    mean_bucket_count_2: float
    mean_bucket_count_4: float
    mean_bucket_count_8: float
    evaluations: dict[str, dict[str, float]]


def triangular_torus_edges(side_length: int) -> tuple[Edge, ...]:
    """Return every undirected edge of a simple triangular torus once."""

    if side_length < 4:
        raise ValueError("side_length must be at least 4")
    edges: set[Edge] = set()
    directions = ((1, 0), (0, 1), (1, -1))
    for y in range(side_length):
        for x in range(side_length):
            vertex = x + side_length * y
            for delta_x, delta_y in directions:
                neighbour_x = (x + delta_x) % side_length
                neighbour_y = (y + delta_y) % side_length
                neighbour = neighbour_x + side_length * neighbour_y
                edges.add(tuple(sorted((vertex, neighbour))))
    result = tuple(sorted(edges))
    expected = 3 * side_length * side_length
    if len(result) != expected:
        raise AssertionError("the torus construction produced duplicate edges")
    return result


def triangular_torus_distance(
    first: int, second: int, side_length: int
) -> int:
    """Exact graph distance for the chosen axial coordinates on the torus."""

    first_x, first_y = first % side_length, first // side_length
    second_x, second_y = second % side_length, second // side_length
    raw_x = first_x - second_x
    raw_y = first_y - second_y
    candidates = []
    for shift_x in (-side_length, 0, side_length):
        for shift_y in (-side_length, 0, side_length):
            delta_x = raw_x + shift_x
            delta_y = raw_y + shift_y
            candidates.append(
                max(abs(delta_x), abs(delta_y), abs(delta_x + delta_y))
            )
    return min(candidates)


class CriticalKruskalForest:
    """Kruskal merge forest retaining every physical cut size."""

    def __init__(
        self,
        side_length: int,
        ranked_edges: Sequence[RankedEdge],
        critical_rank: float = Q_CRITICAL,
    ) -> None:
        self.side_length = side_length
        self.vertex_count = side_length * side_length
        self.edges = triangular_torus_edges(side_length)
        if len(ranked_edges) != len(self.edges):
            raise ValueError("one rank is required for every torus edge")

        vertex_count = self.vertex_count
        self._dsu_parent = list(range(vertex_count))
        self._dsu_size = [1] * vertex_count
        self._tree_root = list(range(vertex_count))
        self.tree_parent = [-1] * (2 * vertex_count)
        self.left_child = [-1] * (2 * vertex_count)
        self.right_child = [-1] * (2 * vertex_count)
        self.merge_rank = [0.0] * (2 * vertex_count)
        self.bucket_size = [0] * (2 * vertex_count)
        self._next_tree_node = vertex_count

        boundary = [dict() for _ in range(vertex_count)]
        for first, second in self.edges:
            boundary[first][second] = boundary[first].get(second, 0) + 1
            boundary[second][first] = boundary[second].get(first, 0) + 1
        self._boundary = boundary

        for rank, first, second in sorted(ranked_edges):
            if rank > critical_rank:
                break
            self._merge(rank, first, second)

        components: dict[int, list[int]] = {}
        for vertex in range(vertex_count):
            root = self._find(vertex)
            components.setdefault(root, []).append(vertex)
        self.components = tuple(tuple(vertices) for vertices in components.values())

    @property
    def internal_nodes(self) -> range:
        """Identifiers of all merger nodes, in creation order."""

        return range(self.vertex_count, self._next_tree_node)

    def connected(self, first: int, second: int) -> bool:
        """Whether two vertices belong to the same final forest tree."""

        return self._find(first) == self._find(second)

    def tree_lca(self, first: int, second: int) -> int:
        """Return the merger node which is the LCA of two connected leaves."""

        if first == second:
            return first
        if not self.connected(first, second):
            raise ValueError("the two leaves belong to different forest trees")
        first_ancestors = set(self._internal_ancestor_chain(first))
        return next(
            node
            for node in self._internal_ancestor_chain(second)
            if node in first_ancestors
        )

    def _find(self, vertex: int) -> int:
        parent = self._dsu_parent
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def _merge(self, rank: float, first: int, second: int) -> None:
        first_root = self._find(first)
        second_root = self._find(second)
        if first_root == second_root:
            return

        if len(self._boundary[first_root]) < len(self._boundary[second_root]):
            first_root, second_root = second_root, first_root

        cut_size = self._boundary[first_root].get(second_root)
        if cut_size is None or cut_size <= 0:
            raise AssertionError("missing physical cut between merging components")

        node = self._next_tree_node
        self._next_tree_node += 1
        first_tree = self._tree_root[first_root]
        second_tree = self._tree_root[second_root]
        self.left_child[node] = first_tree
        self.right_child[node] = second_tree
        self.tree_parent[first_tree] = node
        self.tree_parent[second_tree] = node
        self.merge_rank[node] = rank
        self.bucket_size[node] = cut_size

        del self._boundary[first_root][second_root]
        del self._boundary[second_root][first_root]
        for neighbour, count in tuple(self._boundary[second_root].items()):
            if neighbour == first_root:
                continue
            del self._boundary[neighbour][second_root]
            self._boundary[neighbour][first_root] = (
                self._boundary[neighbour].get(first_root, 0) + count
            )
            self._boundary[first_root][neighbour] = (
                self._boundary[first_root].get(neighbour, 0) + count
            )

        self._boundary[second_root].clear()
        self._dsu_parent[second_root] = first_root
        self._dsu_size[first_root] += self._dsu_size[second_root]
        self._tree_root[first_root] = node

    def _internal_ancestor_chain(self, vertex: int) -> tuple[int, ...]:
        result = []
        node = vertex
        while self.tree_parent[node] != -1:
            node = self.tree_parent[node]
            result.append(node)
        return tuple(result)

    def pair_path_geometry(self, first: int, second: int) -> PairPathGeometry:
        """Return the two arms and LCA bucket exactly once."""

        path_nodes = self.pair_path_nodes(first, second)
        lca = path_nodes[-1]
        component_size = self._dsu_size[self._find(first)]
        return PairPathGeometry(
            distance=triangular_torus_distance(first, second, self.side_length),
            component_size=component_size,
            lca_rank=self.merge_rank[lca],
            path_length=len(path_nodes),
            bucket_sizes=tuple(self.bucket_size[node] for node in path_nodes),
        )

    def pair_path_nodes(self, first: int, second: int) -> tuple[int, ...]:
        """Return both strict arms and their LCA, with no duplicate node."""

        if first == second or not self.connected(first, second):
            raise ValueError("the pair must be distinct and critically connected")
        first_chain = self._internal_ancestor_chain(first)
        second_chain = self._internal_ancestor_chain(second)
        first_positions = {node: position for position, node in enumerate(first_chain)}
        lca = self.tree_lca(first, second)
        first_prefix = first_chain[: first_positions[lca]]
        second_prefix = second_chain[: second_chain.index(lca)]
        return first_prefix + second_prefix + (lca,)

    def sample_far_connected_pair(
        self,
        rng: random.Random,
        distance_fraction: float,
        maximum_attempts: int = 20_000,
    ) -> tuple[int, int]:
        """Sample the exact connected-pair Palm law, then condition on distance."""

        if not 0.0 < distance_fraction < 1.0:
            raise ValueError("distance_fraction must belong to (0, 1)")
        eligible = [component for component in self.components if len(component) >= 2]
        cumulative = []
        total = 0
        for component in eligible:
            total += len(component) * (len(component) - 1)
            cumulative.append(total)
        if total == 0:
            raise RuntimeError("the critical forest has no connected pair")

        minimum_distance = distance_fraction * self.side_length
        for _ in range(maximum_attempts):
            draw = rng.randrange(total)
            component_index = bisect_left(cumulative, draw + 1)
            component = eligible[component_index]
            first, second = rng.sample(component, 2)
            distance = triangular_torus_distance(
                first, second, self.side_length
            )
            if distance >= minimum_distance:
                return first, second
        raise RuntimeError(
            "failed to sample a macroscopically separated connected pair"
        )


def sample_ranked_edges(
    side_length: int, rng: random.Random
) -> tuple[RankedEdge, ...]:
    """Give every physical edge an independent uniform percolation rank."""

    return tuple(
        (rng.random(), first, second)
        for first, second in triangular_torus_edges(side_length)
    )


def clock_time_from_rank(rank: float, p: float) -> float:
    """Invert q_p(t)=p(1-exp(-u_p t))."""

    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2, 1)")
    if not 0.0 <= rank < p:
        raise ValueError("rank must belong to [0, p)")
    coupling = log(p / (1.0 - p))
    return -log(1.0 - rank / p) / coupling


@dataclass(frozen=True)
class RankedPairPathGeometry(PairPathGeometry):
    """Pair path with the clock rank attached to every bucket."""

    bucket_ranks: tuple[float, ...]


def ranked_pair_path_geometry(
    forest: CriticalKruskalForest, first: int, second: int
) -> RankedPairPathGeometry:
    """Return a pair path including all percolation ranks."""

    basic = forest.pair_path_geometry(first, second)
    path_nodes = forest.pair_path_nodes(first, second)
    return RankedPairPathGeometry(
        **asdict(basic),
        bucket_ranks=tuple(forest.merge_rank[node] for node in path_nodes),
    )


def evaluate_ranked_geometry(
    geometry: RankedPairPathGeometry, p: float
) -> PairPathEvaluation:
    """Evaluate exact PATH-FAC attenuation on ranked path data."""

    if Q_CRITICAL >= 2.0 * p - 1.0:
        raise ValueError("the critical cut must occur before time one")
    attenuation = 0.0
    for size, rank in zip(
        geometry.bucket_sizes, geometry.bucket_ranks, strict=True
    ):
        reliability = node_oracle_reliability(
            size, p, clock_time_from_rank(rank, p)
        )
        if reliability == 0.0:
            return PairPathEvaluation(p=p, attenuation=float("inf"), correlation=0.0)
        attenuation -= log(reliability)
    correlation = 0.0 if attenuation == float("inf") else exp(-attenuation)
    return PairPathEvaluation(p=p, attenuation=attenuation, correlation=correlation)


def run_diagnostic(
    side_length: int,
    repetitions: int,
    pairs_per_repetition: int,
    p_values: Sequence[float],
    distance_fraction: float,
    seed: int,
) -> DiagnosticSummary:
    """Run a reproducible finite-volume critical-pair diagnostic."""

    if repetitions <= 0 or pairs_per_repetition <= 0:
        raise ValueError("repetitions and pairs_per_repetition must be positive")
    rng = random.Random(seed)
    geometries: list[RankedPairPathGeometry] = []
    for _ in range(repetitions):
        forest = CriticalKruskalForest(
            side_length, sample_ranked_edges(side_length, rng)
        )
        for _ in range(pairs_per_repetition):
            first, second = forest.sample_far_connected_pair(
                rng, distance_fraction
            )
            geometries.append(ranked_pair_path_geometry(forest, first, second))

    evaluations: dict[str, dict[str, float]] = {}
    for p in p_values:
        values = [evaluate_ranked_geometry(geometry, p) for geometry in geometries]
        evaluations[f"{p:.12g}"] = {
            "mean_attenuation": fmean(value.attenuation for value in values),
            "median_attenuation": median(value.attenuation for value in values),
            "mean_path_fac_correlation": fmean(value.correlation for value in values),
            "mean_same_relation_probability": 0.5
            * (1.0 + fmean(value.correlation for value in values)),
        }

    return DiagnosticSummary(
        side_length=side_length,
        repetitions=repetitions,
        pairs_per_repetition=pairs_per_repetition,
        distance_fraction=distance_fraction,
        pair_count=len(geometries),
        mean_component_size=fmean(item.component_size for item in geometries),
        mean_path_length=fmean(item.path_length for item in geometries),
        median_path_length=median(item.path_length for item in geometries),
        mean_lca_rank_fraction=fmean(
            item.lca_rank / Q_CRITICAL for item in geometries
        ),
        mean_bucket_count_2=fmean(
            item.bounded_bucket_count(2) for item in geometries
        ),
        mean_bucket_count_4=fmean(
            item.bounded_bucket_count(4) for item in geometries
        ),
        mean_bucket_count_8=fmean(
            item.bounded_bucket_count(8) for item in geometries
        ),
        evaluations=evaluations,
    )


def _parse_p_values(raw: str) -> tuple[float, ...]:
    values = tuple(float(value) for value in raw.split(","))
    if not values:
        raise ValueError("at least one p value is required")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=32)
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--pairs", type=int, default=20)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument(
        "--p-values", default="0.794659275831,0.81,0.835805792367,0.86"
    )
    parser.add_argument("--seed", type=int, default=20260718)
    arguments = parser.parse_args()
    summary = run_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        pairs_per_repetition=arguments.pairs,
        p_values=_parse_p_values(arguments.p_values),
        distance_fraction=arguments.distance_fraction,
        seed=arguments.seed,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
