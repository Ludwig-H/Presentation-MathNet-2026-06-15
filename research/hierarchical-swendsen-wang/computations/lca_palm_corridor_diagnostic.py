"""Finite-volume LCA-Palm diagnostics for triangular-torus corridors.

This module keeps two different Palm constructions separate.

``critical_cut_intensity``
    At the deterministic rank ``q_c``, enumerate every adjacent pair of
    components ``A, B``.  Its predictable LCA-Palm weight, up to the common
    rank hazard, is exactly ``m(A, B) * N_rho(A, B)``.  The two components
    are joined by a synthetic merger at ``q_c`` and the two Kruskal arms
    below that merger are recorded.  This is only a fixed-level intensity
    benchmark: the synthetic merger changes the skeleton and is *not* a
    Blackwell-favourable version of the real final corridor.

``final_realized_event``
    Run Kruskal to ``q_1 = 2p - 1`` and enumerate actual merger nodes.  A
    realized node is weighted by ``N_rho`` only: its occurrence law already
    contains the cut-size bias ``m``.  Multiplying by ``m * N_rho`` here
    would double-count the merger race and produce an erroneous ``m^2``
    bias.  Both the actual ranks and the favourable ranks
    ``min(q_v, q_c)`` are recorded on the unchanged corridor skeleton.
    This rank-capping operation, unlike the synthetic benchmark above, is
    the geometry-fixed favourable criticalization relevant to a Blackwell
    comparison.

For a cut of size ``m`` at rank ``q``, the reported geometric charge is

    m * h_p(q)^2,  h_p(q) = (2p - 1 - q) / (1 - q).

Counts based on bounded size or bounded charge are *geometric proxies only*.
No exterior message ``B``, lateral-port screening, or replicated transfer
coefficient is computed.  Consequently the output never labels a cut as
screened and cannot by itself establish hierarchical contraction or a weak-
recovery threshold.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, sqrt
from typing import Callable, Iterable, Literal, Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_component_boundary import closed_categories
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    RankedEdge,
    clock_time_from_rank,
    sample_ranked_edges,
    triangular_torus_distance,
)


PalmMode = Literal["critical_cut_intensity", "final_realized_event"]


@dataclass(frozen=True)
class CorridorCut:
    """One physical Kruskal bucket on a two-arm corridor.

    ``favourable_rank`` is obtained by criticalizing the time only.  The
    bucket, its size, and every other item of the skeleton remain unchanged.
    """

    node: int | None
    is_lca: bool
    bucket_size: int
    actual_rank: float
    favourable_rank: float
    actual_beta: float
    favourable_beta: float
    actual_charge: float
    favourable_charge: float


@dataclass(frozen=True)
class PalmCorridorObservation:
    """One sampled endpoint pair attached to one Palm cut or event.

    The physical buckets represented by ``cuts`` have been checked to be
    pairwise edge-disjoint.  This does not mean that neighbourhoods around
    the cuts are screened from the exterior.
    """

    mode: PalmMode
    repetition: int
    first: int
    second: int
    distance: int
    n_rho: int
    lca_bucket_size: int
    palm_weight: int
    palm_weight_definition: str
    cuts: tuple[CorridorCut, ...]
    physical_buckets_pairwise_disjoint: bool
    screening_computed: bool = False
    screened_cut_count: int | None = None

    def bounded_size_cut_count(self, maximum_size: int) -> int:
        """Count disjoint corridor cuts satisfying ``2 <= m <= maximum``."""

        if maximum_size < 2:
            raise ValueError("maximum_size must be at least 2")
        return sum(
            2 <= cut.bucket_size <= maximum_size for cut in self.cuts
        )

    def bounded_charge_cut_count(
        self, maximum_charge: float, *, favourable: bool
    ) -> int:
        """Count nontrivial cuts whose selected geometric charge is bounded."""

        if maximum_charge < 0.0:
            raise ValueError("maximum_charge must be nonnegative")
        return sum(
            cut.bucket_size >= 2
            and (
                cut.favourable_charge if favourable else cut.actual_charge
            )
            <= maximum_charge
            for cut in self.cuts
        )

    def geometric_proxy_count(
        self,
        maximum_size: int,
        maximum_charge: float,
        *,
        favourable: bool,
    ) -> int:
        """Count cuts passing both geometric filters, without screening."""

        if maximum_size < 2:
            raise ValueError("maximum_size must be at least 2")
        if maximum_charge < 0.0:
            raise ValueError("maximum_charge must be nonnegative")
        return sum(
            2 <= cut.bucket_size <= maximum_size
            and (
                cut.favourable_charge if favourable else cut.actual_charge
            )
            <= maximum_charge
            for cut in self.cuts
        )


@dataclass(frozen=True)
class PalmModeSummary:
    """Weighted empirical Campbell summary for one Palm convention."""

    mode: PalmMode
    weight_convention: str
    observation_count: int
    total_palm_weight: int
    effective_observation_count: float
    repetition_count: int
    repetitions_with_positive_palm_mass: int
    effective_repetition_count: float
    weighted_mean_distance: float
    jackknife_standard_error_distance: float | None
    weighted_mean_corridor_cut_count: float
    jackknife_standard_error_corridor_cut_count: float | None
    weighted_mean_lca_bucket_size: float
    jackknife_standard_error_lca_bucket_size: float | None
    weighted_mean_lca_rank: float
    weighted_postcritical_lca_fraction: float
    weighted_mean_disjoint_bucket_size_two_count: float
    jackknife_standard_error_disjoint_bucket_size_two_count: float | None
    weighted_mean_disjoint_bounded_size_count: float
    jackknife_standard_error_disjoint_bounded_size_count: float | None
    weighted_mean_disjoint_actual_bounded_charge_count: float
    weighted_mean_disjoint_favourable_bounded_charge_count: float
    weighted_mean_disjoint_actual_geometric_proxy_count: float
    jackknife_standard_error_disjoint_actual_geometric_proxy_count: (
        float | None
    )
    weighted_mean_disjoint_favourable_geometric_proxy_count: float
    jackknife_standard_error_disjoint_favourable_geometric_proxy_count: (
        float | None
    )
    weighted_pairwise_disjoint_fraction: float
    interpretation: str
    uncertainty_note: str
    screening_computed: bool
    weighted_mean_screened_cut_count: float | None


@dataclass(frozen=True)
class CampbellWeightAudit:
    """Deterministic comparison of the two equivalent Palm representations."""

    cut_sizes: tuple[int, ...]
    distant_pair_counts: tuple[int, ...]
    cut_intensity_contributions: tuple[int, ...]
    event_rate_factors: tuple[int, ...]
    realized_event_weights: tuple[int, ...]
    event_expected_contributions: tuple[int, ...]
    incorrectly_double_biased_contributions: tuple[int, ...]
    representations_agree: bool
    double_bias_detected: bool


@dataclass(frozen=True)
class ChargeIdentityAudit:
    """Independent algebraic check of the rank-to-charge conversion."""

    checked_rank_count: int
    maximum_margin_error: float
    passed: bool


@dataclass(frozen=True)
class LcaPairPartitionAudit:
    """Check that distant connected ordered pairs partition by realized LCA."""

    checked_environment_count: int
    realized_event_pair_total: int
    connected_pair_total: int
    maximum_environment_difference: int
    passed: bool


@dataclass(frozen=True)
class LcaPalmDiagnostic:
    """Complete reproducible output of the two finite-volume diagnostics."""

    side_length: int
    repetitions: int
    p: float
    critical_rank: float
    final_rank: float
    distance_fraction: float
    maximum_bucket_size: int
    maximum_charge: float
    seed: int
    critical_snapshot_benchmark: PalmModeSummary
    final_corridor: PalmModeSummary
    campbell_weight_audit: CampbellWeightAudit
    charge_identity_audit: ChargeIdentityAudit
    lca_pair_partition_audit: LcaPairPartitionAudit


@dataclass(frozen=True)
class _SnapshotPalmCut:
    """Adjacent components and their fixed-level cut-intensity weight."""

    left_root: int
    right_root: int
    left_vertices: frozenset[int]
    right_vertices: frozenset[int]
    cut_edges: tuple[tuple[int, int], ...]
    distant_pairs: tuple[tuple[int, int], ...]

    @property
    def bucket_size(self) -> int:
        return len(self.cut_edges)

    @property
    def n_rho(self) -> int:
        return 2 * len(self.distant_pairs)

    @property
    def palm_weight(self) -> int:
        return cut_intensity_palm_weight(self.bucket_size, self.n_rho)


def cut_intensity_palm_weight(bucket_size: int, n_rho: int) -> int:
    """Return ``m*N_rho`` for a predictable cut at a fixed level."""

    if bucket_size <= 0:
        raise ValueError("bucket_size must be positive")
    if n_rho < 0:
        raise ValueError("n_rho must be nonnegative")
    return bucket_size * n_rho


def realized_event_palm_weight(n_rho: int) -> int:
    """Return the Palm weight of an already-realized Kruskal merger.

    The merger event itself was selected by a race with rate proportional to
    its physical cut size, so no additional factor ``m`` is allowed here.
    """

    if n_rho < 0:
        raise ValueError("n_rho must be nonnegative")
    return n_rho


def campbell_weight_counter_audit(
    cut_sizes: Sequence[int] = (2, 5),
    distant_pair_counts: Sequence[int] = (3, 7),
) -> CampbellWeightAudit:
    """Check ``m*N = (event rate m)*(event weight N)`` exactly.

    The common time-hazard factor is omitted because it cancels between cuts
    at a fixed level.  The deliberately wrong last contribution corresponds
    to multiplying realized events by ``m*N`` and therefore exposes ``m^2``.
    """

    sizes = tuple(cut_sizes)
    counts = tuple(distant_pair_counts)
    if len(sizes) != len(counts) or not sizes:
        raise ValueError("the two nonempty sequences must have equal length")
    if any(size <= 0 for size in sizes):
        raise ValueError("cut sizes must be positive")
    if any(count < 0 for count in counts):
        raise ValueError("pair counts must be nonnegative")

    cut_contributions = tuple(
        cut_intensity_palm_weight(size, count)
        for size, count in zip(sizes, counts, strict=True)
    )
    event_weights = tuple(realized_event_palm_weight(count) for count in counts)
    event_contributions = tuple(
        rate * weight
        for rate, weight in zip(sizes, event_weights, strict=True)
    )
    double_biased = tuple(
        size * cut_intensity_palm_weight(size, count)
        for size, count in zip(sizes, counts, strict=True)
    )
    return CampbellWeightAudit(
        cut_sizes=sizes,
        distant_pair_counts=counts,
        cut_intensity_contributions=cut_contributions,
        event_rate_factors=sizes,
        realized_event_weights=event_weights,
        event_expected_contributions=event_contributions,
        incorrectly_double_biased_contributions=double_biased,
        representations_agree=cut_contributions == event_contributions,
        double_bias_detected=double_biased != cut_contributions,
    )


def residual_margin_from_rank(p: float, rank: float) -> float:
    """Return ``h_p`` directly in the percolation-rank coordinate."""

    final_rank = 2.0 * p - 1.0
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2, 1)")
    if rank < 0.0 or rank > final_rank + 1e-14:
        raise ValueError("rank must belong to [0, 2p-1]")
    rank = min(rank, final_rank)
    return (final_rank - rank) / (1.0 - rank)


def geometric_charge(bucket_size: int, p: float, rank: float) -> float:
    """Return the proxy ``m*h_p(rank)^2`` for one geometry-fixed cut."""

    if bucket_size <= 0:
        raise ValueError("bucket_size must be positive")
    margin = residual_margin_from_rank(p, rank)
    return bucket_size * margin * margin


def charge_identity_counter_audit(
    p: float, ranks: Iterable[float]
) -> ChargeIdentityAudit:
    """Compare the direct rank formula with the residual closed-edge law."""

    rank_tuple = tuple(ranks)
    if not rank_tuple:
        raise ValueError("at least one rank is required")
    errors = []
    for rank in rank_tuple:
        beta = clock_time_from_rank(rank, p)
        from_categories = closed_categories(p, beta).signed_margin
        direct = residual_margin_from_rank(p, rank)
        errors.append(abs(from_categories - direct))
    maximum_error = max(errors)
    return ChargeIdentityAudit(
        checked_rank_count=len(rank_tuple),
        maximum_margin_error=maximum_error,
        passed=maximum_error <= 2e-12,
    )


def _used_tree_nodes(forest: CriticalKruskalForest) -> tuple[int, ...]:
    return tuple(range(forest.vertex_count)) + tuple(forest.internal_nodes)


def _tree_roots(forest: CriticalKruskalForest) -> tuple[int, ...]:
    return tuple(
        node
        for node in _used_tree_nodes(forest)
        if forest.tree_parent[node] == -1
    )


def _leaf_sets(forest: CriticalKruskalForest) -> dict[int, frozenset[int]]:
    """Reconstruct all child sets without depending on mutable DSU internals."""

    leaves: dict[int, frozenset[int]] = {
        vertex: frozenset((vertex,)) for vertex in range(forest.vertex_count)
    }
    for node in forest.internal_nodes:
        left = forest.left_child[node]
        right = forest.right_child[node]
        if left not in leaves or right not in leaves:
            raise AssertionError("Kruskal nodes are not in topological order")
        leaves[node] = leaves[left] | leaves[right]
    return leaves


def _internal_ancestor_chain(
    forest: CriticalKruskalForest, vertex: int
) -> tuple[int, ...]:
    result = []
    node = vertex
    while forest.tree_parent[node] != -1:
        node = forest.tree_parent[node]
        result.append(node)
    return tuple(result)


def _bucket_edges(
    forest: CriticalKruskalForest,
    node: int,
    leaf_sets: dict[int, frozenset[int]],
) -> tuple[tuple[int, int], ...]:
    left = leaf_sets[forest.left_child[node]]
    right = leaf_sets[forest.right_child[node]]
    edges = tuple(
        edge
        for edge in forest.edges
        if (edge[0] in left and edge[1] in right)
        or (edge[1] in left and edge[0] in right)
    )
    if len(edges) != forest.bucket_size[node]:
        raise AssertionError("reconstructed cut differs from Kruskal bucket")
    return edges


def _far_cross_pairs(
    left: Iterable[int],
    right: Iterable[int],
    side_length: int,
    distance_fraction: float,
) -> tuple[tuple[int, int], ...]:
    minimum_distance = distance_fraction * side_length
    return tuple(
        (first, second)
        for first in sorted(left)
        for second in sorted(right)
        if triangular_torus_distance(first, second, side_length)
        >= minimum_distance
    )


def _sample_oriented_pair(
    pairs: Sequence[tuple[int, int]], rng: random.Random
) -> tuple[int, int]:
    if not pairs:
        raise ValueError("at least one distant cross-pair is required")
    first, second = pairs[rng.randrange(len(pairs))]
    return (first, second) if rng.randrange(2) == 0 else (second, first)


def lca_pair_partition_counts(
    forest: CriticalKruskalForest,
    distance_fraction: float,
) -> tuple[int, int]:
    """Return event-``N_rho`` sum and connected distant ordered-pair count.

    Equality is an exact finite-volume identity: every connected ordered
    pair has one and only one realized Kruskal LCA.
    """

    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0, 1)")
    leaf_sets = _leaf_sets(forest)
    event_pair_total = 0
    for node in forest.internal_nodes:
        event_pair_total += 2 * len(
            _far_cross_pairs(
                leaf_sets[forest.left_child[node]],
                leaf_sets[forest.right_child[node]],
                forest.side_length,
                distance_fraction,
            )
        )

    minimum_distance = distance_fraction * forest.side_length
    connected_pair_total = 0
    for root in _tree_roots(forest):
        vertices = sorted(leaf_sets[root])
        connected_pair_total += sum(
            first != second
            and triangular_torus_distance(
                first, second, forest.side_length
            )
            >= minimum_distance
            for first in vertices
            for second in vertices
        )
    return event_pair_total, connected_pair_total


def _snapshot_palm_cuts(
    forest: CriticalKruskalForest,
    distance_fraction: float,
) -> tuple[_SnapshotPalmCut, ...]:
    """Enumerate fixed-rank adjacent components with their exact ``m*N``."""

    leaf_sets = _leaf_sets(forest)
    roots = _tree_roots(forest)
    membership: dict[int, int] = {}
    for root_index, root in enumerate(roots):
        for vertex in leaf_sets[root]:
            membership[vertex] = root_index
    if len(membership) != forest.vertex_count:
        raise AssertionError("tree roots do not partition the torus")

    boundary_edges: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for edge in forest.edges:
        left_index = membership[edge[0]]
        right_index = membership[edge[1]]
        if left_index == right_index:
            continue
        key = tuple(sorted((left_index, right_index)))
        boundary_edges.setdefault(key, []).append(edge)

    result = []
    for (left_index, right_index), edges in sorted(boundary_edges.items()):
        left_root = roots[left_index]
        right_root = roots[right_index]
        left_vertices = leaf_sets[left_root]
        right_vertices = leaf_sets[right_root]
        pairs = _far_cross_pairs(
            left_vertices,
            right_vertices,
            forest.side_length,
            distance_fraction,
        )
        if not pairs:
            continue
        result.append(
            _SnapshotPalmCut(
                left_root=left_root,
                right_root=right_root,
                left_vertices=left_vertices,
                right_vertices=right_vertices,
                cut_edges=tuple(sorted(edges)),
                distant_pairs=pairs,
            )
        )
    return tuple(result)


def _make_corridor_cut(
    *,
    node: int | None,
    is_lca: bool,
    bucket_size: int,
    actual_rank: float,
    p: float,
    critical_rank: float,
) -> CorridorCut:
    favourable_rank = min(actual_rank, critical_rank)
    return CorridorCut(
        node=node,
        is_lca=is_lca,
        bucket_size=bucket_size,
        actual_rank=actual_rank,
        favourable_rank=favourable_rank,
        actual_beta=clock_time_from_rank(actual_rank, p),
        favourable_beta=clock_time_from_rank(favourable_rank, p),
        actual_charge=geometric_charge(bucket_size, p, actual_rank),
        favourable_charge=geometric_charge(bucket_size, p, favourable_rank),
    )


def _verify_pairwise_disjoint(
    buckets: Sequence[Sequence[tuple[int, int]]],
) -> bool:
    seen: set[tuple[int, int]] = set()
    for bucket in buckets:
        current = set(bucket)
        if len(current) != len(bucket) or seen.intersection(current):
            return False
        seen.update(current)
    return True


def critical_cut_intensity_observations(
    critical_forest: CriticalKruskalForest,
    *,
    repetition: int,
    p: float,
    distance_fraction: float,
    rng: random.Random,
) -> tuple[PalmCorridorObservation, ...]:
    """Build the synthetic fixed-level benchmark with explicit ``m*N`` bias.

    ``critical_forest`` must be the snapshot censored at ``Q_CRITICAL``.
    This benchmark does not preserve the final-corridor skeleton and hence
    does not assert a favourable Blackwell domination.
    """

    forest = critical_forest
    if any(
        forest.merge_rank[node] > Q_CRITICAL
        for node in forest.internal_nodes
    ):
        raise ValueError("critical_forest contains a postcritical merger")
    leaf_sets = _leaf_sets(forest)
    observations = []
    for candidate in _snapshot_palm_cuts(forest, distance_fraction):
        first, second = _sample_oriented_pair(candidate.distant_pairs, rng)
        first_root = (
            candidate.left_root
            if first in candidate.left_vertices
            else candidate.right_root
        )
        second_root = (
            candidate.right_root
            if first_root == candidate.left_root
            else candidate.left_root
        )
        first_chain = _internal_ancestor_chain(forest, first)
        second_chain = _internal_ancestor_chain(forest, second)
        if not first_chain or first_chain[-1] != first_root:
            if first != first_root:
                raise AssertionError("first endpoint arm misses its component root")
        if not second_chain or second_chain[-1] != second_root:
            if second != second_root:
                raise AssertionError("second endpoint arm misses its component root")

        path_nodes = first_chain + second_chain
        bucket_edge_sets = [
            _bucket_edges(forest, node, leaf_sets) for node in path_nodes
        ]
        bucket_edge_sets.append(candidate.cut_edges)
        disjoint = _verify_pairwise_disjoint(bucket_edge_sets)
        if not disjoint:
            raise AssertionError("critical-oracle corridor buckets overlap")

        cuts = tuple(
            _make_corridor_cut(
                node=node,
                is_lca=False,
                bucket_size=forest.bucket_size[node],
                actual_rank=forest.merge_rank[node],
                p=p,
                critical_rank=Q_CRITICAL,
            )
            for node in path_nodes
        ) + (
            _make_corridor_cut(
                node=None,
                is_lca=True,
                bucket_size=candidate.bucket_size,
                actual_rank=Q_CRITICAL,
                p=p,
                critical_rank=Q_CRITICAL,
            ),
        )
        observations.append(
            PalmCorridorObservation(
                mode="critical_cut_intensity",
                repetition=repetition,
                first=first,
                second=second,
                distance=triangular_torus_distance(
                    first, second, forest.side_length
                ),
                n_rho=candidate.n_rho,
                lca_bucket_size=candidate.bucket_size,
                palm_weight=candidate.palm_weight,
                palm_weight_definition="m(A,B) * N_rho(A,B)",
                cuts=cuts,
                physical_buckets_pairwise_disjoint=disjoint,
            )
        )
    return tuple(observations)


def final_realized_event_observations(
    forest: CriticalKruskalForest,
    *,
    repetition: int,
    p: float,
    distance_fraction: float,
    rng: random.Random,
    critical_rank: float = Q_CRITICAL,
) -> tuple[PalmCorridorObservation, ...]:
    """Enumerate actual LCA events, using ``N_rho`` rather than ``m*N_rho``."""

    leaf_sets = _leaf_sets(forest)
    observations = []
    for lca in forest.internal_nodes:
        left_vertices = leaf_sets[forest.left_child[lca]]
        right_vertices = leaf_sets[forest.right_child[lca]]
        pairs = _far_cross_pairs(
            left_vertices,
            right_vertices,
            forest.side_length,
            distance_fraction,
        )
        if not pairs:
            continue
        first, second = _sample_oriented_pair(pairs, rng)
        path_nodes = forest.pair_path_nodes(first, second)
        if path_nodes[-1] != lca:
            raise AssertionError("sampled cross-pair has the wrong LCA")

        bucket_edge_sets = [
            _bucket_edges(forest, node, leaf_sets) for node in path_nodes
        ]
        disjoint = _verify_pairwise_disjoint(bucket_edge_sets)
        if not disjoint:
            raise AssertionError("realized-event corridor buckets overlap")
        cuts = tuple(
            _make_corridor_cut(
                node=node,
                is_lca=node == lca,
                bucket_size=forest.bucket_size[node],
                actual_rank=forest.merge_rank[node],
                p=p,
                critical_rank=critical_rank,
            )
            for node in path_nodes
        )
        n_rho = 2 * len(pairs)
        observations.append(
            PalmCorridorObservation(
                mode="final_realized_event",
                repetition=repetition,
                first=first,
                second=second,
                distance=triangular_torus_distance(
                    first, second, forest.side_length
                ),
                n_rho=n_rho,
                lca_bucket_size=forest.bucket_size[lca],
                palm_weight=realized_event_palm_weight(n_rho),
                palm_weight_definition=(
                    "N_rho(A,B); the realized Kruskal event already carries m"
                ),
                cuts=cuts,
                physical_buckets_pairwise_disjoint=disjoint,
            )
        )
    return tuple(observations)


def _weighted_mean(
    observations: Sequence[PalmCorridorObservation],
    statistic: Callable[[PalmCorridorObservation], float],
) -> float:
    total_weight = sum(item.palm_weight for item in observations)
    if total_weight <= 0:
        raise ValueError("positive total Palm weight is required")
    return fsum(item.palm_weight * statistic(item) for item in observations) / (
        total_weight
    )


def _cluster_jackknife_standard_error(
    observations: Sequence[PalmCorridorObservation],
    statistic: Callable[[PalmCorridorObservation], float],
    repetition_count: int,
) -> float | None:
    """Delete-one-environment jackknife error for a weighted ratio.

    Palm observations from one rank environment are highly dependent.  They
    are therefore removed as one cluster.  Repetitions with zero Palm mass
    are retained as zero-contribution clusters.
    """

    if repetition_count < 2:
        return None
    weights = [0] * repetition_count
    numerators = [0.0] * repetition_count
    for item in observations:
        if not 0 <= item.repetition < repetition_count:
            raise ValueError("observation repetition is outside the run")
        weights[item.repetition] += item.palm_weight
        numerators[item.repetition] += item.palm_weight * statistic(item)
    total_weight = sum(weights)
    total_numerator = fsum(numerators)
    leave_one_out = []
    for weight, numerator in zip(weights, numerators, strict=True):
        remaining_weight = total_weight - weight
        if remaining_weight > 0:
            leave_one_out.append(
                (total_numerator - numerator) / remaining_weight
            )
    if len(leave_one_out) < 2:
        return None
    mean_leave_one_out = fsum(leave_one_out) / len(leave_one_out)
    variance = (len(leave_one_out) - 1.0) / len(leave_one_out) * fsum(
        (value - mean_leave_one_out) ** 2 for value in leave_one_out
    )
    return sqrt(max(0.0, variance))


def summarize_observations(
    observations: Sequence[PalmCorridorObservation],
    *,
    maximum_bucket_size: int,
    maximum_charge: float,
    repetition_count: int | None = None,
) -> PalmModeSummary:
    """Aggregate observations with their correct Campbell weights."""

    if not observations:
        raise ValueError("at least one Palm observation is required")
    modes = {item.mode for item in observations}
    definitions = {item.palm_weight_definition for item in observations}
    if len(modes) != 1 or len(definitions) != 1:
        raise ValueError("observations must share one Palm convention")
    if any(item.screening_computed for item in observations):
        raise ValueError("this geometric diagnostic must not claim screening")
    total_weight = sum(item.palm_weight for item in observations)
    squared_weight_sum = sum(item.palm_weight**2 for item in observations)
    effective_count = total_weight * total_weight / squared_weight_sum
    if repetition_count is None:
        repetition_count = max(item.repetition for item in observations) + 1
    if repetition_count <= 0:
        raise ValueError("repetition_count must be positive")
    repetition_weights = [0] * repetition_count
    for item in observations:
        if not 0 <= item.repetition < repetition_count:
            raise ValueError("observation repetition is outside the run")
        repetition_weights[item.repetition] += item.palm_weight
    positive_repetition_weights = [
        weight for weight in repetition_weights if weight > 0
    ]
    effective_repetition_count = total_weight * total_weight / sum(
        weight * weight for weight in positive_repetition_weights
    )

    distance = lambda item: float(item.distance)
    corridor_cut_count = lambda item: float(len(item.cuts))
    lca_bucket_size = lambda item: float(item.lca_bucket_size)
    bucket_size_two_count = lambda item: float(
        sum(cut.bucket_size == 2 for cut in item.cuts)
    )
    bounded_size_count = lambda item: float(
        item.bounded_size_cut_count(maximum_bucket_size)
    )
    actual_proxy_count = lambda item: float(
        item.geometric_proxy_count(
            maximum_bucket_size,
            maximum_charge,
            favourable=False,
        )
    )
    favourable_proxy_count = lambda item: float(
        item.geometric_proxy_count(
            maximum_bucket_size,
            maximum_charge,
            favourable=True,
        )
    )

    return PalmModeSummary(
        mode=next(iter(modes)),
        weight_convention=next(iter(definitions)),
        observation_count=len(observations),
        total_palm_weight=total_weight,
        effective_observation_count=effective_count,
        repetition_count=repetition_count,
        repetitions_with_positive_palm_mass=len(positive_repetition_weights),
        effective_repetition_count=effective_repetition_count,
        weighted_mean_distance=_weighted_mean(observations, distance),
        jackknife_standard_error_distance=_cluster_jackknife_standard_error(
            observations, distance, repetition_count
        ),
        weighted_mean_corridor_cut_count=_weighted_mean(
            observations, corridor_cut_count
        ),
        jackknife_standard_error_corridor_cut_count=(
            _cluster_jackknife_standard_error(
                observations, corridor_cut_count, repetition_count
            )
        ),
        weighted_mean_lca_bucket_size=_weighted_mean(
            observations, lca_bucket_size
        ),
        jackknife_standard_error_lca_bucket_size=(
            _cluster_jackknife_standard_error(
                observations, lca_bucket_size, repetition_count
            )
        ),
        weighted_mean_lca_rank=_weighted_mean(
            observations,
            lambda item: next(
                cut.actual_rank for cut in item.cuts if cut.is_lca
            ),
        ),
        weighted_postcritical_lca_fraction=_weighted_mean(
            observations,
            lambda item: float(
                next(cut.actual_rank for cut in item.cuts if cut.is_lca)
                > next(cut.favourable_rank for cut in item.cuts if cut.is_lca)
            ),
        ),
        weighted_mean_disjoint_bucket_size_two_count=_weighted_mean(
            observations, bucket_size_two_count
        ),
        jackknife_standard_error_disjoint_bucket_size_two_count=(
            _cluster_jackknife_standard_error(
                observations, bucket_size_two_count, repetition_count
            )
        ),
        weighted_mean_disjoint_bounded_size_count=_weighted_mean(
            observations, bounded_size_count
        ),
        jackknife_standard_error_disjoint_bounded_size_count=(
            _cluster_jackknife_standard_error(
                observations, bounded_size_count, repetition_count
            )
        ),
        weighted_mean_disjoint_actual_bounded_charge_count=_weighted_mean(
            observations,
            lambda item: float(
                item.bounded_charge_cut_count(
                    maximum_charge, favourable=False
                )
            ),
        ),
        weighted_mean_disjoint_favourable_bounded_charge_count=_weighted_mean(
            observations,
            lambda item: float(
                item.bounded_charge_cut_count(
                    maximum_charge, favourable=True
                )
            ),
        ),
        weighted_mean_disjoint_actual_geometric_proxy_count=_weighted_mean(
            observations, actual_proxy_count
        ),
        jackknife_standard_error_disjoint_actual_geometric_proxy_count=(
            _cluster_jackknife_standard_error(
                observations, actual_proxy_count, repetition_count
            )
        ),
        weighted_mean_disjoint_favourable_geometric_proxy_count=_weighted_mean(
            observations, favourable_proxy_count
        ),
        jackknife_standard_error_disjoint_favourable_geometric_proxy_count=(
            _cluster_jackknife_standard_error(
                observations, favourable_proxy_count, repetition_count
            )
        ),
        weighted_pairwise_disjoint_fraction=_weighted_mean(
            observations,
            lambda item: float(item.physical_buckets_pairwise_disjoint),
        ),
        interpretation=(
            (
                "synthetic fixed-level intensity benchmark, not a "
                "Blackwell domination; "
                if next(iter(modes)) == "critical_cut_intensity"
                else "geometry-fixed rank criticalization; "
            )
            + "geometric proxy only: B, lateral screening, and replicated "
            "transfer contraction are not computed"
        ),
        uncertainty_note=(
            "standard errors are delete-one-rank-environment jackknives; "
            "effective observation count is descriptive, not an iid claim"
        ),
        screening_computed=False,
        weighted_mean_screened_cut_count=None,
    )


def run_diagnostic(
    *,
    side_length: int,
    repetitions: int,
    p: float,
    distance_fraction: float,
    maximum_bucket_size: int,
    maximum_charge: float,
    seed: int,
) -> LcaPalmDiagnostic:
    """Run both Palm diagnostics on shared reproducible rank environments."""

    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0, 1)")
    if maximum_bucket_size < 2:
        raise ValueError("maximum_bucket_size must be at least 2")
    if maximum_charge < 0.0:
        raise ValueError("maximum_charge must be nonnegative")
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2, 1)")
    final_rank = 2.0 * p - 1.0
    if final_rank < Q_CRITICAL:
        raise ValueError("p is too small for q_c to occur before censoring")

    rng = random.Random(seed)
    critical_observations: list[PalmCorridorObservation] = []
    final_observations: list[PalmCorridorObservation] = []
    audited_ranks = {0.0, Q_CRITICAL, final_rank}
    realized_event_pair_total = 0
    connected_pair_total = 0
    maximum_pair_partition_difference = 0
    for repetition in range(repetitions):
        ranked_edges: tuple[RankedEdge, ...] = sample_ranked_edges(
            side_length, rng
        )
        critical_forest = CriticalKruskalForest(
            side_length, ranked_edges, critical_rank=Q_CRITICAL
        )
        final_forest = CriticalKruskalForest(
            side_length, ranked_edges, critical_rank=final_rank
        )
        critical_observations.extend(
            critical_cut_intensity_observations(
                critical_forest,
                repetition=repetition,
                p=p,
                distance_fraction=distance_fraction,
                rng=rng,
            )
        )
        final_observations.extend(
            final_realized_event_observations(
                final_forest,
                repetition=repetition,
                p=p,
                distance_fraction=distance_fraction,
                rng=rng,
            )
        )
        event_pairs, connected_pairs = lca_pair_partition_counts(
            final_forest, distance_fraction
        )
        realized_event_pair_total += event_pairs
        connected_pair_total += connected_pairs
        maximum_pair_partition_difference = max(
            maximum_pair_partition_difference,
            abs(event_pairs - connected_pairs),
        )
        audited_ranks.update(
            final_forest.merge_rank[node] for node in final_forest.internal_nodes
        )

    if not critical_observations:
        raise RuntimeError("no distant critical cut-intensity observation")
    if not final_observations:
        raise RuntimeError("no distant realized-event observation")
    if (
        maximum_pair_partition_difference != 0
        or realized_event_pair_total != connected_pair_total
        or sum(item.palm_weight for item in final_observations)
        != realized_event_pair_total
    ):
        raise AssertionError("realized LCA events do not partition connected pairs")

    campbell_audit = campbell_weight_counter_audit()
    if not campbell_audit.representations_agree:
        raise AssertionError("the two Campbell representations disagree")
    if not campbell_audit.double_bias_detected:
        raise AssertionError("the deterministic audit missed the m^2 bias")
    charge_audit = charge_identity_counter_audit(p, sorted(audited_ranks))
    if not charge_audit.passed:
        raise AssertionError("rank and beta charge calculations disagree")

    return LcaPalmDiagnostic(
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        critical_rank=Q_CRITICAL,
        final_rank=final_rank,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_charge=maximum_charge,
        seed=seed,
        critical_snapshot_benchmark=summarize_observations(
            critical_observations,
            maximum_bucket_size=maximum_bucket_size,
            maximum_charge=maximum_charge,
            repetition_count=repetitions,
        ),
        final_corridor=summarize_observations(
            final_observations,
            maximum_bucket_size=maximum_bucket_size,
            maximum_charge=maximum_charge,
            repetition_count=repetitions,
        ),
        campbell_weight_audit=campbell_audit,
        charge_identity_audit=charge_audit,
        lca_pair_partition_audit=LcaPairPartitionAudit(
            checked_environment_count=repetitions,
            realized_event_pair_total=realized_event_pair_total,
            connected_pair_total=connected_pair_total,
            maximum_environment_difference=maximum_pair_partition_difference,
            passed=(
                maximum_pair_partition_difference == 0
                and realized_event_pair_total == connected_pair_total
            ),
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=16)
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-bucket-size", type=int, default=8)
    parser.add_argument("--maximum-charge", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=20260719)
    arguments = parser.parse_args()
    diagnostic = run_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        maximum_bucket_size=arguments.maximum_bucket_size,
        maximum_charge=arguments.maximum_charge,
        seed=arguments.seed,
    )
    print(json.dumps(asdict(diagnostic), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
