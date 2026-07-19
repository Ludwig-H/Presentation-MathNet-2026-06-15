"""Geometry-only T2 signatures on final triangular-torus corridors.

The diagnostic starts from the *real* Kruskal forest censored at
``q_1 = 2p-1``.  A realized merger ``v`` is an LCA-Palm event with weight
``N_rho(v)`` only.  The cut size ``m(v)`` is already present in the event
rate and is never multiplied into that Palm weight a second time.  One
uniform distant cross-pair is sampled at every eligible LCA; its two tree
arms define the corridor whose buckets are inspected.

Every local signature is reconstructed without edge marks.  It records the
physical cut size, merger rank, favourable charge after rank capping, child
sizes, an oriented comb attachment, the physical boundary of the merged
union, and the number of neighbouring rank-components touching that
boundary immediately after the merger.  The latter number is called a
``port_count`` only as a geometric proxy.  It does not establish lateral
screening, control an exterior message, or certify a T2 transfer deficit.

The reconstruction is counter-audited in finite volume.  Physical cut sizes
are compared with the Kruskal buckets; the winning rank is recovered as the
minimum rank across the child cut; external boundaries are checked by the
degree-sum identity; and port partitions obtained from a rank DSU are
compared with the active Kruskal-tree partition.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, sqrt
from typing import Callable, Iterable, Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    RankedEdge,
    sample_ranked_edges,
    triangular_torus_distance,
)
from lca_palm_corridor_diagnostic import (
    geometric_charge,
    lca_pair_partition_counts,
    realized_event_palm_weight,
)

Edge = tuple[int, int]


@dataclass(frozen=True)
class CorridorT2Signature:
    """A mark-free local signature for one bucket on a pair corridor.

    At a strict arm node, ``spine_child_size`` is the child containing the
    relevant endpoint and ``attachment_child_size`` is its off-spine
    sibling.  At the LCA no child is an off-spine attachment, so both fields
    are ``None``.  This avoids falsely calling the smaller LCA arm, or the
    endpoint leaf itself, a comb attachment.
    """

    node: int
    is_lca: bool
    bucket_size: int
    rank: float
    favourable_rank: float
    favourable_charge: float
    left_child_size: int
    right_child_size: int
    small_child_size: int
    large_child_size: int
    small_child_fraction: float
    child_asymmetry: float
    spine_child_size: int | None
    attachment_child_size: int | None
    attachment_fraction: float | None
    attachment_is_smaller_child: bool | None
    external_physical_edge_count: int
    port_count: int
    ports_are_geometric_proxy: bool = True
    screening_computed: bool = False

    def is_base_candidate(
        self,
        *,
        maximum_bucket_size: int,
        maximum_charge: float,
    ) -> bool:
        """Whether the strict-arm bucket passes only the ``m,J`` filters."""

        return (
            not self.is_lca
            and 2 <= self.bucket_size <= maximum_bucket_size
            and self.favourable_charge <= maximum_charge
        )

    def is_small_attachment_candidate(
        self,
        *,
        maximum_bucket_size: int,
        maximum_charge: float,
        maximum_attachment_size: int,
    ) -> bool:
        """Whether the base candidate has a small *oriented* attachment."""

        return (
            self.is_base_candidate(
                maximum_bucket_size=maximum_bucket_size,
                maximum_charge=maximum_charge,
            )
            and self.attachment_child_size is not None
            and self.attachment_child_size <= maximum_attachment_size
        )

    def is_low_port_candidate(
        self,
        *,
        maximum_bucket_size: int,
        maximum_charge: float,
        maximum_ports: int,
    ) -> bool:
        """Whether the base candidate has at most the proxy port cutoff."""

        return (
            self.is_base_candidate(
                maximum_bucket_size=maximum_bucket_size,
                maximum_charge=maximum_charge,
            )
            and self.port_count <= maximum_ports
        )

    def is_t2_proxy_candidate(
        self,
        *,
        maximum_bucket_size: int,
        maximum_charge: float,
        maximum_ports: int,
        maximum_attachment_size: int,
    ) -> bool:
        """Apply all requested geometric filters, without claiming screening."""

        return (
            self.is_small_attachment_candidate(
                maximum_bucket_size=maximum_bucket_size,
                maximum_charge=maximum_charge,
                maximum_attachment_size=maximum_attachment_size,
            )
            and self.port_count <= maximum_ports
        )


@dataclass(frozen=True)
class EventPalmCorridor:
    """One uniformly sampled pair corridor attached to a realized LCA."""

    repetition: int
    lca_node: int
    first: int
    second: int
    distance: int
    n_rho: int
    palm_weight: int
    palm_weight_definition: str
    signatures: tuple[CorridorT2Signature, ...]
    one_uniform_pair_sampled_within_event: bool = True
    screening_computed: bool = False

    def candidate_count(
        self,
        *,
        maximum_bucket_size: int,
        maximum_charge: float,
        maximum_ports: int,
        maximum_attachment_size: int,
    ) -> int:
        """Count full T2 geometry proxies on this sampled corridor."""

        return sum(
            signature.is_t2_proxy_candidate(
                maximum_bucket_size=maximum_bucket_size,
                maximum_charge=maximum_charge,
                maximum_ports=maximum_ports,
                maximum_attachment_size=maximum_attachment_size,
            )
            for signature in self.signatures
        )


@dataclass(frozen=True)
class ReconstructionAudit:
    """Independent identities checked while rebuilding local signatures."""

    checked_internal_node_count: int
    event_mapping_mismatch_count: int
    child_component_mismatch_count: int
    union_component_mismatch_count: int
    bucket_size_mismatch_count: int
    winning_rank_mismatch_count: int
    external_boundary_identity_mismatch_count: int
    port_partition_mismatch_count: int
    maximum_winning_rank_error: float
    passed: bool


@dataclass(frozen=True)
class PalmPartitionAudit:
    """Exact event-``N_rho`` partition identity over rank environments."""

    checked_environment_count: int
    event_palm_weight_total: int
    connected_distant_ordered_pair_total: int
    maximum_environment_difference: int
    passed: bool


@dataclass(frozen=True)
class T2SignatureSummary:
    """Event-Palm weighted counts with delete-one-environment uncertainty."""

    side_length: int
    repetitions: int
    p: float
    critical_rank: float
    final_rank: float
    distance_fraction: float
    maximum_bucket_size: int
    maximum_charge: float
    maximum_ports: int
    maximum_attachment_size: int
    seed: int
    observation_count: int
    total_palm_weight: int
    effective_observation_count: float
    effective_environment_count: float
    weighted_mean_corridor_bucket_count: float
    jackknife_standard_error_corridor_bucket_count: float | None
    weighted_mean_base_candidate_count: float
    jackknife_standard_error_base_candidate_count: float | None
    weighted_mean_small_attachment_candidate_count: float
    jackknife_standard_error_small_attachment_candidate_count: float | None
    weighted_mean_low_port_candidate_count: float
    jackknife_standard_error_low_port_candidate_count: float | None
    weighted_mean_t2_proxy_candidate_count: float
    jackknife_standard_error_t2_proxy_candidate_count: float | None
    weighted_mean_t2_proxy_density: float
    jackknife_standard_error_t2_proxy_density: float | None
    weighted_fraction_with_t2_proxy_candidate: float
    jackknife_standard_error_fraction_with_t2_proxy_candidate: float | None
    weighted_mean_external_edge_count_per_bucket: float
    weighted_mean_port_count_per_bucket: float
    weighted_base_candidate_count_by_bucket_size: dict[str, float]
    weighted_t2_proxy_count_by_port_count: dict[str, float]
    palm_weight_convention: str
    interpretation: str
    uncertainty_note: str
    ports_are_geometric_proxy: bool
    screening_computed: bool
    reconstruction_audit: ReconstructionAudit
    palm_partition_audit: PalmPartitionAudit


@dataclass(frozen=True)
class _NodeGeometry:
    """Reconstructed node data, including leaves used to orient a corridor."""

    node: int
    rank: float
    bucket_size: int
    left_vertices: frozenset[int]
    right_vertices: frozenset[int]
    union_vertices: frozenset[int]
    external_physical_edge_count: int
    port_count: int


class _DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, vertex: int) -> int:
        parent = self.parent
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(self, first: int, second: int) -> bool:
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root == second_root:
            return False
        if self.size[first_root] < self.size[second_root]:
            first_root, second_root = second_root, first_root
        self.parent[second_root] = first_root
        self.size[first_root] += self.size[second_root]
        return True


def reconstruct_leaf_sets(
    forest: CriticalKruskalForest,
) -> dict[int, frozenset[int]]:
    """Rebuild every merger's leaf set from the immutable tree arrays."""

    leaves: dict[int, frozenset[int]] = {
        vertex: frozenset((vertex,)) for vertex in range(forest.vertex_count)
    }
    for node in forest.internal_nodes:
        left = forest.left_child[node]
        right = forest.right_child[node]
        if left not in leaves or right not in leaves:
            raise AssertionError("Kruskal nodes are not topologically ordered")
        if leaves[left].intersection(leaves[right]):
            raise AssertionError("the two Kruskal children overlap")
        leaves[node] = leaves[left] | leaves[right]
    return leaves


def _normalised_edge_ranks(
    forest: CriticalKruskalForest,
    ranked_edges: Sequence[RankedEdge],
) -> tuple[dict[Edge, float], tuple[RankedEdge, ...]]:
    expected_edges = set(forest.edges)
    edge_ranks: dict[Edge, float] = {}
    seen_ranks: set[float] = set()
    normalised = []
    for rank, first, second in ranked_edges:
        edge = tuple(sorted((first, second)))
        if edge not in expected_edges:
            raise ValueError("ranked_edges contains a non-torus edge")
        if edge in edge_ranks:
            raise ValueError("every torus edge must occur exactly once")
        if not 0.0 <= rank <= 1.0:
            raise ValueError("edge ranks must belong to [0,1]")
        if rank in seen_ranks:
            raise ValueError(
                "distinct ranks are required to define an event-time port state"
            )
        seen_ranks.add(rank)
        edge_ranks[edge] = rank
        normalised.append((rank, edge[0], edge[1]))
    if set(edge_ranks) != expected_edges:
        raise ValueError("one rank is required for every torus edge")
    return edge_ranks, tuple(sorted(normalised))


def _component_vertices(
    dsu: _DisjointSet, vertex_count: int, root: int
) -> frozenset[int]:
    canonical_root = dsu.find(root)
    return frozenset(
        vertex for vertex in range(vertex_count) if dsu.find(vertex) == canonical_root
    )


def _active_tree_root(
    forest: CriticalKruskalForest,
    vertex: int,
    rank: float,
) -> int:
    node = vertex
    while forest.tree_parent[node] != -1:
        parent = forest.tree_parent[node]
        if forest.merge_rank[parent] > rank:
            break
        node = parent
    return node


def _canonical_partition(
    vertices: Iterable[int],
    label: Callable[[int], int],
) -> tuple[tuple[int, ...], ...]:
    groups: dict[int, set[int]] = {}
    for vertex in vertices:
        groups.setdefault(label(vertex), set()).add(vertex)
    return tuple(sorted(tuple(sorted(group)) for group in groups.values()))


def reconstruct_node_geometries(
    forest: CriticalKruskalForest,
    ranked_edges: Sequence[RankedEdge],
    *,
    censor_rank: float,
) -> tuple[dict[int, _NodeGeometry], ReconstructionAudit]:
    """Reconstruct all local event states and return independent audits.

    Ports are the distinct rank-components, immediately after the winning
    edge, which contain an outside endpoint of a physical boundary edge.
    The returned audit compares their full partition on boundary vertices,
    not merely the number of parts.
    """

    if not 0.0 <= censor_rank <= 1.0:
        raise ValueError("censor_rank must belong to [0,1]")
    edge_ranks, ordered_edges = _normalised_edge_ranks(forest, ranked_edges)
    leaves = reconstruct_leaf_sets(forest)
    ranks_to_nodes: dict[float, int] = {}
    for node in forest.internal_nodes:
        rank = forest.merge_rank[node]
        if rank > censor_rank + 1e-14:
            raise ValueError("the forest contains a merger after censor_rank")
        if rank in ranks_to_nodes:
            raise ValueError("distinct merger ranks are required")
        ranks_to_nodes[rank] = node

    degrees = [0] * forest.vertex_count
    for first, second in forest.edges:
        degrees[first] += 1
        degrees[second] += 1

    event_mapping_mismatches = 0
    child_component_mismatches = 0
    union_component_mismatches = 0
    bucket_size_mismatches = 0
    winning_rank_mismatches = 0
    boundary_identity_mismatches = 0
    port_partition_mismatches = 0
    maximum_winning_rank_error = 0.0
    geometries: dict[int, _NodeGeometry] = {}
    dsu = _DisjointSet(forest.vertex_count)

    for rank, first, second in ordered_edges:
        if rank > censor_rank:
            break
        node = ranks_to_nodes.get(rank)
        accepted = dsu.find(first) != dsu.find(second)
        if accepted != (node is not None):
            event_mapping_mismatches += 1

        if node is None:
            dsu.union(first, second)
            continue

        left_vertices = leaves[forest.left_child[node]]
        right_vertices = leaves[forest.right_child[node]]
        before_first = _component_vertices(dsu, forest.vertex_count, first)
        before_second = _component_vertices(dsu, forest.vertex_count, second)
        expected_children = {left_vertices, right_vertices}
        actual_children = {before_first, before_second}
        if expected_children != actual_children:
            child_component_mismatches += 1

        child_cut_edges = tuple(
            edge
            for edge in forest.edges
            if (edge[0] in left_vertices and edge[1] in right_vertices)
            or (edge[1] in left_vertices and edge[0] in right_vertices)
        )
        reconstructed_bucket_size = len(child_cut_edges)
        if reconstructed_bucket_size != forest.bucket_size[node]:
            bucket_size_mismatches += 1
        winning_rank = min(edge_ranks[edge] for edge in child_cut_edges)
        winning_error = abs(winning_rank - rank)
        maximum_winning_rank_error = max(maximum_winning_rank_error, winning_error)
        if winning_error > 1e-14:
            winning_rank_mismatches += 1

        dsu.union(first, second)
        union_vertices = left_vertices | right_vertices
        after_component = _component_vertices(dsu, forest.vertex_count, first)
        if after_component != union_vertices:
            union_component_mismatches += 1

        external_edges = tuple(
            edge
            for edge in forest.edges
            if (edge[0] in union_vertices) != (edge[1] in union_vertices)
        )
        internal_edge_count = sum(
            edge[0] in union_vertices and edge[1] in union_vertices
            for edge in forest.edges
        )
        degree_boundary_count = (
            sum(degrees[vertex] for vertex in union_vertices) - 2 * internal_edge_count
        )
        if len(external_edges) != degree_boundary_count:
            boundary_identity_mismatches += 1

        outside_boundary_vertices = {
            edge[1] if edge[0] in union_vertices else edge[0] for edge in external_edges
        }
        dsu_partition = _canonical_partition(outside_boundary_vertices, dsu.find)
        tree_partition = _canonical_partition(
            outside_boundary_vertices,
            lambda vertex: _active_tree_root(forest, vertex, rank),
        )
        if dsu_partition != tree_partition:
            port_partition_mismatches += 1

        geometries[node] = _NodeGeometry(
            node=node,
            rank=rank,
            bucket_size=reconstructed_bucket_size,
            left_vertices=left_vertices,
            right_vertices=right_vertices,
            union_vertices=union_vertices,
            external_physical_edge_count=len(external_edges),
            port_count=len(dsu_partition),
        )

    checked_node_count = len(tuple(forest.internal_nodes))
    if len(geometries) != checked_node_count:
        event_mapping_mismatches += abs(checked_node_count - len(geometries))
    mismatch_total = (
        event_mapping_mismatches
        + child_component_mismatches
        + union_component_mismatches
        + bucket_size_mismatches
        + winning_rank_mismatches
        + boundary_identity_mismatches
        + port_partition_mismatches
    )
    return geometries, ReconstructionAudit(
        checked_internal_node_count=checked_node_count,
        event_mapping_mismatch_count=event_mapping_mismatches,
        child_component_mismatch_count=child_component_mismatches,
        union_component_mismatch_count=union_component_mismatches,
        bucket_size_mismatch_count=bucket_size_mismatches,
        winning_rank_mismatch_count=winning_rank_mismatches,
        external_boundary_identity_mismatch_count=(boundary_identity_mismatches),
        port_partition_mismatch_count=port_partition_mismatches,
        maximum_winning_rank_error=maximum_winning_rank_error,
        passed=mismatch_total == 0,
    )


def _count_and_sample_distant_cross_pair(
    left: Iterable[int],
    right: Iterable[int],
    *,
    side_length: int,
    distance_fraction: float,
    rng: random.Random,
) -> tuple[int, tuple[int, int] | None]:
    """Count exactly and reservoir-sample one eligible unordered pair."""

    minimum_distance = distance_fraction * side_length
    count = 0
    sampled: tuple[int, int] | None = None
    for first in sorted(left):
        for second in sorted(right):
            if triangular_torus_distance(first, second, side_length) < minimum_distance:
                continue
            count += 1
            if rng.randrange(count) == 0:
                sampled = (first, second)
    if sampled is not None and rng.randrange(2):
        sampled = (sampled[1], sampled[0])
    return 2 * count, sampled


def _make_signature(
    geometry: _NodeGeometry,
    *,
    first: int,
    second: int,
    lca: int,
    p: float,
    critical_rank: float,
) -> CorridorT2Signature:
    left_size = len(geometry.left_vertices)
    right_size = len(geometry.right_vertices)
    small_size = min(left_size, right_size)
    large_size = max(left_size, right_size)
    union_size = left_size + right_size

    spine_size: int | None = None
    attachment_size: int | None = None
    if geometry.node != lca:
        endpoint_membership = (
            first in geometry.union_vertices,
            second in geometry.union_vertices,
        )
        if sum(endpoint_membership) != 1:
            raise AssertionError(
                "a strict corridor arm must contain exactly one endpoint"
            )
        endpoint = first if endpoint_membership[0] else second
        if endpoint in geometry.left_vertices:
            spine_size, attachment_size = left_size, right_size
        elif endpoint in geometry.right_vertices:
            spine_size, attachment_size = right_size, left_size
        else:
            raise AssertionError("the arm endpoint is absent from both children")

    favourable_rank = min(geometry.rank, critical_rank)
    return CorridorT2Signature(
        node=geometry.node,
        is_lca=geometry.node == lca,
        bucket_size=geometry.bucket_size,
        rank=geometry.rank,
        favourable_rank=favourable_rank,
        favourable_charge=geometric_charge(geometry.bucket_size, p, favourable_rank),
        left_child_size=left_size,
        right_child_size=right_size,
        small_child_size=small_size,
        large_child_size=large_size,
        small_child_fraction=small_size / union_size,
        child_asymmetry=(large_size - small_size) / union_size,
        spine_child_size=spine_size,
        attachment_child_size=attachment_size,
        attachment_fraction=(
            None if attachment_size is None else attachment_size / union_size
        ),
        attachment_is_smaller_child=(
            None if attachment_size is None else attachment_size == small_size
        ),
        external_physical_edge_count=(geometry.external_physical_edge_count),
        port_count=geometry.port_count,
    )


def event_palm_corridor_observations(
    forest: CriticalKruskalForest,
    ranked_edges: Sequence[RankedEdge],
    *,
    repetition: int,
    p: float,
    distance_fraction: float,
    rng: random.Random,
    critical_rank: float = Q_CRITICAL,
) -> tuple[tuple[EventPalmCorridor, ...], ReconstructionAudit]:
    """Build final corridors with realized-event weight ``N_rho`` only."""

    if repetition < 0:
        raise ValueError("repetition must be nonnegative")
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2,1)")
    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0,1)")
    final_rank = 2.0 * p - 1.0
    if final_rank < critical_rank:
        raise ValueError("critical_rank must not exceed the final rank")
    node_geometries, audit = reconstruct_node_geometries(
        forest, ranked_edges, censor_rank=final_rank
    )
    if not audit.passed:
        return (), audit

    observations = []
    for lca in forest.internal_nodes:
        geometry = node_geometries[lca]
        n_rho, sampled_pair = _count_and_sample_distant_cross_pair(
            geometry.left_vertices,
            geometry.right_vertices,
            side_length=forest.side_length,
            distance_fraction=distance_fraction,
            rng=rng,
        )
        if sampled_pair is None:
            continue
        first, second = sampled_pair
        path_nodes = forest.pair_path_nodes(first, second)
        if path_nodes[-1] != lca:
            raise AssertionError("the sampled pair has the wrong LCA")
        signatures = tuple(
            _make_signature(
                node_geometries[node],
                first=first,
                second=second,
                lca=lca,
                p=p,
                critical_rank=critical_rank,
            )
            for node in path_nodes
        )
        observations.append(
            EventPalmCorridor(
                repetition=repetition,
                lca_node=lca,
                first=first,
                second=second,
                distance=triangular_torus_distance(first, second, forest.side_length),
                n_rho=n_rho,
                palm_weight=realized_event_palm_weight(n_rho),
                palm_weight_definition=(
                    "N_rho only; the realized Kruskal event already carries m"
                ),
                signatures=signatures,
            )
        )
    return tuple(observations), audit


def _weighted_mean(
    observations: Sequence[EventPalmCorridor],
    statistic: Callable[[EventPalmCorridor], float],
) -> float:
    total_weight = sum(item.palm_weight for item in observations)
    if total_weight <= 0:
        raise ValueError("positive Palm mass is required")
    return (
        fsum(item.palm_weight * statistic(item) for item in observations) / total_weight
    )


def _cluster_jackknife_standard_error(
    observations: Sequence[EventPalmCorridor],
    statistic: Callable[[EventPalmCorridor], float],
    repetition_count: int,
) -> float | None:
    """Delete one complete rank environment from a weighted ratio."""

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
            leave_one_out.append((total_numerator - numerator) / remaining_weight)
    if len(leave_one_out) < 2:
        return None
    leave_one_out_mean = fsum(leave_one_out) / len(leave_one_out)
    variance = (
        (len(leave_one_out) - 1.0)
        / len(leave_one_out)
        * fsum((value - leave_one_out_mean) ** 2 for value in leave_one_out)
    )
    return sqrt(max(0.0, variance))


def _aggregate_reconstruction_audits(
    audits: Sequence[ReconstructionAudit],
) -> ReconstructionAudit:
    if not audits:
        raise ValueError("at least one reconstruction audit is required")
    return ReconstructionAudit(
        checked_internal_node_count=sum(
            audit.checked_internal_node_count for audit in audits
        ),
        event_mapping_mismatch_count=sum(
            audit.event_mapping_mismatch_count for audit in audits
        ),
        child_component_mismatch_count=sum(
            audit.child_component_mismatch_count for audit in audits
        ),
        union_component_mismatch_count=sum(
            audit.union_component_mismatch_count for audit in audits
        ),
        bucket_size_mismatch_count=sum(
            audit.bucket_size_mismatch_count for audit in audits
        ),
        winning_rank_mismatch_count=sum(
            audit.winning_rank_mismatch_count for audit in audits
        ),
        external_boundary_identity_mismatch_count=sum(
            audit.external_boundary_identity_mismatch_count for audit in audits
        ),
        port_partition_mismatch_count=sum(
            audit.port_partition_mismatch_count for audit in audits
        ),
        maximum_winning_rank_error=max(
            audit.maximum_winning_rank_error for audit in audits
        ),
        passed=all(audit.passed for audit in audits),
    )


def run_diagnostic(
    *,
    side_length: int,
    repetitions: int,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 6,
    maximum_charge: float = 1.0,
    maximum_ports: int = 6,
    maximum_attachment_size: int = 4,
    seed: int = 20260719,
) -> T2SignatureSummary:
    """Run the event-Palm signature diagnostic on shared rank environments."""

    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if maximum_bucket_size < 2:
        raise ValueError("maximum_bucket_size must be at least two")
    if maximum_charge < 0.0:
        raise ValueError("maximum_charge must be nonnegative")
    if maximum_ports < 0:
        raise ValueError("maximum_ports must be nonnegative")
    if maximum_attachment_size <= 0:
        raise ValueError("maximum_attachment_size must be positive")
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2,1)")
    final_rank = 2.0 * p - 1.0
    if final_rank < Q_CRITICAL:
        raise ValueError("p is too small for favourable criticalization")

    rng = random.Random(seed)
    observations: list[EventPalmCorridor] = []
    audits = []
    event_palm_weight_total = 0
    connected_pair_total = 0
    maximum_environment_difference = 0
    for repetition in range(repetitions):
        ranked_edges = sample_ranked_edges(side_length, rng)
        forest = CriticalKruskalForest(
            side_length, ranked_edges, critical_rank=final_rank
        )
        environment_observations, audit = event_palm_corridor_observations(
            forest,
            ranked_edges,
            repetition=repetition,
            p=p,
            distance_fraction=distance_fraction,
            rng=rng,
        )
        if not audit.passed:
            raise AssertionError("local Kruskal reconstruction audit failed")
        audits.append(audit)
        environment_weight = sum(item.palm_weight for item in environment_observations)
        event_pairs, connected_pairs = lca_pair_partition_counts(
            forest, distance_fraction
        )
        environment_difference = max(
            abs(environment_weight - event_pairs),
            abs(event_pairs - connected_pairs),
        )
        maximum_environment_difference = max(
            maximum_environment_difference, environment_difference
        )
        event_palm_weight_total += environment_weight
        connected_pair_total += connected_pairs
        observations.extend(environment_observations)

    if not observations:
        raise RuntimeError("no distant realized-event Palm observation")
    reconstruction_audit = _aggregate_reconstruction_audits(audits)
    palm_partition_audit = PalmPartitionAudit(
        checked_environment_count=repetitions,
        event_palm_weight_total=event_palm_weight_total,
        connected_distant_ordered_pair_total=connected_pair_total,
        maximum_environment_difference=maximum_environment_difference,
        passed=(
            maximum_environment_difference == 0
            and event_palm_weight_total == connected_pair_total
        ),
    )
    if not palm_partition_audit.passed:
        raise AssertionError("event N_rho weights do not partition pair mass")

    def base_count(item: EventPalmCorridor) -> float:
        return float(
            sum(
                signature.is_base_candidate(
                    maximum_bucket_size=maximum_bucket_size,
                    maximum_charge=maximum_charge,
                )
                for signature in item.signatures
            )
        )

    def attachment_count(item: EventPalmCorridor) -> float:
        return float(
            sum(
                signature.is_small_attachment_candidate(
                    maximum_bucket_size=maximum_bucket_size,
                    maximum_charge=maximum_charge,
                    maximum_attachment_size=maximum_attachment_size,
                )
                for signature in item.signatures
            )
        )

    def low_port_count(item: EventPalmCorridor) -> float:
        return float(
            sum(
                signature.is_low_port_candidate(
                    maximum_bucket_size=maximum_bucket_size,
                    maximum_charge=maximum_charge,
                    maximum_ports=maximum_ports,
                )
                for signature in item.signatures
            )
        )

    candidate_arguments = dict(
        maximum_bucket_size=maximum_bucket_size,
        maximum_charge=maximum_charge,
        maximum_ports=maximum_ports,
        maximum_attachment_size=maximum_attachment_size,
    )

    def t2_count(item: EventPalmCorridor) -> float:
        return float(item.candidate_count(**candidate_arguments))

    def t2_density(item: EventPalmCorridor) -> float:
        return t2_count(item) / len(item.signatures)

    def has_t2(item: EventPalmCorridor) -> float:
        return float(t2_count(item) > 0.0)

    def mean_external_edges(item: EventPalmCorridor) -> float:
        return fsum(
            signature.external_physical_edge_count for signature in item.signatures
        ) / len(item.signatures)

    def mean_ports(item: EventPalmCorridor) -> float:
        return fsum(signature.port_count for signature in item.signatures) / len(
            item.signatures
        )

    total_weight = sum(item.palm_weight for item in observations)
    squared_weight_sum = sum(item.palm_weight**2 for item in observations)
    environment_weights = [0] * repetitions
    for item in observations:
        environment_weights[item.repetition] += item.palm_weight
    positive_environment_weights = [
        weight for weight in environment_weights if weight > 0
    ]
    effective_environment_count = (
        total_weight
        * total_weight
        / sum(weight * weight for weight in positive_environment_weights)
    )

    bucket_histogram = {
        str(size): _weighted_mean(
            observations,
            lambda item, selected_size=size: float(
                sum(
                    signature.bucket_size == selected_size
                    and signature.is_base_candidate(
                        maximum_bucket_size=maximum_bucket_size,
                        maximum_charge=maximum_charge,
                    )
                    for signature in item.signatures
                )
            ),
        )
        for size in range(2, maximum_bucket_size + 1)
    }
    port_histogram = {
        str(port_count): _weighted_mean(
            observations,
            lambda item, selected_port_count=port_count: float(
                sum(
                    signature.port_count == selected_port_count
                    and signature.is_t2_proxy_candidate(**candidate_arguments)
                    for signature in item.signatures
                )
            ),
        )
        for port_count in range(maximum_ports + 1)
    }

    return T2SignatureSummary(
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        critical_rank=Q_CRITICAL,
        final_rank=final_rank,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_charge=maximum_charge,
        maximum_ports=maximum_ports,
        maximum_attachment_size=maximum_attachment_size,
        seed=seed,
        observation_count=len(observations),
        total_palm_weight=total_weight,
        effective_observation_count=(total_weight * total_weight / squared_weight_sum),
        effective_environment_count=effective_environment_count,
        weighted_mean_corridor_bucket_count=_weighted_mean(
            observations, lambda item: float(len(item.signatures))
        ),
        jackknife_standard_error_corridor_bucket_count=(
            _cluster_jackknife_standard_error(
                observations,
                lambda item: float(len(item.signatures)),
                repetitions,
            )
        ),
        weighted_mean_base_candidate_count=_weighted_mean(observations, base_count),
        jackknife_standard_error_base_candidate_count=(
            _cluster_jackknife_standard_error(observations, base_count, repetitions)
        ),
        weighted_mean_small_attachment_candidate_count=_weighted_mean(
            observations, attachment_count
        ),
        jackknife_standard_error_small_attachment_candidate_count=(
            _cluster_jackknife_standard_error(
                observations, attachment_count, repetitions
            )
        ),
        weighted_mean_low_port_candidate_count=_weighted_mean(
            observations, low_port_count
        ),
        jackknife_standard_error_low_port_candidate_count=(
            _cluster_jackknife_standard_error(observations, low_port_count, repetitions)
        ),
        weighted_mean_t2_proxy_candidate_count=_weighted_mean(observations, t2_count),
        jackknife_standard_error_t2_proxy_candidate_count=(
            _cluster_jackknife_standard_error(observations, t2_count, repetitions)
        ),
        weighted_mean_t2_proxy_density=_weighted_mean(observations, t2_density),
        jackknife_standard_error_t2_proxy_density=(
            _cluster_jackknife_standard_error(observations, t2_density, repetitions)
        ),
        weighted_fraction_with_t2_proxy_candidate=_weighted_mean(observations, has_t2),
        jackknife_standard_error_fraction_with_t2_proxy_candidate=(
            _cluster_jackknife_standard_error(observations, has_t2, repetitions)
        ),
        weighted_mean_external_edge_count_per_bucket=_weighted_mean(
            observations, mean_external_edges
        ),
        weighted_mean_port_count_per_bucket=_weighted_mean(observations, mean_ports),
        weighted_base_candidate_count_by_bucket_size=bucket_histogram,
        weighted_t2_proxy_count_by_port_count=port_histogram,
        palm_weight_convention=(
            "realized event weighted by N_rho only; one uniform pair "
            "sampled within each event"
        ),
        interpretation=(
            "mark-free local geometry diagnostic for the real final "
            "corridor; port_count and a small off-spine attachment are T2 "
            "proxies only, not screening or a transfer certificate"
        ),
        uncertainty_note=(
            "standard errors are delete-one-rank-environment jackknives; "
            "LCA observations are not treated as iid"
        ),
        ports_are_geometric_proxy=True,
        screening_computed=False,
        reconstruction_audit=reconstruction_audit,
        palm_partition_audit=palm_partition_audit,
    )


def run_size_series(
    side_lengths: Sequence[int],
    *,
    repetitions: int,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 6,
    maximum_charge: float = 1.0,
    maximum_ports: int = 6,
    maximum_attachment_size: int = 4,
    seed: int = 20260719,
) -> tuple[T2SignatureSummary, ...]:
    """Run modest sizes with deterministic, distinct random streams."""

    sizes = tuple(side_lengths)
    if not sizes:
        raise ValueError("at least one side length is required")
    return tuple(
        run_diagnostic(
            side_length=side_length,
            repetitions=repetitions,
            p=p,
            distance_fraction=distance_fraction,
            maximum_bucket_size=maximum_bucket_size,
            maximum_charge=maximum_charge,
            maximum_ports=maximum_ports,
            maximum_attachment_size=maximum_attachment_size,
            seed=seed + 1_000_003 * index,
        )
        for index, side_length in enumerate(sizes)
    )


def _parse_side_lengths(raw: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in raw.split(","))
    if not values:
        raise ValueError("at least one side length is required")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sides", default="8,12,16")
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-bucket-size", type=int, default=6)
    parser.add_argument("--maximum-charge", type=float, default=1.0)
    parser.add_argument("--maximum-ports", type=int, default=6)
    parser.add_argument("--maximum-attachment-size", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260719)
    arguments = parser.parse_args()
    result = run_size_series(
        _parse_side_lengths(arguments.sides),
        repetitions=arguments.repetitions,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        maximum_bucket_size=arguments.maximum_bucket_size,
        maximum_charge=arguments.maximum_charge,
        maximum_ports=arguments.maximum_ports,
        maximum_attachment_size=arguments.maximum_attachment_size,
        seed=arguments.seed,
    )
    print(json.dumps([asdict(item) for item in result], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
