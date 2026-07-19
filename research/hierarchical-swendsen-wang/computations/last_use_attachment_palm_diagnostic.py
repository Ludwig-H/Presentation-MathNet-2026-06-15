"""Incidence upper bounds on last use for real-rank corridor attachments.

Fix a strict-arm merger ``u`` on a distant-pair corridor.  Its child which
contains the marked endpoint is the ``spine`` and its other child is the
small ``attachment``.  Use the spine orientation as a gauge.  Flipping the
attachment relative to the spine can change a future hierarchical factor
only if that factor's physical bucket contains an edge incident to the
attachment.  Consequently, if ``v_0=u,v_1,...`` is the real ancestor chain,

    last_possible_attachment_use = max{d : bucket(v_d) meets the attachment}.

After the factor at that depth has been processed, every later factor is
structurally invariant under the attachment flip, so that relative
orientation can be summed out at the factor-graph level.  The converse need
not hold: incident contributions can cancel, so this incidence depth is a
conservative upper bound on the last actual or function-level dependence,
not its exact value.  The statement is conditional on retaining the spine
and exterior state.  It does not assert that a proposed multi-update dynamics
is Markov after the projection.

The diagnostic also records ``last_union_use`` with the same upper-bound
convention.  This is the last *strict ancestor* bucket incident to either
original child and is a conservative obstruction to discarding their
common/global orientation.  Depth zero is used by convention when no strict
ancestor sees the union, because the local bucket is already invariant under
a common flip.  This possible use is generally later than the attachment
one.  Merely asking whether one future bucket has incidences from both
children is not a safe criterion: a bucket incident only to the attachment
can still depend on the relative orientation.

All buckets and ranks are those of the genuine Kruskal forest censored at
``q_1=2p-1``.  Realized LCA events carry Palm weight ``N_rho`` only, with one
uniform distant pair sampled inside each event.  Standard errors delete one
complete rank environment.  No criticalization, transfer deficit, mixing,
screening, or weak-recovery conclusion is computed.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, sqrt
from typing import Callable, Sequence

from corridor_t2_signature_diagnostic import (
    CorridorT2Signature,
    EventPalmCorridor,
    event_palm_corridor_observations,
)
from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import RankedEdge, sample_ranked_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from lca_palm_corridor_diagnostic import lca_pair_partition_counts


@dataclass(frozen=True)
class AncestorIncidence:
    """Physical incidences of the two original children in one real bucket."""

    node: int
    depth_from_update: int
    rank: float
    bucket_size: int
    attachment_incident_edge_count: int
    spine_incident_edge_count: int
    union_incident_edge_count: int
    both_children_have_incidence: bool
    attachment_orientation_sensitive: bool
    union_orientation_sensitive: bool


@dataclass(frozen=True)
class AttachmentLastUse:
    """Conservative incidence upper bound for one oriented attachment."""

    node: int
    rank: float
    bucket_size: int
    spine_child: int
    attachment_child: int
    spine_child_size: int
    attachment_child_size: int
    ancestor_chain_length: int
    pair_lca_depth: int
    last_attachment_use_depth: int
    last_attachment_use_node: int
    last_attachment_use_rank: float
    first_attachment_invariant_depth: int | None
    first_attachment_invariant_rank: float | None
    last_union_use_depth: int
    last_union_use_node: int
    last_union_use_rank: float
    first_union_invariant_depth: int | None
    first_union_invariant_rank: float | None
    last_attachment_use_no_later_than_pair_lca: bool
    last_attachment_use_is_final_root: bool
    last_union_use_is_final_root: bool
    strict_ancestor_joint_incidence_count: int
    incidences: tuple[AncestorIncidence, ...]
    incidence_upper_bound_only: bool = True
    mark_values_used: bool = False
    criticalization_used: bool = False

    def attachment_forgettable_within(self, ancestor_window: int) -> bool:
        """Whether no attachment incidence remains after depth ``k``."""

        if ancestor_window < 0:
            raise ValueError("ancestor_window must be nonnegative")
        return self.last_attachment_use_depth <= ancestor_window

    def union_forgettable_within(self, ancestor_window: int) -> bool:
        """Whether no original-union incidence remains after depth ``k``."""

        if ancestor_window < 0:
            raise ValueError("ancestor_window must be nonnegative")
        return self.last_union_use_depth <= ancestor_window


@dataclass(frozen=True)
class PairLastUseObservation:
    """Incidence last-use bounds on one event-Palm pair corridor."""

    repetition: int
    palm_weight: int
    palm_weight_definition: str
    first: int
    second: int
    lca_node: int
    corridor_node_count: int
    candidates: tuple[AttachmentLastUse, ...]


@dataclass(frozen=True)
class LastUseReconstructionAudit:
    """Independent bucket/LCA identities behind the structural criterion."""

    checked_candidate_count: int
    bucket_size_mismatch_count: int
    child_size_mismatch_count: int
    ancestor_chain_mismatch_count: int
    nonincreasing_rank_chain_count: int
    attachment_boundary_partition_mismatch_count: int
    attachment_edge_lca_mismatch_count: int
    last_use_mismatch_count: int
    future_attachment_incidence_mismatch_count: int
    passed: bool


@dataclass(frozen=True)
class PalmPartitionAudit:
    """Finite-volume identity for realized-event ``N_rho`` pair mass."""

    checked_environment_count: int
    event_palm_weight_total: int
    connected_distant_ordered_pair_total: int
    maximum_environment_difference: int
    passed: bool


@dataclass(frozen=True)
class LastUsePalmSummary:
    """Candidate-Palm incidence bounds with environment jackknife errors."""

    side_length: int
    repetitions: int
    p: float
    final_rank: float
    distance_fraction: float
    maximum_bucket_size: int
    maximum_attachment_size: int
    minimum_rank: float
    ancestor_windows: tuple[int, ...]
    seed: int
    pair_observation_count: int
    total_palm_weight: int
    effective_environment_count: float
    weighted_mean_candidate_count_per_pair: float
    jackknife_standard_error_candidate_count_per_pair: float | None
    weighted_candidate_mass: int
    weighted_mean_last_attachment_use_depth: float
    jackknife_standard_error_last_attachment_use_depth: float | None
    weighted_median_last_attachment_use_depth: int
    weighted_ninetieth_percentile_last_attachment_use_depth: int
    weighted_attachment_forgettable_fractions: tuple[float, ...]
    jackknife_standard_errors_attachment_forgettable_fractions: tuple[float | None, ...]
    weighted_union_forgettable_fractions: tuple[float, ...]
    jackknife_standard_errors_union_forgettable_fractions: tuple[float | None, ...]
    weighted_fraction_attachment_used_by_strict_ancestor: float
    jackknife_standard_error_attachment_used_by_strict_ancestor: float | None
    weighted_fraction_attachment_last_use_no_later_than_pair_lca: float
    jackknife_standard_error_attachment_last_use_no_later_than_pair_lca: float | None
    weighted_fraction_attachment_last_use_is_final_root: float
    weighted_fraction_union_last_use_is_final_root: float
    weighted_fraction_with_strict_ancestor_joint_incidence: float
    weighted_mean_last_attachment_use_rank: float
    weighted_mean_last_union_use_depth: float
    palm_weight_convention: str
    exact_criterion: str
    interpretation: str
    uncertainty_note: str
    transfer_deficit_computed: bool
    weak_recovery_claimed: bool
    reconstruction_audit: LastUseReconstructionAudit
    palm_partition_audit: PalmPartitionAudit


@dataclass
class _MutableAudit:
    checked_candidate_count: int = 0
    bucket_size_mismatch_count: int = 0
    child_size_mismatch_count: int = 0
    ancestor_chain_mismatch_count: int = 0
    nonincreasing_rank_chain_count: int = 0
    attachment_boundary_partition_mismatch_count: int = 0
    attachment_edge_lca_mismatch_count: int = 0
    last_use_mismatch_count: int = 0
    future_attachment_incidence_mismatch_count: int = 0

    def freeze(self) -> LastUseReconstructionAudit:
        values = (
            self.bucket_size_mismatch_count,
            self.child_size_mismatch_count,
            self.ancestor_chain_mismatch_count,
            self.nonincreasing_rank_chain_count,
            self.attachment_boundary_partition_mismatch_count,
            self.attachment_edge_lca_mismatch_count,
            self.last_use_mismatch_count,
            self.future_attachment_incidence_mismatch_count,
        )
        return LastUseReconstructionAudit(
            checked_candidate_count=self.checked_candidate_count,
            bucket_size_mismatch_count=self.bucket_size_mismatch_count,
            child_size_mismatch_count=self.child_size_mismatch_count,
            ancestor_chain_mismatch_count=self.ancestor_chain_mismatch_count,
            nonincreasing_rank_chain_count=self.nonincreasing_rank_chain_count,
            attachment_boundary_partition_mismatch_count=(
                self.attachment_boundary_partition_mismatch_count
            ),
            attachment_edge_lca_mismatch_count=(
                self.attachment_edge_lca_mismatch_count
            ),
            last_use_mismatch_count=self.last_use_mismatch_count,
            future_attachment_incidence_mismatch_count=(
                self.future_attachment_incidence_mismatch_count
            ),
            passed=not any(values),
        )


def _mask_size(mask: int) -> int:
    return mask.bit_count()


def _edge_meets_mask(edge: tuple[int, int], mask: int) -> bool:
    return bool((mask >> edge[0]) & 1) != bool((mask >> edge[1]) & 1)


def _oriented_children(
    environment: HierarchicalSweepEnvironment,
    pair: EventPalmCorridor,
    signature: CorridorT2Signature,
) -> tuple[int, int]:
    node = signature.node
    if signature.is_lca or node == pair.lca_node:
        raise ValueError("last-use candidates must be strict-arm nodes")
    if signature not in pair.signatures:
        raise ValueError("signature does not belong to the pair corridor")
    left = environment.forest.left_child[node]
    right = environment.forest.right_child[node]
    union_mask = environment.cluster_mask[node]
    membership = (
        bool(union_mask & (1 << pair.first)),
        bool(union_mask & (1 << pair.second)),
    )
    if sum(membership) != 1:
        raise AssertionError("a strict corridor arm contains one marked endpoint")
    endpoint = pair.first if membership[0] else pair.second
    if environment.cluster_mask[left] & (1 << endpoint):
        return left, right
    if environment.cluster_mask[right] & (1 << endpoint):
        return right, left
    raise AssertionError("the arm endpoint is absent from both children")


def analyze_attachment_last_use(
    environment: HierarchicalSweepEnvironment,
    pair: EventPalmCorridor,
    signature: CorridorT2Signature,
) -> tuple[AttachmentLastUse, LastUseReconstructionAudit]:
    """Return a safe incidence upper bound on the last factor use.

    The structural test uses edge endpoints only.  Edge marks and the
    planted/current spin state cannot turn a zero incidence into a hidden
    dependence, so no Nishimori mark resampling is involved.  They can make
    incident contributions cancel, however; the returned depth is therefore
    not claimed to be the minimal realized or function-level last use.
    """

    audit = _MutableAudit(checked_candidate_count=1)
    spine_child, attachment_child = _oriented_children(environment, pair, signature)
    forest = environment.forest
    node = signature.node
    spine_mask = environment.cluster_mask[spine_child]
    attachment_mask = environment.cluster_mask[attachment_child]
    union_mask = spine_mask | attachment_mask
    if spine_mask & attachment_mask or union_mask != environment.cluster_mask[node]:
        audit.child_size_mismatch_count += 1
    if signature.spine_child_size != _mask_size(
        spine_mask
    ) or signature.attachment_child_size != _mask_size(attachment_mask):
        audit.child_size_mismatch_count += 1

    chain = environment.factor_chains[node]
    if not chain or chain[0] != node or pair.lca_node not in chain:
        audit.ancestor_chain_mismatch_count += 1
        raise AssertionError("the factor chain does not contain the marked LCA")
    ranks = tuple(forest.merge_rank[factor_node] for factor_node in chain)
    if any(first >= second for first, second in zip(ranks, ranks[1:])):
        audit.nonincreasing_rank_chain_count += 1

    incidences = []
    scanned_attachment_edges: set[tuple[int, int]] = set()
    lca_depths: dict[int, int] = {
        factor_node: depth for depth, factor_node in enumerate(chain)
    }
    for depth, factor_node in enumerate(chain):
        bucket = environment.bucket_edges[factor_node]
        if len(bucket) != forest.bucket_size[factor_node]:
            audit.bucket_size_mismatch_count += 1
        attachment_edges = tuple(
            edge for edge in bucket if _edge_meets_mask(edge, attachment_mask)
        )
        spine_edges = tuple(
            edge for edge in bucket if _edge_meets_mask(edge, spine_mask)
        )
        union_edges = tuple(
            edge for edge in bucket if _edge_meets_mask(edge, union_mask)
        )
        scanned_attachment_edges.update(attachment_edges)
        incidences.append(
            AncestorIncidence(
                node=factor_node,
                depth_from_update=depth,
                rank=forest.merge_rank[factor_node],
                bucket_size=len(bucket),
                attachment_incident_edge_count=len(attachment_edges),
                spine_incident_edge_count=len(spine_edges),
                union_incident_edge_count=len(union_edges),
                both_children_have_incidence=bool(attachment_edges and spine_edges),
                attachment_orientation_sensitive=bool(attachment_edges),
                union_orientation_sensitive=bool(union_edges),
            )
        )

    factor_boundary_edges = {
        edge
        for edge in forest.edges
        if forest.connected(*edge) and _edge_meets_mask(edge, attachment_mask)
    }
    if scanned_attachment_edges != factor_boundary_edges:
        audit.attachment_boundary_partition_mismatch_count += 1

    edge_lca_depths = []
    for edge in sorted(factor_boundary_edges):
        edge_lca = forest.tree_lca(*edge)
        depth = lca_depths.get(edge_lca)
        if depth is None or edge not in environment.bucket_edges[edge_lca]:
            audit.attachment_edge_lca_mismatch_count += 1
            continue
        edge_lca_depths.append(depth)

    attachment_use_depths = tuple(
        item.depth_from_update
        for item in incidences
        if item.attachment_orientation_sensitive
    )
    strict_union_use_depths = tuple(
        item.depth_from_update
        for item in incidences
        if item.union_orientation_sensitive
    )
    if not attachment_use_depths or attachment_use_depths[0] != 0:
        audit.last_use_mismatch_count += 1
        raise AssertionError("the local merger bucket must meet the attachment")
    last_attachment_depth = max(attachment_use_depths)
    # A common flip of the two children leaves the local bucket invariant.
    # Depth zero is the natural elimination point if no strict ancestor uses
    # the union orientation.
    last_union_depth = max((0, *strict_union_use_depths))
    if edge_lca_depths and last_attachment_depth != max(edge_lca_depths):
        audit.last_use_mismatch_count += 1
    if any(
        item.attachment_orientation_sensitive
        for item in incidences[last_attachment_depth + 1 :]
    ):
        audit.future_attachment_incidence_mismatch_count += 1

    last_chain_depth = len(chain) - 1
    pair_lca_depth = chain.index(pair.lca_node)
    first_attachment_invariant_depth = (
        last_attachment_depth + 1 if last_attachment_depth < last_chain_depth else None
    )
    first_union_invariant_depth = (
        last_union_depth + 1 if last_union_depth < last_chain_depth else None
    )
    result = AttachmentLastUse(
        node=node,
        rank=forest.merge_rank[node],
        bucket_size=forest.bucket_size[node],
        spine_child=spine_child,
        attachment_child=attachment_child,
        spine_child_size=_mask_size(spine_mask),
        attachment_child_size=_mask_size(attachment_mask),
        ancestor_chain_length=len(chain),
        pair_lca_depth=pair_lca_depth,
        last_attachment_use_depth=last_attachment_depth,
        last_attachment_use_node=chain[last_attachment_depth],
        last_attachment_use_rank=ranks[last_attachment_depth],
        first_attachment_invariant_depth=first_attachment_invariant_depth,
        first_attachment_invariant_rank=(
            None
            if first_attachment_invariant_depth is None
            else ranks[first_attachment_invariant_depth]
        ),
        last_union_use_depth=last_union_depth,
        last_union_use_node=chain[last_union_depth],
        last_union_use_rank=ranks[last_union_depth],
        first_union_invariant_depth=first_union_invariant_depth,
        first_union_invariant_rank=(
            None
            if first_union_invariant_depth is None
            else ranks[first_union_invariant_depth]
        ),
        last_attachment_use_no_later_than_pair_lca=(
            last_attachment_depth <= pair_lca_depth
        ),
        last_attachment_use_is_final_root=(last_attachment_depth == last_chain_depth),
        last_union_use_is_final_root=last_union_depth == last_chain_depth,
        strict_ancestor_joint_incidence_count=sum(
            item.depth_from_update > 0 and item.both_children_have_incidence
            for item in incidences
        ),
        incidences=tuple(incidences),
    )
    return result, audit.freeze()


def analyze_pair_last_use(
    environment: HierarchicalSweepEnvironment,
    pair: EventPalmCorridor,
    *,
    maximum_bucket_size: int,
    maximum_attachment_size: int,
    minimum_rank: float,
) -> tuple[PairLastUseObservation, LastUseReconstructionAudit]:
    """Analyze all selected strict-arm attachments on one Palm corridor."""

    if maximum_bucket_size < 2 or maximum_attachment_size <= 0:
        raise ValueError("candidate cutoffs must be positive")
    if not 0.0 <= minimum_rank <= 2.0 * environment.p - 1.0:
        raise ValueError("minimum_rank lies outside the realized rank window")
    candidates = []
    audits = []
    for signature in pair.signatures:
        if (
            not signature.is_lca
            and 2 <= signature.bucket_size <= maximum_bucket_size
            and signature.attachment_child_size is not None
            and signature.attachment_child_size <= maximum_attachment_size
            and signature.rank >= minimum_rank
        ):
            candidate, audit = analyze_attachment_last_use(environment, pair, signature)
            candidates.append(candidate)
            audits.append(audit)
    return (
        PairLastUseObservation(
            repetition=pair.repetition,
            palm_weight=pair.palm_weight,
            palm_weight_definition=pair.palm_weight_definition,
            first=pair.first,
            second=pair.second,
            lca_node=pair.lca_node,
            corridor_node_count=len(pair.signatures),
            candidates=tuple(candidates),
        ),
        _aggregate_reconstruction_audits(audits),
    )


def _aggregate_reconstruction_audits(
    audits: Sequence[LastUseReconstructionAudit],
) -> LastUseReconstructionAudit:
    fields = (
        "checked_candidate_count",
        "bucket_size_mismatch_count",
        "child_size_mismatch_count",
        "ancestor_chain_mismatch_count",
        "nonincreasing_rank_chain_count",
        "attachment_boundary_partition_mismatch_count",
        "attachment_edge_lca_mismatch_count",
        "last_use_mismatch_count",
        "future_attachment_incidence_mismatch_count",
    )
    totals = {name: sum(getattr(audit, name) for audit in audits) for name in fields}
    mismatch_total = sum(
        value for name, value in totals.items() if name != "checked_candidate_count"
    )
    return LastUseReconstructionAudit(**totals, passed=mismatch_total == 0)


def _ratio(
    observations: Sequence[PairLastUseObservation],
    numerator: Callable[[PairLastUseObservation], float],
    denominator: Callable[[PairLastUseObservation], float],
) -> float:
    numerator_total = fsum(item.palm_weight * numerator(item) for item in observations)
    denominator_total = fsum(
        item.palm_weight * denominator(item) for item in observations
    )
    if denominator_total <= 0.0:
        raise ValueError("a positive weighted denominator is required")
    return numerator_total / denominator_total


def _cluster_jackknife_ratio(
    observations: Sequence[PairLastUseObservation],
    numerator: Callable[[PairLastUseObservation], float],
    denominator: Callable[[PairLastUseObservation], float],
    repetition_count: int,
) -> float | None:
    if repetition_count < 2:
        return None
    numerators = [0.0] * repetition_count
    denominators = [0.0] * repetition_count
    for item in observations:
        numerators[item.repetition] += item.palm_weight * numerator(item)
        denominators[item.repetition] += item.palm_weight * denominator(item)
    total_numerator = fsum(numerators)
    total_denominator = fsum(denominators)
    leave_one_out = []
    for numerator_value, denominator_value in zip(
        numerators, denominators, strict=True
    ):
        remaining_denominator = total_denominator - denominator_value
        if remaining_denominator > 0.0:
            leave_one_out.append(
                (total_numerator - numerator_value) / remaining_denominator
            )
    if len(leave_one_out) < 2:
        return None
    mean = fsum(leave_one_out) / len(leave_one_out)
    variance = (
        (len(leave_one_out) - 1.0)
        / len(leave_one_out)
        * fsum((value - mean) ** 2 for value in leave_one_out)
    )
    return sqrt(max(0.0, variance))


def _weighted_depth_quantile(
    observations: Sequence[PairLastUseObservation], quantile: float
) -> int:
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must belong to [0,1]")
    masses: dict[int, int] = {}
    for item in observations:
        for candidate in item.candidates:
            depth = candidate.last_attachment_use_depth
            masses[depth] = masses.get(depth, 0) + item.palm_weight
    total = sum(masses.values())
    if total <= 0:
        raise ValueError("positive candidate mass is required")
    target = quantile * total
    cumulative = 0
    for depth, mass in sorted(masses.items()):
        cumulative += mass
        if cumulative >= target:
            return depth
    return max(masses)


def run_diagnostic(
    *,
    side_length: int,
    repetitions: int,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 8,
    maximum_attachment_size: int = 4,
    minimum_rank: float = Q_CRITICAL,
    ancestor_windows: Sequence[int] = (0, 1, 2, 4, 8),
    seed: int = 20260725,
) -> LastUsePalmSummary:
    """Run the genuine event-Palm last-use diagnostic."""

    windows = tuple(int(value) for value in ancestor_windows)
    if side_length < 4 or repetitions <= 0:
        raise ValueError("side_length and repetitions must be admissible")
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2,1)")
    final_rank = 2.0 * p - 1.0
    if not 0.0 <= minimum_rank <= final_rank:
        raise ValueError("minimum_rank lies outside the realized rank window")
    if maximum_bucket_size < 2 or maximum_attachment_size <= 0:
        raise ValueError("candidate cutoffs must be positive")
    if (
        not windows
        or any(value < 0 for value in windows)
        or tuple(sorted(set(windows))) != windows
    ):
        raise ValueError("ancestor_windows must be strictly increasing and nonnegative")

    rng = random.Random(seed)
    observations: list[PairLastUseObservation] = []
    reconstruction_audits = []
    environment_weights = [0] * repetitions
    event_weight_total = 0
    connected_pair_total = 0
    maximum_environment_difference = 0
    for repetition in range(repetitions):
        ranked_edges: tuple[RankedEdge, ...] = sample_ranked_edges(side_length, rng)
        environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
        pairs, geometry_audit = event_palm_corridor_observations(
            environment.forest,
            ranked_edges,
            repetition=repetition,
            p=p,
            distance_fraction=distance_fraction,
            rng=rng,
        )
        if not geometry_audit.passed:
            raise AssertionError("the existing real-corridor reconstruction failed")
        environment_weight = sum(pair.palm_weight for pair in pairs)
        event_pairs, connected_pairs = lca_pair_partition_counts(
            environment.forest, distance_fraction
        )
        difference = max(
            abs(environment_weight - event_pairs),
            abs(event_pairs - connected_pairs),
        )
        maximum_environment_difference = max(maximum_environment_difference, difference)
        environment_weights[repetition] = environment_weight
        event_weight_total += environment_weight
        connected_pair_total += connected_pairs
        for pair in pairs:
            observation, audit = analyze_pair_last_use(
                environment,
                pair,
                maximum_bucket_size=maximum_bucket_size,
                maximum_attachment_size=maximum_attachment_size,
                minimum_rank=minimum_rank,
            )
            observations.append(observation)
            reconstruction_audits.append(audit)

    reconstruction_audit = _aggregate_reconstruction_audits(reconstruction_audits)
    if not reconstruction_audit.passed:
        raise AssertionError("last-use reconstruction audit failed")
    palm_partition_audit = PalmPartitionAudit(
        checked_environment_count=repetitions,
        event_palm_weight_total=event_weight_total,
        connected_distant_ordered_pair_total=connected_pair_total,
        maximum_environment_difference=maximum_environment_difference,
        passed=(
            maximum_environment_difference == 0
            and event_weight_total == connected_pair_total
        ),
    )
    if not palm_partition_audit.passed:
        raise AssertionError("event N_rho weights do not partition pair mass")
    if not observations:
        raise RuntimeError("no distant realized-event Palm observation")
    candidate_mass = sum(
        item.palm_weight * len(item.candidates) for item in observations
    )
    if candidate_mass <= 0:
        raise RuntimeError("no selected small attachment in the Palm sample")

    pair_unit = lambda item: 1.0
    candidate_count = lambda item: float(len(item.candidates))
    depth_sum = lambda item: float(
        sum(candidate.last_attachment_use_depth for candidate in item.candidates)
    )
    union_depth_sum = lambda item: float(
        sum(candidate.last_union_use_depth for candidate in item.candidates)
    )
    rank_sum = lambda item: fsum(
        candidate.last_attachment_use_rank for candidate in item.candidates
    )

    attachment_fractions = []
    attachment_errors = []
    union_fractions = []
    union_errors = []
    for window in windows:
        attachment_count = lambda item, selected=window: float(
            sum(
                candidate.attachment_forgettable_within(selected)
                for candidate in item.candidates
            )
        )
        union_count = lambda item, selected=window: float(
            sum(
                candidate.union_forgettable_within(selected)
                for candidate in item.candidates
            )
        )
        attachment_fractions.append(
            _ratio(observations, attachment_count, candidate_count)
        )
        attachment_errors.append(
            _cluster_jackknife_ratio(
                observations,
                attachment_count,
                candidate_count,
                repetitions,
            )
        )
        union_fractions.append(_ratio(observations, union_count, candidate_count))
        union_errors.append(
            _cluster_jackknife_ratio(
                observations, union_count, candidate_count, repetitions
            )
        )

    strict_use = lambda item: float(
        sum(candidate.last_attachment_use_depth > 0 for candidate in item.candidates)
    )
    before_lca = lambda item: float(
        sum(
            candidate.last_attachment_use_no_later_than_pair_lca
            for candidate in item.candidates
        )
    )
    attachment_at_root = lambda item: float(
        sum(
            candidate.last_attachment_use_is_final_root for candidate in item.candidates
        )
    )
    union_at_root = lambda item: float(
        sum(candidate.last_union_use_is_final_root for candidate in item.candidates)
    )
    joint_incidence = lambda item: float(
        sum(
            candidate.strict_ancestor_joint_incidence_count > 0
            for candidate in item.candidates
        )
    )

    positive_environment_weights = [value for value in environment_weights if value]
    effective_environment_count = (
        event_weight_total
        * event_weight_total
        / sum(value * value for value in positive_environment_weights)
    )
    return LastUsePalmSummary(
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        final_rank=final_rank,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_attachment_size=maximum_attachment_size,
        minimum_rank=minimum_rank,
        ancestor_windows=windows,
        seed=seed,
        pair_observation_count=len(observations),
        total_palm_weight=event_weight_total,
        effective_environment_count=effective_environment_count,
        weighted_mean_candidate_count_per_pair=_ratio(
            observations, candidate_count, pair_unit
        ),
        jackknife_standard_error_candidate_count_per_pair=(
            _cluster_jackknife_ratio(
                observations, candidate_count, pair_unit, repetitions
            )
        ),
        weighted_candidate_mass=candidate_mass,
        weighted_mean_last_attachment_use_depth=_ratio(
            observations, depth_sum, candidate_count
        ),
        jackknife_standard_error_last_attachment_use_depth=(
            _cluster_jackknife_ratio(
                observations, depth_sum, candidate_count, repetitions
            )
        ),
        weighted_median_last_attachment_use_depth=_weighted_depth_quantile(
            observations, 0.5
        ),
        weighted_ninetieth_percentile_last_attachment_use_depth=(
            _weighted_depth_quantile(observations, 0.9)
        ),
        weighted_attachment_forgettable_fractions=tuple(attachment_fractions),
        jackknife_standard_errors_attachment_forgettable_fractions=tuple(
            attachment_errors
        ),
        weighted_union_forgettable_fractions=tuple(union_fractions),
        jackknife_standard_errors_union_forgettable_fractions=tuple(union_errors),
        weighted_fraction_attachment_used_by_strict_ancestor=_ratio(
            observations, strict_use, candidate_count
        ),
        jackknife_standard_error_attachment_used_by_strict_ancestor=(
            _cluster_jackknife_ratio(
                observations, strict_use, candidate_count, repetitions
            )
        ),
        weighted_fraction_attachment_last_use_no_later_than_pair_lca=_ratio(
            observations, before_lca, candidate_count
        ),
        jackknife_standard_error_attachment_last_use_no_later_than_pair_lca=(
            _cluster_jackknife_ratio(
                observations, before_lca, candidate_count, repetitions
            )
        ),
        weighted_fraction_attachment_last_use_is_final_root=_ratio(
            observations, attachment_at_root, candidate_count
        ),
        weighted_fraction_union_last_use_is_final_root=_ratio(
            observations, union_at_root, candidate_count
        ),
        weighted_fraction_with_strict_ancestor_joint_incidence=_ratio(
            observations, joint_incidence, candidate_count
        ),
        weighted_mean_last_attachment_use_rank=_ratio(
            observations, rank_sum, candidate_count
        ),
        weighted_mean_last_union_use_depth=_ratio(
            observations, union_depth_sum, candidate_count
        ),
        palm_weight_convention=(
            "realized LCA event weighted by N_rho only; one uniform distant "
            "pair sampled within the event"
        ),
        exact_criterion=(
            "zero later incidence exactly certifies invariance in the spine "
            "gauge; the last incident bucket is only a conservative upper "
            "bound because incident contributions may cancel"
        ),
        interpretation=(
            "structural incidence upper bound only; a uniformly shallow bound "
            "may motivate a multi-update block, but no minimal functional last "
            "use, no projected Markov closure, and no T2 deficit, screening, "
            "mixing, or weak-recovery bound is certified"
        ),
        uncertainty_note=(
            "ratios use candidate Palm mass and standard errors are "
            "delete-one-rank-environment jackknives"
        ),
        transfer_deficit_computed=False,
        weak_recovery_claimed=False,
        reconstruction_audit=reconstruction_audit,
        palm_partition_audit=palm_partition_audit,
    )


def run_size_series(
    side_lengths: Sequence[int],
    *,
    repetitions: int,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 8,
    maximum_attachment_size: int = 4,
    minimum_rank: float = Q_CRITICAL,
    ancestor_windows: Sequence[int] = (0, 1, 2, 4, 8),
    seed: int = 20260725,
) -> tuple[LastUsePalmSummary, ...]:
    """Run modest sizes with independent deterministic random streams."""

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
            maximum_attachment_size=maximum_attachment_size,
            minimum_rank=minimum_rank,
            ancestor_windows=ancestor_windows,
            seed=seed + index * 1_000_003,
        )
        for index, side_length in enumerate(sizes)
    )


def _parse_integer_tuple(raw: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in raw.split(","))
    if not values:
        raise ValueError("at least one integer is required")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sides", default="8,12,16")
    parser.add_argument("--repetitions", type=int, default=8)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-bucket-size", type=int, default=8)
    parser.add_argument("--maximum-attachment-size", type=int, default=4)
    parser.add_argument("--minimum-rank", type=float, default=Q_CRITICAL)
    parser.add_argument("--ancestor-windows", default="0,1,2,4,8")
    parser.add_argument("--seed", type=int, default=20260725)
    arguments = parser.parse_args()
    result = run_size_series(
        _parse_integer_tuple(arguments.sides),
        repetitions=arguments.repetitions,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        maximum_bucket_size=arguments.maximum_bucket_size,
        maximum_attachment_size=arguments.maximum_attachment_size,
        minimum_rank=arguments.minimum_rank,
        ancestor_windows=_parse_integer_tuple(arguments.ancestor_windows),
        seed=arguments.seed,
    )
    print(json.dumps([asdict(item) for item in result], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
