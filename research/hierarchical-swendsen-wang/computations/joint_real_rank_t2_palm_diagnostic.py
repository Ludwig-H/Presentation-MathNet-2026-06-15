"""Joint real-rank geometry and ancestor-message Palm diagnostic.

The primary filter in this module is evaluated at the *real* Kruskal rank.
For every final realized LCA event, one uniformly sampled distant cross-pair
is taken from :func:`event_palm_corridor_observations`.  The event receives
Palm weight ``N_rho`` only.  On that exact pair corridor, a candidate is a
strict-arm merger with

* physical bucket size ``2 <= m <= M``;
* oriented off-spine attachment size at most ``a0``; and
* real-rank charge ``m*h_p(q_v)^2 <= J``.

The exact strict-ancestor message ``B_actual`` is then obtained from
:func:`analyze_node_polarization` on the same environment and node.  Counts
are reported before and after intersecting with bounded-message windows.

For comparison only, the module also reports the older charge obtained by
replacing ``q_v`` with ``min(q_v, q_c)``.  Every such field is explicitly
named ``criticalized_charge_proxy``.  No Blackwell domination is assumed.
The global ``port_count`` is likewise a mark-free geometric statistic only:
it is not a screening test.  The diagnostic computes neither a T2 transfer
deficit nor a weak-recovery bound.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, isfinite, sqrt
from typing import Callable, Sequence

from ancestral_polarization_palm_diagnostic import analyze_node_polarization
from corridor_t2_signature_diagnostic import (
    EventPalmCorridor,
    PalmPartitionAudit,
    ReconstructionAudit,
    event_palm_corridor_observations,
)
from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import RankedEdge, sample_ranked_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from lca_palm_corridor_diagnostic import (
    geometric_charge,
    lca_pair_partition_counts,
)


@dataclass(frozen=True)
class JointRealRankNode:
    """One oriented small-attachment node with its genuine ancestor field."""

    node: int
    bucket_size: int
    actual_rank: float
    criticalized_rank_proxy: float
    attachment_child_size: int
    port_count: int
    actual_charge: float
    criticalized_charge_proxy: float
    passes_actual_charge: bool
    passes_criticalized_charge_proxy: bool
    actual_external_message: float
    actual_identity_error: float | None
    port_count_is_geometric_proxy: bool = True
    screening_computed: bool = False


@dataclass(frozen=True)
class JointRealRankPairObservation:
    """Counts on one pair sampled inside one realized LCA event."""

    repetition: int
    palm_weight: int
    first: int
    second: int
    corridor_node_count: int
    oriented_small_attachment_count: int
    actual_charge_candidate_count: int
    actual_charge_and_bounded_message_counts: tuple[int, ...]
    criticalized_charge_proxy_candidate_count: int
    criticalized_charge_proxy_and_bounded_message_counts: tuple[int, ...]
    message_thresholds: tuple[float, ...]
    global_mean_port_count: float
    nodes: tuple[JointRealRankNode, ...]
    palm_weight_definition: str
    port_count_is_geometric_proxy: bool = True
    screening_computed: bool = False
    blackwell_domination_assumed: bool = False
    transfer_deficit_computed: bool = False
    weak_recovery_claimed: bool = False


@dataclass(frozen=True)
class JointConsistencyAudit:
    """Counter-audits linking geometry, ranks, factors, and nested filters."""

    checked_environment_count: int
    checked_pair_count: int
    checked_small_attachment_node_count: int
    node_identifier_mismatch_count: int
    bucket_size_mismatch_count: int
    actual_rank_mismatch_count: int
    criticalized_rank_mismatch_count: int
    actual_charge_mismatch_count: int
    criticalized_charge_proxy_mismatch_count: int
    nonfinite_noninfinite_message_count: int
    threshold_nesting_mismatch_count: int
    proxy_subset_mismatch_count: int
    maximum_charge_reconstruction_error: float
    maximum_finite_log_odds_identity_error: float
    passed: bool


@dataclass(frozen=True)
class JointRealRankPalmSummary:
    """Event-Palm means with delete-one-rank-environment errors."""

    side_length: int
    repetitions: int
    p: float
    critical_rank: float
    final_rank: float
    distance_fraction: float
    maximum_bucket_size: int
    maximum_attachment_size: int
    maximum_charge: float
    message_thresholds: tuple[float, ...]
    seed: int
    pair_observation_count: int
    total_palm_weight: int
    effective_environment_count: float
    weighted_mean_corridor_node_count: float
    weighted_mean_oriented_small_attachment_count: float
    jackknife_standard_error_oriented_small_attachment_count: float | None
    weighted_mean_actual_charge_candidate_count: float
    jackknife_standard_error_actual_charge_candidate_count: float | None
    weighted_mean_actual_charge_and_bounded_message_counts: tuple[float, ...]
    jackknife_standard_errors_actual_charge_and_bounded_message_counts: tuple[
        float | None, ...
    ]
    weighted_mean_criticalized_charge_proxy_candidate_count: float
    jackknife_standard_error_criticalized_charge_proxy_candidate_count: float | None
    weighted_mean_criticalized_charge_proxy_and_bounded_message_counts: tuple[
        float, ...
    ]
    jackknife_standard_errors_criticalized_charge_proxy_and_bounded_message_counts: (
        tuple[float | None, ...]
    )
    weighted_mean_global_port_count_per_corridor_node: float
    jackknife_standard_error_global_port_count_per_corridor_node: float | None
    palm_weight_convention: str
    uncertainty_note: str
    interpretation: str
    port_count_is_geometric_proxy: bool
    screening_computed: bool
    blackwell_domination_assumed: bool
    transfer_deficit_computed: bool
    weak_recovery_claimed: bool
    reconstruction_audit: ReconstructionAudit
    palm_partition_audit: PalmPartitionAudit
    joint_consistency_audit: JointConsistencyAudit


def _validated_thresholds(values: Sequence[float]) -> tuple[float, ...]:
    thresholds = tuple(float(value) for value in values)
    if not thresholds or any(
        not isfinite(value) or value < 0.0 for value in thresholds
    ):
        raise ValueError("message thresholds must be finite and nonnegative")
    if thresholds != tuple(sorted(thresholds)) or len(set(thresholds)) != len(
        thresholds
    ):
        raise ValueError("message thresholds must be strictly increasing")
    return thresholds


def analyze_joint_pair(
    environment: HierarchicalSweepEnvironment,
    pair: EventPalmCorridor,
    *,
    maximum_bucket_size: int,
    maximum_attachment_size: int,
    maximum_charge: float,
    message_thresholds: Sequence[float] = (1.0, 2.0, 4.0),
) -> JointRealRankPairObservation:
    """Evaluate real-rank filters and exact messages on one shared corridor."""

    if maximum_bucket_size < 2:
        raise ValueError("maximum_bucket_size must be at least two")
    if maximum_attachment_size <= 0:
        raise ValueError("maximum_attachment_size must be positive")
    if maximum_charge < 0.0:
        raise ValueError("maximum_charge must be nonnegative")
    thresholds = _validated_thresholds(message_thresholds)
    if pair.palm_weight != pair.n_rho:
        raise ValueError("the pair must use realized-event N_rho Palm weight")

    nodes = []
    for signature in pair.signatures:
        is_oriented_small_attachment = (
            not signature.is_lca
            and 2 <= signature.bucket_size <= maximum_bucket_size
            and signature.attachment_child_size is not None
            and signature.attachment_child_size <= maximum_attachment_size
        )
        if not is_oriented_small_attachment:
            continue

        polarization = analyze_node_polarization(
            environment,
            signature.node,
            reference_state=0,
        )
        actual_charge = geometric_charge(
            signature.bucket_size,
            environment.p,
            signature.rank,
        )
        criticalized_rank = min(signature.rank, Q_CRITICAL)
        criticalized_charge_proxy = geometric_charge(
            signature.bucket_size,
            environment.p,
            criticalized_rank,
        )
        nodes.append(
            JointRealRankNode(
                node=signature.node,
                bucket_size=signature.bucket_size,
                actual_rank=signature.rank,
                criticalized_rank_proxy=criticalized_rank,
                attachment_child_size=signature.attachment_child_size,
                port_count=signature.port_count,
                actual_charge=actual_charge,
                criticalized_charge_proxy=criticalized_charge_proxy,
                passes_actual_charge=actual_charge <= maximum_charge,
                passes_criticalized_charge_proxy=(
                    criticalized_charge_proxy <= maximum_charge
                ),
                actual_external_message=polarization.actual_external_message,
                actual_identity_error=polarization.actual_identity_error,
            )
        )

    def bounded_count(*, proxy: bool, threshold: float) -> int:
        return sum(
            (
                node.passes_criticalized_charge_proxy
                if proxy
                else node.passes_actual_charge
            )
            and abs(node.actual_external_message) <= threshold
            for node in nodes
        )

    return JointRealRankPairObservation(
        repetition=pair.repetition,
        palm_weight=pair.palm_weight,
        first=pair.first,
        second=pair.second,
        corridor_node_count=len(pair.signatures),
        oriented_small_attachment_count=len(nodes),
        actual_charge_candidate_count=sum(node.passes_actual_charge for node in nodes),
        actual_charge_and_bounded_message_counts=tuple(
            bounded_count(proxy=False, threshold=threshold) for threshold in thresholds
        ),
        criticalized_charge_proxy_candidate_count=sum(
            node.passes_criticalized_charge_proxy for node in nodes
        ),
        criticalized_charge_proxy_and_bounded_message_counts=tuple(
            bounded_count(proxy=True, threshold=threshold) for threshold in thresholds
        ),
        message_thresholds=thresholds,
        global_mean_port_count=(
            fsum(signature.port_count for signature in pair.signatures)
            / len(pair.signatures)
        ),
        nodes=tuple(nodes),
        palm_weight_definition=pair.palm_weight_definition,
    )


def _weighted_mean(
    observations: Sequence[JointRealRankPairObservation],
    statistic: Callable[[JointRealRankPairObservation], float],
) -> float:
    total_weight = sum(item.palm_weight for item in observations)
    if total_weight <= 0:
        raise ValueError("positive Palm mass is required")
    return (
        fsum(item.palm_weight * statistic(item) for item in observations) / total_weight
    )


def _cluster_jackknife_standard_error(
    observations: Sequence[JointRealRankPairObservation],
    statistic: Callable[[JointRealRankPairObservation], float],
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
    mean = fsum(leave_one_out) / len(leave_one_out)
    variance = (
        (len(leave_one_out) - 1.0)
        / len(leave_one_out)
        * fsum((value - mean) ** 2 for value in leave_one_out)
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


def _joint_consistency_audit(
    observations: Sequence[JointRealRankPairObservation],
    environments: Sequence[HierarchicalSweepEnvironment],
    *,
    maximum_charge: float,
) -> JointConsistencyAudit:
    counters = {
        "node": 0,
        "bucket": 0,
        "rank": 0,
        "critical_rank": 0,
        "actual_charge": 0,
        "proxy_charge": 0,
        "bad_message": 0,
        "threshold": 0,
        "proxy_subset": 0,
    }
    maximum_charge_error = 0.0
    maximum_identity_error = 0.0
    checked_nodes = 0
    for item in observations:
        environment = environments[item.repetition]
        thresholds = item.message_thresholds
        actual_counts = item.actual_charge_and_bounded_message_counts
        proxy_counts = item.criticalized_charge_proxy_and_bounded_message_counts
        if any(
            first > second for first, second in zip(actual_counts, actual_counts[1:])
        ):
            counters["threshold"] += 1
        if any(first > second for first, second in zip(proxy_counts, proxy_counts[1:])):
            counters["threshold"] += 1
        if any(proxy > actual for proxy, actual in zip(proxy_counts, actual_counts)):
            counters["proxy_subset"] += 1
        if (
            item.criticalized_charge_proxy_candidate_count
            > item.actual_charge_candidate_count
        ):
            counters["proxy_subset"] += 1
        if actual_counts and actual_counts[-1] > item.actual_charge_candidate_count:
            counters["threshold"] += 1
        if (
            proxy_counts
            and proxy_counts[-1] > item.criticalized_charge_proxy_candidate_count
        ):
            counters["threshold"] += 1

        for node in item.nodes:
            checked_nodes += 1
            if node.node not in environment.forest.internal_nodes:
                counters["node"] += 1
                continue
            expected_bucket = environment.forest.bucket_size[node.node]
            expected_rank = environment.forest.merge_rank[node.node]
            if node.bucket_size != expected_bucket:
                counters["bucket"] += 1
            if abs(node.actual_rank - expected_rank) > 1e-14:
                counters["rank"] += 1
            expected_critical_rank = min(expected_rank, Q_CRITICAL)
            if abs(node.criticalized_rank_proxy - expected_critical_rank) > 1e-14:
                counters["critical_rank"] += 1
            expected_actual_charge = geometric_charge(
                expected_bucket, environment.p, expected_rank
            )
            expected_proxy_charge = geometric_charge(
                expected_bucket, environment.p, expected_critical_rank
            )
            actual_error = abs(node.actual_charge - expected_actual_charge)
            proxy_error = abs(node.criticalized_charge_proxy - expected_proxy_charge)
            maximum_charge_error = max(maximum_charge_error, actual_error, proxy_error)
            if actual_error > 1e-14 or (
                node.passes_actual_charge != (expected_actual_charge <= maximum_charge)
            ):
                counters["actual_charge"] += 1
            if proxy_error > 1e-14 or (
                node.passes_criticalized_charge_proxy
                != (expected_proxy_charge <= maximum_charge)
            ):
                counters["proxy_charge"] += 1
            message = node.actual_external_message
            if not isfinite(message) and message not in (
                float("inf"),
                float("-inf"),
            ):
                counters["bad_message"] += 1
            if node.actual_identity_error is not None:
                maximum_identity_error = max(
                    maximum_identity_error, node.actual_identity_error
                )

    mismatch_total = sum(counters.values())
    return JointConsistencyAudit(
        checked_environment_count=len(environments),
        checked_pair_count=len(observations),
        checked_small_attachment_node_count=checked_nodes,
        node_identifier_mismatch_count=counters["node"],
        bucket_size_mismatch_count=counters["bucket"],
        actual_rank_mismatch_count=counters["rank"],
        criticalized_rank_mismatch_count=counters["critical_rank"],
        actual_charge_mismatch_count=counters["actual_charge"],
        criticalized_charge_proxy_mismatch_count=counters["proxy_charge"],
        nonfinite_noninfinite_message_count=counters["bad_message"],
        threshold_nesting_mismatch_count=counters["threshold"],
        proxy_subset_mismatch_count=counters["proxy_subset"],
        maximum_charge_reconstruction_error=maximum_charge_error,
        maximum_finite_log_odds_identity_error=maximum_identity_error,
        passed=mismatch_total == 0 and maximum_identity_error < 5e-12,
    )


def summarize_joint_diagnostic(
    observations: Sequence[JointRealRankPairObservation],
    environments: Sequence[HierarchicalSweepEnvironment],
    reconstruction_audits: Sequence[ReconstructionAudit],
    palm_partition_audit: PalmPartitionAudit,
    *,
    side_length: int,
    repetitions: int,
    p: float,
    distance_fraction: float,
    maximum_bucket_size: int,
    maximum_attachment_size: int,
    maximum_charge: float,
    seed: int,
) -> JointRealRankPalmSummary:
    """Aggregate one shared run without treating LCA observations as iid."""

    if not observations:
        raise ValueError("at least one pair observation is required")
    if len(environments) != repetitions:
        raise ValueError("one environment is required per repetition")
    thresholds = observations[0].message_thresholds
    if any(item.message_thresholds != thresholds for item in observations):
        raise ValueError("all observations must use identical thresholds")

    reconstruction_audit = _aggregate_reconstruction_audits(reconstruction_audits)
    consistency_audit = _joint_consistency_audit(
        observations, environments, maximum_charge=maximum_charge
    )
    if not reconstruction_audit.passed:
        raise AssertionError("Kruskal reconstruction audit failed")
    if not palm_partition_audit.passed:
        raise AssertionError("event N_rho Palm partition audit failed")
    if not consistency_audit.passed:
        raise AssertionError("joint geometry/message consistency audit failed")

    total_weight = sum(item.palm_weight for item in observations)
    environment_weights = [0] * repetitions
    for item in observations:
        environment_weights[item.repetition] += item.palm_weight
    effective_environment_count = (
        total_weight
        * total_weight
        / sum(weight * weight for weight in environment_weights if weight > 0)
    )

    def indexed_count(
        attribute: str, index: int
    ) -> Callable[[JointRealRankPairObservation], float]:
        return lambda item: float(getattr(item, attribute)[index])

    actual_attribute = "actual_charge_and_bounded_message_counts"
    proxy_attribute = "criticalized_charge_proxy_and_bounded_message_counts"
    actual_means = tuple(
        _weighted_mean(observations, indexed_count(actual_attribute, index))
        for index in range(len(thresholds))
    )
    actual_errors = tuple(
        _cluster_jackknife_standard_error(
            observations,
            indexed_count(actual_attribute, index),
            repetitions,
        )
        for index in range(len(thresholds))
    )
    proxy_means = tuple(
        _weighted_mean(observations, indexed_count(proxy_attribute, index))
        for index in range(len(thresholds))
    )
    proxy_errors = tuple(
        _cluster_jackknife_standard_error(
            observations,
            indexed_count(proxy_attribute, index),
            repetitions,
        )
        for index in range(len(thresholds))
    )

    small_count = lambda item: float(item.oriented_small_attachment_count)
    actual_count = lambda item: float(item.actual_charge_candidate_count)
    proxy_count = lambda item: float(item.criticalized_charge_proxy_candidate_count)
    global_ports = lambda item: item.global_mean_port_count
    return JointRealRankPalmSummary(
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        critical_rank=Q_CRITICAL,
        final_rank=2.0 * p - 1.0,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_attachment_size=maximum_attachment_size,
        maximum_charge=maximum_charge,
        message_thresholds=thresholds,
        seed=seed,
        pair_observation_count=len(observations),
        total_palm_weight=total_weight,
        effective_environment_count=effective_environment_count,
        weighted_mean_corridor_node_count=_weighted_mean(
            observations, lambda item: float(item.corridor_node_count)
        ),
        weighted_mean_oriented_small_attachment_count=_weighted_mean(
            observations, small_count
        ),
        jackknife_standard_error_oriented_small_attachment_count=(
            _cluster_jackknife_standard_error(observations, small_count, repetitions)
        ),
        weighted_mean_actual_charge_candidate_count=_weighted_mean(
            observations, actual_count
        ),
        jackknife_standard_error_actual_charge_candidate_count=(
            _cluster_jackknife_standard_error(observations, actual_count, repetitions)
        ),
        weighted_mean_actual_charge_and_bounded_message_counts=actual_means,
        jackknife_standard_errors_actual_charge_and_bounded_message_counts=(
            actual_errors
        ),
        weighted_mean_criticalized_charge_proxy_candidate_count=(
            _weighted_mean(observations, proxy_count)
        ),
        jackknife_standard_error_criticalized_charge_proxy_candidate_count=(
            _cluster_jackknife_standard_error(observations, proxy_count, repetitions)
        ),
        weighted_mean_criticalized_charge_proxy_and_bounded_message_counts=(
            proxy_means
        ),
        jackknife_standard_errors_criticalized_charge_proxy_and_bounded_message_counts=(
            proxy_errors
        ),
        weighted_mean_global_port_count_per_corridor_node=_weighted_mean(
            observations, global_ports
        ),
        jackknife_standard_error_global_port_count_per_corridor_node=(
            _cluster_jackknife_standard_error(observations, global_ports, repetitions)
        ),
        palm_weight_convention=(
            "realized LCA event weighted by N_rho only; one uniform distant "
            "cross-pair sampled within each event"
        ),
        uncertainty_note=(
            "standard errors delete one complete rank environment; LCA "
            "observations within an environment are not treated as iid"
        ),
        interpretation=(
            "primary counts use the real-rank charge and the actual strict-"
            "ancestor message; criticalized charges and port_count are "
            "geometric proxies only, with no domination, screening, T2 "
            "deficit, or weak-recovery conclusion"
        ),
        port_count_is_geometric_proxy=True,
        screening_computed=False,
        blackwell_domination_assumed=False,
        transfer_deficit_computed=False,
        weak_recovery_claimed=False,
        reconstruction_audit=reconstruction_audit,
        palm_partition_audit=palm_partition_audit,
        joint_consistency_audit=consistency_audit,
    )


def run_diagnostic(
    *,
    side_length: int,
    repetitions: int,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 8,
    maximum_attachment_size: int = 4,
    maximum_charge: float = 1.0,
    message_thresholds: Sequence[float] = (1.0, 2.0, 4.0),
    seed: int = 20260723,
) -> JointRealRankPalmSummary:
    """Run the joint diagnostic on shared geometry/message environments."""

    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2,1)")
    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0,1)")
    if maximum_bucket_size < 2:
        raise ValueError("maximum_bucket_size must be at least two")
    if maximum_attachment_size <= 0:
        raise ValueError("maximum_attachment_size must be positive")
    if maximum_charge < 0.0:
        raise ValueError("maximum_charge must be nonnegative")
    thresholds = _validated_thresholds(message_thresholds)
    final_rank = 2.0 * p - 1.0
    if final_rank <= Q_CRITICAL:
        raise ValueError("the final rank must strictly exceed q_c")

    rng = random.Random(seed)
    environments = []
    observations: list[JointRealRankPairObservation] = []
    reconstruction_audits = []
    event_palm_weight_total = 0
    connected_pair_total = 0
    maximum_environment_difference = 0
    for repetition in range(repetitions):
        ranked_edges: tuple[RankedEdge, ...] = sample_ranked_edges(side_length, rng)
        environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
        environments.append(environment)
        pairs, audit = event_palm_corridor_observations(
            environment.forest,
            ranked_edges,
            repetition=repetition,
            p=p,
            distance_fraction=distance_fraction,
            rng=rng,
        )
        reconstruction_audits.append(audit)
        if not audit.passed:
            raise AssertionError("local Kruskal reconstruction audit failed")

        environment_weight = sum(pair.palm_weight for pair in pairs)
        event_pairs, connected_pairs = lca_pair_partition_counts(
            environment.forest, distance_fraction
        )
        maximum_environment_difference = max(
            maximum_environment_difference,
            abs(environment_weight - event_pairs),
            abs(event_pairs - connected_pairs),
        )
        event_palm_weight_total += environment_weight
        connected_pair_total += connected_pairs
        observations.extend(
            analyze_joint_pair(
                environment,
                pair,
                maximum_bucket_size=maximum_bucket_size,
                maximum_attachment_size=maximum_attachment_size,
                maximum_charge=maximum_charge,
                message_thresholds=thresholds,
            )
            for pair in pairs
        )

    if not observations:
        raise RuntimeError("no distant realized-event Palm observation")
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
    return summarize_joint_diagnostic(
        observations,
        environments,
        reconstruction_audits,
        palm_partition_audit,
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_attachment_size=maximum_attachment_size,
        maximum_charge=maximum_charge,
        seed=seed,
    )


def run_size_series(
    side_lengths: Sequence[int],
    repetitions_by_side: Sequence[int],
    *,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 8,
    maximum_attachment_size: int = 4,
    maximum_charge: float = 1.0,
    message_thresholds: Sequence[float] = (1.0, 2.0, 4.0),
    seed: int = 20260723,
) -> tuple[JointRealRankPalmSummary, ...]:
    """Run several sizes with explicit efforts and independent streams."""

    sizes = tuple(side_lengths)
    repetitions = tuple(repetitions_by_side)
    if not sizes or len(sizes) != len(repetitions):
        raise ValueError("side lengths and repetition counts must align")
    return tuple(
        run_diagnostic(
            side_length=side_length,
            repetitions=repetition_count,
            p=p,
            distance_fraction=distance_fraction,
            maximum_bucket_size=maximum_bucket_size,
            maximum_attachment_size=maximum_attachment_size,
            maximum_charge=maximum_charge,
            message_thresholds=message_thresholds,
            seed=seed + 1_000_003 * index,
        )
        for index, (side_length, repetition_count) in enumerate(
            zip(sizes, repetitions, strict=True)
        )
    )


def _parse_ints(raw: str) -> tuple[int, ...]:
    return tuple(int(value) for value in raw.split(","))


def _parse_floats(raw: str) -> tuple[float, ...]:
    return tuple(float(value) for value in raw.split(","))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sides", default="8,12,16")
    parser.add_argument("--repetitions", default="24,10,5")
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-bucket-size", type=int, default=8)
    parser.add_argument("--maximum-attachment-size", type=int, default=4)
    parser.add_argument("--maximum-charge", type=float, default=1.0)
    parser.add_argument("--message-thresholds", default="1,2,4")
    parser.add_argument("--seed", type=int, default=20260723)
    arguments = parser.parse_args()
    result = run_size_series(
        _parse_ints(arguments.sides),
        _parse_ints(arguments.repetitions),
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        maximum_bucket_size=arguments.maximum_bucket_size,
        maximum_attachment_size=arguments.maximum_attachment_size,
        maximum_charge=arguments.maximum_charge,
        message_thresholds=_parse_floats(arguments.message_thresholds),
        seed=arguments.seed,
    )
    print(json.dumps([asdict(item) for item in result], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
