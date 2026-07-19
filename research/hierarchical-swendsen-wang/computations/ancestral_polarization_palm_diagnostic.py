"""Palm diagnostic for ancestral polarization along the real corridor.

The geometric proxy ``m*h_p(q)^2`` does not reveal whether a small merger
bucket is already pinned by its strict ancestors.  This module measures that
missing quantity in the genuine finite-volume joint Nishimori experiment.

For a realized Kruskal node ``u`` and the four flips ``(a,b)`` of its two
children, it evaluates separately

    local_u(a,b) = log F_u(Lambda_u^{ab}),
    external_u(a,b) = sum_{v strict ancestor of u} log F_v(Lambda_v^{ab}),

where ``F_v(x)=x*exp((1-beta_v)*x)``.  The external parity log-odds is

    B_u = LSE(external_00, external_11)
          - LSE(external_10, external_01).

The planted reference state is a sample from the joint posterior/dendrogram
law before gauging.  Averaging these messages over ranked-edge environments
therefore gives a Nishimori diagnostic, not an adversarial boundary scan.

Pairs are represented through already-realized LCA events.  Such events are
weighted by ``N_rho`` only: the Kruskal race has already supplied the cut-size
bias.  Late ranks can also be replaced by ``min(q_v,q_c)`` on the unchanged
skeleton.  The resulting ``favourable_message`` is an evaluation of the
criticalized factors on the original joint sample; it is not a claim that
the sample itself has been redrawn from the criticalized conditional law.

This module measures ancestor fields only.  A bounded field is not a proof of
lateral screening, transfer contraction, or a weak-recovery threshold.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import exp, fsum, isfinite, log, sqrt, tanh
from typing import Callable, Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    RankedEdge,
    clock_time_from_rank,
    sample_ranked_edges,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from lca_palm_corridor_diagnostic import (
    PalmCorridorObservation,
    final_realized_event_observations,
)


@dataclass(frozen=True)
class NodePolarization:
    """Exact local and ancestral log-odds at one realized corridor node."""

    node: int
    bucket_size: int
    actual_rank: float
    favourable_rank: float
    strict_ancestor_count: int
    actual_external_message: float
    favourable_external_message: float
    actual_local_log_odds: float
    favourable_local_log_odds: float
    actual_total_log_odds: float
    favourable_total_log_odds: float
    actual_identity_error: float | None
    favourable_identity_error: float | None

    @property
    def actual_local_reliability(self) -> float:
        return _log_odds_reliability(self.actual_local_log_odds)

    @property
    def favourable_local_reliability(self) -> float:
        return _log_odds_reliability(self.favourable_local_log_odds)

    @property
    def actual_total_reliability(self) -> float:
        return _log_odds_reliability(self.actual_total_log_odds)

    @property
    def favourable_total_reliability(self) -> float:
        return _log_odds_reliability(self.favourable_total_log_odds)


@dataclass(frozen=True)
class PairPolarizationObservation:
    """Candidate counts and messages for one LCA-Palm endpoint pair."""

    repetition: int
    palm_weight: int
    first: int
    second: int
    corridor_node_count: int
    candidate_count: int
    message_thresholds: tuple[float, ...]
    actual_bounded_message_counts: tuple[int, ...]
    favourable_bounded_message_counts: tuple[int, ...]
    actual_infinite_message_count: int
    favourable_infinite_message_count: int
    actual_local_reliability_sum: float
    favourable_local_reliability_sum: float
    actual_total_reliability_sum: float
    favourable_total_reliability_sum: float
    candidates: tuple[NodePolarization, ...]


@dataclass(frozen=True)
class PolarizationPalmSummary:
    """Weighted finite-volume summary with environment-cluster errors."""

    side_length: int
    repetitions: int
    p: float
    critical_rank: float
    final_rank: float
    distance_fraction: float
    maximum_bucket_size: int
    maximum_charge: float
    seed: int
    pair_observation_count: int
    total_palm_weight: int
    effective_repetition_count: float
    weighted_mean_corridor_node_count: float
    weighted_mean_candidate_count: float
    jackknife_standard_error_candidate_count: float | None
    message_thresholds: tuple[float, ...]
    weighted_mean_actual_bounded_message_counts: tuple[float, ...]
    jackknife_standard_errors_actual_bounded_message_counts: tuple[float | None, ...]
    weighted_mean_favourable_bounded_message_counts: tuple[float, ...]
    jackknife_standard_errors_favourable_bounded_message_counts: tuple[
        float | None, ...
    ]
    actual_infinite_message_fraction_among_candidates: float
    favourable_infinite_message_fraction_among_candidates: float
    actual_finite_message_median: float | None
    favourable_finite_message_median: float | None
    actual_finite_message_ninetieth_percentile: float | None
    favourable_finite_message_ninetieth_percentile: float | None
    actual_mean_local_reliability_among_candidates: float
    favourable_mean_local_reliability_among_candidates: float
    actual_mean_total_reliability_among_candidates: float
    favourable_mean_total_reliability_among_candidates: float
    maximum_finite_log_odds_identity_error: float
    interpretation: str


def _logsumexp_pair(first: float, second: float) -> float:
    if first == float("-inf"):
        return second
    if second == float("-inf"):
        return first
    maximum = max(first, second)
    return maximum + log(exp(first - maximum) + exp(second - maximum))


def parity_log_odds(log_weights: Sequence[float]) -> float:
    """Return log((w_00+w_11)/(w_10+w_01)) robustly."""

    if len(log_weights) != 4:
        raise ValueError("four orbit weights are required")
    same = _logsumexp_pair(log_weights[0], log_weights[3])
    opposite = _logsumexp_pair(log_weights[1], log_weights[2])
    if same == float("-inf") and opposite == float("-inf"):
        raise ValueError("the complete four-state orbit has zero mass")
    if opposite == float("-inf"):
        return float("inf")
    if same == float("-inf"):
        return float("-inf")
    return same - opposite


def _log_odds_reliability(value: float) -> float:
    if value in (float("inf"), float("-inf")):
        return 1.0
    if not isfinite(value):
        raise ValueError("a log-odds value must be finite or signed infinity")
    return tanh(value / 2.0) ** 2


def _factor_log_weight(
    environment: HierarchicalSweepEnvironment,
    factor_node: int,
    state: int,
    *,
    favourable: bool,
) -> float:
    satisfied_count = environment._satisfied_count(factor_node, state)
    if satisfied_count == 0:
        return float("-inf")
    rank = environment.forest.merge_rank[factor_node]
    if favourable:
        rank = min(rank, Q_CRITICAL)
    beta = clock_time_from_rank(rank, environment.p)
    return log(satisfied_count) + (
        (1.0 - beta) * environment.coupling * satisfied_count
    )


def _orbit_log_weights(
    environment: HierarchicalSweepEnvironment,
    update_node: int,
    factor_nodes: Sequence[int],
    *,
    favourable: bool,
    reference_state: int,
) -> tuple[float, ...]:
    values = []
    for mask in environment.proposal_masks(update_node):
        state = reference_state ^ mask
        value = 0.0
        for factor_node in factor_nodes:
            factor = _factor_log_weight(
                environment,
                factor_node,
                state,
                favourable=favourable,
            )
            if factor == float("-inf"):
                value = factor
                break
            value += factor
        values.append(value)
    return tuple(values)


def _identity_error(total: float, external: float, local: float) -> float | None:
    if all(isfinite(value) for value in (total, external, local)):
        return abs(total - external - local)
    return None


def analyze_node_polarization(
    environment: HierarchicalSweepEnvironment,
    node: int,
    *,
    reference_state: int = 0,
) -> NodePolarization:
    """Evaluate the exact node/ancestor decomposition in both time schemes."""

    if node not in environment.forest.internal_nodes:
        raise ValueError("node must be a realized internal Kruskal node")
    chain = environment.factor_chains[node]
    if not chain or chain[0] != node:
        raise AssertionError("the factor chain must start at the update node")

    results: dict[str, tuple[float, float, float]] = {}
    for name, favourable in (("actual", False), ("favourable", True)):
        local_weights = _orbit_log_weights(
            environment,
            node,
            (node,),
            favourable=favourable,
            reference_state=reference_state,
        )
        external_weights = _orbit_log_weights(
            environment,
            node,
            chain[1:],
            favourable=favourable,
            reference_state=reference_state,
        )
        full_weights = tuple(
            local + external
            for local, external in zip(local_weights, external_weights, strict=True)
        )
        local = parity_log_odds(local_weights)
        external = parity_log_odds(external_weights)
        total = parity_log_odds(full_weights)
        results[name] = external, local, total

    actual_external, actual_local, actual_total = results["actual"]
    favourable_external, favourable_local, favourable_total = results["favourable"]
    actual_rank = environment.forest.merge_rank[node]
    return NodePolarization(
        node=node,
        bucket_size=environment.forest.bucket_size[node],
        actual_rank=actual_rank,
        favourable_rank=min(actual_rank, Q_CRITICAL),
        strict_ancestor_count=len(chain) - 1,
        actual_external_message=actual_external,
        favourable_external_message=favourable_external,
        actual_local_log_odds=actual_local,
        favourable_local_log_odds=favourable_local,
        actual_total_log_odds=actual_total,
        favourable_total_log_odds=favourable_total,
        actual_identity_error=_identity_error(
            actual_total, actual_external, actual_local
        ),
        favourable_identity_error=_identity_error(
            favourable_total, favourable_external, favourable_local
        ),
    )


def analyze_pair_polarization(
    environment: HierarchicalSweepEnvironment,
    pair: PalmCorridorObservation,
    *,
    maximum_bucket_size: int,
    maximum_charge: float,
    message_thresholds: Sequence[float],
) -> PairPolarizationObservation:
    """Attach exact ancestor messages to the geometric corridor candidates."""

    thresholds = tuple(float(value) for value in message_thresholds)
    if not thresholds or any(
        not isfinite(value) or value < 0.0 for value in thresholds
    ):
        raise ValueError("message thresholds must be finite and nonnegative")
    if tuple(sorted(thresholds)) != thresholds or len(set(thresholds)) != len(
        thresholds
    ):
        raise ValueError("message thresholds must be strictly increasing")
    if maximum_bucket_size < 2 or maximum_charge < 0.0:
        raise ValueError("invalid geometric candidate filter")
    if pair.mode != "final_realized_event":
        raise ValueError("ancestor messages require a realized final corridor")

    candidates = []
    for cut in pair.cuts:
        if cut.node is None:
            raise AssertionError("a final-corridor cut must have a tree node")
        if (
            2 <= cut.bucket_size <= maximum_bucket_size
            and cut.favourable_charge <= maximum_charge
        ):
            candidates.append(analyze_node_polarization(environment, cut.node))

    actual_counts = tuple(
        sum(
            abs(candidate.actual_external_message) <= threshold
            for candidate in candidates
        )
        for threshold in thresholds
    )
    favourable_counts = tuple(
        sum(
            abs(candidate.favourable_external_message) <= threshold
            for candidate in candidates
        )
        for threshold in thresholds
    )
    return PairPolarizationObservation(
        repetition=pair.repetition,
        palm_weight=pair.palm_weight,
        first=pair.first,
        second=pair.second,
        corridor_node_count=len(pair.cuts),
        candidate_count=len(candidates),
        message_thresholds=thresholds,
        actual_bounded_message_counts=actual_counts,
        favourable_bounded_message_counts=favourable_counts,
        actual_infinite_message_count=sum(
            not isfinite(candidate.actual_external_message) for candidate in candidates
        ),
        favourable_infinite_message_count=sum(
            not isfinite(candidate.favourable_external_message)
            for candidate in candidates
        ),
        actual_local_reliability_sum=fsum(
            candidate.actual_local_reliability for candidate in candidates
        ),
        favourable_local_reliability_sum=fsum(
            candidate.favourable_local_reliability for candidate in candidates
        ),
        actual_total_reliability_sum=fsum(
            candidate.actual_total_reliability for candidate in candidates
        ),
        favourable_total_reliability_sum=fsum(
            candidate.favourable_total_reliability for candidate in candidates
        ),
        candidates=tuple(candidates),
    )


def _weighted_mean(
    observations: Sequence[PairPolarizationObservation],
    statistic: Callable[[PairPolarizationObservation], float],
) -> float:
    total_weight = sum(item.palm_weight for item in observations)
    if total_weight <= 0:
        raise ValueError("positive Palm mass is required")
    return (
        fsum(item.palm_weight * statistic(item) for item in observations) / total_weight
    )


def _cluster_jackknife_standard_error(
    observations: Sequence[PairPolarizationObservation],
    statistic: Callable[[PairPolarizationObservation], float],
    repetition_count: int,
) -> float | None:
    if repetition_count < 2:
        return None
    weights = [0] * repetition_count
    numerators = [0.0] * repetition_count
    for item in observations:
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


def _weighted_quantile(
    values_and_weights: Sequence[tuple[float, int]], quantile: float
) -> float | None:
    finite = sorted(
        (value, weight)
        for value, weight in values_and_weights
        if isfinite(value) and weight > 0
    )
    if not finite:
        return None
    total = sum(weight for _, weight in finite)
    target = quantile * total
    cumulative = 0
    for value, weight in finite:
        cumulative += weight
        if cumulative >= target:
            return value
    return finite[-1][0]


def summarize_polarization(
    observations: Sequence[PairPolarizationObservation],
    *,
    side_length: int,
    repetitions: int,
    p: float,
    distance_fraction: float,
    maximum_bucket_size: int,
    maximum_charge: float,
    seed: int,
) -> PolarizationPalmSummary:
    if not observations:
        raise ValueError("at least one pair observation is required")
    thresholds = observations[0].message_thresholds
    if any(item.message_thresholds != thresholds for item in observations):
        raise ValueError("all observations must use the same thresholds")
    total_weight = sum(item.palm_weight for item in observations)
    repetition_weights = [0] * repetitions
    for item in observations:
        repetition_weights[item.repetition] += item.palm_weight
    effective_repetitions = (
        total_weight
        * total_weight
        / sum(weight * weight for weight in repetition_weights if weight > 0)
    )
    candidate_mass = fsum(
        item.palm_weight * item.candidate_count for item in observations
    )
    if candidate_mass <= 0:
        raise RuntimeError("the diagnostic found no geometric candidate")

    actual_messages = []
    favourable_messages = []
    maximum_error = 0.0
    for item in observations:
        for candidate in item.candidates:
            actual_messages.append(
                (abs(candidate.actual_external_message), item.palm_weight)
            )
            favourable_messages.append(
                (abs(candidate.favourable_external_message), item.palm_weight)
            )
            for error in (
                candidate.actual_identity_error,
                candidate.favourable_identity_error,
            ):
                if error is not None:
                    maximum_error = max(maximum_error, error)

    def indexed_count(
        attribute: str, index: int
    ) -> Callable[[PairPolarizationObservation], float]:
        return lambda item: float(getattr(item, attribute)[index])

    actual_means = tuple(
        _weighted_mean(
            observations,
            indexed_count("actual_bounded_message_counts", index),
        )
        for index in range(len(thresholds))
    )
    actual_errors = tuple(
        _cluster_jackknife_standard_error(
            observations,
            indexed_count("actual_bounded_message_counts", index),
            repetitions,
        )
        for index in range(len(thresholds))
    )
    favourable_means = tuple(
        _weighted_mean(
            observations,
            indexed_count("favourable_bounded_message_counts", index),
        )
        for index in range(len(thresholds))
    )
    favourable_errors = tuple(
        _cluster_jackknife_standard_error(
            observations,
            indexed_count("favourable_bounded_message_counts", index),
            repetitions,
        )
        for index in range(len(thresholds))
    )

    def total(attribute: str) -> float:
        return fsum(
            item.palm_weight * getattr(item, attribute) for item in observations
        )

    return PolarizationPalmSummary(
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        critical_rank=Q_CRITICAL,
        final_rank=2.0 * p - 1.0,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_charge=maximum_charge,
        seed=seed,
        pair_observation_count=len(observations),
        total_palm_weight=total_weight,
        effective_repetition_count=effective_repetitions,
        weighted_mean_corridor_node_count=_weighted_mean(
            observations, lambda item: float(item.corridor_node_count)
        ),
        weighted_mean_candidate_count=_weighted_mean(
            observations, lambda item: float(item.candidate_count)
        ),
        jackknife_standard_error_candidate_count=(
            _cluster_jackknife_standard_error(
                observations,
                lambda item: float(item.candidate_count),
                repetitions,
            )
        ),
        message_thresholds=thresholds,
        weighted_mean_actual_bounded_message_counts=actual_means,
        jackknife_standard_errors_actual_bounded_message_counts=actual_errors,
        weighted_mean_favourable_bounded_message_counts=favourable_means,
        jackknife_standard_errors_favourable_bounded_message_counts=(favourable_errors),
        actual_infinite_message_fraction_among_candidates=(
            total("actual_infinite_message_count") / candidate_mass
        ),
        favourable_infinite_message_fraction_among_candidates=(
            total("favourable_infinite_message_count") / candidate_mass
        ),
        actual_finite_message_median=_weighted_quantile(actual_messages, 0.5),
        favourable_finite_message_median=_weighted_quantile(favourable_messages, 0.5),
        actual_finite_message_ninetieth_percentile=_weighted_quantile(
            actual_messages, 0.9
        ),
        favourable_finite_message_ninetieth_percentile=_weighted_quantile(
            favourable_messages, 0.9
        ),
        actual_mean_local_reliability_among_candidates=(
            total("actual_local_reliability_sum") / candidate_mass
        ),
        favourable_mean_local_reliability_among_candidates=(
            total("favourable_local_reliability_sum") / candidate_mass
        ),
        actual_mean_total_reliability_among_candidates=(
            total("actual_total_reliability_sum") / candidate_mass
        ),
        favourable_mean_total_reliability_among_candidates=(
            total("favourable_total_reliability_sum") / candidate_mass
        ),
        maximum_finite_log_odds_identity_error=maximum_error,
        interpretation=(
            "joint-Nishimori ancestor-message diagnostic; bounded B is not "
            "lateral screening and the favourable state law is not resampled"
        ),
    )


def run_diagnostic(
    *,
    side_length: int,
    repetitions: int,
    p: float,
    distance_fraction: float,
    maximum_bucket_size: int,
    maximum_charge: float,
    message_thresholds: Sequence[float],
    seed: int,
) -> PolarizationPalmSummary:
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0,1)")
    if Q_CRITICAL >= 2.0 * p - 1.0:
        raise ValueError("the critical rank must precede final censoring")
    rng = random.Random(seed)
    observations: list[PairPolarizationObservation] = []
    for repetition in range(repetitions):
        ranked_edges: tuple[RankedEdge, ...] = sample_ranked_edges(side_length, rng)
        environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
        pairs = final_realized_event_observations(
            environment.forest,
            repetition=repetition,
            p=p,
            distance_fraction=distance_fraction,
            rng=rng,
        )
        observations.extend(
            analyze_pair_polarization(
                environment,
                pair,
                maximum_bucket_size=maximum_bucket_size,
                maximum_charge=maximum_charge,
                message_thresholds=message_thresholds,
            )
            for pair in pairs
        )
    return summarize_polarization(
        observations,
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        distance_fraction=distance_fraction,
        maximum_bucket_size=maximum_bucket_size,
        maximum_charge=maximum_charge,
        seed=seed,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=8)
    parser.add_argument("--repetitions", type=int, default=30)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-bucket-size", type=int, default=8)
    parser.add_argument("--maximum-charge", type=float, default=1.0)
    parser.add_argument("--message-thresholds", default="1,2,4")
    parser.add_argument("--seed", type=int, default=20260722)
    arguments = parser.parse_args()
    summary = run_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        maximum_bucket_size=arguments.maximum_bucket_size,
        maximum_charge=arguments.maximum_charge,
        message_thresholds=tuple(
            float(value) for value in arguments.message_thresholds.split(",")
        ),
        seed=arguments.seed,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
