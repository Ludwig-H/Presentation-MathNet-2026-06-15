"""Finite-volume audit of critical, late, and root-separated pairs.

Pairs are classified by their minimax Kruskal rank in the final forest:

* ``early``: the LCA is below a fixed left critical window;
* ``critical``: the LCA belongs to that left window;
* ``late``: the LCA is above the critical rank but before time one;
* ``separate``: the vertices are in different final roots.

One pair per nonempty class and environment is sampled uniformly.  Its
contribution is weighted by the exact number of pairs in that class, so the
ratio estimates the annealed law of a uniform far pair conditional on the
class.  The same sweep outputs are reused for all selected classes.
Critical-minus-late errors use a paired delete-one-environment jackknife.

This is a diagnostic of HF-S2, not a proof of favorable domination.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, sqrt
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    sample_ranked_edges,
    triangular_torus_distance,
)
from joint_hierarchical_sweep import (
    HierarchicalSweepEnvironment,
    pair_relation,
    unbiased_rademacher_mean_square,
)


PAIR_CLASSES = ("early", "critical", "late", "separate")


@dataclass(frozen=True)
class PairClassEstimate:
    """Weighted conditional estimates for one pair class."""

    estimated_pair_mass: float
    contributing_environments: int
    annealed_correlation: float
    annealed_correlation_jackknife_se: float
    conditional_second_moment: float
    conditional_second_moment_jackknife_se: float


@dataclass(frozen=True)
class FavorabilityContrast:
    """Paired critical-minus-late second-moment contrast."""

    critical_minus_late_second_moment: float
    critical_minus_late_jackknife_se: float


@dataclass(frozen=True)
class FavorabilitySummary:
    """Comparison of pair classes for systematic hierarchical sweeps."""

    side_length: int
    repetitions: int
    sweeps_per_environment: int
    passes_per_sweep: int
    distance_fraction: float
    critical_window_width: float
    results: dict[str, dict[str, dict[str, PairClassEstimate]]]
    contrasts: dict[str, dict[str, FavorabilityContrast]]


def classify_pair(
    forest: CriticalKruskalForest,
    first: int,
    second: int,
    critical_window_width: float,
) -> str:
    """Classify one distinct pair from its final Kruskal forest."""

    if first == second:
        raise ValueError("the pair must be distinct")
    if not 0.0 < critical_window_width < Q_CRITICAL:
        raise ValueError("critical_window_width must belong to (0, q_c)")
    if not forest.connected(first, second):
        return "separate"
    rank = forest.merge_rank[forest.tree_lca(first, second)]
    if rank < Q_CRITICAL - critical_window_width:
        return "early"
    if rank <= Q_CRITICAL:
        return "critical"
    return "late"


def far_pair_classes(
    forest: CriticalKruskalForest,
    distance_fraction: float,
    critical_window_width: float,
) -> dict[str, list[tuple[int, int]]]:
    """Enumerate unordered far pairs and partition them exactly."""

    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0, 1)")
    classes = {name: [] for name in PAIR_CLASSES}
    minimum_distance = distance_fraction * forest.side_length
    for first in range(forest.vertex_count):
        for second in range(first + 1, forest.vertex_count):
            if triangular_torus_distance(
                first, second, forest.side_length
            ) < minimum_distance:
                continue
            name = classify_pair(
                forest, first, second, critical_window_width
            )
            classes[name].append((first, second))
    return classes


def _weighted_ratio(records: Sequence[tuple[float, float]]) -> float:
    denominator = fsum(weight for weight, _ in records)
    if denominator <= 0.0:
        return float("nan")
    return fsum(weight * value for weight, value in records) / denominator


def _jackknife_ratio_standard_error(
    records: Sequence[tuple[float, float]],
) -> float:
    """Delete-one-environment jackknife error of a weighted ratio."""

    if len(records) < 2:
        return 0.0
    total_weight = fsum(weight for weight, _ in records)
    total_value = fsum(weight * value for weight, value in records)
    leave_one = []
    for weight, value in records:
        denominator = total_weight - weight
        if denominator > 0.0:
            leave_one.append((total_value - weight * value) / denominator)
    if len(leave_one) < 2:
        return 0.0
    mean = fsum(leave_one) / len(leave_one)
    variance = (len(leave_one) - 1.0) / len(leave_one) * fsum(
        (value - mean) ** 2 for value in leave_one
    )
    return sqrt(variance)


def _jackknife_ratio_difference_standard_error(
    first: Sequence[tuple[float, float]],
    second: Sequence[tuple[float, float]],
) -> float:
    """Paired delete-one error for the difference of two weighted ratios."""

    if len(first) != len(second):
        raise ValueError("paired records must have the same length")
    if len(first) < 2:
        return 0.0
    first_weight = fsum(weight for weight, _ in first)
    first_value = fsum(weight * value for weight, value in first)
    second_weight = fsum(weight for weight, _ in second)
    second_value = fsum(weight * value for weight, value in second)
    leave_one = []
    for (weight_a, value_a), (weight_b, value_b) in zip(
        first, second, strict=True
    ):
        denominator_a = first_weight - weight_a
        denominator_b = second_weight - weight_b
        if denominator_a > 0.0 and denominator_b > 0.0:
            leave_one.append(
                (first_value - weight_a * value_a) / denominator_a
                - (second_value - weight_b * value_b) / denominator_b
            )
    if len(leave_one) < 2:
        return 0.0
    mean = fsum(leave_one) / len(leave_one)
    variance = (len(leave_one) - 1.0) / len(leave_one) * fsum(
        (value - mean) ** 2 for value in leave_one
    )
    return sqrt(variance)


def run_favorability_diagnostic(
    side_length: int,
    repetitions: int,
    sweeps_per_environment: int,
    passes_per_sweep: int,
    p_values: Sequence[float],
    distance_fraction: float,
    critical_window_width: float,
    seed: int,
) -> FavorabilitySummary:
    """Estimate the HF-S2 profile with exact pair-count reweighting."""

    if repetitions <= 0 or sweeps_per_environment < 2:
        raise ValueError("positive repetitions and at least two sweeps are required")
    if passes_per_sweep <= 0:
        raise ValueError("passes_per_sweep must be positive")
    geometry_rng = random.Random(seed)
    stored: dict[
        tuple[float, str, str], list[tuple[float, float, float]]
    ] = {
        (p, order, name): []
        for p in p_values
        for order in ("bottom-up", "top-down")
        for name in PAIR_CLASSES
    }
    far_pair_count: int | None = None

    for repetition in range(repetitions):
        ranked_edges = sample_ranked_edges(side_length, geometry_rng)
        for p_index, p in enumerate(p_values):
            environment = HierarchicalSweepEnvironment(
                side_length, ranked_edges, p
            )
            classes = far_pair_classes(
                environment.forest,
                distance_fraction,
                critical_window_width,
            )
            current_far_count = sum(len(pairs) for pairs in classes.values())
            if far_pair_count is None:
                far_pair_count = current_far_count
            elif far_pair_count != current_far_count:
                raise AssertionError("the number of geometric far pairs changed")

            chosen: dict[str, tuple[int, int]] = {}
            for class_index, name in enumerate(PAIR_CLASSES):
                pairs = classes[name]
                if not pairs:
                    continue
                pair_seed = (
                    seed
                    ^ ((repetition + 1) * 0x9E3779B97F4A7C15)
                    ^ ((p_index + 1) * 0xBF58476D1CE4E5B9)
                    ^ ((class_index + 1) * 0x94D049BB133111EB)
                ) & ((1 << 64) - 1)
                chosen[name] = random.Random(pair_seed).choice(pairs)

            for order_index, order in enumerate(("bottom-up", "top-down")):
                sweep_seed = (
                    seed
                    ^ ((repetition + 1) * 0xD2B74407B1CE6E93)
                    ^ ((p_index + 1) * 0xCA5A826395121157)
                    ^ (order_index * 0x9E3779B185EBCA87)
                ) & ((1 << 64) - 1)
                sweep_rng = random.Random(sweep_seed)
                states = tuple(
                    environment.sweep(
                        order, sweep_rng, passes_per_sweep
                    )
                    for _ in range(sweeps_per_environment)
                )
                for name in PAIR_CLASSES:
                    if name not in chosen:
                        stored[(p, order, name)].append((0.0, 0.0, 0.0))
                        continue
                    first, second = chosen[name]
                    samples = tuple(
                        pair_relation(state, first, second)
                        for state in states
                    )
                    correlation = fsum(samples) / len(samples)
                    second_moment = unbiased_rademacher_mean_square(samples)
                    stored[(p, order, name)].append(
                        (float(len(classes[name])), correlation, second_moment)
                    )

    if far_pair_count is None or far_pair_count == 0:
        raise RuntimeError("no far pair was available")

    results: dict[str, dict[str, dict[str, PairClassEstimate]]] = {}
    contrasts: dict[str, dict[str, FavorabilityContrast]] = {}
    for p in p_values:
        p_result: dict[str, dict[str, PairClassEstimate]] = {}
        p_contrasts: dict[str, FavorabilityContrast] = {}
        for order in ("bottom-up", "top-down"):
            order_result: dict[str, PairClassEstimate] = {}
            for name in PAIR_CLASSES:
                values = stored[(p, order, name)]
                correlation_records = [
                    (weight, correlation)
                    for weight, correlation, _ in values
                ]
                second_records = [
                    (weight, second)
                    for weight, _, second in values
                ]
                total_class_weight = fsum(
                    weight for weight, _, _ in values
                )
                order_result[name] = PairClassEstimate(
                    estimated_pair_mass=total_class_weight
                    / (repetitions * far_pair_count),
                    contributing_environments=sum(
                        weight > 0.0 for weight, _, _ in values
                    ),
                    annealed_correlation=_weighted_ratio(correlation_records),
                    annealed_correlation_jackknife_se=(
                        _jackknife_ratio_standard_error(
                            correlation_records
                        )
                    ),
                    conditional_second_moment=_weighted_ratio(second_records),
                    conditional_second_moment_jackknife_se=(
                        _jackknife_ratio_standard_error(
                            second_records
                        )
                    ),
                )
            p_result[order] = order_result
            critical_records = [
                (weight, second)
                for weight, _, second in stored[(p, order, "critical")]
            ]
            late_records = [
                (weight, second)
                for weight, _, second in stored[(p, order, "late")]
            ]
            p_contrasts[order] = FavorabilityContrast(
                critical_minus_late_second_moment=(
                    _weighted_ratio(critical_records)
                    - _weighted_ratio(late_records)
                ),
                critical_minus_late_jackknife_se=(
                    _jackknife_ratio_difference_standard_error(
                        critical_records, late_records
                    )
                ),
            )
        results[f"{p:.12g}"] = p_result
        contrasts[f"{p:.12g}"] = p_contrasts

    return FavorabilitySummary(
        side_length=side_length,
        repetitions=repetitions,
        sweeps_per_environment=sweeps_per_environment,
        passes_per_sweep=passes_per_sweep,
        distance_fraction=distance_fraction,
        critical_window_width=critical_window_width,
        results=results,
        contrasts=contrasts,
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
    parser.add_argument("--critical-window", type=float, default=0.05)
    parser.add_argument("--p-values", default="0.8")
    parser.add_argument("--seed", type=int, default=20260719)
    arguments = parser.parse_args()
    summary = run_favorability_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        sweeps_per_environment=arguments.sweeps,
        passes_per_sweep=arguments.passes,
        p_values=_parse_p_values(arguments.p_values),
        distance_fraction=arguments.distance_fraction,
        critical_window_width=arguments.critical_window,
        seed=arguments.seed,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
