"""Exact odd-sector L2 dissipation on tiny real-rank corridors.

This module evaluates the global projection pivot which survives the local
Feynman--Kac no-go.  Conditional on one genuine Nishimori/Kruskal hierarchy,
it enumerates the full posterior ``pi_D`` and a uniformly sampled distant
pair.  If the pair lies in one final tree, it integrates successively larger
sets of corridor orientations and records the exact Pythagorean losses

    ||M_{k-1}||_2^2 - ||M_k||_2^2 = ||M_{k-1}-M_k||_2^2.

Every collapsed update is implemented as a conditional expectation over the
linear span of the corresponding cluster-flip masks.  No boundary-state
compression, criticalization, Feynman--Kac envelope, or independent replica
environment is used.  Disconnected pairs have persistence zero after the
exact final-root recolouring.

The enumeration has ``2^(L^2)`` states and is intentionally restricted to
``L <= 4``.  Its output is a finite diagnostic, not an asymptotic or
weak-recovery result.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import exp, fsum, sqrt
from statistics import fmean
from typing import Sequence

from critical_pair_path_geometry import (
    sample_ranked_edges,
    triangular_torus_distance,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


@dataclass(frozen=True)
class NestedProjectionStep:
    """One enlargement of the collapsed orientation set."""

    depth: int
    added_nodes: tuple[int, ...]
    added_ranks: tuple[float, ...]
    generated_flip_rank: int
    norm_before: float
    norm_after: float
    absolute_loss: float
    difference_norm_square: float
    pythagorean_error: float
    relative_loss: float
    logarithmic_loss: float


@dataclass(frozen=True)
class PairProjectionDiagnostic:
    """Exact projection data for one environment and one distant pair."""

    repetition: int
    first: int
    second: int
    distance: int
    connected_in_final_forest: bool
    path_node_count: int
    posterior_state_count: int
    positive_posterior_state_count: int
    posterior_mean_square: float
    final_collapsed_persistence: float
    cumulative_logarithmic_loss: float | None
    positive_relative_loss_count: int
    effective_absolute_loss_step_count: float
    pythagorean_error: float
    posterior_below_collapsed: bool
    steps: tuple[NestedProjectionStep, ...]


@dataclass(frozen=True)
class DepthProjectionSummary:
    """Mean loss among connected corridors which reach one depth."""

    depth: int
    observation_count: int
    mean_absolute_loss: float
    mean_relative_loss: float
    mean_logarithmic_loss: float
    positive_relative_loss_fraction: float
    minimum_relative_loss: float
    tenth_percentile_relative_loss: float
    median_relative_loss: float
    ninetieth_percentile_relative_loss: float


@dataclass(frozen=True)
class NestedProjectionSummary:
    """Environment-level means for the tiny-volume exact diagnostic."""

    side_length: int
    repetitions: int
    p: float
    distance_fraction: float
    seed: int
    connected_environment_count: int
    connected_fraction: float
    mean_unconditional_collapsed_persistence: float
    standard_error_unconditional_collapsed_persistence: float | None
    mean_connected_collapsed_persistence: float | None
    standard_error_connected_collapsed_persistence: float | None
    mean_connected_posterior_square: float | None
    mean_connected_path_node_count: float | None
    mean_connected_positive_loss_count: float | None
    mean_connected_effective_loss_step_count: float | None
    mean_connected_cumulative_logarithmic_loss: float | None
    mean_connected_pre_lca_logarithmic_loss: float | None
    mean_connected_lca_logarithmic_loss: float | None
    mean_connected_lca_absolute_loss_share: float | None
    positive_dissipation_environment_count: int
    mean_positive_dissipation_largest_absolute_loss_share: float | None
    second_pre_lca_observation_count: int
    second_pre_lca_positive_loss_fraction: float | None
    second_pre_lca_loss_above_one_percent_fraction: float | None
    second_pre_lca_energy_weighted_relative_loss: float | None
    mean_second_pre_lca_logarithmic_loss: float | None
    median_second_pre_lca_relative_loss: float | None
    ninetieth_percentile_second_pre_lca_relative_loss: float | None
    pre_lca_depth_profile: tuple[DepthProjectionSummary, ...]
    maximum_pythagorean_error: float
    every_posterior_below_collapsed: bool
    exact_full_posterior_enumeration: bool
    ranks_are_realized: bool
    boundary_compression_used: bool
    weak_recovery_claimed: bool
    interpretation: str


def _sample_uniform_distant_pair(
    side_length: int,
    distance_fraction: float,
    rng: random.Random,
) -> tuple[int, int]:
    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0,1)")
    vertex_count = side_length * side_length
    minimum_distance = distance_fraction * side_length
    for _ in range(100_000):
        first = rng.randrange(vertex_count)
        second = rng.randrange(vertex_count - 1)
        if second >= first:
            second += 1
        if triangular_torus_distance(first, second, side_length) >= minimum_distance:
            return first, second
    raise RuntimeError("unable to sample a distant pair")


def posterior_weights(
    environment: HierarchicalSweepEnvironment,
) -> tuple[float, ...]:
    """Enumerate and normalize the exact conditional law ``pi_D``."""

    if environment.vertex_count > 16:
        raise ValueError("full posterior enumeration is restricted to L <= 4")
    log_weights = tuple(
        environment.full_log_weight(state)
        for state in range(1 << environment.vertex_count)
    )
    maximum = max(log_weights)
    if maximum == float("-inf"):
        raise AssertionError("the conditional posterior has zero total mass")
    raw = tuple(
        0.0 if value == float("-inf") else exp(value - maximum) for value in log_weights
    )
    normalizer = fsum(raw)
    return tuple(value / normalizer for value in raw)


def pair_character(state_count: int, first: int, second: int) -> tuple[float, ...]:
    """Return the spin--spin character as a function on bit states."""

    if state_count <= 0 or state_count & (state_count - 1):
        raise ValueError("state_count must be a positive power of two")
    vertex_count = state_count.bit_length() - 1
    if first == second or not (
        0 <= first < vertex_count and 0 <= second < vertex_count
    ):
        raise ValueError("pair endpoints are distinct vertices of the state space")
    return tuple(
        1.0 if (((state >> first) ^ (state >> second)) & 1) == 0 else -1.0
        for state in range(state_count)
    )


def weighted_norm_square(values: Sequence[float], weights: Sequence[float]) -> float:
    if len(values) != len(weights):
        raise ValueError("values and weights must have the same length")
    return fsum(weight * value * value for value, weight in zip(values, weights))


def _xor_basis(generators: Sequence[int]) -> tuple[int, ...]:
    pivots: dict[int, int] = {}
    for generator in generators:
        value = generator
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return tuple(pivots[pivot] for pivot in sorted(pivots, reverse=True))


def _quotient_key(state: int, basis: Sequence[int]) -> int:
    value = state
    for vector in basis:
        pivot = vector.bit_length() - 1
        if (value >> pivot) & 1:
            value ^= vector
    return value


def collapsed_projection(
    values: Sequence[float],
    weights: Sequence[float],
    generators: Sequence[int],
) -> tuple[float, ...]:
    """Condition on cosets of the flip group generated by ``generators``."""

    if len(values) != len(weights):
        raise ValueError("values and weights must have the same length")
    basis = _xor_basis(generators)
    mass: dict[int, float] = {}
    signed_mass: dict[int, float] = {}
    for state, (value, weight) in enumerate(zip(values, weights, strict=True)):
        key = _quotient_key(state, basis)
        mass[key] = mass.get(key, 0.0) + weight
        signed_mass[key] = signed_mass.get(key, 0.0) + weight * value
    means = {
        key: (0.0 if denominator == 0.0 else signed_mass[key] / denominator)
        for key, denominator in mass.items()
    }
    return tuple(means[_quotient_key(state, basis)] for state in range(len(values)))


def _corridor_depth_packets(
    environment: HierarchicalSweepEnvironment,
    first: int,
    second: int,
) -> tuple[tuple[int, ...], ...]:
    forest = environment.forest
    lca = forest.tree_lca(first, second)

    def chain_to_lca(vertex: int) -> tuple[int, ...]:
        result = []
        node = vertex
        while node != lca:
            node = forest.tree_parent[node]
            if node == -1:
                raise AssertionError("the LCA is missing from an endpoint chain")
            result.append(node)
        return tuple(result)

    first_chain = chain_to_lca(first)
    second_chain = chain_to_lca(second)
    packets = []
    strict_depth_count = max(len(first_chain), len(second_chain)) - 1
    for depth in range(max(0, strict_depth_count)):
        nodes = []
        if depth < len(first_chain) - 1:
            nodes.append(first_chain[depth])
        if depth < len(second_chain) - 1:
            nodes.append(second_chain[depth])
        if nodes:
            packets.append(tuple(nodes))
    packets.append((lca,))
    flattened = tuple(node for packet in packets for node in packet)
    expected = set(forest.pair_path_nodes(first, second))
    if len(flattened) != len(set(flattened)) or set(flattened) != expected:
        raise AssertionError("depth packets do not partition the pair corridor")
    return tuple(packets)


def _node_generators(
    environment: HierarchicalSweepEnvironment, nodes: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        mask for node in nodes for mask in environment.proposal_masks(node) if mask
    )


def analyze_pair(
    environment: HierarchicalSweepEnvironment,
    first: int,
    second: int,
    repetition: int,
) -> PairProjectionDiagnostic:
    """Compute the nested collapsed martingale for one fixed pair."""

    distance = triangular_torus_distance(first, second, environment.side_length)
    state_count = 1 << environment.vertex_count
    if not environment.forest.connected(first, second):
        return PairProjectionDiagnostic(
            repetition=repetition,
            first=first,
            second=second,
            distance=distance,
            connected_in_final_forest=False,
            path_node_count=0,
            posterior_state_count=state_count,
            positive_posterior_state_count=0,
            posterior_mean_square=0.0,
            final_collapsed_persistence=0.0,
            cumulative_logarithmic_loss=None,
            positive_relative_loss_count=1,
            effective_absolute_loss_step_count=1.0,
            pythagorean_error=0.0,
            posterior_below_collapsed=True,
            steps=(),
        )

    weights = posterior_weights(environment)
    values = pair_character(state_count, first, second)
    posterior_mean = fsum(
        weight * value for value, weight in zip(values, weights, strict=True)
    )
    posterior_square = posterior_mean * posterior_mean
    packets = _corridor_depth_packets(environment, first, second)
    generators: list[int] = []
    steps = []
    current = values
    initial_norm = weighted_norm_square(current, weights)
    previous_norm = initial_norm
    absolute_loss_sum = 0.0
    for depth, packet in enumerate(packets):
        generators.extend(_node_generators(environment, packet))
        projected = collapsed_projection(current, weights, generators)
        norm_after = weighted_norm_square(projected, weights)
        absolute_loss = previous_norm - norm_after
        if absolute_loss < 0.0 and absolute_loss > -1e-12:
            absolute_loss = 0.0
        if absolute_loss < 0.0:
            raise AssertionError("a conditional expectation increased the L2 norm")
        difference_norm_square = weighted_norm_square(
            tuple(
                before - after for before, after in zip(current, projected, strict=True)
            ),
            weights,
        )
        pythagorean_step_error = abs(absolute_loss - difference_norm_square)
        relative_loss = 0.0 if previous_norm == 0.0 else absolute_loss / previous_norm
        if norm_after == 0.0:
            logarithmic_loss = float("inf")
        elif previous_norm == 0.0:
            logarithmic_loss = 0.0
        else:
            from math import log

            logarithmic_loss = -log(norm_after / previous_norm)
        basis_rank = len(_xor_basis(generators))
        steps.append(
            NestedProjectionStep(
                depth=depth,
                added_nodes=packet,
                added_ranks=tuple(
                    environment.forest.merge_rank[node] for node in packet
                ),
                generated_flip_rank=basis_rank,
                norm_before=previous_norm,
                norm_after=norm_after,
                absolute_loss=absolute_loss,
                difference_norm_square=difference_norm_square,
                pythagorean_error=pythagorean_step_error,
                relative_loss=relative_loss,
                logarithmic_loss=logarithmic_loss,
            )
        )
        current = projected
        previous_norm = norm_after
        absolute_loss_sum += absolute_loss

    final_norm = previous_norm
    losses = tuple(step.absolute_loss for step in steps)
    squared_loss_sum = fsum(loss * loss for loss in losses)
    effective_count = (
        0.0
        if squared_loss_sum == 0.0
        else absolute_loss_sum * absolute_loss_sum / squared_loss_sum
    )
    finite_relative_losses = tuple(
        step.relative_loss for step in steps if step.relative_loss < 1.0
    )
    if any(step.relative_loss >= 1.0 for step in steps):
        cumulative_logarithmic_loss = float("inf")
    else:
        from math import log

        cumulative_logarithmic_loss = fsum(
            -log(1.0 - loss) for loss in finite_relative_losses
        )
    telescoping_error = abs(initial_norm - final_norm - absolute_loss_sum)
    pythagorean_error = max(
        (telescoping_error, *(step.pythagorean_error for step in steps))
    )
    return PairProjectionDiagnostic(
        repetition=repetition,
        first=first,
        second=second,
        distance=distance,
        connected_in_final_forest=True,
        path_node_count=sum(len(packet) for packet in packets),
        posterior_state_count=state_count,
        positive_posterior_state_count=sum(weight > 0.0 for weight in weights),
        posterior_mean_square=posterior_square,
        final_collapsed_persistence=final_norm,
        cumulative_logarithmic_loss=cumulative_logarithmic_loss,
        positive_relative_loss_count=sum(step.absolute_loss > 1e-14 for step in steps),
        effective_absolute_loss_step_count=effective_count,
        pythagorean_error=pythagorean_error,
        posterior_below_collapsed=posterior_square <= final_norm + 1e-12,
        steps=tuple(steps),
    )


def _standard_error(values: Sequence[float]) -> float | None:
    if len(values) < 2:
        return None
    mean = fmean(values)
    return sqrt(
        fsum((value - mean) ** 2 for value in values)
        / (len(values) * (len(values) - 1))
    )


def summarize(
    diagnostics: Sequence[PairProjectionDiagnostic],
    side_length: int,
    p: float,
    distance_fraction: float,
    seed: int,
) -> NestedProjectionSummary:
    if not diagnostics:
        raise ValueError("at least one diagnostic is required")
    connected = tuple(item for item in diagnostics if item.connected_in_final_forest)
    unconditional = tuple(item.final_collapsed_persistence for item in diagnostics)

    def connected_mean(attribute: str) -> float | None:
        if not connected:
            return None
        return fmean(float(getattr(item, attribute)) for item in connected)

    def mean_step_statistic(statistic) -> float | None:
        if not connected:
            return None
        return fmean(statistic(item) for item in connected)

    def pre_lca_logarithmic_loss(item: PairProjectionDiagnostic) -> float:
        return fsum(step.logarithmic_loss for step in item.steps[:-1])

    def lca_logarithmic_loss(item: PairProjectionDiagnostic) -> float:
        return item.steps[-1].logarithmic_loss

    def lca_absolute_loss_share(item: PairProjectionDiagnostic) -> float:
        total = fsum(step.absolute_loss for step in item.steps)
        return 0.0 if total == 0.0 else item.steps[-1].absolute_loss / total

    connected_persistence = tuple(
        item.final_collapsed_persistence for item in connected
    )

    def empirical_quantile(values: Sequence[float], probability: float) -> float:
        if not values:
            raise ValueError("an empirical quantile requires data")
        if not 0.0 <= probability <= 1.0:
            raise ValueError("probability must belong to [0,1]")
        ordered = sorted(values)
        return ordered[int(probability * (len(ordered) - 1))]

    maximum_depth_count = max(
        (max(0, len(item.steps) - 1) for item in connected), default=0
    )
    pre_lca_depth_profile = []
    for depth in range(maximum_depth_count):
        depth_steps = tuple(
            item.steps[depth] for item in connected if depth < len(item.steps) - 1
        )
        relative_losses = tuple(step.relative_loss for step in depth_steps)
        pre_lca_depth_profile.append(
            DepthProjectionSummary(
                depth=depth,
                observation_count=len(depth_steps),
                mean_absolute_loss=fmean(step.absolute_loss for step in depth_steps),
                mean_relative_loss=fmean(step.relative_loss for step in depth_steps),
                mean_logarithmic_loss=fmean(
                    step.logarithmic_loss for step in depth_steps
                ),
                positive_relative_loss_fraction=(
                    sum(loss > 1e-14 for loss in relative_losses) / len(relative_losses)
                ),
                minimum_relative_loss=min(relative_losses),
                tenth_percentile_relative_loss=empirical_quantile(relative_losses, 0.1),
                median_relative_loss=empirical_quantile(relative_losses, 0.5),
                ninetieth_percentile_relative_loss=empirical_quantile(
                    relative_losses, 0.9
                ),
            )
        )
    positive_dissipation = tuple(
        item
        for item in connected
        if fsum(step.absolute_loss for step in item.steps) > 1e-14
    )
    largest_absolute_loss_shares = tuple(
        max(step.absolute_loss for step in item.steps)
        / fsum(step.absolute_loss for step in item.steps)
        for item in positive_dissipation
    )
    second_pre_lca_steps = tuple(
        item.steps[1] for item in connected if len(item.steps) >= 3
    )
    second_pre_lca_relative_losses = tuple(
        step.relative_loss for step in second_pre_lca_steps
    )
    second_pre_lca_energy_before = fsum(
        step.norm_before for step in second_pre_lca_steps
    )
    return NestedProjectionSummary(
        side_length=side_length,
        repetitions=len(diagnostics),
        p=p,
        distance_fraction=distance_fraction,
        seed=seed,
        connected_environment_count=len(connected),
        connected_fraction=len(connected) / len(diagnostics),
        mean_unconditional_collapsed_persistence=fmean(unconditional),
        standard_error_unconditional_collapsed_persistence=_standard_error(
            unconditional
        ),
        mean_connected_collapsed_persistence=(
            None if not connected_persistence else fmean(connected_persistence)
        ),
        standard_error_connected_collapsed_persistence=_standard_error(
            connected_persistence
        ),
        mean_connected_posterior_square=connected_mean("posterior_mean_square"),
        mean_connected_path_node_count=connected_mean("path_node_count"),
        mean_connected_positive_loss_count=connected_mean(
            "positive_relative_loss_count"
        ),
        mean_connected_effective_loss_step_count=connected_mean(
            "effective_absolute_loss_step_count"
        ),
        mean_connected_cumulative_logarithmic_loss=connected_mean(
            "cumulative_logarithmic_loss"
        ),
        mean_connected_pre_lca_logarithmic_loss=mean_step_statistic(
            pre_lca_logarithmic_loss
        ),
        mean_connected_lca_logarithmic_loss=mean_step_statistic(lca_logarithmic_loss),
        mean_connected_lca_absolute_loss_share=mean_step_statistic(
            lca_absolute_loss_share
        ),
        positive_dissipation_environment_count=len(positive_dissipation),
        mean_positive_dissipation_largest_absolute_loss_share=(
            None
            if not largest_absolute_loss_shares
            else fmean(largest_absolute_loss_shares)
        ),
        second_pre_lca_observation_count=len(second_pre_lca_steps),
        second_pre_lca_positive_loss_fraction=(
            None
            if not second_pre_lca_steps
            else sum(loss > 1e-14 for loss in second_pre_lca_relative_losses)
            / len(second_pre_lca_relative_losses)
        ),
        second_pre_lca_loss_above_one_percent_fraction=(
            None
            if not second_pre_lca_steps
            else sum(loss > 0.01 for loss in second_pre_lca_relative_losses)
            / len(second_pre_lca_relative_losses)
        ),
        second_pre_lca_energy_weighted_relative_loss=(
            None
            if second_pre_lca_energy_before == 0.0
            else fsum(step.absolute_loss for step in second_pre_lca_steps)
            / second_pre_lca_energy_before
        ),
        mean_second_pre_lca_logarithmic_loss=(
            None
            if not second_pre_lca_steps
            else fmean(step.logarithmic_loss for step in second_pre_lca_steps)
        ),
        median_second_pre_lca_relative_loss=(
            None
            if not second_pre_lca_steps
            else empirical_quantile(second_pre_lca_relative_losses, 0.5)
        ),
        ninetieth_percentile_second_pre_lca_relative_loss=(
            None
            if not second_pre_lca_steps
            else empirical_quantile(second_pre_lca_relative_losses, 0.9)
        ),
        pre_lca_depth_profile=tuple(pre_lca_depth_profile),
        maximum_pythagorean_error=max(item.pythagorean_error for item in diagnostics),
        every_posterior_below_collapsed=all(
            item.posterior_below_collapsed for item in diagnostics
        ),
        exact_full_posterior_enumeration=True,
        ranks_are_realized=True,
        boundary_compression_used=False,
        weak_recovery_claimed=False,
        interpretation=(
            "finite exact target-specific L2 diagnostic; constants prevent a global "
            "operator contraction and no thermodynamic conclusion is claimed"
        ),
    )


def run_diagnostic(
    side_length: int = 4,
    repetitions: int = 12,
    p: float = 0.805,
    distance_fraction: float = 0.5,
    seed: int = 20260726,
) -> tuple[NestedProjectionSummary, tuple[PairProjectionDiagnostic, ...]]:
    if side_length != 4:
        raise ValueError("the triangular torus diagnostic requires side_length=4")
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    rng = random.Random(seed)
    diagnostics = []
    for repetition in range(repetitions):
        ranked_edges = sample_ranked_edges(side_length, rng)
        environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
        first, second = _sample_uniform_distant_pair(
            side_length, distance_fraction, rng
        )
        diagnostics.append(analyze_pair(environment, first, second, repetition))
    result = tuple(diagnostics)
    return (
        summarize(result, side_length, p, distance_fraction, seed),
        result,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=4)
    parser.add_argument("--repetitions", type=int, default=12)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument("--include-pairs", action="store_true")
    arguments = parser.parse_args()
    summary, diagnostics = run_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        seed=arguments.seed,
    )
    payload: dict[str, object] = {"summary": asdict(summary)}
    if arguments.include_pairs:
        payload["pairs"] = [asdict(item) for item in diagnostics]
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
