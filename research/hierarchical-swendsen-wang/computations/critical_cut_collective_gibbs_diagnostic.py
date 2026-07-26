"""Exact collective critical-cut Gibbs persistence on tiny tori.

Fix one realization of the uniform edge ranks.  The forest at the
triangular-lattice critical rank ``q_c`` partitions the vertices into
critical blocks ``C_a``.  Inside each final Kruskal root, this module
enumerates all joint block orientations

    z = (z_a) in {-1,+1}^k

and assigns to the corresponding spin state the exact conditional weight

    exp(environment.full_log_weight(state)).

The resulting matrix ``M_ab = E[z_a z_b]`` gives the collective heat-bath
persistence

    G_L^c = n^-2 sum_R sum_{a,b in R} |C_a| |C_b| M_ab^2.

Correlations between distinct final roots are set to zero, as required by
their independent fair global recolourings.  The diagonal contribution is

    S_L^c = n^-2 sum_a |C_a|^2.

Enumeration is exponential in the number of blocks.  ``maximum_block_count``
therefore excludes an environment if either its total number of critical
blocks or the number in one final root exceeds the cap.  Any such exclusion
is reported as selection bias.  This is a reproducible finite-volume
diagnostic, not an asymptotic estimate or a weak-recovery theorem.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from dataclasses import asdict, dataclass
from math import exp, fsum, sqrt
from statistics import fmean
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    RankedEdge,
    sample_ranked_edges,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


BASELINE_P = 0.809439
P_SW = (1.0 + Q_CRITICAL) / 2.0
ROOT_CAP_SKIP = "maximum_root_block_count_exceeds_maximum_block_count"
TOTAL_CAP_SKIP = "total_critical_block_count_exceeds_maximum_block_count"


@dataclass(frozen=True)
class CriticalCutPartition:
    """Critical blocks grouped by their final Kruskal component."""

    blocks: tuple[tuple[int, ...], ...]
    block_masks: tuple[int, ...]
    block_sizes: tuple[int, ...]
    final_components: tuple[tuple[int, ...], ...]
    final_component_block_indices: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class RootCollectiveGibbsDiagnostic:
    """Exact orientation correlation matrix inside one final root."""

    final_component_index: int
    representative_vertex: int
    final_component_size: int
    block_indices: tuple[int, ...]
    block_sizes: tuple[int, ...]
    orientation_state_count: int
    positive_orientation_state_count: int
    correlation_matrix: tuple[tuple[float, ...], ...]
    unnormalized_persistence: float
    diagonal_unnormalized_persistence: float
    maximum_diagonal_error: float
    maximum_symmetry_error: float
    maximum_global_flip_log_weight_error: float


@dataclass(frozen=True)
class EnvironmentCollectiveGibbsDiagnostic:
    """Exact result, or an explicit cap exclusion, for one environment."""

    repetition: int
    total_critical_block_count: int
    final_root_count: int
    maximum_root_block_count: int
    used: bool
    skip_reason: str | None
    critical_diagonal_persistence: float
    collective_persistence: float | None
    off_diagonal_persistence: float | None
    maximum_decomposition_error: float | None
    maximum_diagonal_error: float | None
    maximum_symmetry_error: float | None
    maximum_global_flip_log_weight_error: float | None
    roots: tuple[RootCollectiveGibbsDiagnostic, ...]


@dataclass(frozen=True)
class ScalarEstimate:
    """Environment-clustered sample mean and standard error."""

    mean: float
    standard_error: float | None


@dataclass(frozen=True)
class CriticalCutCollectiveGibbsSummary:
    """Monte Carlo summary with cap selection kept fully visible."""

    side_length: int
    vertex_count: int
    p: float
    critical_rank: float
    final_rank: float
    repetitions: int
    seed: int
    maximum_block_count: int
    used_environment_count: int
    skipped_environment_count: int
    skipped_for_root_block_cap_count: int
    skipped_for_total_block_cap_count: int
    all_sample_critical_diagonal_persistence: ScalarEstimate
    used_sample_collective_persistence: ScalarEstimate | None
    used_sample_critical_diagonal_persistence: ScalarEstimate | None
    used_sample_off_diagonal_persistence: ScalarEstimate | None
    used_sample_normalized_off_diagonal_persistence: float | None
    unconditional_normalized_off_diagonal_persistence: float | None
    maximum_decomposition_error: float | None
    maximum_diagonal_error: float | None
    maximum_symmetry_error: float | None
    maximum_global_flip_log_weight_error: float | None
    exact_orientation_enumeration_on_used_environments: bool
    maximum_block_count_selection_bias: bool
    selection_bias_warning: str
    weak_recovery_claimed: bool
    interpretation: str


def _validate_inputs(
    side_length: int,
    p: float,
    maximum_block_count: int,
) -> None:
    if side_length < 4:
        raise ValueError("side_length must be at least 4")
    if not P_SW < p < 1.0:
        raise ValueError(f"p must satisfy {P_SW} < p < 1")
    if maximum_block_count <= 0:
        raise ValueError("maximum_block_count must be positive")


def _component_membership(
    vertex_count: int,
    components: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    membership = [-1] * vertex_count
    for component_index, component in enumerate(components):
        for vertex in component:
            if membership[vertex] != -1:
                raise AssertionError("components overlap")
            membership[vertex] = component_index
    if any(index < 0 for index in membership):
        raise AssertionError("components do not cover every vertex")
    return tuple(membership)


def critical_cut_partition(
    environment: HierarchicalSweepEnvironment,
    ranked_edges: Sequence[RankedEdge],
    critical_rank: float = Q_CRITICAL,
) -> CriticalCutPartition:
    """Build the critical partition and group it by final roots."""

    if not 0.0 <= critical_rank < 2.0 * environment.p - 1.0:
        raise ValueError("critical_rank must lie below the final rank")
    critical_forest = CriticalKruskalForest(
        environment.side_length,
        ranked_edges,
        critical_rank=critical_rank,
    )
    blocks = tuple(tuple(component) for component in critical_forest.components)
    block_masks = tuple(
        sum(1 << vertex for vertex in block)
        for block in blocks
    )
    block_sizes = tuple(len(block) for block in blocks)
    final_components = tuple(
        tuple(component) for component in environment.forest.components
    )
    final_membership = _component_membership(
        environment.vertex_count,
        final_components,
    )
    grouped: list[list[int]] = [[] for _ in final_components]
    for block_index, block in enumerate(blocks):
        component_indices = {final_membership[vertex] for vertex in block}
        if len(component_indices) != 1:
            raise AssertionError("a critical block crosses final roots")
        grouped[next(iter(component_indices))].append(block_index)
    if sum(len(indices) for indices in grouped) != len(blocks):
        raise AssertionError("some critical blocks were not assigned")
    for component_index, block_indices in enumerate(grouped):
        covered = {
            vertex
            for block_index in block_indices
            for vertex in blocks[block_index]
        }
        if covered != set(final_components[component_index]):
            raise AssertionError("critical blocks do not partition a final root")
    return CriticalCutPartition(
        blocks=blocks,
        block_masks=block_masks,
        block_sizes=block_sizes,
        final_components=final_components,
        final_component_block_indices=tuple(
            tuple(indices) for indices in grouped
        ),
    )


def _normalized_weights(log_weights: Sequence[float]) -> tuple[float, ...]:
    maximum = max(log_weights)
    if maximum == float("-inf"):
        raise AssertionError("the collective orientation orbit has zero mass")
    raw = tuple(
        0.0 if value == float("-inf") else exp(value - maximum)
        for value in log_weights
    )
    normalizer = fsum(raw)
    if normalizer <= 0.0:
        raise AssertionError("the collective orientation normalizer vanished")
    return tuple(value / normalizer for value in raw)


def _global_flip_log_weight_error(log_weights: Sequence[float]) -> float:
    full_mask = len(log_weights) - 1
    errors = []
    for configuration, first in enumerate(log_weights):
        second = log_weights[configuration ^ full_mask]
        if first == float("-inf") and second == float("-inf"):
            errors.append(0.0)
        elif first == float("-inf") or second == float("-inf"):
            errors.append(float("inf"))
        else:
            errors.append(abs(first - second))
    return max(errors, default=0.0)


def analyze_final_root(
    environment: HierarchicalSweepEnvironment,
    partition: CriticalCutPartition,
    final_component_index: int,
) -> RootCollectiveGibbsDiagnostic:
    """Enumerate the exact joint orientation heat bath in one final root."""

    if not 0 <= final_component_index < len(partition.final_components):
        raise ValueError("final_component_index is out of range")
    block_indices = partition.final_component_block_indices[
        final_component_index
    ]
    if not block_indices:
        raise AssertionError("a final root contains no critical block")
    block_masks = tuple(partition.block_masks[index] for index in block_indices)
    block_sizes = tuple(partition.block_sizes[index] for index in block_indices)
    block_count = len(block_indices)
    orientation_state_count = 1 << block_count

    states = []
    signs = []
    for configuration in range(orientation_state_count):
        state = 0
        orientation = []
        for local_index, block_mask in enumerate(block_masks):
            flipped = bool((configuration >> local_index) & 1)
            if flipped:
                state ^= block_mask
            orientation.append(-1 if flipped else 1)
        states.append(state)
        signs.append(tuple(orientation))
    log_weights = tuple(environment.full_log_weight(state) for state in states)
    probabilities = _normalized_weights(log_weights)

    correlations = []
    for first in range(block_count):
        correlations.append(
            tuple(
                fsum(
                    probability
                    * orientation[first]
                    * orientation[second]
                    for probability, orientation in zip(
                        probabilities,
                        signs,
                        strict=True,
                    )
                )
                for second in range(block_count)
            )
        )
    matrix = tuple(correlations)
    diagonal = fsum(size * size for size in block_sizes)
    persistence = fsum(
        block_sizes[first]
        * block_sizes[second]
        * matrix[first][second]
        * matrix[first][second]
        for first in range(block_count)
        for second in range(block_count)
    )
    maximum_diagonal_error = max(
        abs(matrix[index][index] - 1.0)
        for index in range(block_count)
    )
    maximum_symmetry_error = max(
        abs(matrix[first][second] - matrix[second][first])
        for first in range(block_count)
        for second in range(block_count)
    )
    component = partition.final_components[final_component_index]
    if sum(block_sizes) != len(component):
        raise AssertionError("block sizes disagree with final component size")
    if persistence < diagonal - 1e-10:
        raise AssertionError("off-diagonal squared persistence is negative")
    if persistence > len(component) ** 2 + 1e-10:
        raise AssertionError("root persistence exceeds its trivial bound")
    return RootCollectiveGibbsDiagnostic(
        final_component_index=final_component_index,
        representative_vertex=min(component),
        final_component_size=len(component),
        block_indices=block_indices,
        block_sizes=block_sizes,
        orientation_state_count=orientation_state_count,
        positive_orientation_state_count=sum(
            probability > 0.0 for probability in probabilities
        ),
        correlation_matrix=matrix,
        unnormalized_persistence=persistence,
        diagonal_unnormalized_persistence=diagonal,
        maximum_diagonal_error=maximum_diagonal_error,
        maximum_symmetry_error=maximum_symmetry_error,
        maximum_global_flip_log_weight_error=(
            _global_flip_log_weight_error(log_weights)
        ),
    )


def diagnose_environment(
    side_length: int,
    ranked_edges: Sequence[RankedEdge],
    p: float,
    maximum_block_count: int,
    repetition: int = 0,
) -> EnvironmentCollectiveGibbsDiagnostic:
    """Compute ``G_L^c`` exactly, unless the declared cap excludes it."""

    _validate_inputs(side_length, p, maximum_block_count)
    environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
    partition = critical_cut_partition(environment, ranked_edges)
    total_block_count = len(partition.blocks)
    root_block_counts = tuple(
        len(indices) for indices in partition.final_component_block_indices
    )
    maximum_root_block_count = max(root_block_counts)
    vertex_count = environment.vertex_count
    diagonal = fsum(size * size for size in partition.block_sizes) / (
        vertex_count * vertex_count
    )

    if maximum_root_block_count > maximum_block_count:
        skip_reason = ROOT_CAP_SKIP
    elif total_block_count > maximum_block_count:
        skip_reason = TOTAL_CAP_SKIP
    else:
        skip_reason = None
    if skip_reason is not None:
        return EnvironmentCollectiveGibbsDiagnostic(
            repetition=repetition,
            total_critical_block_count=total_block_count,
            final_root_count=len(partition.final_components),
            maximum_root_block_count=maximum_root_block_count,
            used=False,
            skip_reason=skip_reason,
            critical_diagonal_persistence=diagonal,
            collective_persistence=None,
            off_diagonal_persistence=None,
            maximum_decomposition_error=None,
            maximum_diagonal_error=None,
            maximum_symmetry_error=None,
            maximum_global_flip_log_weight_error=None,
            roots=(),
        )

    roots = tuple(
        analyze_final_root(environment, partition, component_index)
        for component_index in range(len(partition.final_components))
    )
    collective = fsum(
        root.unnormalized_persistence for root in roots
    ) / (vertex_count * vertex_count)
    diagonal_from_roots = fsum(
        root.diagonal_unnormalized_persistence for root in roots
    ) / (vertex_count * vertex_count)
    decomposition_error = abs(diagonal - diagonal_from_roots)
    off_diagonal = collective - diagonal
    if off_diagonal < 0.0 and off_diagonal > -1e-12:
        off_diagonal = 0.0
    if off_diagonal < 0.0:
        raise AssertionError("collective persistence is below its diagonal")
    return EnvironmentCollectiveGibbsDiagnostic(
        repetition=repetition,
        total_critical_block_count=total_block_count,
        final_root_count=len(partition.final_components),
        maximum_root_block_count=maximum_root_block_count,
        used=True,
        skip_reason=None,
        critical_diagonal_persistence=diagonal,
        collective_persistence=collective,
        off_diagonal_persistence=off_diagonal,
        maximum_decomposition_error=decomposition_error,
        maximum_diagonal_error=max(root.maximum_diagonal_error for root in roots),
        maximum_symmetry_error=max(root.maximum_symmetry_error for root in roots),
        maximum_global_flip_log_weight_error=max(
            root.maximum_global_flip_log_weight_error for root in roots
        ),
        roots=roots,
    )


def scalar_estimate(values: Sequence[float]) -> ScalarEstimate:
    """Return a mean with an environment-clustered standard error."""

    if not values:
        raise ValueError("at least one value is required")
    mean = fmean(values)
    if len(values) < 2:
        standard_error = None
    else:
        standard_error = sqrt(
            fsum((value - mean) ** 2 for value in values)
            / (len(values) * (len(values) - 1))
        )
    return ScalarEstimate(mean=mean, standard_error=standard_error)


def summarize(
    diagnostics: Sequence[EnvironmentCollectiveGibbsDiagnostic],
    side_length: int,
    p: float,
    maximum_block_count: int,
    seed: int,
) -> CriticalCutCollectiveGibbsSummary:
    """Aggregate environments without hiding cap-induced conditioning."""

    if not diagnostics:
        raise ValueError("at least one diagnostic is required")
    used = tuple(item for item in diagnostics if item.used)
    skipped = tuple(item for item in diagnostics if not item.used)
    skipped_root = sum(item.skip_reason == ROOT_CAP_SKIP for item in skipped)
    skipped_total = sum(item.skip_reason == TOTAL_CAP_SKIP for item in skipped)
    if skipped_root + skipped_total != len(skipped):
        raise AssertionError("an excluded environment has an unknown reason")

    all_diagonal = scalar_estimate(
        tuple(item.critical_diagonal_persistence for item in diagnostics)
    )
    if used:
        collective = scalar_estimate(
            tuple(
                float(item.collective_persistence)
                for item in used
                if item.collective_persistence is not None
            )
        )
        used_diagonal = scalar_estimate(
            tuple(item.critical_diagonal_persistence for item in used)
        )
        off_diagonal = scalar_estimate(
            tuple(
                float(item.off_diagonal_persistence)
                for item in used
                if item.off_diagonal_persistence is not None
            )
        )
        denominator = 1.0 - used_diagonal.mean
        normalized_used = (
            None
            if denominator <= 1e-15
            else (collective.mean - used_diagonal.mean) / denominator
        )
        maximum_decomposition_error = max(
            float(item.maximum_decomposition_error)
            for item in used
            if item.maximum_decomposition_error is not None
        )
        maximum_diagonal_error = max(
            float(item.maximum_diagonal_error)
            for item in used
            if item.maximum_diagonal_error is not None
        )
        maximum_symmetry_error = max(
            float(item.maximum_symmetry_error)
            for item in used
            if item.maximum_symmetry_error is not None
        )
        maximum_global_flip_error = max(
            float(item.maximum_global_flip_log_weight_error)
            for item in used
            if item.maximum_global_flip_log_weight_error is not None
        )
    else:
        collective = None
        used_diagonal = None
        off_diagonal = None
        normalized_used = None
        maximum_decomposition_error = None
        maximum_diagonal_error = None
        maximum_symmetry_error = None
        maximum_global_flip_error = None

    selection_bias = bool(skipped)
    if selection_bias:
        warning = (
            "WARNING: maximum_block_count excluded "
            f"{len(skipped)}/{len(diagnostics)} environments. Gibbs means and "
            "the used-sample normalized persistence are conditional on passing "
            "the cap and are selection-biased; they do not estimate the "
            "unconditional E[G_L^c]."
        )
        unconditional_normalized = None
    else:
        warning = (
            "No maximum_block_count exclusion occurred; the sampled-environment "
            "means are unconditional with respect to the cap."
        )
        unconditional_normalized = normalized_used

    return CriticalCutCollectiveGibbsSummary(
        side_length=side_length,
        vertex_count=side_length * side_length,
        p=p,
        critical_rank=Q_CRITICAL,
        final_rank=2.0 * p - 1.0,
        repetitions=len(diagnostics),
        seed=seed,
        maximum_block_count=maximum_block_count,
        used_environment_count=len(used),
        skipped_environment_count=len(skipped),
        skipped_for_root_block_cap_count=skipped_root,
        skipped_for_total_block_cap_count=skipped_total,
        all_sample_critical_diagonal_persistence=all_diagonal,
        used_sample_collective_persistence=collective,
        used_sample_critical_diagonal_persistence=used_diagonal,
        used_sample_off_diagonal_persistence=off_diagonal,
        used_sample_normalized_off_diagonal_persistence=normalized_used,
        unconditional_normalized_off_diagonal_persistence=(
            unconditional_normalized
        ),
        maximum_decomposition_error=maximum_decomposition_error,
        maximum_diagonal_error=maximum_diagonal_error,
        maximum_symmetry_error=maximum_symmetry_error,
        maximum_global_flip_log_weight_error=maximum_global_flip_error,
        exact_orientation_enumeration_on_used_environments=True,
        maximum_block_count_selection_bias=selection_bias,
        selection_bias_warning=warning,
        weak_recovery_claimed=False,
        interpretation=(
            "finite exact collective critical-block heat-bath diagnostic; "
            "no thermodynamic extrapolation or weak-recovery threshold is claimed"
        ),
    )


def run_diagnostic(
    side_length: int = 4,
    repetitions: int = 8,
    p: float = BASELINE_P,
    maximum_block_count: int = 16,
    seed: int = 20260726,
) -> tuple[
    CriticalCutCollectiveGibbsSummary,
    tuple[EnvironmentCollectiveGibbsDiagnostic, ...],
]:
    """Sample rank environments and evaluate the exact capped diagnostic."""

    _validate_inputs(side_length, p, maximum_block_count)
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    rng = random.Random(seed)
    diagnostics = tuple(
        diagnose_environment(
            side_length=side_length,
            ranked_edges=sample_ranked_edges(side_length, rng),
            p=p,
            maximum_block_count=maximum_block_count,
            repetition=repetition,
        )
        for repetition in range(repetitions)
    )
    return (
        summarize(
            diagnostics,
            side_length=side_length,
            p=p,
            maximum_block_count=maximum_block_count,
            seed=seed,
        ),
        diagnostics,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=4)
    parser.add_argument("--repetitions", type=int, default=8)
    parser.add_argument("--p", type=float, default=BASELINE_P)
    parser.add_argument("--maximum-block-count", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument("--include-environments", action="store_true")
    arguments = parser.parse_args()
    summary, diagnostics = run_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        p=arguments.p,
        maximum_block_count=arguments.maximum_block_count,
        seed=arguments.seed,
    )
    payload: dict[str, object] = {"summary": asdict(summary)}
    if arguments.include_environments:
        payload["environments"] = [asdict(item) for item in diagnostics]
    print(json.dumps(payload, indent=2, sort_keys=True))
    if summary.maximum_block_count_selection_bias:
        print(summary.selection_bias_warning, file=sys.stderr)


if __name__ == "__main__":
    main()
