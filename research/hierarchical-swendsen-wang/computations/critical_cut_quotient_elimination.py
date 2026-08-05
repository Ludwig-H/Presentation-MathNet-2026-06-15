"""Exact critical-block correlations by binary factor elimination.

The exhaustive diagnostic in ``critical_cut_collective_gibbs_diagnostic``
enumerates all ``2^k`` orientations of the ``k`` critical blocks in a final
Kruskal root.  This module represents the very same conditional law as a
factor graph on those block orientations and eliminates variables with a
deterministic min-fill order.

For one hierarchy factor ``u``, let ``E_u`` be its edge bucket and let
``N_u(z)`` be the number of satisfied edges after orienting critical blocks
by ``z``.  Up to an orientation-independent constant, its factor is

    psi_u(z) = N_u(z) exp((1-beta_u) J N_u(z)).

Only critical blocks incident to an inter-block edge of ``E_u`` enter its
scope.  Factors entirely inside one critical block cancel from normalized
correlations.  The construction is exact conditional on the fixed all-plus
internal representative used by the exhaustive diagnostic; it does not
integrate the internal spin degrees of freedom of the full conditional
Gibbs law.

Pair correlations are obtained from four constrained partition functions.
Their cost is exponential in the induced factor-graph width, not directly
in the total number of critical blocks.  Every reported width is an observed
finite-instance complexity diagnostic, not an asymptotic treewidth bound.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import exp, fsum, log
from statistics import fmean
from time import perf_counter
from typing import Mapping, Sequence

from critical_cut_collective_gibbs_diagnostic import (
    BASELINE_P,
    P_SW,
    CriticalCutPartition,
    critical_cut_partition,
)
from critical_pair_path_geometry import (
    clock_time_from_rank,
    sample_ranked_edges,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


@dataclass(frozen=True)
class BinaryFactor:
    """A nonnegative table indexed by binary assignments to ``variables``."""

    variables: tuple[int, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        if tuple(sorted(set(self.variables))) != self.variables:
            raise ValueError("factor variables must be sorted and distinct")
        if len(self.values) != 1 << len(self.variables):
            raise ValueError("factor table size disagrees with its scope")
        if any(value < 0.0 for value in self.values):
            raise ValueError("factor values must be nonnegative")


@dataclass(frozen=True)
class EliminationAudit:
    """One constrained partition function and its elimination complexity."""

    fixed_assignments: tuple[tuple[int, int], ...]
    log_partition: float
    elimination_order: tuple[int, ...]
    induced_width: int
    maximum_intermediate_scope_size: int


@dataclass(frozen=True)
class RootQuotientEliminationDiagnostic:
    """Exact correlation matrix and factor-width audits for one final root."""

    final_component_index: int
    block_indices: tuple[int, ...]
    block_sizes: tuple[int, ...]
    block_count: int
    orientation_state_count: int
    nonconstant_factor_count: int
    maximum_initial_factor_scope_size: int
    scaled_unconstrained_log_partition: float
    unconstrained_induced_width: int
    maximum_pair_induced_width: int
    maximum_pair_intermediate_scope_size: int
    constrained_partition_evaluation_count: int
    correlation_matrix: tuple[tuple[float, ...], ...]
    maximum_diagonal_error: float
    maximum_symmetry_error: float
    maximum_global_flip_log_partition_error: float
    exact_factor_elimination: bool
    conditional_on_fixed_internal_representative: bool
    factor_log_normalization_constants_omitted: bool
    scaled_log_partition_definition: str


@dataclass(frozen=True)
class EnvironmentQuotientEliminationDiagnostic:
    """Exact quotient result and complexity audit for one rank environment."""

    repetition: int
    total_critical_block_count: int
    final_root_count: int
    maximum_root_block_count: int
    maximum_initial_factor_scope_size: int
    maximum_unconstrained_induced_width: int
    maximum_pair_induced_width: int
    maximum_pair_intermediate_scope_size: int
    maximum_orientation_state_count_avoided: int
    critical_diagonal_persistence: float
    collective_persistence: float
    off_diagonal_persistence: float
    roots: tuple[RootQuotientEliminationDiagnostic, ...]


@dataclass(frozen=True)
class QuotientEliminationSummary:
    """Uncapped finite-sample summary of exact quotient eliminations."""

    side_length: int
    vertex_count: int
    p: float
    repetitions: int
    seed: int
    mean_total_critical_block_count: float
    maximum_total_critical_block_count: int
    mean_maximum_root_block_count: float
    maximum_root_block_count: int
    mean_maximum_unconstrained_induced_width: float
    maximum_unconstrained_induced_width: int
    mean_maximum_pair_induced_width: float
    maximum_pair_induced_width: int
    maximum_initial_factor_scope_size: int
    maximum_pair_intermediate_scope_size: int
    maximum_orientation_state_count_avoided: int
    mean_critical_diagonal_persistence: float
    mean_collective_persistence: float
    mean_off_diagonal_persistence: float
    cap_exclusion_count: int
    selection_bias_from_complexity_cap: bool
    every_environment_eliminated_exactly: bool
    conditional_on_fixed_internal_representative: bool
    factor_log_normalization_constants_omitted: bool
    weak_recovery_claimed: bool
    interpretation: str


def _normalized_factor(
    variables: Sequence[int],
    log_values: Sequence[float],
) -> BinaryFactor:
    """Exponentiate a log table after removing a common scale."""

    if not log_values:
        raise ValueError("a factor table cannot be empty")
    maximum = max(log_values)
    if maximum == float("-inf"):
        raise AssertionError("a hierarchy factor vanishes on every assignment")
    values = tuple(
        0.0 if value == float("-inf") else exp(value - maximum)
        for value in log_values
    )
    return BinaryFactor(tuple(variables), values)


def critical_root_factors(
    environment: HierarchicalSweepEnvironment,
    partition: CriticalCutPartition,
    final_component_index: int,
) -> tuple[BinaryFactor, ...]:
    """Tabulate the exact nonconstant quotient factors in one final root."""

    if not 0 <= final_component_index < len(partition.final_components):
        raise ValueError("final_component_index is out of range")
    block_indices = partition.final_component_block_indices[
        final_component_index
    ]
    local_index = {
        block_index: index
        for index, block_index in enumerate(block_indices)
    }
    vertex_block = {}
    for block_index, block in enumerate(partition.blocks):
        for vertex in block:
            if vertex in vertex_block:
                raise AssertionError("critical blocks overlap")
            vertex_block[vertex] = block_index
    if len(vertex_block) != environment.vertex_count:
        raise AssertionError("critical blocks do not cover the environment")

    factors = []
    for factor_node in environment.forest.internal_nodes:
        edges = environment.bucket_edges[factor_node]
        scope_global = {
            block_index
            for first, second in edges
            for block_index in (
                vertex_block[first],
                vertex_block[second],
            )
            if vertex_block[first] != vertex_block[second]
        }
        if not scope_global:
            continue
        if not scope_global.issubset(local_index):
            # This factor belongs to a different final root.  A Kruskal
            # bucket cannot meet two distinct final roots.
            if scope_global.isdisjoint(local_index):
                continue
            raise AssertionError("a hierarchy factor crosses final roots")
        scope = tuple(sorted(local_index[index] for index in scope_global))
        scope_position = {
            block_indices[local]: position
            for position, local in enumerate(scope)
        }
        log_values = []
        for assignment in range(1 << len(scope)):
            satisfied_count = 0
            for first, second in edges:
                first_block = vertex_block[first]
                second_block = vertex_block[second]
                if first_block == second_block:
                    relative_parity = 0
                else:
                    first_bit = (
                        assignment >> scope_position[first_block]
                    ) & 1
                    second_bit = (
                        assignment >> scope_position[second_block]
                    ) & 1
                    relative_parity = first_bit ^ second_bit
                truth_satisfied = (
                    environment.edge_rank[tuple(sorted((first, second)))]
                    <= environment.p
                )
                satisfied_count += (
                    truth_satisfied != bool(relative_parity)
                )
            if satisfied_count == 0:
                log_values.append(float("-inf"))
                continue
            time = clock_time_from_rank(
                environment.forest.merge_rank[factor_node],
                environment.p,
            )
            log_values.append(
                log(satisfied_count)
                + (1.0 - time)
                * environment.coupling
                * satisfied_count
            )
        factors.append(_normalized_factor(scope, log_values))
    return tuple(factors)


def _restrict_factor(
    factor: BinaryFactor,
    fixed_assignments: Mapping[int, int],
) -> BinaryFactor:
    """Restrict a table and remove the fixed variables from its scope."""

    remaining = tuple(
        variable
        for variable in factor.variables
        if variable not in fixed_assignments
    )
    positions = {variable: index for index, variable in enumerate(remaining)}
    values = []
    for assignment in range(1 << len(remaining)):
        original_assignment = 0
        for original_position, variable in enumerate(factor.variables):
            if variable in fixed_assignments:
                bit = fixed_assignments[variable]
            else:
                bit = (assignment >> positions[variable]) & 1
            original_assignment |= bit << original_position
        values.append(factor.values[original_assignment])
    return BinaryFactor(remaining, tuple(values))


def _min_fill_order(
    factors: Sequence[BinaryFactor],
    variables: Sequence[int],
) -> tuple[tuple[int, ...], int]:
    """Return a deterministic min-fill order and its induced width."""

    graph = {variable: set() for variable in variables}
    for factor in factors:
        scope = tuple(
            variable for variable in factor.variables if variable in graph
        )
        for variable in scope:
            graph[variable].update(set(scope) - {variable})
    order = []
    induced_width = 0
    while graph:
        def score(variable: int) -> tuple[int, int, int]:
            neighbours = tuple(sorted(graph[variable]))
            missing_edges = sum(
                second not in graph[first]
                for index, first in enumerate(neighbours)
                for second in neighbours[index + 1 :]
            )
            return missing_edges, len(neighbours), variable

        variable = min(graph, key=score)
        neighbours = set(graph[variable])
        induced_width = max(induced_width, len(neighbours))
        for neighbour in neighbours:
            graph[neighbour].update(neighbours - {neighbour})
            graph[neighbour].discard(variable)
        del graph[variable]
        order.append(variable)
    return tuple(order), induced_width


def _factor_value_from_union_assignment(
    factor: BinaryFactor,
    union_positions: Mapping[int, int],
    union_assignment: int,
) -> float:
    local_assignment = 0
    for local_position, variable in enumerate(factor.variables):
        bit = (union_assignment >> union_positions[variable]) & 1
        local_assignment |= bit << local_position
    return factor.values[local_assignment]


def eliminate_log_partition(
    factors: Sequence[BinaryFactor],
    variable_count: int,
    fixed_assignments: Mapping[int, int] | None = None,
) -> EliminationAudit:
    """Compute one constrained partition function by exact variable elimination."""

    if variable_count < 0:
        raise ValueError("variable_count must be nonnegative")
    fixed = dict(fixed_assignments or {})
    if any(
        not 0 <= variable < variable_count
        for variable in fixed
    ):
        raise ValueError("a fixed variable is out of range")
    if any(bit not in (0, 1) for bit in fixed.values()):
        raise ValueError("fixed assignments must be binary")
    if any(
        variable < 0 or variable >= variable_count
        for factor in factors
        for variable in factor.variables
    ):
        raise ValueError("a factor variable is out of range")

    active = [_restrict_factor(factor, fixed) for factor in factors]
    free_variables = tuple(
        variable
        for variable in range(variable_count)
        if variable not in fixed
    )
    order, induced_width = _min_fill_order(active, free_variables)
    log_scale = 0.0
    maximum_scope = 0

    scalars = [factor for factor in active if not factor.variables]
    active = [factor for factor in active if factor.variables]
    for factor in scalars:
        value = factor.values[0]
        if value == 0.0:
            return EliminationAudit(
                fixed_assignments=tuple(sorted(fixed.items())),
                log_partition=float("-inf"),
                elimination_order=order,
                induced_width=induced_width,
                maximum_intermediate_scope_size=maximum_scope,
            )
        log_scale += log(value)

    for variable in order:
        bucket = [
            factor for factor in active if variable in factor.variables
        ]
        if not bucket:
            log_scale += log(2.0)
            continue
        active = [
            factor for factor in active if variable not in factor.variables
        ]
        union = tuple(sorted({
            item
            for factor in bucket
            for item in factor.variables
        }))
        maximum_scope = max(maximum_scope, len(union))
        remaining = tuple(item for item in union if item != variable)
        union_positions = {
            item: index for index, item in enumerate(union)
        }
        variable_position = union_positions[variable]
        remaining_positions = {
            item: index for index, item in enumerate(remaining)
        }
        new_values = []
        for assignment in range(1 << len(remaining)):
            total = 0.0
            for variable_bit in (0, 1):
                union_assignment = variable_bit << variable_position
                for item in remaining:
                    bit = (
                        assignment >> remaining_positions[item]
                    ) & 1
                    union_assignment |= bit << union_positions[item]
                product = 1.0
                for factor in bucket:
                    product *= _factor_value_from_union_assignment(
                        factor,
                        union_positions,
                        union_assignment,
                    )
                total += product
            new_values.append(total)
        maximum = max(new_values)
        if maximum == 0.0:
            return EliminationAudit(
                fixed_assignments=tuple(sorted(fixed.items())),
                log_partition=float("-inf"),
                elimination_order=order,
                induced_width=induced_width,
                maximum_intermediate_scope_size=maximum_scope,
            )
        log_scale += log(maximum)
        normalized = tuple(value / maximum for value in new_values)
        if remaining:
            active.append(BinaryFactor(remaining, normalized))
        else:
            log_scale += log(normalized[0])

    if active:
        raise AssertionError("variable elimination left nonconstant factors")
    return EliminationAudit(
        fixed_assignments=tuple(sorted(fixed.items())),
        log_partition=log_scale,
        elimination_order=order,
        induced_width=induced_width,
        maximum_intermediate_scope_size=maximum_scope,
    )


def _correlation_from_log_partitions(
    log_partitions: Sequence[float],
) -> float:
    """Return ``E[(-1)^(x+y)]`` from the four pair partition functions."""

    if len(log_partitions) != 4:
        raise ValueError("four pair partition functions are required")
    maximum = max(log_partitions)
    if maximum == float("-inf"):
        raise AssertionError("all pair assignments have zero mass")
    weights = tuple(
        0.0 if value == float("-inf") else exp(value - maximum)
        for value in log_partitions
    )
    denominator = fsum(weights)
    return (
        weights[0] - weights[1] - weights[2] + weights[3]
    ) / denominator


def analyze_root_by_elimination(
    environment: HierarchicalSweepEnvironment,
    partition: CriticalCutPartition,
    final_component_index: int,
) -> RootQuotientEliminationDiagnostic:
    """Compute the exact quotient correlation matrix without ``2^k`` enumeration."""

    block_indices = partition.final_component_block_indices[
        final_component_index
    ]
    block_sizes = tuple(
        partition.block_sizes[index] for index in block_indices
    )
    block_count = len(block_indices)
    if block_count <= 0:
        raise AssertionError("a final root contains no critical block")
    factors = critical_root_factors(
        environment,
        partition,
        final_component_index,
    )
    unconstrained = eliminate_log_partition(factors, block_count)
    matrix = [[0.0] * block_count for _ in range(block_count)]
    for index in range(block_count):
        matrix[index][index] = 1.0
    audits = []
    flip_errors = []
    for first in range(block_count):
        for second in range(first + 1, block_count):
            pair_audits = tuple(
                eliminate_log_partition(
                    factors,
                    block_count,
                    {
                        first: first_bit,
                        second: second_bit,
                    },
                )
                for first_bit, second_bit in (
                    (0, 0),
                    (0, 1),
                    (1, 0),
                    (1, 1),
                )
            )
            logs = tuple(audit.log_partition for audit in pair_audits)
            value = _correlation_from_log_partitions(logs)
            matrix[first][second] = value
            matrix[second][first] = value
            audits.extend(pair_audits)
            for left, right in ((logs[0], logs[3]), (logs[1], logs[2])):
                if left == float("-inf") and right == float("-inf"):
                    flip_errors.append(0.0)
                elif left == float("-inf") or right == float("-inf"):
                    flip_errors.append(float("inf"))
                else:
                    flip_errors.append(abs(left - right))
    correlations = tuple(tuple(row) for row in matrix)
    return RootQuotientEliminationDiagnostic(
        final_component_index=final_component_index,
        block_indices=block_indices,
        block_sizes=block_sizes,
        block_count=block_count,
        orientation_state_count=1 << block_count,
        nonconstant_factor_count=len(factors),
        maximum_initial_factor_scope_size=max(
            (len(factor.variables) for factor in factors),
            default=0,
        ),
        scaled_unconstrained_log_partition=unconstrained.log_partition,
        unconstrained_induced_width=unconstrained.induced_width,
        maximum_pair_induced_width=max(
            (audit.induced_width for audit in audits),
            default=0,
        ),
        maximum_pair_intermediate_scope_size=max(
            (audit.maximum_intermediate_scope_size for audit in audits),
            default=0,
        ),
        constrained_partition_evaluation_count=len(audits),
        correlation_matrix=correlations,
        maximum_diagonal_error=max(
            abs(correlations[index][index] - 1.0)
            for index in range(block_count)
        ),
        maximum_symmetry_error=max(
            abs(correlations[first][second] - correlations[second][first])
            for first in range(block_count)
            for second in range(block_count)
        ),
        maximum_global_flip_log_partition_error=max(
            flip_errors,
            default=0.0,
        ),
        exact_factor_elimination=True,
        conditional_on_fixed_internal_representative=True,
        factor_log_normalization_constants_omitted=True,
        scaled_log_partition_definition=(
            "log partition of the nonconstant quotient factors after each "
            "factor table is divided by its maximum; orientation-constant "
            "hierarchy factors and the removed per-factor maxima are omitted"
        ),
    )


def diagnose_environment_by_elimination(
    side_length: int,
    p: float,
    ranked_edges: Sequence[tuple[float, int, int]],
    repetition: int = 0,
) -> EnvironmentQuotientEliminationDiagnostic:
    """Eliminate every final root exactly, without a block-count cap."""

    environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
    partition = critical_cut_partition(environment, ranked_edges)
    roots = tuple(
        analyze_root_by_elimination(
            environment,
            partition,
            component_index,
        )
        for component_index in range(len(partition.final_components))
    )
    vertex_count = environment.vertex_count
    normalization = vertex_count * vertex_count
    diagonal = fsum(
        size * size for size in partition.block_sizes
    ) / normalization
    collective = fsum(
        root.block_sizes[first]
        * root.block_sizes[second]
        * root.correlation_matrix[first][second] ** 2
        for root in roots
        for first in range(root.block_count)
        for second in range(root.block_count)
    ) / normalization
    off_diagonal = collective - diagonal
    if off_diagonal < 0.0 and off_diagonal > -1e-12:
        off_diagonal = 0.0
    if off_diagonal < 0.0:
        raise AssertionError("collective persistence is below its diagonal")
    return EnvironmentQuotientEliminationDiagnostic(
        repetition=repetition,
        total_critical_block_count=len(partition.blocks),
        final_root_count=len(partition.final_components),
        maximum_root_block_count=max(root.block_count for root in roots),
        maximum_initial_factor_scope_size=max(
            root.maximum_initial_factor_scope_size for root in roots
        ),
        maximum_unconstrained_induced_width=max(
            root.unconstrained_induced_width for root in roots
        ),
        maximum_pair_induced_width=max(
            root.maximum_pair_induced_width for root in roots
        ),
        maximum_pair_intermediate_scope_size=max(
            root.maximum_pair_intermediate_scope_size for root in roots
        ),
        maximum_orientation_state_count_avoided=max(
            root.orientation_state_count for root in roots
        ),
        critical_diagonal_persistence=diagonal,
        collective_persistence=collective,
        off_diagonal_persistence=off_diagonal,
        roots=roots,
    )


def summarize_eliminations(
    diagnostics: Sequence[EnvironmentQuotientEliminationDiagnostic],
    side_length: int,
    p: float,
    seed: int,
) -> QuotientEliminationSummary:
    """Aggregate an uncapped collection of exact environment diagnostics."""

    if not diagnostics:
        raise ValueError("at least one environment diagnostic is required")
    return QuotientEliminationSummary(
        side_length=side_length,
        vertex_count=side_length * side_length,
        p=p,
        repetitions=len(diagnostics),
        seed=seed,
        mean_total_critical_block_count=fmean(
            item.total_critical_block_count for item in diagnostics
        ),
        maximum_total_critical_block_count=max(
            item.total_critical_block_count for item in diagnostics
        ),
        mean_maximum_root_block_count=fmean(
            item.maximum_root_block_count for item in diagnostics
        ),
        maximum_root_block_count=max(
            item.maximum_root_block_count for item in diagnostics
        ),
        mean_maximum_unconstrained_induced_width=fmean(
            item.maximum_unconstrained_induced_width
            for item in diagnostics
        ),
        maximum_unconstrained_induced_width=max(
            item.maximum_unconstrained_induced_width
            for item in diagnostics
        ),
        mean_maximum_pair_induced_width=fmean(
            item.maximum_pair_induced_width for item in diagnostics
        ),
        maximum_pair_induced_width=max(
            item.maximum_pair_induced_width for item in diagnostics
        ),
        maximum_initial_factor_scope_size=max(
            item.maximum_initial_factor_scope_size for item in diagnostics
        ),
        maximum_pair_intermediate_scope_size=max(
            item.maximum_pair_intermediate_scope_size
            for item in diagnostics
        ),
        maximum_orientation_state_count_avoided=max(
            item.maximum_orientation_state_count_avoided
            for item in diagnostics
        ),
        mean_critical_diagonal_persistence=fmean(
            item.critical_diagonal_persistence for item in diagnostics
        ),
        mean_collective_persistence=fmean(
            item.collective_persistence for item in diagnostics
        ),
        mean_off_diagonal_persistence=fmean(
            item.off_diagonal_persistence for item in diagnostics
        ),
        cap_exclusion_count=0,
        selection_bias_from_complexity_cap=False,
        every_environment_eliminated_exactly=True,
        conditional_on_fixed_internal_representative=True,
        factor_log_normalization_constants_omitted=True,
        weak_recovery_claimed=False,
        interpretation=(
            "uncapped exact factor elimination for critical-block orientations "
            "conditional on the fixed all-plus internal representative; "
            "observed min-fill widths are finite-instance complexity data, "
            "with no full conditional-Gibbs or asymptotic recovery claim"
        ),
    )


def run_elimination_diagnostic(
    side_length: int = 6,
    repetitions: int = 8,
    p: float = BASELINE_P,
    seed: int = 20260726,
) -> tuple[
    QuotientEliminationSummary,
    tuple[EnvironmentQuotientEliminationDiagnostic, ...],
]:
    """Sample environments and run the uncapped exact quotient diagnostic."""

    if side_length < 4:
        raise ValueError("side_length must be at least 4")
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if not P_SW < p < 1.0:
        raise ValueError(f"p must satisfy {P_SW} < p < 1")
    rng = random.Random(seed)
    diagnostics = tuple(
        diagnose_environment_by_elimination(
            side_length=side_length,
            p=p,
            ranked_edges=sample_ranked_edges(side_length, rng),
            repetition=repetition,
        )
        for repetition in range(repetitions)
    )
    return (
        summarize_eliminations(
            diagnostics,
            side_length=side_length,
            p=p,
            seed=seed,
        ),
        diagnostics,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=6)
    parser.add_argument("--repetitions", type=int, default=8)
    parser.add_argument("--p", type=float, default=BASELINE_P)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument("--include-environments", action="store_true")
    arguments = parser.parse_args()
    started = perf_counter()
    summary, diagnostics = run_elimination_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        p=arguments.p,
        seed=arguments.seed,
    )
    payload: dict[str, object] = {
        "summary": asdict(summary),
        "benchmark_wall_clock_seconds": perf_counter() - started,
        "benchmark_timing_is_machine_dependent": True,
    }
    if arguments.include_environments:
        payload["environments"] = [asdict(item) for item in diagnostics]
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
