"""Exact tiny-volume diagnostic for two independent Gibbs hierarchies.

The observation is a signed triangular torus generated under the all-plus
truth.  Its posterior ``mu_O`` is enumerated on every spin state.  Two
posterior reference spins are then drawn independently.  Conditional on
each pair ``(O, sigma)``, an independent clock hierarchy is sampled:
satisfied edges receive a uniform rank in ``[0,p]`` and unsatisfied edges a
rank in ``(p,1)``.  In the gauge of ``sigma`` this is precisely the ranked
environment expected by :class:`HierarchicalSweepEnvironment`.

For each hierarchy, the full conditional Gibbs law ``pi_{O,D}`` is
enumerated through ``full_log_weight``.  No postcritical factor is removed.
The two physical pair-correlation matrices are multiplied entrywise.  The
diagnostic reports their normalized total, its exact decomposition over
intersections of final roots, the contribution of the intersection of the
two largest roots, and the diagonal mass of the common refinement of the
two critical cuts.  On the intersection of the two largest roots, it also
splits the signed pair-correlation product into pairs lying in the same
common critical cell and pairs lying in distinct cells.  The latter term is
kept signed; no absolute-value envelope is substituted for it.

The direct posterior matrix ``c_ij = mu_O(sigma_i sigma_j)`` and
``n^-2 sum_ij c_ij^2`` are also computed.  The independent-hierarchy
quantity is an unbiased Monte Carlo estimator of this direct target after
averaging over the posterior references and hierarchies; equality is not
expected for one sampled replica pair.

All Gibbs sums are exact at ``L=4``.  Observation, posterior-reference, and
hierarchy draws remain Monte Carlo.  This is a finite diagnostic, not an
asymptotic double-giant reduction or a weak-recovery result.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import exp, fsum, log, sqrt
from statistics import fmean
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    RankedEdge,
    triangular_torus_edges,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment


BASELINE_P = 0.809439
DEFAULT_MAXIMUM_STATE_COUNT = 1 << 16
DEFAULT_MAXIMUM_OBSERVATION_COUNT = 4
DEFAULT_MAXIMUM_REPLICA_PAIRS = 8
AUDIT_TOLERANCE = 5e-11


@dataclass(frozen=True)
class ScalarEstimate:
    """Sample mean and raw standard error over Monte Carlo draws."""

    mean: float
    standard_error: float | None


@dataclass(frozen=True)
class ReplicaHierarchyDiagnostic:
    """One exact conditional Gibbs enumeration in a sampled hierarchy."""

    replica_index: int
    posterior_spin_seed: int
    hierarchy_seed: int
    reference_state: int
    reference_spins: tuple[int, ...]
    exact_gibbs_state_count: int
    positive_gibbs_state_count: int
    final_component_sizes: tuple[int, ...]
    largest_final_component_index: int
    largest_final_component_vertices: tuple[int, ...]
    critical_component_sizes: tuple[int, ...]
    pair_correlation_matrix: tuple[tuple[float, ...], ...]
    maximum_matrix_diagonal_error: float
    maximum_matrix_symmetry_error: float
    maximum_cross_final_root_correlation: float
    conditional_ranks_match_observation_in_reference_gauge: bool


@dataclass(frozen=True)
class RootIntersectionContribution:
    """Overlap energy supported by one pair of final roots."""

    first_root_index: int
    second_root_index: int
    vertices: tuple[int, ...]
    size: int
    normalized_pair_product_contribution: float


@dataclass(frozen=True)
class ReplicaPairDiagnostic:
    """Two independent hierarchies and their exact replicated decomposition."""

    observation_index: int
    replica_pair_index: int
    first: ReplicaHierarchyDiagnostic
    second: ReplicaHierarchyDiagnostic
    shared_hierarchy_used: bool
    independent_hierarchy_seeds: bool
    normalized_pair_correlation_product: float
    root_intersection_contributions: tuple[RootIntersectionContribution, ...]
    root_intersection_contribution_sum: float
    root_intersection_decomposition_error: float
    separated_in_at_least_one_hierarchy_pair_count: int
    maximum_separated_pair_product_absolute_value: float
    separated_pair_products_zero_within_tolerance: bool
    largest_root_intersection_vertices: tuple[int, ...]
    largest_root_intersection_size: int
    largest_root_intersection_fraction: float
    largest_root_intersection_contribution: float
    common_critical_refinement_cell_sizes: tuple[int, ...]
    common_critical_refinement_diagonal_mass: float
    largest_root_intersection_refinement_cell_sizes: tuple[int, ...]
    largest_root_intersection_refinement_diagonal_mass: float
    largest_root_intersection_critical_cell_diagonal_contribution: float
    largest_root_intersection_critical_cell_off_diagonal_signed_contribution: float
    largest_root_intersection_critical_cell_contribution_sum: float
    largest_root_intersection_critical_cell_decomposition_error: float
    exact_root_intersection_decomposition_passed: bool
    exact_largest_root_intersection_critical_cell_decomposition_passed: bool


@dataclass(frozen=True)
class ObservationDiagnostic:
    """One signed observation, its exact posterior target, and sampled lifts."""

    observation_index: int
    observation_seed: int
    observation_signs: tuple[int, ...]
    positive_observation_edge_count: int
    negative_observation_edge_count: int
    exact_posterior_state_count: int
    direct_posterior_pair_correlation_matrix: tuple[tuple[float, ...], ...]
    direct_posterior_pair_persistence: float
    replica_pairs: tuple[ReplicaPairDiagnostic, ...]
    replicated_hierarchy_estimate: ScalarEstimate
    replicated_estimate_minus_direct_target: float
    largest_root_intersection_contribution_estimate: ScalarEstimate
    largest_root_intersection_critical_cell_diagonal_contribution_estimate: ScalarEstimate
    largest_root_intersection_critical_cell_off_diagonal_signed_contribution_estimate: ScalarEstimate
    maximum_largest_root_intersection_critical_cell_decomposition_error: float
    every_largest_root_intersection_critical_cell_decomposition_passed: bool
    equality_with_direct_target_expected_per_sample: bool
    independent_hierarchy_estimator_targets_direct_persistence: bool


@dataclass(frozen=True)
class DoubleGiantReplicatedGibbsSummary:
    """Finite exact-enumeration summary with explicit Monte Carlo scope."""

    side_length: int
    vertex_count: int
    edge_count: int
    p: float
    critical_rank: float
    final_rank: float
    observation_count: int
    replica_pairs_per_observation: int
    total_replica_pair_count: int
    seed: int
    maximum_state_count: int
    maximum_observation_count: int
    maximum_replica_pairs_per_observation: int
    direct_posterior_pair_persistence: ScalarEstimate
    independent_hierarchy_pair_product_estimate: ScalarEstimate
    mean_replicated_estimate_minus_matched_direct_target: float
    matched_independent_minus_direct_estimate: ScalarEstimate
    largest_root_intersection_contribution_estimate: ScalarEstimate
    largest_root_intersection_critical_cell_diagonal_contribution_estimate: ScalarEstimate
    largest_root_intersection_critical_cell_off_diagonal_signed_contribution_estimate: ScalarEstimate
    maximum_root_intersection_decomposition_error: float
    maximum_largest_root_intersection_critical_cell_decomposition_error: float
    maximum_separated_pair_product_absolute_value: float
    every_root_intersection_decomposition_passed: bool
    every_largest_root_intersection_critical_cell_decomposition_passed: bool
    every_separated_pair_product_zero_within_tolerance: bool
    exact_mu_o_enumeration: bool
    exact_pi_o_d_enumeration: bool
    posterior_reference_draws_are_monte_carlo: bool
    hierarchy_draws_are_monte_carlo: bool
    summary_standard_errors_clustered_by_observation: bool
    two_hierarchies_independent_conditional_on_observation_and_references: bool
    shared_hierarchy_used: bool
    independent_hierarchy_estimator_targets_q_direct: bool
    weak_recovery_claimed: bool
    asymptotic_double_giant_reduction_claimed: bool
    interpretation: str


@dataclass(frozen=True)
class _ReplicaWork:
    """Public diagnostic plus partitions needed for the paired audit."""

    diagnostic: ReplicaHierarchyDiagnostic
    final_components: tuple[tuple[int, ...], ...]
    critical_components: tuple[tuple[int, ...], ...]


def _derived_seed(seed: int, *coordinates: int) -> int:
    """Derive stable independent-looking 64-bit streams."""

    value = seed & ((1 << 64) - 1)
    for coordinate in coordinates:
        value ^= (coordinate + 1) * 0x9E3779B97F4A7C15
        value &= (1 << 64) - 1
        value ^= value >> 30
        value = (value * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)
        value ^= value >> 27
        value = (value * 0x94D049BB133111EB) & ((1 << 64) - 1)
        value ^= value >> 31
    return value


def _validate_inputs(
    side_length: int,
    p: float,
    observation_count: int,
    replica_pairs_per_observation: int,
    maximum_state_count: int,
    maximum_observation_count: int,
    maximum_replica_pairs_per_observation: int,
) -> int:
    if side_length != 4:
        raise ValueError("the exact diagnostic is restricted to side_length=4")
    if not (1.0 + Q_CRITICAL) / 2.0 < p < 1.0:
        raise ValueError("p must put the final forest strictly above q_c")
    if observation_count <= 0 or replica_pairs_per_observation <= 0:
        raise ValueError("observation and replica-pair counts must be positive")
    if maximum_state_count <= 0:
        raise ValueError("maximum_state_count must be positive")
    if maximum_observation_count <= 0:
        raise ValueError("maximum_observation_count must be positive")
    if maximum_replica_pairs_per_observation <= 0:
        raise ValueError(
            "maximum_replica_pairs_per_observation must be positive"
        )
    if observation_count > maximum_observation_count:
        raise ValueError("observation_count exceeds maximum_observation_count")
    if replica_pairs_per_observation > maximum_replica_pairs_per_observation:
        raise ValueError(
            "replica_pairs_per_observation exceeds its explicit cap"
        )
    state_count = 1 << (side_length * side_length)
    if state_count > maximum_state_count:
        raise ValueError("exact state enumeration exceeds maximum_state_count")
    return state_count


def state_spins(state: int, vertex_count: int) -> tuple[int, ...]:
    """Return the physical spin vector represented by one bit state."""

    if not 0 <= state < 1 << vertex_count:
        raise ValueError("state is outside the spin space")
    return tuple(-1 if (state >> vertex) & 1 else 1 for vertex in range(vertex_count))


def edge_relation(state: int, first: int, second: int) -> int:
    """Return ``sigma_first * sigma_second`` for a bit state."""

    return -1 if ((state >> first) ^ (state >> second)) & 1 else 1


def generate_signed_observation(
    side_length: int,
    p: float,
    rng: random.Random,
) -> tuple[int, ...]:
    """Generate edge signs under the deterministic all-plus truth."""

    return tuple(
        1 if rng.random() < p else -1
        for _ in triangular_torus_edges(side_length)
    )


def enumerate_observation_posterior(
    side_length: int,
    observation_signs: Sequence[int],
    p: float,
    maximum_state_count: int = DEFAULT_MAXIMUM_STATE_COUNT,
) -> tuple[float, ...]:
    """Enumerate and normalize the exact posterior ``mu_O``."""

    edges = triangular_torus_edges(side_length)
    if len(observation_signs) != len(edges):
        raise ValueError("one observation sign is required for every edge")
    if any(sign not in (-1, 1) for sign in observation_signs):
        raise ValueError("observation signs must belong to {-1,+1}")
    vertex_count = side_length * side_length
    state_count = 1 << vertex_count
    if state_count > maximum_state_count:
        raise ValueError("exact state enumeration exceeds maximum_state_count")
    coupling = log(p / (1.0 - p))
    log_weights = []
    for state in range(state_count):
        satisfied_count = sum(
            sign * edge_relation(state, first, second) == 1
            for sign, (first, second) in zip(
                observation_signs,
                edges,
                strict=True,
            )
        )
        log_weights.append(coupling * satisfied_count)
    return normalize_log_weights(log_weights)


def normalize_log_weights(log_weights: Sequence[float]) -> tuple[float, ...]:
    """Normalize a finite vector of possibly vanishing log weights."""

    if not log_weights:
        raise ValueError("at least one log weight is required")
    maximum = max(log_weights)
    if maximum == float("-inf"):
        raise AssertionError("all states have zero conditional mass")
    raw = tuple(
        0.0 if value == float("-inf") else exp(value - maximum)
        for value in log_weights
    )
    normalizer = fsum(raw)
    if normalizer <= 0.0:
        raise AssertionError("the conditional normalizer vanished")
    return tuple(value / normalizer for value in raw)


def sample_categorical(weights: Sequence[float], rng: random.Random) -> int:
    """Draw an index from an explicitly normalized probability vector."""

    if not weights or min(weights) < 0.0:
        raise ValueError("weights must be a nonempty nonnegative sequence")
    if abs(fsum(weights) - 1.0) > 2e-12:
        raise ValueError("weights must sum to one")
    draw = rng.random()
    cumulative = 0.0
    for index, weight in enumerate(weights):
        cumulative += weight
        if draw <= cumulative:
            return index
    return len(weights) - 1


def walsh_transform(probabilities: Sequence[float]) -> tuple[float, ...]:
    """Return every Walsh coefficient of a probability vector."""

    state_count = len(probabilities)
    if state_count <= 0 or state_count & (state_count - 1):
        raise ValueError("the probability vector length must be a power of two")
    values = list(probabilities)
    width = 1
    while width < state_count:
        step = 2 * width
        for start in range(0, state_count, step):
            for offset in range(width):
                first = values[start + offset]
                second = values[start + offset + width]
                values[start + offset] = first + second
                values[start + offset + width] = first - second
        width = step
    return tuple(values)


def physical_pair_correlation_matrix(
    probabilities: Sequence[float],
    reference_state: int,
) -> tuple[tuple[float, ...], ...]:
    """Compute all physical pair correlations by an exact Walsh transform."""

    state_count = len(probabilities)
    if state_count <= 0 or state_count & (state_count - 1):
        raise ValueError("the state count must be a positive power of two")
    vertex_count = state_count.bit_length() - 1
    if not 0 <= reference_state < state_count:
        raise ValueError("reference_state is outside the spin space")
    coefficients = walsh_transform(probabilities)
    matrix = [[0.0] * vertex_count for _ in range(vertex_count)]
    for first in range(vertex_count):
        matrix[first][first] = coefficients[0]
        for second in range(first + 1, vertex_count):
            gauge_correlation = coefficients[(1 << first) | (1 << second)]
            physical = (
                edge_relation(reference_state, first, second)
                * gauge_correlation
            )
            matrix[first][second] = physical
            matrix[second][first] = physical
    return tuple(tuple(row) for row in matrix)


def pair_persistence(matrix: Sequence[Sequence[float]]) -> float:
    """Return ``n^-2`` times the squared Frobenius norm."""

    vertex_count = len(matrix)
    if vertex_count <= 0 or any(len(row) != vertex_count for row in matrix):
        raise ValueError("a nonempty square matrix is required")
    return fsum(value * value for row in matrix for value in row) / (
        vertex_count * vertex_count
    )


def sample_conditional_ranked_edges(
    side_length: int,
    observation_signs: Sequence[int],
    reference_state: int,
    p: float,
    rng: random.Random,
) -> tuple[RankedEdge, ...]:
    """Sample clocks conditional on ``(O, reference_state)``.

    Satisfied edges have a finite exponential clock, represented by a
    uniform percolation rank in ``[0,p]``.  Unsatisfied edges have no finite
    clock and receive an irrelevant rank strictly above ``p``.
    """

    edges = triangular_torus_edges(side_length)
    if len(observation_signs) != len(edges):
        raise ValueError("one observation sign is required for every edge")
    ranked_edges = []
    for sign, (first, second) in zip(observation_signs, edges, strict=True):
        satisfied = sign * edge_relation(reference_state, first, second) == 1
        if satisfied:
            rank = p * rng.random()
        else:
            rank = p + (1.0 - p) * rng.random()
            if rank <= p:
                from math import nextafter

                rank = nextafter(p, 1.0)
        ranked_edges.append((rank, first, second))
    return tuple(ranked_edges)


def conditional_rank_audit(
    observation_signs: Sequence[int],
    reference_state: int,
    ranked_edges: Sequence[RankedEdge],
    p: float,
) -> bool:
    """Check that rank satisfaction equals the observation in the new gauge."""

    return all(
        (rank <= p)
        == (sign * edge_relation(reference_state, first, second) == 1)
        for sign, (rank, first, second) in zip(
            observation_signs,
            ranked_edges,
            strict=True,
        )
    )


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
    if any(value < 0 for value in membership):
        raise AssertionError("components do not cover every vertex")
    return tuple(membership)


def _largest_component_index(components: Sequence[Sequence[int]]) -> int:
    if not components:
        raise ValueError("at least one component is required")
    return max(
        range(len(components)),
        key=lambda index: (len(components[index]), -min(components[index])),
    )


def _matrix_audits(
    matrix: Sequence[Sequence[float]],
    components: Sequence[Sequence[int]],
) -> tuple[float, float, float]:
    vertex_count = len(matrix)
    membership = _component_membership(vertex_count, components)
    diagonal_error = max(
        abs(matrix[index][index] - 1.0) for index in range(vertex_count)
    )
    symmetry_error = max(
        abs(matrix[first][second] - matrix[second][first])
        for first in range(vertex_count)
        for second in range(vertex_count)
    )
    cross_root = max(
        (
            abs(matrix[first][second])
            for first in range(vertex_count)
            for second in range(vertex_count)
            if membership[first] != membership[second]
        ),
        default=0.0,
    )
    return diagonal_error, symmetry_error, cross_root


def analyze_replica(
    side_length: int,
    observation_signs: Sequence[int],
    observation_posterior: Sequence[float],
    p: float,
    replica_index: int,
    posterior_spin_seed: int,
    hierarchy_seed: int,
    maximum_state_count: int,
) -> _ReplicaWork:
    """Draw one reference and enumerate its exact hierarchy Gibbs law."""

    vertex_count = side_length * side_length
    reference_state = sample_categorical(
        observation_posterior,
        random.Random(posterior_spin_seed),
    )
    ranked_edges = sample_conditional_ranked_edges(
        side_length,
        observation_signs,
        reference_state,
        p,
        random.Random(hierarchy_seed),
    )
    environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
    state_count = 1 << vertex_count
    if state_count > maximum_state_count:
        raise ValueError("exact hierarchy enumeration exceeds maximum_state_count")
    probabilities = normalize_log_weights(
        tuple(
            environment.full_log_weight(state)
            for state in range(state_count)
        )
    )
    matrix = physical_pair_correlation_matrix(probabilities, reference_state)
    final_components = tuple(
        tuple(component) for component in environment.forest.components
    )
    critical_forest = CriticalKruskalForest(
        side_length,
        ranked_edges,
        critical_rank=Q_CRITICAL,
    )
    critical_components = tuple(
        tuple(component) for component in critical_forest.components
    )
    largest_index = _largest_component_index(final_components)
    diagonal_error, symmetry_error, cross_root = _matrix_audits(
        matrix,
        final_components,
    )
    diagnostic = ReplicaHierarchyDiagnostic(
        replica_index=replica_index,
        posterior_spin_seed=posterior_spin_seed,
        hierarchy_seed=hierarchy_seed,
        reference_state=reference_state,
        reference_spins=state_spins(reference_state, vertex_count),
        exact_gibbs_state_count=state_count,
        positive_gibbs_state_count=sum(value > 0.0 for value in probabilities),
        final_component_sizes=tuple(
            len(component) for component in final_components
        ),
        largest_final_component_index=largest_index,
        largest_final_component_vertices=final_components[largest_index],
        critical_component_sizes=tuple(
            len(component) for component in critical_components
        ),
        pair_correlation_matrix=matrix,
        maximum_matrix_diagonal_error=diagonal_error,
        maximum_matrix_symmetry_error=symmetry_error,
        maximum_cross_final_root_correlation=cross_root,
        conditional_ranks_match_observation_in_reference_gauge=(
            conditional_rank_audit(
                observation_signs,
                reference_state,
                ranked_edges,
                p,
            )
        ),
    )
    return _ReplicaWork(
        diagnostic=diagnostic,
        final_components=final_components,
        critical_components=critical_components,
    )


def _common_refinement_cells(
    vertex_count: int,
    first_components: Sequence[Sequence[int]],
    second_components: Sequence[Sequence[int]],
    allowed_vertices: set[int] | None = None,
) -> tuple[tuple[int, ...], ...]:
    first_membership = _component_membership(vertex_count, first_components)
    second_membership = _component_membership(vertex_count, second_components)
    cells: dict[tuple[int, int], list[int]] = {}
    for vertex in range(vertex_count):
        if allowed_vertices is not None and vertex not in allowed_vertices:
            continue
        key = (first_membership[vertex], second_membership[vertex])
        cells.setdefault(key, []).append(vertex)
    return tuple(
        tuple(vertices)
        for _, vertices in sorted(cells.items())
        if vertices
    )


def analyze_replica_pair(
    observation_index: int,
    replica_pair_index: int,
    first: _ReplicaWork,
    second: _ReplicaWork,
    tolerance: float = AUDIT_TOLERANCE,
) -> ReplicaPairDiagnostic:
    """Compute the exact root-intersection decomposition for two hierarchies."""

    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    first_matrix = first.diagnostic.pair_correlation_matrix
    second_matrix = second.diagnostic.pair_correlation_matrix
    vertex_count = len(first_matrix)
    if len(second_matrix) != vertex_count:
        raise ValueError("replicas have different vertex counts")
    product = tuple(
        tuple(
            first_matrix[first_vertex][second_vertex]
            * second_matrix[first_vertex][second_vertex]
            for second_vertex in range(vertex_count)
        )
        for first_vertex in range(vertex_count)
    )
    normalization = vertex_count * vertex_count
    total = fsum(value for row in product for value in row) / normalization

    first_membership = _component_membership(
        vertex_count,
        first.final_components,
    )
    second_membership = _component_membership(
        vertex_count,
        second.final_components,
    )
    separated_values = tuple(
        abs(product[first_vertex][second_vertex])
        for first_vertex in range(vertex_count)
        for second_vertex in range(vertex_count)
        if (
            first_membership[first_vertex] != first_membership[second_vertex]
            or second_membership[first_vertex]
            != second_membership[second_vertex]
        )
    )

    root_cells = _common_refinement_cells(
        vertex_count,
        first.final_components,
        second.final_components,
    )
    contributions = []
    for vertices in root_cells:
        first_root = first_membership[vertices[0]]
        second_root = second_membership[vertices[0]]
        contribution = fsum(
            product[first_vertex][second_vertex]
            for first_vertex in vertices
            for second_vertex in vertices
        ) / normalization
        if contribution < 0.0 and contribution > -tolerance:
            contribution = 0.0
        contributions.append(
            RootIntersectionContribution(
                first_root_index=first_root,
                second_root_index=second_root,
                vertices=vertices,
                size=len(vertices),
                normalized_pair_product_contribution=contribution,
            )
        )
    contribution_sum = fsum(
        item.normalized_pair_product_contribution for item in contributions
    )
    decomposition_error = abs(total - contribution_sum)

    first_largest = first.diagnostic.largest_final_component_index
    second_largest = second.diagnostic.largest_final_component_index
    largest_contribution = next(
        (
            item
            for item in contributions
            if item.first_root_index == first_largest
            and item.second_root_index == second_largest
        ),
        None,
    )
    largest_vertices = (
        largest_contribution.vertices
        if largest_contribution is not None
        else ()
    )
    largest_vertex_set = set(largest_vertices)

    critical_cells = _common_refinement_cells(
        vertex_count,
        first.critical_components,
        second.critical_components,
    )
    critical_sizes = tuple(len(cell) for cell in critical_cells)
    critical_diagonal = fsum(size * size for size in critical_sizes) / normalization
    largest_critical_cells = _common_refinement_cells(
        vertex_count,
        first.critical_components,
        second.critical_components,
        allowed_vertices=largest_vertex_set,
    )
    largest_critical_sizes = tuple(len(cell) for cell in largest_critical_cells)
    largest_critical_diagonal = fsum(
        size * size for size in largest_critical_sizes
    ) / normalization
    if sum(largest_critical_sizes) != len(largest_vertices):
        raise AssertionError(
            "common critical cells do not partition the largest-root intersection"
        )
    largest_critical_cell_diagonal_contribution = fsum(
        product[first_vertex][second_vertex]
        for cell in largest_critical_cells
        for first_vertex in cell
        for second_vertex in cell
    ) / normalization
    largest_critical_cell_off_diagonal_signed_contribution = fsum(
        product[first_vertex][second_vertex]
        for first_cell_index, first_cell in enumerate(largest_critical_cells)
        for second_cell_index, second_cell in enumerate(largest_critical_cells)
        if first_cell_index != second_cell_index
        for first_vertex in first_cell
        for second_vertex in second_cell
    ) / normalization
    largest_critical_cell_contribution_sum = (
        largest_critical_cell_diagonal_contribution
        + largest_critical_cell_off_diagonal_signed_contribution
    )
    largest_total_contribution = (
        largest_contribution.normalized_pair_product_contribution
        if largest_contribution is not None
        else 0.0
    )
    largest_critical_cell_decomposition_error = abs(
        largest_total_contribution - largest_critical_cell_contribution_sum
    )
    if (
        abs(largest_critical_cell_diagonal_contribution)
        > largest_critical_diagonal + tolerance
    ):
        raise AssertionError(
            "signed same-cell contribution exceeds its geometric mass bound"
        )
    maximum_separated = max(separated_values, default=0.0)
    passed = decomposition_error <= tolerance
    largest_critical_cell_passed = (
        largest_critical_cell_decomposition_error <= tolerance
    )
    return ReplicaPairDiagnostic(
        observation_index=observation_index,
        replica_pair_index=replica_pair_index,
        first=first.diagnostic,
        second=second.diagnostic,
        shared_hierarchy_used=False,
        independent_hierarchy_seeds=(
            first.diagnostic.hierarchy_seed
            != second.diagnostic.hierarchy_seed
        ),
        normalized_pair_correlation_product=total,
        root_intersection_contributions=tuple(contributions),
        root_intersection_contribution_sum=contribution_sum,
        root_intersection_decomposition_error=decomposition_error,
        separated_in_at_least_one_hierarchy_pair_count=len(separated_values),
        maximum_separated_pair_product_absolute_value=maximum_separated,
        separated_pair_products_zero_within_tolerance=(
            maximum_separated <= tolerance
        ),
        largest_root_intersection_vertices=largest_vertices,
        largest_root_intersection_size=len(largest_vertices),
        largest_root_intersection_fraction=(
            len(largest_vertices) / vertex_count
        ),
        largest_root_intersection_contribution=largest_total_contribution,
        common_critical_refinement_cell_sizes=critical_sizes,
        common_critical_refinement_diagonal_mass=critical_diagonal,
        largest_root_intersection_refinement_cell_sizes=(
            largest_critical_sizes
        ),
        largest_root_intersection_refinement_diagonal_mass=(
            largest_critical_diagonal
        ),
        largest_root_intersection_critical_cell_diagonal_contribution=(
            largest_critical_cell_diagonal_contribution
        ),
        largest_root_intersection_critical_cell_off_diagonal_signed_contribution=(
            largest_critical_cell_off_diagonal_signed_contribution
        ),
        largest_root_intersection_critical_cell_contribution_sum=(
            largest_critical_cell_contribution_sum
        ),
        largest_root_intersection_critical_cell_decomposition_error=(
            largest_critical_cell_decomposition_error
        ),
        exact_root_intersection_decomposition_passed=passed,
        exact_largest_root_intersection_critical_cell_decomposition_passed=(
            largest_critical_cell_passed
        ),
    )


def scalar_estimate(values: Sequence[float]) -> ScalarEstimate:
    """Return a sample mean and raw standard error."""

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


def analyze_observation(
    side_length: int,
    p: float,
    observation_index: int,
    observation_seed: int,
    replica_pairs_per_observation: int,
    seed: int,
    maximum_state_count: int,
) -> ObservationDiagnostic:
    """Enumerate one posterior and sample its independent hierarchy pairs."""

    observation = generate_signed_observation(
        side_length,
        p,
        random.Random(observation_seed),
    )
    posterior = enumerate_observation_posterior(
        side_length,
        observation,
        p,
        maximum_state_count,
    )
    direct_matrix = physical_pair_correlation_matrix(posterior, 0)
    direct_persistence = pair_persistence(direct_matrix)
    pairs = []
    for pair_index in range(replica_pairs_per_observation):
        work = []
        for replica_index in range(2):
            posterior_seed = _derived_seed(
                seed,
                101,
                observation_index,
                pair_index,
                replica_index,
            )
            hierarchy_seed = _derived_seed(
                seed,
                211,
                observation_index,
                pair_index,
                replica_index,
            )
            work.append(
                analyze_replica(
                    side_length=side_length,
                    observation_signs=observation,
                    observation_posterior=posterior,
                    p=p,
                    replica_index=replica_index,
                    posterior_spin_seed=posterior_seed,
                    hierarchy_seed=hierarchy_seed,
                    maximum_state_count=maximum_state_count,
                )
            )
        pairs.append(
            analyze_replica_pair(
                observation_index,
                pair_index,
                work[0],
                work[1],
            )
        )
    pair_values = tuple(
        item.normalized_pair_correlation_product for item in pairs
    )
    largest_root_values = tuple(
        item.largest_root_intersection_contribution for item in pairs
    )
    largest_root_diagonal_values = tuple(
        item.largest_root_intersection_critical_cell_diagonal_contribution
        for item in pairs
    )
    largest_root_off_diagonal_values = tuple(
        item.largest_root_intersection_critical_cell_off_diagonal_signed_contribution
        for item in pairs
    )
    estimate = scalar_estimate(pair_values)
    return ObservationDiagnostic(
        observation_index=observation_index,
        observation_seed=observation_seed,
        observation_signs=observation,
        positive_observation_edge_count=sum(sign == 1 for sign in observation),
        negative_observation_edge_count=sum(sign == -1 for sign in observation),
        exact_posterior_state_count=len(posterior),
        direct_posterior_pair_correlation_matrix=direct_matrix,
        direct_posterior_pair_persistence=direct_persistence,
        replica_pairs=tuple(pairs),
        replicated_hierarchy_estimate=estimate,
        replicated_estimate_minus_direct_target=(
            estimate.mean - direct_persistence
        ),
        largest_root_intersection_contribution_estimate=scalar_estimate(
            largest_root_values
        ),
        largest_root_intersection_critical_cell_diagonal_contribution_estimate=(
            scalar_estimate(largest_root_diagonal_values)
        ),
        largest_root_intersection_critical_cell_off_diagonal_signed_contribution_estimate=(
            scalar_estimate(largest_root_off_diagonal_values)
        ),
        maximum_largest_root_intersection_critical_cell_decomposition_error=max(
            item.largest_root_intersection_critical_cell_decomposition_error
            for item in pairs
        ),
        every_largest_root_intersection_critical_cell_decomposition_passed=all(
            item.exact_largest_root_intersection_critical_cell_decomposition_passed
            for item in pairs
        ),
        equality_with_direct_target_expected_per_sample=False,
        independent_hierarchy_estimator_targets_direct_persistence=True,
    )


def run_diagnostic(
    side_length: int = 4,
    p: float = BASELINE_P,
    observation_count: int = 1,
    replica_pairs_per_observation: int = 1,
    seed: int = 20260726,
    maximum_state_count: int = DEFAULT_MAXIMUM_STATE_COUNT,
    maximum_observation_count: int = DEFAULT_MAXIMUM_OBSERVATION_COUNT,
    maximum_replica_pairs_per_observation: int = DEFAULT_MAXIMUM_REPLICA_PAIRS,
) -> tuple[
    DoubleGiantReplicatedGibbsSummary,
    tuple[ObservationDiagnostic, ...],
]:
    """Run the seeded exact-enumeration, Monte Carlo hierarchy diagnostic."""

    _validate_inputs(
        side_length,
        p,
        observation_count,
        replica_pairs_per_observation,
        maximum_state_count,
        maximum_observation_count,
        maximum_replica_pairs_per_observation,
    )
    observations = tuple(
        analyze_observation(
            side_length=side_length,
            p=p,
            observation_index=observation_index,
            observation_seed=_derived_seed(seed, 17, observation_index),
            replica_pairs_per_observation=replica_pairs_per_observation,
            seed=seed,
            maximum_state_count=maximum_state_count,
        )
        for observation_index in range(observation_count)
    )
    pairs = tuple(
        pair
        for observation in observations
        for pair in observation.replica_pairs
    )
    direct_values = tuple(
        observation.direct_posterior_pair_persistence
        for observation in observations
    )
    observation_pair_means = tuple(
        observation.replicated_hierarchy_estimate.mean
        for observation in observations
    )
    observation_largest_root_means = tuple(
        observation.largest_root_intersection_contribution_estimate.mean
        for observation in observations
    )
    observation_largest_root_diagonal_means = tuple(
        observation
        .largest_root_intersection_critical_cell_diagonal_contribution_estimate
        .mean
        for observation in observations
    )
    observation_largest_root_off_diagonal_means = tuple(
        observation
        .largest_root_intersection_critical_cell_off_diagonal_signed_contribution_estimate
        .mean
        for observation in observations
    )
    matched_differences = tuple(
        pair.normalized_pair_correlation_product
        - observations[pair.observation_index].direct_posterior_pair_persistence
        for pair in pairs
    )
    observation_matched_differences = tuple(
        observation.replicated_hierarchy_estimate.mean
        - observation.direct_posterior_pair_persistence
        for observation in observations
    )
    summary = DoubleGiantReplicatedGibbsSummary(
        side_length=side_length,
        vertex_count=side_length * side_length,
        edge_count=len(triangular_torus_edges(side_length)),
        p=p,
        critical_rank=Q_CRITICAL,
        final_rank=2.0 * p - 1.0,
        observation_count=observation_count,
        replica_pairs_per_observation=replica_pairs_per_observation,
        total_replica_pair_count=len(pairs),
        seed=seed,
        maximum_state_count=maximum_state_count,
        maximum_observation_count=maximum_observation_count,
        maximum_replica_pairs_per_observation=(
            maximum_replica_pairs_per_observation
        ),
        direct_posterior_pair_persistence=scalar_estimate(direct_values),
        independent_hierarchy_pair_product_estimate=scalar_estimate(
            observation_pair_means
        ),
        mean_replicated_estimate_minus_matched_direct_target=fmean(
            matched_differences
        ),
        matched_independent_minus_direct_estimate=scalar_estimate(
            observation_matched_differences
        ),
        largest_root_intersection_contribution_estimate=scalar_estimate(
            observation_largest_root_means
        ),
        largest_root_intersection_critical_cell_diagonal_contribution_estimate=(
            scalar_estimate(observation_largest_root_diagonal_means)
        ),
        largest_root_intersection_critical_cell_off_diagonal_signed_contribution_estimate=(
            scalar_estimate(observation_largest_root_off_diagonal_means)
        ),
        maximum_root_intersection_decomposition_error=max(
            pair.root_intersection_decomposition_error for pair in pairs
        ),
        maximum_largest_root_intersection_critical_cell_decomposition_error=max(
            pair.largest_root_intersection_critical_cell_decomposition_error
            for pair in pairs
        ),
        maximum_separated_pair_product_absolute_value=max(
            pair.maximum_separated_pair_product_absolute_value
            for pair in pairs
        ),
        every_root_intersection_decomposition_passed=all(
            pair.exact_root_intersection_decomposition_passed for pair in pairs
        ),
        every_largest_root_intersection_critical_cell_decomposition_passed=all(
            pair.exact_largest_root_intersection_critical_cell_decomposition_passed
            for pair in pairs
        ),
        every_separated_pair_product_zero_within_tolerance=all(
            pair.separated_pair_products_zero_within_tolerance for pair in pairs
        ),
        exact_mu_o_enumeration=True,
        exact_pi_o_d_enumeration=True,
        posterior_reference_draws_are_monte_carlo=True,
        hierarchy_draws_are_monte_carlo=True,
        summary_standard_errors_clustered_by_observation=True,
        two_hierarchies_independent_conditional_on_observation_and_references=True,
        shared_hierarchy_used=False,
        independent_hierarchy_estimator_targets_q_direct=True,
        weak_recovery_claimed=False,
        asymptotic_double_giant_reduction_claimed=False,
        interpretation=(
            "finite L=4 exact posterior and conditional-Gibbs enumerations; "
            "independent posterior-reference and hierarchy draws estimate Q_direct "
            "only after Monte Carlo averaging; summary standard errors use "
            "independent observations as clusters, with no asymptotic claim"
        ),
    )
    return summary, observations


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=4)
    parser.add_argument("--p", type=float, default=BASELINE_P)
    parser.add_argument("--observations", type=int, default=1)
    parser.add_argument("--replica-pairs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument(
        "--maximum-state-count",
        type=int,
        default=DEFAULT_MAXIMUM_STATE_COUNT,
    )
    parser.add_argument(
        "--maximum-observations",
        type=int,
        default=DEFAULT_MAXIMUM_OBSERVATION_COUNT,
    )
    parser.add_argument(
        "--maximum-replica-pairs",
        type=int,
        default=DEFAULT_MAXIMUM_REPLICA_PAIRS,
    )
    arguments = parser.parse_args()
    summary, observations = run_diagnostic(
        side_length=arguments.side,
        p=arguments.p,
        observation_count=arguments.observations,
        replica_pairs_per_observation=arguments.replica_pairs,
        seed=arguments.seed,
        maximum_state_count=arguments.maximum_state_count,
        maximum_observation_count=arguments.maximum_observations,
        maximum_replica_pairs_per_observation=(
            arguments.maximum_replica_pairs
        ),
    )
    print(
        json.dumps(
            {
                "summary": asdict(summary),
                "observations": [asdict(item) for item in observations],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
