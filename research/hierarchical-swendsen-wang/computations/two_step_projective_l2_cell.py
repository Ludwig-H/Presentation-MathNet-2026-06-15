"""Exact two-step projective L2 cell from one realized hierarchy.

The active cell consists of two consecutive strict-arm merger nodes.  At the
lower node, the endpoint branch and one off-spine attachment are integrated.
At its parent, the only new flip direction is the ancestral sibling.  These
are three disjoint cluster-orientation bits.  The sibling contains neither
endpoint, so projecting the original pair character along that direction
alone has exactly zero loss.  A nonzero second loss therefore measures the
effect of enlarging the *already collapsed* first block.

No exterior field is chosen by hand.  The full posterior ``pi_D`` is
enumerated on a triangular torus of side four.  For every positive outside
coset, the strict-ancestor factors induce an eight-entry projective exterior
potential on the three active bits.  Thus the finite family of potentials
reported here is attained by an actual Nishimori/Kruskal environment.

The default example was selected after an exploratory finite scan because it
has two positive Pythagorean losses.  It is a witness, not an unbiased Palm
sample and not evidence that such cells have positive asymptotic density.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, isfinite, log
from statistics import median
from typing import Sequence

from critical_pair_path_geometry import (
    clock_time_from_rank,
    sample_ranked_edges,
    triangular_torus_distance,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from nested_projection_l2_diagnostic import (
    _quotient_key,
    _sample_uniform_distant_pair,
    _xor_basis,
    collapsed_projection,
    pair_character,
    posterior_weights,
    weighted_norm_square,
)

DEFAULT_SEED = 20260726
DEFAULT_REPETITION = 60
DEFAULT_P = 0.805
DEFAULT_PAIR = (0, 13)
DEFAULT_LOWER_NODE = 21


@dataclass(frozen=True)
class TwoStepCellGeometry:
    """Three active cluster bits in two consecutive strict-arm updates."""

    side_length: int
    p: float
    pair_first: int
    pair_second: int
    pair_distance: int
    endpoint_below_cell: int
    pair_lca: int
    lower_node: int
    upper_node: int
    spine_child: int
    attachment_child: int
    ancestral_sibling: int
    spine_size: int
    attachment_size: int
    ancestral_sibling_size: int
    lower_bucket_size: int
    upper_bucket_size: int
    lower_rank: float
    upper_rank: float
    lower_beta: float
    upper_beta: float
    first_flip_rank: int
    cumulative_flip_rank: int
    upper_is_strictly_before_lca: bool
    ancestral_sibling_contains_no_endpoint: bool


@dataclass(frozen=True)
class AttainableExteriorPotential:
    """One exact projective potential conditional on an outside coset."""

    outside_coset_key: int
    posterior_mass: float
    positive_posterior_orientation_count: int
    strictly_positive_exterior: bool
    exterior_projective_log_weights: tuple[float, ...]
    cell_projective_log_weights: tuple[float, ...]
    full_projective_log_weights: tuple[float, ...]
    incoming_energy: float
    second_absolute_loss: float
    second_relative_loss: float | None
    exterior_projective_span: float | None
    projective_factorization_error: float


@dataclass(frozen=True)
class TwoStepProjectiveL2Diagnostic:
    """Exact losses and the attained exterior-potential family."""

    geometry: TwoStepCellGeometry
    initial_norm_square: float
    after_first_norm_square: float
    after_second_norm_square: float
    first_absolute_loss: float
    first_relative_loss: float
    second_absolute_loss: float
    second_relative_loss: float
    first_difference_norm_square: float
    second_difference_norm_square: float
    maximum_pythagorean_error: float
    direct_ancestral_sibling_loss_on_pair_character: float
    attainable_potential_count: int
    strictly_positive_exterior_potential_count: int
    boundary_exterior_potential_count: int
    strictly_positive_exterior_posterior_mass: float
    strictly_positive_exterior_incoming_energy: float
    strictly_positive_exterior_second_loss: float
    strictly_positive_exterior_energy_ratio: float | None
    minimum_strictly_positive_exterior_second_relative_loss: float | None
    median_strictly_positive_exterior_second_relative_loss: float | None
    minimum_all_attainable_second_relative_loss: float
    fraction_incoming_energy_on_positive_second_loss: float
    maximum_projective_factorization_error: float
    orbit_second_loss_audit_error: float
    potential_orbits: tuple[AttainableExteriorPotential, ...]
    full_posterior_enumerated: bool
    exterior_potentials_are_attained: bool
    second_projection_applied_to_first_output: bool
    witness_selected_after_exploratory_scan: bool
    unbiased_palm_sample: bool
    weak_recovery_claimed: bool
    interpretation: str


def _factor_log_weight(
    environment: HierarchicalSweepEnvironment,
    factor_node: int,
    state: int,
) -> float:
    count = 0
    for first, second in environment.bucket_edges[factor_node]:
        edge = tuple(sorted((first, second)))
        truth_satisfied = environment.edge_rank[edge] <= environment.p
        relative_parity = ((state >> first) ^ (state >> second)) & 1
        count += truth_satisfied != bool(relative_parity)
    if count == 0:
        return float("-inf")
    beta = clock_time_from_rank(
        environment.forest.merge_rank[factor_node], environment.p
    )
    return log(count) + (1.0 - beta) * environment.coupling * count


def _sum_log_factors(
    environment: HierarchicalSweepEnvironment,
    factor_nodes: Sequence[int],
    state: int,
) -> float:
    values = tuple(
        _factor_log_weight(environment, node, state) for node in factor_nodes
    )
    if any(value == float("-inf") for value in values):
        return float("-inf")
    return fsum(values)


def _projectivize(log_weights: Sequence[float]) -> tuple[float, ...]:
    maximum = max(log_weights)
    if maximum == float("-inf"):
        raise ValueError("a projective potential cannot vanish identically")
    return tuple(
        value if value == float("-inf") else value - maximum for value in log_weights
    )


def _projective_error(first: Sequence[float], second: Sequence[float]) -> float:
    if len(first) != len(second):
        raise ValueError("projective vectors must have the same dimension")
    error = 0.0
    for left, right in zip(first, second, strict=True):
        if left == float("-inf") or right == float("-inf"):
            if left != right:
                return float("inf")
        else:
            error = max(error, abs(left - right))
    return error


def _difference_norm_square(
    first: Sequence[float],
    second: Sequence[float],
    weights: Sequence[float],
) -> float:
    return weighted_norm_square(
        tuple(left - right for left, right in zip(first, second, strict=True)),
        weights,
    )


def _nonnegative_norm_loss(before: float, after: float) -> float:
    loss = before - after
    if -1e-12 < loss < 0.0:
        return 0.0
    if loss < 0.0:
        raise AssertionError("a conditional projection increased the L2 norm")
    return loss


def _active_orbit_masks(generators: Sequence[int]) -> tuple[int, ...]:
    if len(generators) != 3 or len(_xor_basis(generators)) != 3:
        raise ValueError("the two-step cell requires three independent bits")
    return tuple(
        (generators[0] if bits & 1 else 0)
        ^ (generators[1] if bits & 2 else 0)
        ^ (generators[2] if bits & 4 else 0)
        for bits in range(8)
    )


def analyze_two_step_cell(
    environment: HierarchicalSweepEnvironment,
    first: int,
    second: int,
    lower_node: int,
    *,
    witness_selected_after_exploratory_scan: bool = False,
) -> TwoStepProjectiveL2Diagnostic:
    """Analyze two consecutive strict-arm collapsed projections exactly."""

    forest = environment.forest
    if first == second or not forest.connected(first, second):
        raise ValueError("the marked pair must be distinct and connected")
    if lower_node not in forest.internal_nodes:
        raise ValueError("lower_node must be an internal merger node")
    pair_lca = forest.tree_lca(first, second)
    endpoint_membership = tuple(
        bool(environment.cluster_mask[lower_node] & (1 << endpoint))
        for endpoint in (first, second)
    )
    if sum(endpoint_membership) != 1:
        raise ValueError("lower_node must lie on one strict arm of the pair")
    endpoint = first if endpoint_membership[0] else second

    upper_node = forest.tree_parent[lower_node]
    if upper_node == -1 or upper_node == pair_lca:
        raise ValueError("the parent update must also be strictly before the LCA")
    if upper_node not in forest.internal_nodes:
        raise AssertionError("the parent of a merger must be a merger")
    upper_membership = tuple(
        bool(environment.cluster_mask[upper_node] & (1 << marked))
        for marked in (first, second)
    )
    if sum(upper_membership) != 1:
        raise AssertionError("the upper node is not on the same strict arm")

    lower_left = forest.left_child[lower_node]
    lower_right = forest.right_child[lower_node]
    if environment.cluster_mask[lower_left] & (1 << endpoint):
        spine_child, attachment_child = lower_left, lower_right
    else:
        spine_child, attachment_child = lower_right, lower_left

    upper_left = forest.left_child[upper_node]
    upper_right = forest.right_child[upper_node]
    if lower_node == upper_left:
        ancestral_sibling = upper_right
    elif lower_node == upper_right:
        ancestral_sibling = upper_left
    else:
        raise AssertionError("the two merger nodes are not consecutive")

    sibling_mask = environment.cluster_mask[ancestral_sibling]
    sibling_contains_no_endpoint = not any(
        sibling_mask & (1 << marked) for marked in (first, second)
    )
    if not sibling_contains_no_endpoint:
        raise AssertionError("a strict-arm sibling contains no marked endpoint")

    spine_mask = environment.cluster_mask[spine_child]
    attachment_mask = environment.cluster_mask[attachment_child]
    first_generators = (spine_mask, attachment_mask)
    cumulative_generators = first_generators + (sibling_mask,)
    first_rank = len(_xor_basis(first_generators))
    cumulative_rank = len(_xor_basis(cumulative_generators))
    if first_rank != 2 or cumulative_rank != 3:
        raise AssertionError("the three disjoint clusters must give three bits")

    lower_proposal_basis = _xor_basis(
        tuple(mask for mask in environment.proposal_masks(lower_node) if mask)
    )
    cumulative_proposal_basis = _xor_basis(
        tuple(mask for mask in environment.proposal_masks(lower_node) if mask)
        + tuple(mask for mask in environment.proposal_masks(upper_node) if mask)
    )
    if lower_proposal_basis != _xor_basis(first_generators):
        raise AssertionError("the first group is not the lower update orbit")
    if cumulative_proposal_basis != _xor_basis(cumulative_generators):
        raise AssertionError("the second group is not the cumulative update orbit")

    weights = posterior_weights(environment)
    pair_values = pair_character(len(weights), first, second)
    after_first = collapsed_projection(pair_values, weights, first_generators)
    after_second = collapsed_projection(after_first, weights, cumulative_generators)
    initial_norm = weighted_norm_square(pair_values, weights)
    first_norm = weighted_norm_square(after_first, weights)
    second_norm = weighted_norm_square(after_second, weights)
    first_loss = _nonnegative_norm_loss(initial_norm, first_norm)
    second_loss = _nonnegative_norm_loss(first_norm, second_norm)
    first_difference = _difference_norm_square(pair_values, after_first, weights)
    second_difference = _difference_norm_square(after_first, after_second, weights)

    sibling_projection = collapsed_projection(pair_values, weights, (sibling_mask,))
    direct_sibling_loss = _difference_norm_square(
        pair_values, sibling_projection, weights
    )

    active_masks = _active_orbit_masks(cumulative_generators)
    cumulative_basis = _xor_basis(cumulative_generators)
    orbit_keys = sorted(
        {
            _quotient_key(state, cumulative_basis)
            for state, weight in enumerate(weights)
            if weight > 0.0
        }
    )
    strict_ancestor_nodes = environment.factor_chains[upper_node][1:]
    if environment.factor_chains[lower_node][:2] != (lower_node, upper_node):
        raise AssertionError("the local factor chain is not consecutive")

    potential_orbits = []
    for key in orbit_keys:
        states = tuple(key ^ mask for mask in active_masks)
        orbit_weights = tuple(weights[state] for state in states)
        posterior_mass = fsum(orbit_weights)
        if posterior_mass <= 0.0:
            raise AssertionError("a retained orbit must have positive mass")

        lower_logs = tuple(
            _factor_log_weight(environment, lower_node, state) for state in states
        )
        upper_logs = tuple(
            _factor_log_weight(environment, upper_node, state) for state in states
        )
        exterior_logs = tuple(
            _sum_log_factors(environment, strict_ancestor_nodes, state)
            for state in states
        )
        cell_logs = tuple(
            (
                float("-inf")
                if lower == float("-inf") or upper == float("-inf")
                else lower + upper
            )
            for lower, upper in zip(lower_logs, upper_logs, strict=True)
        )
        combined_logs = tuple(
            (
                float("-inf")
                if cell == float("-inf") or exterior == float("-inf")
                else cell + exterior
            )
            for cell, exterior in zip(cell_logs, exterior_logs, strict=True)
        )
        full_logs = tuple(environment.full_log_weight(state) for state in states)
        exterior_projective = _projectivize(exterior_logs)
        cell_projective = _projectivize(cell_logs)
        combined_projective = _projectivize(combined_logs)
        full_projective = _projectivize(full_logs)
        factorization_error = _projective_error(combined_projective, full_projective)
        if factorization_error > 5e-12:
            raise AssertionError("the two-node cell missed a varying factor")

        incoming_energy = fsum(
            weight * after_first[state] * after_first[state]
            for state, weight in zip(states, orbit_weights, strict=True)
        )
        orbit_second_loss = fsum(
            weight * (after_first[state] - after_second[state]) ** 2
            for state, weight in zip(states, orbit_weights, strict=True)
        )
        relative_loss = (
            None if incoming_energy == 0.0 else orbit_second_loss / incoming_energy
        )
        strictly_positive = all(isfinite(value) for value in exterior_projective)
        exterior_span = (
            max(exterior_projective) - min(exterior_projective)
            if strictly_positive
            else None
        )
        potential_orbits.append(
            AttainableExteriorPotential(
                outside_coset_key=key,
                posterior_mass=posterior_mass,
                positive_posterior_orientation_count=sum(
                    weight > 0.0 for weight in orbit_weights
                ),
                strictly_positive_exterior=strictly_positive,
                exterior_projective_log_weights=exterior_projective,
                cell_projective_log_weights=cell_projective,
                full_projective_log_weights=full_projective,
                incoming_energy=incoming_energy,
                second_absolute_loss=orbit_second_loss,
                second_relative_loss=relative_loss,
                exterior_projective_span=exterior_span,
                projective_factorization_error=factorization_error,
            )
        )

    positive_exterior = tuple(
        orbit for orbit in potential_orbits if orbit.strictly_positive_exterior
    )
    boundary_exterior = tuple(
        orbit for orbit in potential_orbits if not orbit.strictly_positive_exterior
    )
    positive_relative_losses = tuple(
        orbit.second_relative_loss
        for orbit in positive_exterior
        if orbit.second_relative_loss is not None
    )
    all_relative_losses = tuple(
        orbit.second_relative_loss
        for orbit in potential_orbits
        if orbit.second_relative_loss is not None
    )
    if not all_relative_losses:
        raise AssertionError("the attained family has no incoming odd energy")

    positive_mass = fsum(orbit.posterior_mass for orbit in positive_exterior)
    positive_energy = fsum(orbit.incoming_energy for orbit in positive_exterior)
    positive_loss = fsum(orbit.second_absolute_loss for orbit in positive_exterior)
    total_orbit_loss = fsum(orbit.second_absolute_loss for orbit in potential_orbits)
    total_orbit_energy = fsum(orbit.incoming_energy for orbit in potential_orbits)
    energy_on_loss = fsum(
        orbit.incoming_energy
        for orbit in potential_orbits
        if orbit.second_absolute_loss > 1e-14
    )
    orbit_loss_error = abs(total_orbit_loss - second_difference)
    if orbit_loss_error > 5e-12:
        raise AssertionError("orbit losses do not reconstruct the global loss")

    geometry = TwoStepCellGeometry(
        side_length=environment.side_length,
        p=environment.p,
        pair_first=first,
        pair_second=second,
        pair_distance=triangular_torus_distance(first, second, environment.side_length),
        endpoint_below_cell=endpoint,
        pair_lca=pair_lca,
        lower_node=lower_node,
        upper_node=upper_node,
        spine_child=spine_child,
        attachment_child=attachment_child,
        ancestral_sibling=ancestral_sibling,
        spine_size=spine_mask.bit_count(),
        attachment_size=attachment_mask.bit_count(),
        ancestral_sibling_size=sibling_mask.bit_count(),
        lower_bucket_size=forest.bucket_size[lower_node],
        upper_bucket_size=forest.bucket_size[upper_node],
        lower_rank=forest.merge_rank[lower_node],
        upper_rank=forest.merge_rank[upper_node],
        lower_beta=clock_time_from_rank(forest.merge_rank[lower_node], environment.p),
        upper_beta=clock_time_from_rank(forest.merge_rank[upper_node], environment.p),
        first_flip_rank=first_rank,
        cumulative_flip_rank=cumulative_rank,
        upper_is_strictly_before_lca=upper_node != pair_lca,
        ancestral_sibling_contains_no_endpoint=sibling_contains_no_endpoint,
    )
    return TwoStepProjectiveL2Diagnostic(
        geometry=geometry,
        initial_norm_square=initial_norm,
        after_first_norm_square=first_norm,
        after_second_norm_square=second_norm,
        first_absolute_loss=first_loss,
        first_relative_loss=(0.0 if initial_norm == 0.0 else first_loss / initial_norm),
        second_absolute_loss=second_loss,
        second_relative_loss=(0.0 if first_norm == 0.0 else second_loss / first_norm),
        first_difference_norm_square=first_difference,
        second_difference_norm_square=second_difference,
        maximum_pythagorean_error=max(
            abs(first_loss - first_difference),
            abs(second_loss - second_difference),
        ),
        direct_ancestral_sibling_loss_on_pair_character=direct_sibling_loss,
        attainable_potential_count=len(potential_orbits),
        strictly_positive_exterior_potential_count=len(positive_exterior),
        boundary_exterior_potential_count=len(boundary_exterior),
        strictly_positive_exterior_posterior_mass=positive_mass,
        strictly_positive_exterior_incoming_energy=positive_energy,
        strictly_positive_exterior_second_loss=positive_loss,
        strictly_positive_exterior_energy_ratio=(
            None if positive_energy == 0.0 else positive_loss / positive_energy
        ),
        minimum_strictly_positive_exterior_second_relative_loss=(
            None if not positive_relative_losses else min(positive_relative_losses)
        ),
        median_strictly_positive_exterior_second_relative_loss=(
            None if not positive_relative_losses else median(positive_relative_losses)
        ),
        minimum_all_attainable_second_relative_loss=min(all_relative_losses),
        fraction_incoming_energy_on_positive_second_loss=(
            energy_on_loss / total_orbit_energy
        ),
        maximum_projective_factorization_error=max(
            orbit.projective_factorization_error for orbit in potential_orbits
        ),
        orbit_second_loss_audit_error=orbit_loss_error,
        potential_orbits=tuple(potential_orbits),
        full_posterior_enumerated=True,
        exterior_potentials_are_attained=True,
        second_projection_applied_to_first_output=True,
        witness_selected_after_exploratory_scan=(
            witness_selected_after_exploratory_scan
        ),
        unbiased_palm_sample=False,
        weak_recovery_claimed=False,
        interpretation=(
            "exact selected finite witness: the second collapsed enlargement "
            "dissipates the function produced by the first one, even though "
            "the new sibling flip leaves the original pair character fixed; "
            "the positive interior-potential margin is only for this finite "
            "attained family and gives neither abundance nor an asymptotic bound"
        ),
    )


def build_preselected_witness() -> TwoStepProjectiveL2Diagnostic:
    """Reconstruct the declared real-rank witness from the exploratory scan."""

    rng = random.Random(DEFAULT_SEED)
    environment = None
    pair = None
    for _ in range(DEFAULT_REPETITION + 1):
        environment = HierarchicalSweepEnvironment(
            4, sample_ranked_edges(4, rng), DEFAULT_P
        )
        pair = _sample_uniform_distant_pair(4, 0.25, rng)
    if environment is None or pair != DEFAULT_PAIR:
        raise AssertionError("the declared witness is not reproducible")
    return analyze_two_step_cell(
        environment,
        *pair,
        DEFAULT_LOWER_NODE,
        witness_selected_after_exploratory_scan=True,
    )


def _json_safe(value):
    if isinstance(value, float) and not isfinite(value):
        return "+infinity" if value > 0.0 else "-infinity"
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-potentials", action="store_true")
    arguments = parser.parse_args()
    result = build_preselected_witness()
    payload = asdict(result)
    if not arguments.include_potentials:
        payload.pop("potential_orbits")
    print(json.dumps(_json_safe(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
