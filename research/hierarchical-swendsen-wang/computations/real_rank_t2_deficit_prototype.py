"""One-step real-rank T2 deficit audit on a genuine corridor attachment.

This module answers a deliberately narrow question.  Fix one final
Kruskal/Nishimori environment, one distant pair, and one strict-arm merger
whose off-spine child is a small attachment.  The genuine four-state
hierarchical heat bath at that merger flips the spine child and the
attachment child independently.  Every factor at the merger and at all its
strict ancestors is evaluated at its *realized* rank.

For two posterior replicas with independent heat-bath draws in the same
environment, let ``s_k`` be the orientation chosen for the spine child in
replica ``k``.  Since exactly one endpoint of the marked pair lies below the
merger, the replicated twist is

    epsilon = s_1 * s_2.

The module constructs two positive/signed transfer pairs.

``faithful``
    The target state retains both complete output spin states.  It therefore
    retains the exterior boundary exactly.  The twist is measurable from
    the target state, hence ``|U| = K`` on every positive transition and the
    one-step Feynman--Kac deficit is identically zero.  This is a useful
    no-go certificate: a single fully resolved hierarchical update cannot
    supply contraction.

``projected``
    The target retains only the within-merger relative orientations
    ``delta_k = s_k*a_k`` and drops the two global spine orientations.  The
    resulting cancellation can give ``|U| < K``.  It is an exact
    coarse-graining of this one heat bath, but it is *not* a composable T2
    transfer: the dropped orientations change the strict-ancestor potential
    seen by the next update.  Until a port gauge and the induced transition
    of the full projective exterior potential are specified, its deficit is
    only a diagnostic and proves no decorrelation bound.

The exterior state stored here is substantially richer than a scalar field:
for each replica it records the four strict-ancestor factor vectors, their
real ranks, and the full projective four-orientation potential.  No
criticalization, resampling, lateral independence, or independent replica
environment is used.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from itertools import product
from math import copysign, exp, fsum, isfinite, log
from typing import Sequence

from corridor_t2_signature_diagnostic import (
    CorridorT2Signature,
    EventPalmCorridor,
    event_palm_corridor_observations,
)
from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    RankedEdge,
    clock_time_from_rank,
    sample_ranked_edges,
)
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from twisted_feynman_kac_composition import (
    DominatedTransfer,
    dominated_transfer,
    finite_horizon_doob_certificate,
)

Orientation = tuple[int, int]
Matrix = tuple[tuple[float, ...], ...]

LEFT_RIGHT_ORIENTATIONS: tuple[Orientation, ...] = (
    (1, 1),
    (-1, 1),
    (1, -1),
    (-1, -1),
)
REPLICATED_RELATIVE_STATES: tuple[tuple[int, int], ...] = tuple(
    product((-1, 1), repeat=2)
)


@dataclass(frozen=True)
class FactorOrbitState:
    """One real-rank Lambda factor on the four-state update orbit."""

    node: int
    rank: float
    beta: float
    bucket_size: int
    bucket_edges: tuple[tuple[int, int], ...]
    satisfied_counts: tuple[int, ...]
    log_weights: tuple[float, ...]


@dataclass(frozen=True)
class ReplicaExteriorState:
    """Exact local and exterior potential for one replica before the update."""

    reference_state: int
    orbit_masks: tuple[int, ...]
    local_factor: FactorOrbitState
    strict_ancestor_factors: tuple[FactorOrbitState, ...]
    exterior_log_weights: tuple[float, ...]
    exterior_projective_log_weights: tuple[float, ...]
    full_log_weights: tuple[float, ...]
    full_projective_log_weights: tuple[float, ...]
    heat_bath_probabilities: tuple[float, ...]
    heat_bath_audit_error: float


@dataclass(frozen=True)
class AttachmentCellGeometry:
    """Marked strict-arm geometry defining the local T2 update."""

    side_length: int
    p: float
    pair_first: int
    pair_second: int
    pair_lca: int
    pair_palm_weight: int
    node: int
    actual_rank: float
    beta: float
    bucket_size: int
    spine_child: int
    attachment_child: int
    spine_child_size: int
    attachment_child_size: int
    strict_ancestor_nodes: tuple[int, ...]
    is_postcritical: bool
    all_factor_ranks_are_realized: bool = True
    winning_edge_identity_used: bool = False
    criticalization_used: bool = False


@dataclass(frozen=True)
class FaithfulReplicatedTarget:
    """Full output boundary for two replicas and its replicated twist."""

    first_orbit_index: int
    second_orbit_index: int
    first_output_state: int
    second_output_state: int
    first_spine_orientation: int
    second_spine_orientation: int
    first_relative_orientation: int
    second_relative_orientation: int
    replicated_twist: int


@dataclass(frozen=True)
class LocalDeficitSummary:
    """A dominated transfer and its one-step Feynman--Kac diagnostics."""

    transfer: DominatedTransfer
    transition_deficits: tuple[float, ...]
    feynman_kac_envelope: float
    logarithmic_attenuation: float
    finite_horizon_doob_envelope_error: float
    finite_horizon_doob_inequality_holds: bool


@dataclass(frozen=True)
class RealRankT2DeficitPrototype:
    """Faithful no-go and nonclosed projection on one genuine attachment."""

    geometry: AttachmentCellGeometry
    replica_exterior_states: tuple[ReplicaExteriorState, ReplicaExteriorState]
    faithful_targets: tuple[FaithfulReplicatedTarget, ...]
    faithful: LocalDeficitSummary
    projected_relative_states: tuple[tuple[int, int], ...]
    projected: LocalDeficitSummary
    shared_nishimori_environment: bool
    replica_draws_conditionally_independent: bool
    full_boundary_retained_by_faithful_transfer: bool
    projected_boundary_is_markov_closed: bool
    composable_t2_deficit_certified: bool
    palm_abundance_estimated: bool
    weak_recovery_claimed: bool
    interpretation: str


def _component_size(mask: int) -> int:
    return mask.bit_count()


def _factor_satisfied_count(
    environment: HierarchicalSweepEnvironment,
    factor_node: int,
    state: int,
) -> int:
    """Recompute a Lambda count without calling the sweep implementation."""

    count = 0
    for first, second in environment.bucket_edges[factor_node]:
        edge = tuple(sorted((first, second)))
        truth_satisfied = environment.edge_rank[edge] <= environment.p
        relative_parity = ((state >> first) ^ (state >> second)) & 1
        count += truth_satisfied != bool(relative_parity)
    return count


def _factor_orbit_state(
    environment: HierarchicalSweepEnvironment,
    factor_node: int,
    orbit_states: Sequence[int],
) -> FactorOrbitState:
    rank = environment.forest.merge_rank[factor_node]
    beta = clock_time_from_rank(rank, environment.p)
    counts = tuple(
        _factor_satisfied_count(environment, factor_node, state)
        for state in orbit_states
    )
    log_weights = tuple(
        (
            float("-inf")
            if count == 0
            else log(count) + (1.0 - beta) * environment.coupling * count
        )
        for count in counts
    )
    return FactorOrbitState(
        node=factor_node,
        rank=rank,
        beta=beta,
        bucket_size=environment.forest.bucket_size[factor_node],
        bucket_edges=environment.bucket_edges[factor_node],
        satisfied_counts=counts,
        log_weights=log_weights,
    )


def _sum_log_vectors(vectors: Sequence[Sequence[float]]) -> tuple[float, ...]:
    if not vectors:
        return (0.0, 0.0, 0.0, 0.0)
    if any(len(vector) != 4 for vector in vectors):
        raise ValueError("every orbit log vector must have four entries")
    return tuple(
        (
            float("-inf")
            if any(vector[index] == float("-inf") for vector in vectors)
            else fsum(vector[index] for vector in vectors)
        )
        for index in range(4)
    )


def _projectivize(log_weights: Sequence[float]) -> tuple[float, ...]:
    if len(log_weights) != 4:
        raise ValueError("four orbit weights are required")
    maximum = max(log_weights)
    if maximum == float("-inf"):
        raise ValueError("an orbit potential cannot vanish identically")
    return tuple(
        value if value == float("-inf") else value - maximum for value in log_weights
    )


def _softmax(log_weights: Sequence[float]) -> tuple[float, ...]:
    projective = _projectivize(log_weights)
    weights = tuple(
        0.0 if value == float("-inf") else exp(value) for value in projective
    )
    normalizer = fsum(weights)
    return tuple(value / normalizer for value in weights)


def _replica_exterior_state(
    environment: HierarchicalSweepEnvironment,
    node: int,
    reference_state: int,
) -> ReplicaExteriorState:
    if reference_state < 0 or reference_state >= 1 << environment.vertex_count:
        raise ValueError("reference_state is outside the finite spin space")
    orbit_masks = environment.proposal_masks(node)
    if len(orbit_masks) != 4:
        raise ValueError("the T2 prototype requires an internal four-state orbit")
    orbit_states = tuple(reference_state ^ mask for mask in orbit_masks)
    factor_chain = environment.factor_chains[node]
    if not factor_chain or factor_chain[0] != node:
        raise AssertionError("the factor chain must start at the updated node")

    local = _factor_orbit_state(environment, node, orbit_states)
    ancestors = tuple(
        _factor_orbit_state(environment, factor_node, orbit_states)
        for factor_node in factor_chain[1:]
    )
    exterior = _sum_log_vectors(tuple(factor.log_weights for factor in ancestors))
    full = tuple(
        (
            float("-inf")
            if local_weight == float("-inf") or exterior_weight == float("-inf")
            else local_weight + exterior_weight
        )
        for local_weight, exterior_weight in zip(
            local.log_weights, exterior, strict=True
        )
    )
    probabilities = _softmax(full)

    direct_probabilities = tuple(
        probability
        for _, probability in environment.proposal_probabilities(node, reference_state)
    )
    audit_error = max(
        abs(first - second)
        for first, second in zip(probabilities, direct_probabilities, strict=True)
    )
    if audit_error > 3e-12:
        raise AssertionError("factor decomposition disagrees with the exact heat bath")

    return ReplicaExteriorState(
        reference_state=reference_state,
        orbit_masks=orbit_masks,
        local_factor=local,
        strict_ancestor_factors=ancestors,
        exterior_log_weights=exterior,
        exterior_projective_log_weights=_projectivize(exterior),
        full_log_weights=full,
        full_projective_log_weights=_projectivize(full),
        heat_bath_probabilities=probabilities,
        heat_bath_audit_error=audit_error,
    )


def _spine_attachment_geometry(
    environment: HierarchicalSweepEnvironment,
    pair: EventPalmCorridor,
    signature: CorridorT2Signature,
) -> tuple[int, int]:
    if signature.is_lca or signature.node == pair.lca_node:
        raise ValueError("a small attachment must be a strict-arm node")
    if signature.node not in environment.forest.internal_nodes:
        raise ValueError("the signature node is absent from the environment")
    if signature not in pair.signatures:
        raise ValueError("the signature does not belong to the marked corridor")

    node = signature.node
    left = environment.forest.left_child[node]
    right = environment.forest.right_child[node]
    union_mask = environment.cluster_mask[node]
    endpoint_membership = (
        bool(union_mask & (1 << pair.first)),
        bool(union_mask & (1 << pair.second)),
    )
    if sum(endpoint_membership) != 1:
        raise AssertionError("a strict arm contains exactly one marked endpoint")
    endpoint = pair.first if endpoint_membership[0] else pair.second
    if environment.cluster_mask[left] & (1 << endpoint):
        spine, attachment = left, right
    elif environment.cluster_mask[right] & (1 << endpoint):
        spine, attachment = right, left
    else:
        raise AssertionError("the marked endpoint is missing from both children")

    if signature.spine_child_size != _component_size(environment.cluster_mask[spine]):
        raise AssertionError("the signature has the wrong spine child size")
    if signature.attachment_child_size != _component_size(
        environment.cluster_mask[attachment]
    ):
        raise AssertionError("the signature has the wrong attachment child size")
    return spine, attachment


def _orientations_for_spine(
    environment: HierarchicalSweepEnvironment,
    node: int,
    spine_child: int,
) -> tuple[Orientation, ...]:
    left = environment.forest.left_child[node]
    right = environment.forest.right_child[node]
    if spine_child == left:
        return LEFT_RIGHT_ORIENTATIONS
    if spine_child == right:
        return tuple(
            (right_orientation, left_orientation)
            for left_orientation, right_orientation in LEFT_RIGHT_ORIENTATIONS
        )
    raise ValueError("spine_child must be one of the merger children")


def _repair_roundoff_domination(mass: float, twisted: float) -> float:
    if abs(twisted) <= mass:
        return twisted
    if abs(twisted) - mass <= 5e-15:
        return copysign(mass, twisted)
    raise AssertionError("coarse-graining violated |U| <= K")


def _deficit_summary(
    mass_row: Sequence[float], twisted_row: Sequence[float]
) -> LocalDeficitSummary:
    transfer = dominated_transfer((tuple(mass_row),), (tuple(twisted_row),))
    deficits = []
    for mass, twisted in zip(mass_row, twisted_row, strict=True):
        if mass == 0.0 or twisted == 0.0:
            deficits.append(float("inf"))
        else:
            deficits.append(-log(abs(twisted) / mass))
    envelope = fsum(abs(value) for value in twisted_row)
    if envelope > fsum(mass_row) + 5e-15:
        raise AssertionError("the Feynman--Kac envelope exceeds total mass")
    attenuation = float("inf") if envelope == 0.0 else -log(envelope)
    terminal = (1.0,) * len(mass_row)
    certificate = finite_horizon_doob_certificate((transfer,), terminal, terminal)
    certificate_envelope = float(certificate.feynman_kac_envelope[0])
    return LocalDeficitSummary(
        transfer=transfer,
        transition_deficits=tuple(deficits),
        feynman_kac_envelope=envelope,
        logarithmic_attenuation=attenuation,
        finite_horizon_doob_envelope_error=abs(certificate_envelope - envelope),
        finite_horizon_doob_inequality_holds=certificate.inequality_holds,
    )


def build_real_rank_t2_deficit_prototype(
    environment: HierarchicalSweepEnvironment,
    pair: EventPalmCorridor,
    signature: CorridorT2Signature,
    *,
    first_reference_state: int = 0,
    second_reference_state: int = 0,
) -> RealRankT2DeficitPrototype:
    """Build faithful and projected replicated transfers for one attachment.

    The two heat-bath draws are conditionally independent, which is the exact
    replicated construction for a second moment.  Their environment,
    geometry, physical edge marks, merger ranks, and ancestor factors are
    shared.
    """

    spine_child, attachment_child = _spine_attachment_geometry(
        environment, pair, signature
    )
    node = signature.node
    factor_chain = environment.factor_chains[node]
    rank = environment.forest.merge_rank[node]
    beta = clock_time_from_rank(rank, environment.p)
    geometry = AttachmentCellGeometry(
        side_length=environment.side_length,
        p=environment.p,
        pair_first=pair.first,
        pair_second=pair.second,
        pair_lca=pair.lca_node,
        pair_palm_weight=pair.palm_weight,
        node=node,
        actual_rank=rank,
        beta=beta,
        bucket_size=environment.forest.bucket_size[node],
        spine_child=spine_child,
        attachment_child=attachment_child,
        spine_child_size=_component_size(environment.cluster_mask[spine_child]),
        attachment_child_size=_component_size(
            environment.cluster_mask[attachment_child]
        ),
        strict_ancestor_nodes=tuple(factor_chain[1:]),
        is_postcritical=rank >= Q_CRITICAL,
    )

    first_boundary = _replica_exterior_state(environment, node, first_reference_state)
    second_boundary = _replica_exterior_state(environment, node, second_reference_state)
    orientations = _orientations_for_spine(environment, node, spine_child)

    faithful_targets = []
    faithful_mass = []
    faithful_twisted = []
    projected_mass = {state: [] for state in REPLICATED_RELATIVE_STATES}
    projected_twisted = {state: [] for state in REPLICATED_RELATIVE_STATES}
    for first_index, second_index in product(range(4), repeat=2):
        first_spine, first_attachment = orientations[first_index]
        second_spine, second_attachment = orientations[second_index]
        first_relative = first_spine * first_attachment
        second_relative = second_spine * second_attachment
        epsilon = first_spine * second_spine
        probability = (
            first_boundary.heat_bath_probabilities[first_index]
            * second_boundary.heat_bath_probabilities[second_index]
        )
        target = FaithfulReplicatedTarget(
            first_orbit_index=first_index,
            second_orbit_index=second_index,
            first_output_state=(
                first_reference_state ^ first_boundary.orbit_masks[first_index]
            ),
            second_output_state=(
                second_reference_state ^ second_boundary.orbit_masks[second_index]
            ),
            first_spine_orientation=first_spine,
            second_spine_orientation=second_spine,
            first_relative_orientation=first_relative,
            second_relative_orientation=second_relative,
            replicated_twist=epsilon,
        )
        faithful_targets.append(target)
        faithful_mass.append(probability)
        faithful_twisted.append(epsilon * probability)
        relative_state = (first_relative, second_relative)
        projected_mass[relative_state].append(probability)
        projected_twisted[relative_state].append(epsilon * probability)

    total_mass = fsum(faithful_mass)
    if abs(total_mass - 1.0) > 5e-14:
        raise AssertionError("the replicated heat bath does not have unit mass")
    faithful_twisted = [
        _repair_roundoff_domination(mass, twisted)
        for mass, twisted in zip(faithful_mass, faithful_twisted, strict=True)
    ]
    faithful_summary = _deficit_summary(faithful_mass, faithful_twisted)

    coarse_mass = tuple(
        fsum(projected_mass[state]) for state in REPLICATED_RELATIVE_STATES
    )
    coarse_twisted = tuple(
        _repair_roundoff_domination(
            mass,
            fsum(projected_twisted[state]),
        )
        for state, mass in zip(REPLICATED_RELATIVE_STATES, coarse_mass, strict=True)
    )
    projected_summary = _deficit_summary(coarse_mass, coarse_twisted)

    return RealRankT2DeficitPrototype(
        geometry=geometry,
        replica_exterior_states=(first_boundary, second_boundary),
        faithful_targets=tuple(faithful_targets),
        faithful=faithful_summary,
        projected_relative_states=REPLICATED_RELATIVE_STATES,
        projected=projected_summary,
        shared_nishimori_environment=True,
        replica_draws_conditionally_independent=True,
        full_boundary_retained_by_faithful_transfer=True,
        projected_boundary_is_markov_closed=False,
        composable_t2_deficit_certified=False,
        palm_abundance_estimated=False,
        weak_recovery_claimed=False,
        interpretation=(
            "The faithful one-update transfer has zero deficit because the "
            "complete output boundary reveals the twist.  The projected "
            "attenuation is an exact one-step coarse-graining but is not "
            "composable: it drops global orientations needed to update the "
            "strict-ancestor projective potential.  A nontrivial certified "
            "T2 deficit requires either a multi-update block that eliminates "
            "those orientations after their last use, or an exact port-gauge "
            "transition carrying the full exterior potential."
        ),
    )


def find_sample_prototype(
    *,
    side_length: int = 8,
    p: float = 0.805,
    distance_fraction: float = 0.25,
    maximum_bucket_size: int = 8,
    maximum_attachment_size: int = 4,
    minimum_rank: float = Q_CRITICAL,
    seed: int = 20260724,
    maximum_environments: int = 100,
) -> RealRankT2DeficitPrototype:
    """Return the first deterministic small-attachment example found.

    This is a search utility, not a Palm estimator: stopping at the first
    candidate changes the sampling law.  ``pair_palm_weight`` is retained
    only to audit that the candidate came from a genuine event-Palm corridor.
    By default the search is restricted to genuinely postcritical ranks.
    """

    if side_length < 4:
        raise ValueError("side_length must be at least four")
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2,1)")
    if maximum_bucket_size < 2 or maximum_attachment_size <= 0:
        raise ValueError("the bucket and attachment cutoffs must be positive")
    if not 0.0 <= minimum_rank <= 2.0 * p - 1.0:
        raise ValueError("minimum_rank must lie in the realized rank window")
    if maximum_environments <= 0:
        raise ValueError("maximum_environments must be positive")

    rng = random.Random(seed)
    for repetition in range(maximum_environments):
        ranked_edges: tuple[RankedEdge, ...] = sample_ranked_edges(side_length, rng)
        environment = HierarchicalSweepEnvironment(side_length, ranked_edges, p)
        pairs, audit = event_palm_corridor_observations(
            environment.forest,
            ranked_edges,
            repetition=repetition,
            p=p,
            distance_fraction=distance_fraction,
            rng=rng,
        )
        if not audit.passed:
            raise AssertionError("the Kruskal reconstruction audit failed")
        for pair in pairs:
            for signature in pair.signatures:
                if (
                    not signature.is_lca
                    and 2 <= signature.bucket_size <= maximum_bucket_size
                    and signature.attachment_child_size is not None
                    and signature.attachment_child_size <= maximum_attachment_size
                    and signature.rank >= minimum_rank
                ):
                    return build_real_rank_t2_deficit_prototype(
                        environment,
                        pair,
                        signature,
                    )
    raise RuntimeError("no small-attachment corridor cell was found")


def _json_number(value: float) -> float | str:
    if isfinite(value):
        return value
    return "+infinity" if value > 0 else "-infinity"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=8)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-bucket-size", type=int, default=8)
    parser.add_argument("--maximum-attachment-size", type=int, default=4)
    parser.add_argument("--minimum-rank", type=float, default=Q_CRITICAL)
    parser.add_argument("--seed", type=int, default=20260724)
    arguments = parser.parse_args()
    result = find_sample_prototype(
        side_length=arguments.side,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        maximum_bucket_size=arguments.maximum_bucket_size,
        maximum_attachment_size=arguments.maximum_attachment_size,
        minimum_rank=arguments.minimum_rank,
        seed=arguments.seed,
    )
    payload = {
        "geometry": asdict(result.geometry),
        "external_projective_log_weights": [
            [_json_number(value) for value in state.exterior_projective_log_weights]
            for state in result.replica_exterior_states
        ],
        "faithful_feynman_kac_envelope": (result.faithful.feynman_kac_envelope),
        "faithful_logarithmic_attenuation": _json_number(
            result.faithful.logarithmic_attenuation
        ),
        "projected_feynman_kac_envelope": (result.projected.feynman_kac_envelope),
        "projected_logarithmic_attenuation": _json_number(
            result.projected.logarithmic_attenuation
        ),
        "projected_transition_deficits": [
            _json_number(value) for value in result.projected.transition_deficits
        ],
        "projected_boundary_is_markov_closed": (
            result.projected_boundary_is_markov_closed
        ),
        "composable_t2_deficit_certified": (result.composable_t2_deficit_certified),
        "interpretation": result.interpretation,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
