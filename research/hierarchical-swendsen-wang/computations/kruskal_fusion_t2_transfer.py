"""Exact finite T2-Kruskal transfer and a criticalization counterexample.

This module isolates the smallest quotient cell that simultaneously contains

* an exact merger factor for two target children ``C1,C2`` in an abstract
  quotient cell;
* an exact ancestral merger factor with one lateral attachment;
* a latent winning edge which is marginalized, never revealed;
* the exact factor ``Lambda exp((1-beta)Lambda)`` at both mergers;
* a parametric outside potential ``B`` and attachment coupling ``J``;
* two posterior replicas using the same observed environment.

The latent state is ``x=(x1,x2)``.  At the target merger all physical cut
edges carry the relation ``chi(x)=x1*x2``.  At the ancestral merger, edges
incident to the two target children carry the two distinct relations ``x1``
and ``x2`` to a fixed lateral port.  The latter is precisely the multiport
feature that a scalar bucket count loses.

For a homogeneous bucket observed at time ``beta``, the dendrogram reveals
the merger time but not the winning edge.  Conditional on a latent state,
one latent winner is uniform and satisfied; all other residual satisfactions
are independent with probability ``s_p(beta)``.  Summing the winner gives a
channel whose posterior is exactly

    prior(x) * product_v Lambda_v(x) exp((1-beta_v)Lambda_v(x)).

Two independent calculations are implemented.  The ``edgewise`` route
enumerates every observed sign and every possible latent winner.  The
``grouped`` route eliminates edge identities analytically with binomial
sums and retains only satisfied counts in the relation groups.  Their mass
transfers, shared-replica laws, and chi-square reliabilities are compared in
the tests.

The principal falsification result is finite and exact up to floating-point
evaluation of elementary functions.  At p=0.805, with target time beta_c,
outside field B=4, attachment coupling J=3, and a two-edge ancestral bucket
carrying relations x1,x2, moving the ancestral time from beta_c to 0.8
*increases* the normalized posterior second moment of chi=x1*x2.  Therefore
criticalization is not a Blackwell domination on the complete multiport
state for arbitrary polarized priors.  This does not rule out an annealed
comparison under the actual corridor-message law.

Scope: one fixed finite quotient cell.  There is no Palm law, abundance
estimate, thermodynamic limit, or weak-recovery bound in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import comb, exp, fsum, isfinite, log, sqrt
from typing import Sequence

from critical_band_thresholds import beta_critical

State = tuple[int, int]
EdgewiseEnvironment = tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]
GroupedEnvironment = tuple[tuple[int, ...], tuple[int, ...]]
Matrix = tuple[tuple[float, ...], ...]

STATES: tuple[State, ...] = tuple(product((-1, 1), repeat=2))
REPLICATED_STATES: tuple[tuple[State, State], ...] = tuple(product(STATES, repeat=2))


def coupling(p: float) -> float:
    """Return the homogeneous Nishimori rate log(p/(1-p))."""

    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2,1)")
    return log(p / (1.0 - p))


def closed_satisfaction_probability(p: float, beta: float) -> float:
    """Return P(edge satisfied | its clock is still closed at beta)."""

    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must belong to [0,1]")
    rate = coupling(p)
    numerator = p * exp(-rate * beta)
    return numerator / (1.0 - p + numerator)


@dataclass(frozen=True)
class T2Geometry:
    """The target bucket followed by one multiport ancestral bucket."""

    target_size: int
    attachment_groups: tuple[int, int]
    beta_target: float
    beta_attachment: float
    label: str = "custom"

    def __post_init__(self) -> None:
        if self.target_size < 2:
            raise ValueError("target_size must be at least two")
        if len(self.attachment_groups) != 2 or min(self.attachment_groups) <= 0:
            raise ValueError("attachment_groups must contain two positive sizes")
        if not (0.0 <= self.beta_target <= self.beta_attachment <= 1.0):
            raise ValueError("merger times must be ordered in [0,1]")

    @property
    def bucket_groups(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return (self.target_size,), self.attachment_groups


def p805_geometry(
    attachment_groups: tuple[int, int] = (1, 1),
    *,
    attachment_beta: float | None = None,
) -> T2Geometry:
    """Return the p=161/200 cell with its target merger at beta_c."""

    critical = beta_critical(161.0 / 200.0)
    beta = critical if attachment_beta is None else float(attachment_beta)
    suffix = "critical" if beta == critical else f"beta-{beta:g}"
    return T2Geometry(
        target_size=2,
        attachment_groups=attachment_groups,
        beta_target=critical,
        beta_attachment=beta,
        label=f"T2-{attachment_groups[0]}x{attachment_groups[1]}-{suffix}",
    )


def prior_distribution(
    outside_field: float = 0.0, attachment_coupling: float = 0.0
) -> tuple[float, ...]:
    """Return pi(x) proportional to exp(B(x1+x2)+J*x1*x2)."""

    if not isfinite(outside_field) or not isfinite(attachment_coupling):
        raise ValueError("outside potentials must be finite")
    weights = tuple(
        exp(
            outside_field * (state[0] + state[1])
            + attachment_coupling * state[0] * state[1]
        )
        for state in STATES
    )
    normalizer = fsum(weights)
    return tuple(weight / normalizer for weight in weights)


def edgewise_merger_channel_probability(
    observed_groups: Sequence[Sequence[int]],
    relations: Sequence[int],
    p: float,
    beta: float,
) -> float:
    """Sum the latent winning-edge identity for one observed bucket.

    ``relations[g]`` is the candidate true relation shared by group ``g``.
    The observed signs are physical signs, not satisfaction indicators.
    """

    groups = tuple(tuple(group) for group in observed_groups)
    candidate_relations = tuple(relations)
    if not groups or len(groups) != len(candidate_relations):
        raise ValueError("one relation is required for each nonempty group")
    if any(not group for group in groups):
        raise ValueError("edgewise groups must be nonempty")
    if any(sign not in (-1, 1) for group in groups for sign in group):
        raise ValueError("observed signs must belong to {-1,+1}")
    if any(relation not in (-1, 1) for relation in candidate_relations):
        raise ValueError("relations must belong to {-1,+1}")

    residual = tuple(
        sign * candidate_relations[group_index]
        for group_index, group in enumerate(groups)
        for sign in group
    )
    satisfaction = closed_satisfaction_probability(p, beta)
    total = len(residual)
    probability = 0.0
    for winner, winner_satisfaction in enumerate(residual):
        if winner_satisfaction != 1:
            continue
        conditional = 1.0 / total
        for index, value in enumerate(residual):
            if index == winner:
                continue
            conditional *= satisfaction if value == 1 else 1.0 - satisfaction
        probability += conditional
    return probability


def _binomial_mass(size: int, count: int, probability: float) -> float:
    if not 0 <= count <= size:
        return 0.0
    return (
        comb(size, count) * probability**count * (1.0 - probability) ** (size - count)
    )


def grouped_merger_count_distribution(
    group_sizes: Sequence[int], p: float, beta: float
) -> dict[tuple[int, ...], float]:
    """Eliminate winner and edge identities by grouped binomial sums."""

    sizes = tuple(int(size) for size in group_sizes)
    if not sizes or min(sizes) <= 0:
        raise ValueError("group sizes must be positive")
    satisfaction = closed_satisfaction_probability(p, beta)
    total = sum(sizes)
    result: dict[tuple[int, ...], float] = {}
    for counts in product(*(range(size + 1) for size in sizes)):
        mass = 0.0
        for winner_group, winner_group_size in enumerate(sizes):
            conditional = winner_group_size / total
            for group, (size, count) in enumerate(zip(sizes, counts, strict=True)):
                if group == winner_group:
                    conditional *= _binomial_mass(size - 1, count - 1, satisfaction)
                else:
                    conditional *= _binomial_mass(size, count, satisfaction)
            mass += conditional
        if mass > 0.0:
            result[counts] = mass
    return result


def _all_edgewise_environments(
    geometry: T2Geometry,
) -> tuple[EdgewiseEnvironment, ...]:
    target_signs = tuple(product((-1, 1), repeat=geometry.target_size))
    first_size, second_size = geometry.attachment_groups
    first_signs = tuple(product((-1, 1), repeat=first_size))
    second_signs = tuple(product((-1, 1), repeat=second_size))
    return tuple(
        (target, (first, second))
        for target, first, second in product(target_signs, first_signs, second_signs)
    )


def _all_grouped_environments(
    geometry: T2Geometry,
) -> tuple[GroupedEnvironment, ...]:
    first_size, second_size = geometry.attachment_groups
    return tuple(
        ((target_count,), (first_count, second_count))
        for target_count, first_count, second_count in product(
            range(geometry.target_size + 1),
            range(first_size + 1),
            range(second_size + 1),
        )
    )


def edgewise_cell_channel_probability(
    environment: EdgewiseEnvironment,
    state: State,
    geometry: T2Geometry,
    p: float,
) -> float:
    """Return P(environment | state, fixed unmarked merger times)."""

    if state not in STATES:
        raise ValueError("state must belong to {-1,+1}^2")
    if len(environment) != 2:
        raise ValueError("a T2 environment contains two buckets")
    child_character = state[0] * state[1]
    target = edgewise_merger_channel_probability(
        (environment[0],),
        (child_character,),
        p,
        geometry.beta_target,
    )
    attachment = edgewise_merger_channel_probability(
        environment[1], state, p, geometry.beta_attachment
    )
    return target * attachment


def _observed_count_from_residual(size: int, residual_count: int, relation: int) -> int:
    return residual_count if relation == 1 else size - residual_count


def grouped_cell_channel_probability(
    environment: GroupedEnvironment,
    state: State,
    geometry: T2Geometry,
    p: float,
) -> float:
    """Return the same channel using independently grouped elimination."""

    if state not in STATES:
        raise ValueError("state must belong to {-1,+1}^2")
    if len(environment) != 2 or len(environment[0]) != 1 or len(environment[1]) != 2:
        raise ValueError("grouped environment has the wrong shape")
    target_count = environment[0][0]
    attachment_counts = environment[1]
    if not 0 <= target_count <= geometry.target_size or any(
        not 0 <= count <= size
        for count, size in zip(
            attachment_counts, geometry.attachment_groups, strict=True
        )
    ):
        raise ValueError("grouped counts do not match the geometry")

    target_law = grouped_merger_count_distribution(
        (geometry.target_size,), p, geometry.beta_target
    )
    attachment_law = grouped_merger_count_distribution(
        geometry.attachment_groups, p, geometry.beta_attachment
    )
    character = state[0] * state[1]
    target_residual = _observed_count_from_residual(
        geometry.target_size, target_count, character
    )
    attachment_residual = tuple(
        _observed_count_from_residual(size, count, relation)
        for size, count, relation in zip(
            geometry.attachment_groups,
            attachment_counts,
            state,
            strict=True,
        )
    )
    return target_law.get((target_residual,), 0.0) * attachment_law.get(
        attachment_residual, 0.0
    )


def _log_fusion_factor(rate: float, beta: float) -> float:
    if rate < 0.0:
        raise ValueError("a merger rate cannot be negative")
    if rate == 0.0:
        return float("-inf")
    return log(rate) + (1.0 - beta) * rate


def _normalize_log_weights(values: Sequence[float]) -> tuple[float, ...]:
    maximum = max(values)
    if maximum == float("-inf"):
        raise ValueError("the observed environment is impossible")
    unnormalized = tuple(
        0.0 if value == float("-inf") else exp(value - maximum) for value in values
    )
    normalizer = fsum(unnormalized)
    return tuple(value / normalizer for value in unnormalized)


def _prior_log_weight(
    state: State, outside_field: float, attachment_coupling: float
) -> float:
    return (
        outside_field * (state[0] + state[1])
        + attachment_coupling * state[0] * state[1]
    )


def edgewise_heat_bath(
    environment: EdgewiseEnvironment,
    geometry: T2Geometry,
    p: float,
    outside_field: float = 0.0,
    attachment_coupling: float = 0.0,
) -> tuple[float, ...]:
    """Exact posterior heat bath by direct edge-sign enumeration.

    This route evaluates the two ``Lambda exp`` factors directly; it does not
    call the channel-probability formula used by the Bayes counter-audit.
    """

    if not isfinite(outside_field) or not isfinite(attachment_coupling):
        raise ValueError("outside potentials must be finite")
    if (
        len(environment) != 2
        or len(environment[0]) != geometry.target_size
        or len(environment[1]) != 2
        or tuple(map(len, environment[1])) != geometry.attachment_groups
    ):
        raise ValueError("edgewise environment does not match the geometry")
    if any(
        sign not in (-1, 1)
        for bucket in (environment[0], *environment[1])
        for sign in bucket
    ):
        raise ValueError("observed signs must belong to {-1,+1}")

    rate = coupling(p)
    log_weights = []
    for state in STATES:
        character = state[0] * state[1]
        target_rate = rate * sum(sign * character == 1 for sign in environment[0])
        attachment_rate = rate * fsum(
            sum(sign * relation == 1 for sign in group)
            for group, relation in zip(environment[1], state, strict=True)
        )
        log_weights.append(
            _prior_log_weight(state, outside_field, attachment_coupling)
            + _log_fusion_factor(target_rate, geometry.beta_target)
            + _log_fusion_factor(attachment_rate, geometry.beta_attachment)
        )
    return _normalize_log_weights(log_weights)


def grouped_heat_bath(
    environment: GroupedEnvironment,
    geometry: T2Geometry,
    p: float,
    outside_field: float = 0.0,
    attachment_coupling: float = 0.0,
) -> tuple[float, ...]:
    """Same posterior after analytic elimination of all edge identities."""

    if not isfinite(outside_field) or not isfinite(attachment_coupling):
        raise ValueError("outside potentials must be finite")
    if len(environment) != 2 or len(environment[0]) != 1 or len(environment[1]) != 2:
        raise ValueError("grouped environment has the wrong shape")
    target_count = environment[0][0]
    attachment_counts = environment[1]
    if not 0 <= target_count <= geometry.target_size or any(
        not 0 <= count <= size
        for count, size in zip(
            attachment_counts, geometry.attachment_groups, strict=True
        )
    ):
        raise ValueError("grouped counts do not match the geometry")

    rate = coupling(p)
    log_weights = []
    for state in STATES:
        character = state[0] * state[1]
        target_satisfied = (
            target_count if character == 1 else geometry.target_size - target_count
        )
        attachment_satisfied = fsum(
            count if relation == 1 else size - count
            for count, size, relation in zip(
                attachment_counts,
                geometry.attachment_groups,
                state,
                strict=True,
            )
        )
        log_weights.append(
            _prior_log_weight(state, outside_field, attachment_coupling)
            + _log_fusion_factor(rate * target_satisfied, geometry.beta_target)
            + _log_fusion_factor(rate * attachment_satisfied, geometry.beta_attachment)
        )
    return _normalize_log_weights(log_weights)


def environment_counts(
    environment: EdgewiseEnvironment,
) -> GroupedEnvironment:
    """Map observed edge signs to plus counts in each relation group."""

    return (
        (sum(sign == 1 for sign in environment[0]),),
        tuple(sum(sign == 1 for sign in group) for group in environment[1]),
    )


def _observable_statistics(
    prior: Sequence[float], observable: Sequence[float]
) -> tuple[float, float, tuple[float, ...]]:
    mean = fsum(mass * value for mass, value in zip(prior, observable, strict=True))
    variance = fsum(
        mass * (value - mean) ** 2
        for mass, value in zip(prior, observable, strict=True)
    )
    if variance <= 0.0:
        raise ValueError("the observable has zero prior variance")
    normalized = tuple((value - mean) / sqrt(variance) for value in observable)
    return mean, variance, normalized


@dataclass(frozen=True)
class T2Transfer:
    """Mass transfer, shared-replica law, and the chi-square sector."""

    p: float
    geometry: T2Geometry
    outside_field: float
    attachment_coupling: float
    prior: tuple[float, ...]
    mass: Matrix
    shared_replicated_law: Matrix
    full_replicated_heat_bath: Matrix
    raw_chi_second_moment: float
    chi_square_reliability: float
    prior_chi_mean: float
    mass_row_sum_error: float
    detailed_balance_error: float
    shared_law_mass_error: float
    route: str
    scope_label: str = "finite T2-Kruskal cell; no Palm or weak-recovery claim"


def _build_transfer(
    geometry: T2Geometry,
    p: float,
    outside_field: float,
    attachment_coupling: float,
    route: str,
) -> T2Transfer:
    prior = prior_distribution(outside_field, attachment_coupling)
    chi = tuple(float(state[0] * state[1]) for state in STATES)
    prior_mean, _, normalized_chi = _observable_statistics(prior, chi)

    if route == "edgewise":
        environments: Sequence[EdgewiseEnvironment | GroupedEnvironment] = (
            _all_edgewise_environments(geometry)
        )
        channel = edgewise_cell_channel_probability
        heat_bath = edgewise_heat_bath
    elif route == "grouped":
        environments = _all_grouped_environments(geometry)
        channel = grouped_cell_channel_probability
        heat_bath = grouped_heat_bath
    else:
        raise ValueError("route must be 'edgewise' or 'grouped'")

    mass = [[0.0 for _ in STATES] for _ in STATES]
    shared = [[0.0 for _ in STATES] for _ in STATES]
    raw_second_moment = 0.0
    chi_square = 0.0
    output_mass_sum = 0.0

    for environment in environments:
        channel_row = tuple(
            channel(environment, state, geometry, p) for state in STATES
        )
        output_mass = fsum(
            prior_mass * likelihood
            for prior_mass, likelihood in zip(prior, channel_row, strict=True)
        )
        if output_mass == 0.0:
            continue
        posterior = heat_bath(
            environment,
            geometry,
            p,
            outside_field,
            attachment_coupling,
        )
        bayes_posterior = tuple(
            prior_mass * likelihood / output_mass
            for prior_mass, likelihood in zip(prior, channel_row, strict=True)
        )
        if (
            max(
                abs(first - second)
                for first, second in zip(posterior, bayes_posterior, strict=True)
            )
            > 3e-12
        ):
            raise AssertionError("Lambda-factor heat bath disagrees with Bayes channel")

        output_mass_sum += output_mass
        for latent_index, likelihood in enumerate(channel_row):
            for output_index, posterior_mass in enumerate(posterior):
                mass[latent_index][output_index] += likelihood * posterior_mass
        for first in range(4):
            for second in range(4):
                shared[first][second] += (
                    output_mass * posterior[first] * posterior[second]
                )
        posterior_chi = fsum(
            posterior_mass * value
            for posterior_mass, value in zip(posterior, chi, strict=True)
        )
        posterior_normalized_chi = fsum(
            posterior_mass * value
            for posterior_mass, value in zip(posterior, normalized_chi, strict=True)
        )
        raw_second_moment += output_mass * posterior_chi**2
        chi_square += output_mass * posterior_normalized_chi**2

    if abs(output_mass_sum - 1.0) > 3e-12:
        raise AssertionError("the finite observation channel is not normalized")
    mass_matrix = tuple(tuple(row) for row in mass)
    shared_matrix = tuple(tuple(row) for row in shared)
    flattened_shared = tuple(
        shared_matrix[STATES.index(first)][STATES.index(second)]
        for first, second in REPLICATED_STATES
    )
    full_replicated = tuple(flattened_shared for _ in REPLICATED_STATES)
    detailed_balance_error = max(
        abs(
            prior[first] * mass_matrix[first][second]
            - prior[second] * mass_matrix[second][first]
        )
        for first in range(4)
        for second in range(4)
    )
    return T2Transfer(
        p=p,
        geometry=geometry,
        outside_field=outside_field,
        attachment_coupling=attachment_coupling,
        prior=prior,
        mass=mass_matrix,
        shared_replicated_law=shared_matrix,
        full_replicated_heat_bath=full_replicated,
        raw_chi_second_moment=raw_second_moment,
        chi_square_reliability=chi_square,
        prior_chi_mean=prior_mean,
        mass_row_sum_error=max(abs(fsum(row) - 1.0) for row in mass_matrix),
        detailed_balance_error=detailed_balance_error,
        shared_law_mass_error=abs(
            fsum(value for row in shared_matrix for value in row) - 1.0
        ),
        route=route,
    )


@lru_cache(maxsize=None)
def edgewise_t2_transfer(
    geometry: T2Geometry,
    p: float = 0.805,
    outside_field: float = 0.0,
    attachment_coupling: float = 0.0,
) -> T2Transfer:
    """Build the transfer by full edge-sign/winner enumeration."""

    return _build_transfer(
        geometry,
        p,
        outside_field,
        attachment_coupling,
        "edgewise",
    )


@lru_cache(maxsize=None)
def grouped_t2_transfer(
    geometry: T2Geometry,
    p: float = 0.805,
    outside_field: float = 0.0,
    attachment_coupling: float = 0.0,
) -> T2Transfer:
    """Build the transfer by independent grouped-count elimination."""

    return _build_transfer(
        geometry,
        p,
        outside_field,
        attachment_coupling,
        "grouped",
    )


@dataclass(frozen=True)
class CriticalLateComparison:
    outside_field: float
    attachment_coupling: float
    attachment_groups: tuple[int, int]
    critical_beta: float
    late_beta: float
    critical_reliability: float
    late_reliability: float
    late_minus_critical: float

    @property
    def refutes_uniform_critical_domination(self) -> bool:
        return self.late_minus_critical > 0.0


def compare_critical_and_late(
    p: float = 0.805,
    outside_field: float = 4.0,
    attachment_coupling: float = 3.0,
    late_beta: float = 0.8,
    attachment_groups: tuple[int, int] = (1, 1),
    *,
    route: str = "grouped",
) -> CriticalLateComparison:
    """Compare chi-square reliability without assuming Blackwell order."""

    critical = beta_critical(p)
    if not critical <= late_beta <= 1.0:
        raise ValueError("late_beta must belong to [beta_c,1]")
    critical_geometry = T2Geometry(
        2,
        attachment_groups,
        critical,
        critical,
        "criticalized multiport ancestor",
    )
    late_geometry = T2Geometry(
        2,
        attachment_groups,
        critical,
        late_beta,
        "late multiport ancestor",
    )
    builder = grouped_t2_transfer if route == "grouped" else edgewise_t2_transfer
    if route not in ("grouped", "edgewise"):
        raise ValueError("route must be 'edgewise' or 'grouped'")
    critical_value = builder(
        critical_geometry,
        p,
        outside_field,
        attachment_coupling,
    ).chi_square_reliability
    late_value = builder(
        late_geometry,
        p,
        outside_field,
        attachment_coupling,
    ).chi_square_reliability
    return CriticalLateComparison(
        outside_field=outside_field,
        attachment_coupling=attachment_coupling,
        attachment_groups=attachment_groups,
        critical_beta=critical,
        late_beta=late_beta,
        critical_reliability=critical_value,
        late_reliability=late_value,
        late_minus_critical=late_value - critical_value,
    )


@dataclass(frozen=True)
class P805ScanRow:
    outside_field: float
    attachment_coupling: float
    attachment_groups: tuple[int, int]
    critical_reliability: float
    late_reliability: float
    late_minus_critical: float


def p805_scan(
    fields: Sequence[float] = (0.0, 2.0, 4.0),
    attachment_couplings: Sequence[float] = (0.0, 1.0, 3.0),
    attachment_groups: Sequence[tuple[int, int]] = ((1, 1),),
    late_beta: float = 0.8,
) -> tuple[P805ScanRow, ...]:
    """Scan fields and attachment potentials in the fixed finite cell."""

    if not fields or not attachment_couplings or not attachment_groups:
        raise ValueError("scan parameter families must be nonempty")
    rows = []
    for groups, field, attachment in product(
        attachment_groups, fields, attachment_couplings
    ):
        comparison = compare_critical_and_late(
            p=161.0 / 200.0,
            outside_field=float(field),
            attachment_coupling=float(attachment),
            late_beta=late_beta,
            attachment_groups=groups,
        )
        rows.append(
            P805ScanRow(
                outside_field=float(field),
                attachment_coupling=float(attachment),
                attachment_groups=groups,
                critical_reliability=comparison.critical_reliability,
                late_reliability=comparison.late_reliability,
                late_minus_critical=comparison.late_minus_critical,
            )
        )
    return tuple(rows)


def main() -> None:
    print("p=0.805, target beta=beta_c, late ancestor beta=0.8")
    print("groups    B      J      critical       late          late-critical")
    for row in p805_scan():
        print(
            f"{str(row.attachment_groups):9s} "
            f"{row.outside_field:5.1f} {row.attachment_coupling:6.1f} "
            f"{row.critical_reliability:13.9f} "
            f"{row.late_reliability:13.9f} "
            f"{row.late_minus_critical:+14.9f}"
        )
    counterexample = compare_critical_and_late()
    print(
        "counterexample: late > critical = "
        f"{counterexample.refutes_uniform_critical_domination}"
    )
    print("scope: finite cell only; no Palm law and no weak-recovery bound")


if __name__ == "__main__":
    main()
