"""Exact collapsed certificate on a chain of triangular cactus blocks.

The graph is a chain of triangles.  Consecutive triangles share exactly one
articulation vertex, so the event that the two endpoints are connected at a
percolation rank ``q`` factorizes over the blocks.

For one block, the unmarked Kruskal hierarchy has two possible information
types for the endpoint relation:

* the direct articulation edge merges first and reveals the relation;
* a side edge merges first, after which the completion bucket has size two
  and is exactly an erasure experiment of reliability ``(p-r)/(1-r)`` at
  its completion rank ``r``.

The module derives the resulting replicated transfer in closed form and
counter-audits the size-two part by a direct enumeration of all spins and
all residual marks on chains of two or three triangles.  It is exact for the
stated cactus model.  It is not a theorem for the triangular lattice, where
overlapping cycles create a nontrivial boundary state.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import exp, fsum, log, tanh
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL, inverse_time


Matrix2 = tuple[tuple[float, float], tuple[float, float]]


def _validate_p_q(p: float, q: float) -> None:
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2, 1)")
    maximum = 2.0 * p - 1.0
    if not 0.0 < q <= maximum + 1e-15:
        raise ValueError("q must belong to (0, 2p-1]")


def residual_satisfaction_probability(p: float, rank: float) -> float:
    """Probability that a nonwinning edge is satisfied given rank > ``rank``."""

    _validate_p_q(p, rank)
    return (p - rank) / (1.0 - rank)


def triangle_connection_probability(rank: float) -> float:
    """Probability that the two articulation vertices connect by ``rank``."""

    if not 0.0 <= rank <= 1.0:
        raise ValueError("rank must belong to [0, 1]")
    return rank + rank * rank - rank**3


def triangle_connection_density(rank: float) -> float:
    """Density of the articulation-pair merger rank."""

    if not 0.0 <= rank < 1.0:
        raise ValueError("rank must belong to [0, 1)")
    return 1.0 + 2.0 * rank - 3.0 * rank * rank


def direct_first_connection_mass(rank: float) -> float:
    """Mass of connected histories whose direct edge is the first merger."""

    if not 0.0 <= rank <= 1.0:
        raise ValueError("rank must belong to [0, 1]")
    return rank - rank * rank + rank**3 / 3.0


def path_first_connection_mass(rank: float) -> float:
    """Mass of connected histories whose first merger is a side edge."""

    if not 0.0 <= rank <= 1.0:
        raise ValueError("rank must belong to [0, 1]")
    return 2.0 * rank * rank - 4.0 * rank**3 / 3.0


def direct_first_connection_density(rank: float) -> float:
    """Merger-rank density carried by direct-first histories."""

    if not 0.0 <= rank < 1.0:
        raise ValueError("rank must belong to [0, 1)")
    return (1.0 - rank) ** 2


def path_first_connection_density(rank: float) -> float:
    """Merger-rank density carried by path-first histories."""

    return path_first_completion_density(rank)


def path_first_completion_density(rank: float) -> float:
    """Unconditional density of a path-first completion at the given rank."""

    if not 0.0 <= rank <= 1.0:
        raise ValueError("rank must belong to [0, 1]")
    return 4.0 * rank * (1.0 - rank)


@dataclass(frozen=True)
class TriangleConnectionDecomposition:
    """Exact masses entering one connected cactus-block transfer."""

    connection_mass: float
    direct_first_mass: float
    path_first_mass: float
    informative_path_mass: float
    replicated_reliability: float

    @property
    def direct_first_probability(self) -> float:
        return self.direct_first_mass / self.connection_mass


def connected_triangle_decomposition(
    p: float, rank: float
) -> TriangleConnectionDecomposition:
    """Return the exact critical-connection decomposition of one triangle."""

    _validate_p_q(p, rank)
    connection = triangle_connection_probability(rank)
    direct = direct_first_connection_mass(rank)
    path = path_first_connection_mass(rank)
    informative_path = 2.0 * p * rank * rank - 4.0 * rank**3 / 3.0
    reliability = (direct + informative_path) / connection
    if abs(direct + path - connection) > 1e-12:
        raise AssertionError("the Kruskal history masses do not partition connection")
    return TriangleConnectionDecomposition(
        connection_mass=connection,
        direct_first_mass=direct,
        path_first_mass=path,
        informative_path_mass=informative_path,
        replicated_reliability=reliability,
    )


def connected_triangle_reliability(p: float, rank: float) -> float:
    """Closed-form replicated reliability conditional on local connection."""

    _validate_p_q(p, rank)
    return (
        1.0 + (2.0 * p - 1.0) * rank - rank * rank
    ) / (1.0 + rank - rank * rank)


def merger_flux_triangle_reliability(p: float, rank: float) -> float:
    """Reliability under the Palm density fixing the merger rank at ``rank``."""

    _validate_p_q(p, rank)
    return (
        1.0 + (4.0 * p - 2.0) * rank - 3.0 * rank * rank
    ) / (1.0 + 2.0 * rank - 3.0 * rank * rank)


def screened_connected_triangle_reliability(
    p: float, rank: float, message_bound: float
) -> float:
    """Worst local persistence when the erasure message satisfies ``|B|<=b``."""

    if message_bound < 0.0:
        raise ValueError("message_bound must be nonnegative")
    neutral = connected_triangle_reliability(p, rank)
    return neutral + (1.0 - neutral) * tanh(message_bound / 2.0) ** 2


def connected_triangle_reliability_by_quadrature(
    p: float, rank: float, subdivisions: int = 1000
) -> float:
    """Independently integrate the path-first histories by Simpson's rule."""

    _validate_p_q(p, rank)
    if subdivisions <= 0 or subdivisions % 2:
        raise ValueError("subdivisions must be a positive even integer")
    step = rank / subdivisions

    def integrand(value: float) -> float:
        return path_first_completion_density(value) * (
            (p - value) / (1.0 - value)
        )

    values = [integrand(index * step) for index in range(subdivisions + 1)]
    integral = step / 3.0 * (
        values[0]
        + values[-1]
        + 4.0 * fsum(values[1:-1:2])
        + 2.0 * fsum(values[2:-1:2])
    )
    return (
        direct_first_connection_mass(rank) + integral
    ) / triangle_connection_probability(rank)


def bsc_kernel(reliability: float) -> Matrix2:
    """Return the two-state convolution kernel with the given reliability."""

    if not 0.0 <= reliability <= 1.0:
        raise ValueError("reliability must belong to [0, 1]")
    same = (1.0 + reliability) / 2.0
    different = (1.0 - reliability) / 2.0
    return ((same, different), (different, same))


def multiply_two_by_two(first: Matrix2, second: Matrix2) -> Matrix2:
    """Multiply two 2-by-2 matrices without an external dependency."""

    return tuple(
        tuple(
            fsum(first[row][middle] * second[middle][column] for middle in (0, 1))
            for column in (0, 1)
        )
        for row in (0, 1)
    )  # type: ignore[return-value]


def connected_triangle_replicated_kernel(p: float, rank: float) -> Matrix2:
    """Replicated parity transfer of one connected cactus triangle."""

    return bsc_kernel(connected_triangle_reliability(p, rank))


def critical_to_late_replicated_garbling(
    p: float, early_rank: float, late_rank: float
) -> Matrix2:
    """Explicit BSC degrading an earlier connected block into a later one."""

    _validate_p_q(p, early_rank)
    _validate_p_q(p, late_rank)
    if early_rank > late_rank:
        raise ValueError("require early_rank <= late_rank")
    early = connected_triangle_reliability(p, early_rank)
    late = connected_triangle_reliability(p, late_rank)
    return bsc_kernel(late / early)


def critical_to_late_merger_flux_garbling(
    p: float, early_rank: float, late_rank: float
) -> Matrix2:
    """Explicit BSC degrading an earlier pivotal block into a later one."""

    _validate_p_q(p, early_rank)
    _validate_p_q(p, late_rank)
    if early_rank > late_rank:
        raise ValueError("require early_rank <= late_rank")
    early = merger_flux_triangle_reliability(p, early_rank)
    late = merger_flux_triangle_reliability(p, late_rank)
    return bsc_kernel(late / early)


def connected_cactus_second_moment(
    p: float, rank: float, block_count: int
) -> float:
    """Conditional endpoint second moment on a connected triangle chain."""

    if block_count < 0:
        raise ValueError("block_count must be nonnegative")
    coefficient = connected_triangle_reliability(p, rank)
    return coefficient**block_count


def connected_cactus_conformity_probability(
    p: float, rank: float, block_count: int
) -> float:
    """Probability that recolouring preserves the reference endpoint relation."""

    return 0.5 * (1.0 + connected_cactus_second_moment(p, rank, block_count))


def lca_rank_cactus_second_moment(
    p: float, rank: float, block_count: int
) -> float:
    """Endpoint second moment under the Palm density ``LCA rank = rank``.

    Exactly one triangle realizes the maximal merger rank.  Its transfer is
    the merger-flux coefficient; every other triangle is conditioned only to
    have connected its articulation pair before that rank.
    """

    if block_count <= 0:
        raise ValueError("block_count must be positive")
    return merger_flux_triangle_reliability(p, rank) * (
        connected_triangle_reliability(p, rank) ** (block_count - 1)
    )


def lca_only_cactus_second_moment(
    p: float, rank: float, block_count: int
) -> float:
    """Second moment when only the two children of the global LCA are updated.

    The unique pivotal triangle is resampled.  Every nonpivotal relation on
    the two descendant arms remains fixed and therefore contributes one.
    """

    if block_count <= 0:
        raise ValueError("block_count must be positive")
    return merger_flux_triangle_reliability(p, rank)


def lca_only_cactus_conformity_probability(
    p: float, rank: float, block_count: int
) -> float:
    """Reference-relation conformity after only the global LCA heat bath."""

    return 0.5 * (1.0 + lca_only_cactus_second_moment(p, rank, block_count))


def full_to_lca_cactus_second_moment_ratio(
    p: float, rank: float, block_count: int
) -> float:
    """Exact gain from resampling all descendant blocks instead of only LCA."""

    if block_count <= 0:
        raise ValueError("block_count must be positive")
    return connected_triangle_reliability(p, rank) ** (block_count - 1)


def lca_rank_cactus_conformity_probability(
    p: float, rank: float, block_count: int
) -> float:
    """Conformity probability under the Palm density fixing the LCA rank."""

    return 0.5 * (1.0 + lca_rank_cactus_second_moment(p, rank, block_count))


def unconditional_cactus_second_moment(
    p: float, rank: float, block_count: int
) -> float:
    """Include the exact probability that the two endpoints are connected."""

    if block_count < 0:
        raise ValueError("block_count must be nonnegative")
    decomposition = connected_triangle_decomposition(p, rank)
    informative_mass = (
        decomposition.connection_mass * decomposition.replicated_reliability
    )
    return informative_mass**block_count


def _hierarchical_factor(p: float, rank: float, satisfied: int) -> float:
    if satisfied == 0:
        return 0.0
    coupling = log(p / (1.0 - p))
    level = inverse_time(p, rank)
    rate = coupling * satisfied
    return rate * exp((1.0 - level) * rate)


def _fixed_mark_cases(
    p: float, completion_rank: float
) -> tuple[tuple[tuple[int, int], float], ...]:
    satisfaction = residual_satisfaction_probability(p, completion_rank)
    return (
        ((1, 1), satisfaction),
        ((1, -1), (1.0 - satisfaction) / 2.0),
        ((-1, 1), (1.0 - satisfaction) / 2.0),
    )


def _validate_completion_ranks(p: float, ranks: Sequence[float]) -> None:
    if not ranks:
        raise ValueError("at least one completion rank is required")
    previous = -1.0
    for rank in ranks:
        _validate_p_q(p, rank)
        if rank <= previous:
            raise ValueError("completion ranks must be strictly increasing")
        previous = rank


def fixed_path_first_cactus_direct_second_moment(
    p: float,
    completion_ranks: Sequence[float],
    satellite_ranks: Sequence[float] | None = None,
) -> float:
    """Enumerate every spin and residual mark for a fixed path-first hierarchy.

    Triangle ``r`` has articulation vertices ``a[r-1], a[r]`` and tip
    ``b[r]``.  The edge ``(b[r], a[r])`` merges first in a size-one bucket.
    The prefix then merges with that pair through the two-edge bucket
    ``(a[r-1], b[r]), (a[r-1], a[r])``.
    """

    _validate_completion_ranks(p, completion_ranks)
    block_count = len(completion_ranks)
    if satellite_ranks is None:
        satellite_ranks = tuple(rank / 2.0 for rank in completion_ranks)
    if len(satellite_ranks) != block_count:
        raise ValueError("one satellite rank is required per block")
    for satellite, completion in zip(
        satellite_ranks, completion_ranks, strict=True
    ):
        if not 0.0 < satellite < completion:
            raise ValueError("each satellite rank must precede its completion")

    vertex_count = 2 * block_count + 1
    mark_families = tuple(
        _fixed_mark_cases(p, rank) for rank in completion_ranks
    )
    answer = 0.0
    for marked_blocks in product(*mark_families):
        observations = tuple(case for case, _ in marked_blocks)
        environment_mass = 1.0
        for _, mass in marked_blocks:
            environment_mass *= mass

        denominator = 0.0
        numerator = 0.0
        for state in range(1 << (vertex_count - 1)):
            def spin(vertex: int) -> int:
                if vertex == 0:
                    return 1
                return -1 if state & (1 << (vertex - 1)) else 1

            weight = 1.0
            for index, (completion, satellite, marks) in enumerate(
                zip(
                    completion_ranks,
                    satellite_ranks,
                    observations,
                    strict=True,
                )
            ):
                left = index
                right = index + 1
                tip = block_count + 1 + index
                satellite_satisfied = int(spin(tip) * spin(right) == 1)
                weight *= _hierarchical_factor(
                    p, satellite, satellite_satisfied
                )
                if weight == 0.0:
                    break
                first_mark, second_mark = marks
                parent_satisfied = int(
                    first_mark * spin(left) * spin(tip) == 1
                ) + int(second_mark * spin(left) * spin(right) == 1)
                weight *= _hierarchical_factor(
                    p, completion, parent_satisfied
                )
                if weight == 0.0:
                    break
            denominator += weight
            numerator += weight * spin(block_count)
        if denominator <= 0.0:
            raise AssertionError("a conditioned hierarchy has zero mass")
        correlation = numerator / denominator
        answer += environment_mass * correlation * correlation
    return answer


def fixed_path_first_cactus_transfer_second_moment(
    p: float, completion_ranks: Sequence[float]
) -> float:
    """Evaluate the same fixed hierarchy as a product of local transfers."""

    _validate_completion_ranks(p, completion_ranks)
    result = 1.0
    for rank in completion_ranks:
        result *= residual_satisfaction_probability(p, rank)
    return result


def _certified_q_triangle_interval() -> tuple[Fraction, Fraction]:
    scale = 10**15
    lower_q = Fraction(347296355333860, scale)
    upper_q = Fraction(347296355333861, scale)
    polynomial = lambda value: value**3 - 3 * value + 1
    if not polynomial(lower_q) > 0 > polynomial(upper_q):
        raise AssertionError("the interval does not enclose q_triangle")
    return lower_q, upper_q


def p_eight_critical_cactus_interval(
    block_count: int,
) -> tuple[Fraction, Fraction]:
    """Rational enclosure under conditioning on critical connection."""

    if block_count < 0:
        raise ValueError("block_count must be nonnegative")
    lower_q, upper_q = _certified_q_triangle_interval()
    p = Fraction(4, 5)

    def coefficient(rank: Fraction) -> Fraction:
        return (
            1 + (2 * p - 1) * rank - rank * rank
        ) / (1 + rank - rank * rank)

    # The coefficient is strictly decreasing in q.
    return coefficient(upper_q) ** block_count, coefficient(lower_q) ** block_count


def p_eight_critical_lca_interval(
    block_count: int,
) -> tuple[Fraction, Fraction]:
    """Rational enclosure under the Palm density fixing the LCA at beta_c."""

    if block_count <= 0:
        raise ValueError("block_count must be positive")
    lower_q, upper_q = _certified_q_triangle_interval()
    p = Fraction(4, 5)

    def connected(rank: Fraction) -> Fraction:
        return (
            1 + (2 * p - 1) * rank - rank * rank
        ) / (1 + rank - rank * rank)

    def flux(rank: Fraction) -> Fraction:
        return (
            1 + (4 * p - 2) * rank - 3 * rank * rank
        ) / (1 + 2 * rank - 3 * rank * rank)

    lower = flux(upper_q) * connected(upper_q) ** (block_count - 1)
    upper = flux(lower_q) * connected(lower_q) ** (block_count - 1)
    return lower, upper


def main() -> None:
    p = 4.0 / 5.0
    rank = Q_CRITICAL
    decomposition = connected_triangle_decomposition(p, rank)
    print(
        f"p={p:.12g} q_critical={rank:.12f} "
        f"beta_critical={inverse_time(p, rank):.12f}"
    )
    print(
        "one block: connection="
        f"{decomposition.connection_mass:.12f} "
        f"direct-first|connected={decomposition.direct_first_probability:.12f} "
        f"connected_reliability={decomposition.replicated_reliability:.12f} "
        f"flux_reliability={merger_flux_triangle_reliability(p, rank):.12f}"
    )
    for count in (2, 3, 5, 10, 20, 40):
        lca_only = lca_only_cactus_second_moment(p, rank, count)
        second = connected_cactus_second_moment(p, rank, count)
        conformity = 0.5 * (1.0 + second)
        lca_second = lca_rank_cactus_second_moment(p, rank, count)
        lca_conformity = 0.5 * (1.0 + lca_second)
        print(
            f"blocks={count:2d} lca_only={lca_only:.12g} "
            f"full_over_lca={lca_second/lca_only:.12g} "
            f"second_moment={second:.12g} "
            f"conformity={conformity:.12g} "
            f"lca_second={lca_second:.12g} "
            f"lca_conformity={lca_conformity:.12g}"
        )

    strict_ranks = (rank - 0.002, rank - 0.001, rank)
    direct = fixed_path_first_cactus_direct_second_moment(p, strict_ranks)
    transfer = fixed_path_first_cactus_transfer_second_moment(p, strict_ranks)
    print(
        f"three path-first blocks: direct={direct:.12f} "
        f"transfer={transfer:.12f} gap={direct-transfer:.3g}"
    )


if __name__ == "__main__":
    main()
