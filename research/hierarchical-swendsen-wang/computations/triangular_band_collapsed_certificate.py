"""Replicated collapsed transfer for one width-two triangular-strip cell.

The cell has two left ports ``L0,L1``, two right ports ``R0,R1``, and the
four edges

    L0--R0,  L1--R0,  L1--R1,  R0--R1.

Thus ``L1,R0,R1`` form a triangle.  All four edges are conditioned to be
closed at a percolation rank ``q``.  In the planted gauge their residual
signs are independent, with satisfaction probability

    s = (p-q)/(1-q).

For a fixed residual-sign environment, the two right spins are refreshed by
their exact joint heat bath conditional on the two left spins.  Two replicas
use independent heat-bath randomness but the *same* sign environment.  The
resulting 16-state kernel is reduced under the two independent global-spin
symmetries.  Its trivial Fourier block is the four-state mass transfer, and
its ``chi tensor chi`` block propagates the annealed posterior second moment.

This is an exact ``E1+`` sector test for this fixed, all-closed, depth-one
cell.  It is not the ``E2/T2`` transfer and not a certificate for the critical
Palm corridor: the winning edge, open-skeleton partitions, ancestral Lambda
messages, pivotal status, and compatibility between adjacent Kruskal cells
are absent.  In particular, the strict neutral local contraction below must
not be quoted as a weak-recovery bound.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import exp, fsum, isfinite, log
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL


SpinPair = tuple[int, int]
Environment = tuple[int, int, int, int]
Matrix = tuple[tuple[float, ...], ...]
RationalInterval = tuple[Fraction, Fraction]
IntervalMatrix = tuple[tuple[RationalInterval, ...], ...]

BOUNDARY_SPINS: tuple[SpinPair, ...] = tuple(
    product((-1, 1), repeat=2)
)
ENVIRONMENTS: tuple[Environment, ...] = tuple(
    product((-1, 1), repeat=4)
)
REPLICATED_BOUNDARIES: tuple[tuple[SpinPair, SpinPair], ...] = tuple(
    product(BOUNDARY_SPINS, repeat=2)
)
PARITY_PAIRS: tuple[tuple[int, int], ...] = tuple(
    product((-1, 1), repeat=2)
)


def _validate_p_rank(p: float, rank: float) -> None:
    if not 0.5 < p < 1.0:
        raise ValueError("p must belong to (1/2, 1)")
    if not 0.0 <= rank < 2.0 * p - 1.0:
        raise ValueError("rank must belong to [0, 2p-1)")


def _validate_spin_pair(spins: SpinPair) -> None:
    if len(spins) != 2 or any(spin not in (-1, 1) for spin in spins):
        raise ValueError("a boundary state must contain two spins in {-1,+1}")


def _validate_environment(environment: Environment) -> None:
    if len(environment) != 4 or any(
        sign not in (-1, 1) for sign in environment
    ):
        raise ValueError("an environment must contain four signs in {-1,+1}")


def residual_satisfaction_probability(p: float, rank: float) -> float:
    """Return the edge satisfaction probability conditional on rank > q."""

    _validate_p_rank(p, rank)
    return (p - rank) / (1.0 - rank)


def cell_edge_relations(left: SpinPair, right: SpinPair) -> Environment:
    """Return candidate relations in the documented four-edge order."""

    _validate_spin_pair(left)
    _validate_spin_pair(right)
    return (
        left[0] * right[0],
        left[1] * right[0],
        left[1] * right[1],
        right[0] * right[1],
    )


def residual_environment_probability(
    environment: Environment, p: float, rank: float
) -> float:
    """Probability of a residual-sign environment in the planted gauge."""

    _validate_environment(environment)
    satisfaction = residual_satisfaction_probability(p, rank)
    return _sign_word_probability(environment, satisfaction)


def _sign_word_probability(
    signs: Sequence[int], satisfaction: float
) -> float:
    probability = 1.0
    for sign in signs:
        probability *= satisfaction if sign == 1 else 1.0 - satisfaction
    return probability


def cell_likelihood(
    left: SpinPair,
    right: SpinPair,
    environment: Environment,
    p: float,
    rank: float,
) -> float:
    """Residual BSC likelihood of one right boundary configuration."""

    _validate_environment(environment)
    satisfaction = residual_satisfaction_probability(p, rank)
    relations = cell_edge_relations(left, right)
    agreements = tuple(
        observed * candidate
        for observed, candidate in zip(
            environment, relations, strict=True
        )
    )
    return _sign_word_probability(agreements, satisfaction)


@lru_cache(maxsize=None)
def collapsed_cell_kernel(
    environment: Environment, p: float, rank: float
) -> Matrix:
    """Joint right-port heat bath, with rows and columns in BOUNDARY_SPINS."""

    _validate_environment(environment)
    _validate_p_rank(p, rank)
    rows = []
    for left in BOUNDARY_SPINS:
        weights = tuple(
            cell_likelihood(left, right, environment, p, rank)
            for right in BOUNDARY_SPINS
        )
        normalizer = fsum(weights)
        if normalizer <= 0.0:
            raise AssertionError("the residual cell heat bath has zero mass")
        rows.append(tuple(weight / normalizer for weight in weights))
    return tuple(rows)


@lru_cache(maxsize=None)
def polarized_collapsed_cell_kernel(
    environment: Environment,
    p: float,
    rank: float,
    boundary_field: float,
) -> Matrix:
    """Heat bath after adding a common field to both right ports.

    ``boundary_field`` is the log-odds field per right spin: a configuration
    ``right`` receives the extra factor

        exp(boundary_field * (right[0] + right[1]) / 2).

    A nonzero field breaks the global-spin symmetry, so this kernel cannot be
    inserted into the four-state Fourier block without retaining the field as
    part of the boundary state.
    """

    _validate_environment(environment)
    _validate_p_rank(p, rank)
    if not isfinite(boundary_field):
        raise ValueError("boundary_field must be finite")
    rows = []
    for left in BOUNDARY_SPINS:
        log_weights = tuple(
            log(cell_likelihood(left, right, environment, p, rank))
            + boundary_field * (right[0] + right[1]) / 2.0
            for right in BOUNDARY_SPINS
        )
        maximum = max(log_weights)
        weights = tuple(exp(value - maximum) for value in log_weights)
        normalizer = fsum(weights)
        rows.append(tuple(weight / normalizer for weight in weights))
    return tuple(rows)


@lru_cache(maxsize=None)
def shared_environment_replicated_kernel(p: float, rank: float) -> Matrix:
    """Return E_Z[K_Z tensor K_Z] on explicit pairs of spin boundaries."""

    _validate_p_rank(p, rank)
    kernels = {
        environment: collapsed_cell_kernel(environment, p, rank)
        for environment in ENVIRONMENTS
    }
    environment_masses = {
        environment: residual_environment_probability(environment, p, rank)
        for environment in ENVIRONMENTS
    }
    rows = []
    for left_first, left_second in REPLICATED_BOUNDARIES:
        first_index = BOUNDARY_SPINS.index(left_first)
        second_index = BOUNDARY_SPINS.index(left_second)
        row = []
        for right_first, right_second in REPLICATED_BOUNDARIES:
            right_first_index = BOUNDARY_SPINS.index(right_first)
            right_second_index = BOUNDARY_SPINS.index(right_second)
            row.append(
                fsum(
                    environment_masses[environment]
                    * kernels[environment][first_index][right_first_index]
                    * kernels[environment][second_index][right_second_index]
                    for environment in ENVIRONMENTS
                )
            )
        rows.append(tuple(row))
    return tuple(rows)


@lru_cache(maxsize=None)
def independent_environment_replicated_kernel(p: float, rank: float) -> Matrix:
    """Counterfactual E[K_Z] tensor E[K_Z] using two environments.

    This kernel computes a square of an annealed first moment.  It is exposed
    only as a counter-audit and must not replace the shared-environment
    kernel in a posterior-second-moment calculation.
    """

    _validate_p_rank(p, rank)
    averaged_rows = []
    for left in BOUNDARY_SPINS:
        left_index = BOUNDARY_SPINS.index(left)
        averaged_rows.append(
            tuple(
                fsum(
                    residual_environment_probability(environment, p, rank)
                    * collapsed_cell_kernel(environment, p, rank)[left_index][
                        BOUNDARY_SPINS.index(right)
                    ]
                    for environment in ENVIRONMENTS
                )
                for right in BOUNDARY_SPINS
            )
        )
    return tuple(
        tuple(
            averaged_rows[BOUNDARY_SPINS.index(left_first)][
                BOUNDARY_SPINS.index(right_first)
            ]
            * averaged_rows[BOUNDARY_SPINS.index(left_second)][
                BOUNDARY_SPINS.index(right_second)
            ]
            for right_first, right_second in REPLICATED_BOUNDARIES
        )
        for left_first, left_second in REPLICATED_BOUNDARIES
    )


def _boundary_from_orientation(orientation: int, parity: int) -> SpinPair:
    if orientation not in (-1, 1) or parity not in (-1, 1):
        raise ValueError("orientation and parity must belong to {-1,+1}")
    return orientation, orientation * parity


def replicated_fourier_sector(
    kernel: Matrix, first_character: int, second_character: int
) -> Matrix:
    """Reduce a 16-state equivariant kernel to one four-state Fourier block.

    Each character is 0 (trivial) or 1 (global-spin sign).  Rows and columns
    are indexed by ``PARITY_PAIRS``.  The input orbit representative has top
    orientation +1 in each replica.
    """

    if first_character not in (0, 1) or second_character not in (0, 1):
        raise ValueError("Fourier characters must belong to {0,1}")
    if len(kernel) != 16 or any(len(row) != 16 for row in kernel):
        raise ValueError("the replicated kernel must be 16 by 16")

    rows = []
    for left_first_parity, left_second_parity in PARITY_PAIRS:
        left_state = (
            _boundary_from_orientation(1, left_first_parity),
            _boundary_from_orientation(1, left_second_parity),
        )
        left_index = REPLICATED_BOUNDARIES.index(left_state)
        row = []
        for right_first_parity, right_second_parity in PARITY_PAIRS:
            coefficient = 0.0
            for first_orientation, second_orientation in product(
                (-1, 1), repeat=2
            ):
                right_state = (
                    _boundary_from_orientation(
                        first_orientation, right_first_parity
                    ),
                    _boundary_from_orientation(
                        second_orientation, right_second_parity
                    ),
                )
                sign = (
                    first_orientation**first_character
                    * second_orientation**second_character
                )
                coefficient += sign * kernel[left_index][
                    REPLICATED_BOUNDARIES.index(right_state)
                ]
            row.append(coefficient)
        rows.append(tuple(row))
    return tuple(rows)


def weighted_absolute_contraction(
    matrix: Matrix, weights: Sequence[float]
) -> float:
    """Return max_x sum_y |T(x,y)| w(y)/w(x)."""

    if not matrix or any(len(row) != len(matrix) for row in matrix):
        raise ValueError("matrix must be non-empty and square")
    if len(weights) != len(matrix) or min(weights) <= 0.0:
        raise ValueError("one strictly positive weight is required per state")
    return max(
        fsum(abs(entry) * weights[column] for column, entry in enumerate(row))
        / weights[row_index]
        for row_index, row in enumerate(matrix)
    )


def _apply_matrix(matrix: Matrix, vector: Sequence[float]) -> tuple[float, ...]:
    if len(matrix) != len(vector):
        raise ValueError("matrix and vector dimensions do not match")
    return tuple(
        fsum(entry * value for entry, value in zip(row, vector, strict=True))
        for row in matrix
    )


@dataclass(frozen=True)
class TriangularCellTransfer:
    """Mass and chi-tensor-chi blocks of the finite replicated transfer."""

    p: float
    rank: float
    residual_satisfaction: float
    mass: Matrix
    chi_tensor_chi: Matrix
    independent_chi_tensor_chi: Matrix
    mass_row_sum_error: float
    uniform_weight_contraction: float
    worst_boundary_state: tuple[int, int]
    scope_label: str
    missing_features: tuple[str, ...]

    @property
    def has_strict_uniform_contraction(self) -> bool:
        return self.uniform_weight_contraction < 1.0


def triangular_cell_transfer(
    p: float = 0.805, rank: float = Q_CRITICAL
) -> TriangularCellTransfer:
    """Construct and diagnose the smallest explicit triangular-strip block."""

    shared = shared_environment_replicated_kernel(p, rank)
    independent = independent_environment_replicated_kernel(p, rank)
    mass = replicated_fourier_sector(shared, 0, 0)
    twisted = replicated_fourier_sector(shared, 1, 1)
    independent_twisted = replicated_fourier_sector(independent, 1, 1)
    row_absolute_sums = tuple(fsum(abs(value) for value in row) for row in twisted)
    worst_index = max(range(4), key=row_absolute_sums.__getitem__)
    return TriangularCellTransfer(
        p=p,
        rank=rank,
        residual_satisfaction=residual_satisfaction_probability(p, rank),
        mass=mass,
        chi_tensor_chi=twisted,
        independent_chi_tensor_chi=independent_twisted,
        mass_row_sum_error=max(abs(fsum(row) - 1.0) for row in mass),
        uniform_weight_contraction=max(row_absolute_sums),
        worst_boundary_state=PARITY_PAIRS[worst_index],
        scope_label="E1+ neutral all-closed cell sector test",
        missing_features=(
            "winning edge and pivotal status",
            "open-skeleton boundary partition",
            "ancestral Lambda messages",
            "Kruskal compatibility between cells",
            "critical Palm law",
        ),
    )


def twisted_chain_second_moment(
    p: float,
    depth: int,
    initial_relative_parity: int = 1,
    rank: float = Q_CRITICAL,
) -> float:
    """Apply the chi-tensor-chi block along iid copies of the fixed cell."""

    if depth < 0:
        raise ValueError("depth must be nonnegative")
    if initial_relative_parity not in (-1, 1):
        raise ValueError("initial_relative_parity must belong to {-1,+1}")
    transfer = triangular_cell_transfer(p, rank).chi_tensor_chi
    terminal = (1.0,) * 4
    for _ in range(depth):
        terminal = _apply_matrix(transfer, terminal)
    initial = (initial_relative_parity, initial_relative_parity)
    return terminal[PARITY_PAIRS.index(initial)]


def polarized_cell_second_moment(
    p: float,
    boundary_field: float,
    initial_left: SpinPair = (1, 1),
    rank: float = Q_CRITICAL,
) -> float:
    """Exact one-cell second moment with an external polarizing field.

    The environment is shared by the replicas.  Conditional on an
    environment, ``correlation`` is the expected top-port transport sign of
    one heat bath, and this function averages ``correlation**2``.

    Since every residual likelihood is positive, as ``boundary_field`` tends
    to ``+infinity`` the right boundary converges uniformly to ``(+1,+1)``.
    Consequently this second moment tends to one for either initial top spin.
    Hence no absolute contraction coefficient below one can hold uniformly
    over unbounded exterior boundary potentials.  A full strip certificate
    must retain the boundary polarization and either control its tails or use
    a centered/annealed contraction adapted to its stationary law.
    """

    _validate_spin_pair(initial_left)
    if not isfinite(boundary_field):
        raise ValueError("boundary_field must be finite")
    left_index = BOUNDARY_SPINS.index(initial_left)
    answer = 0.0
    for environment in ENVIRONMENTS:
        kernel = polarized_collapsed_cell_kernel(
            environment, p, rank, boundary_field
        )
        correlation = fsum(
            initial_left[0] * right[0] * kernel[left_index][right_index]
            for right_index, right in enumerate(BOUNDARY_SPINS)
        )
        answer += (
            residual_environment_probability(environment, p, rank)
            * correlation
            * correlation
        )
    return answer


@dataclass(frozen=True)
class BoundaryPolarizationDiagnostic:
    """Finite values plus the exact limiting no-go for unbounded fields."""

    fields: tuple[float, ...]
    second_moments: tuple[float, ...]
    limiting_second_moment: float = 1.0

    @property
    def has_uniform_absolute_contraction(self) -> bool:
        return False


def boundary_polarization_diagnostic(
    p: float = 0.805,
    fields: Sequence[float] = (0.0, 1.0, 2.0, 4.0, 8.0, 12.0),
    rank: float = Q_CRITICAL,
) -> BoundaryPolarizationDiagnostic:
    """Evaluate the approach to the point-mass boundary obstruction."""

    values = tuple(float(field) for field in fields)
    if not values or any(not isfinite(field) or field < 0.0 for field in values):
        raise ValueError("fields must be a non-empty family of finite values >= 0")
    return BoundaryPolarizationDiagnostic(
        fields=values,
        second_moments=tuple(
            polarized_cell_second_moment(p, field, rank=rank)
            for field in values
        ),
    )


def _interval_add(
    first: RationalInterval, second: RationalInterval
) -> RationalInterval:
    return first[0] + second[0], first[1] + second[1]


def _interval_negate(value: RationalInterval) -> RationalInterval:
    return -value[1], -value[0]


def _interval_multiply(
    first: RationalInterval, second: RationalInterval
) -> RationalInterval:
    products = tuple(
        first[index] * second[other]
        for index in (0, 1)
        for other in (0, 1)
    )
    return min(products), max(products)


def _interval_divide(
    numerator: RationalInterval, denominator: RationalInterval
) -> RationalInterval:
    if denominator[0] <= 0.0:
        raise ZeroDivisionError("the interval denominator is not positive")
    reciprocal = Fraction(1, denominator[1]), Fraction(1, denominator[0])
    return _interval_multiply(numerator, reciprocal)


def _interval_sum(values: Sequence[RationalInterval]) -> RationalInterval:
    answer = Fraction(0), Fraction(0)
    for value in values:
        answer = _interval_add(answer, value)
    return answer


def _interval_product(values: Sequence[RationalInterval]) -> RationalInterval:
    answer = Fraction(1), Fraction(1)
    for value in values:
        answer = _interval_multiply(answer, value)
    return answer


def _interval_absolute_upper(value: RationalInterval) -> Fraction:
    return max(abs(value[0]), abs(value[1]))


def strict_decimal_upper_bound(
    value: Fraction, decimal_places: int
) -> str:
    """Return a fixed-point decimal strictly larger than ``value``.

    The conversion uses integer arithmetic and always advances by one unit
    in the last displayed place.  It therefore remains an upper bound rather
    than a nearest-decimal rendering which could round downward.
    """

    if value < 0:
        raise ValueError("the certified coefficient must be nonnegative")
    if decimal_places < 0:
        raise ValueError("decimal_places must be nonnegative")
    scale = 10**decimal_places
    scaled_upper = value.numerator * scale // value.denominator + 1
    if decimal_places == 0:
        return str(scaled_upper)
    integer_part, fractional_part = divmod(scaled_upper, scale)
    return f"{integer_part}.{fractional_part:0{decimal_places}d}"


@dataclass(frozen=True)
class P805UniformContractionCertificate:
    """Rational interval certificate for w=1 in the neutral E1+ cell only."""

    p: Fraction
    rank: RationalInterval
    residual_satisfaction: RationalInterval
    twisted_entries: IntervalMatrix
    row_sum_upper_bounds: tuple[Fraction, ...]
    contraction_upper_bound: Fraction

    @property
    def is_strict(self) -> bool:
        return self.contraction_upper_bound < 1


@lru_cache(maxsize=1)
def p_805_uniform_contraction_certificate(
) -> P805UniformContractionCertificate:
    """Prove by rational intervals that ||U_chi-chi||_infinity < 0.3.

    The triangular bond threshold is enclosed by consecutive decimal
    rationals.  The sign change of ``q^3-3q+1`` and monotonicity on ``(0,1)``
    certify that the desired algebraic root lies in this interval.  Every
    subsequent operation, including all heat-bath normalizers, uses outward
    rational interval arithmetic.
    """

    scale = 10**15
    rank = (
        Fraction(347296355333860, scale),
        Fraction(347296355333861, scale),
    )
    polynomial = lambda value: value**3 - 3 * value + 1
    if not 0 < rank[0] < rank[1] < 1:
        raise AssertionError("the critical-rank enclosure must lie in (0,1)")
    if not polynomial(rank[0]) > 0 > polynomial(rank[1]):
        raise AssertionError("the interval does not enclose q_triangle")

    p = Fraction(161, 200)
    satisfaction = (
        (p - rank[1]) / (1 - rank[1]),
        (p - rank[0]) / (1 - rank[0]),
    )
    failure = 1 - satisfaction[1], 1 - satisfaction[0]

    def sign_probability(signs: Sequence[int]) -> RationalInterval:
        return _interval_product(
            tuple(satisfaction if sign == 1 else failure for sign in signs)
        )

    environment_masses = {
        environment: sign_probability(environment)
        for environment in ENVIRONMENTS
    }
    kernels: dict[tuple[Environment, SpinPair], tuple[RationalInterval, ...]] = {}
    for environment in ENVIRONMENTS:
        for left in BOUNDARY_SPINS:
            weights = []
            for right in BOUNDARY_SPINS:
                relations = cell_edge_relations(left, right)
                agreements = tuple(
                    observed * candidate
                    for observed, candidate in zip(
                        environment, relations, strict=True
                    )
                )
                weights.append(sign_probability(agreements))
            normalizer = _interval_sum(tuple(weights))
            kernels[environment, left] = tuple(
                _interval_divide(weight, normalizer) for weight in weights
            )

    def shared_entry(
        left_first: SpinPair,
        left_second: SpinPair,
        right_first: SpinPair,
        right_second: SpinPair,
    ) -> RationalInterval:
        right_first_index = BOUNDARY_SPINS.index(right_first)
        right_second_index = BOUNDARY_SPINS.index(right_second)
        return _interval_sum(
            tuple(
                _interval_multiply(
                    environment_masses[environment],
                    _interval_multiply(
                        kernels[environment, left_first][right_first_index],
                        kernels[environment, left_second][right_second_index],
                    ),
                )
                for environment in ENVIRONMENTS
            )
        )

    twisted_rows = []
    for left_first_parity, left_second_parity in PARITY_PAIRS:
        left_first = _boundary_from_orientation(1, left_first_parity)
        left_second = _boundary_from_orientation(1, left_second_parity)
        row = []
        for right_first_parity, right_second_parity in PARITY_PAIRS:
            terms = []
            for first_orientation, second_orientation in product(
                (-1, 1), repeat=2
            ):
                entry = shared_entry(
                    left_first,
                    left_second,
                    _boundary_from_orientation(
                        first_orientation, right_first_parity
                    ),
                    _boundary_from_orientation(
                        second_orientation, right_second_parity
                    ),
                )
                terms.append(
                    entry
                    if first_orientation * second_orientation == 1
                    else _interval_negate(entry)
                )
            row.append(_interval_sum(tuple(terms)))
        twisted_rows.append(tuple(row))
    twisted = tuple(twisted_rows)
    row_bounds = tuple(
        sum((_interval_absolute_upper(entry) for entry in row), Fraction(0))
        for row in twisted
    )
    upper_bound = max(row_bounds)
    if upper_bound >= Fraction(3, 10):
        raise AssertionError("the rational interval audit lost the 0.3 margin")
    return P805UniformContractionCertificate(
        p=p,
        rank=rank,
        residual_satisfaction=satisfaction,
        twisted_entries=twisted,
        row_sum_upper_bounds=row_bounds,
        contraction_upper_bound=upper_bound,
    )


def main() -> None:
    transfer = triangular_cell_transfer()
    certificate = p_805_uniform_contraction_certificate()
    print(
        f"p={transfer.p:.12g} q={transfer.rank:.12f} "
        f"s={transfer.residual_satisfaction:.12f}"
    )
    print(f"scope={transfer.scope_label}")
    print(f"mass row-sum error={transfer.mass_row_sum_error:.3g}")
    print(
        "shared chi-x-chi uniform coefficient="
        f"{transfer.uniform_weight_contraction:.12f} "
        f"at parity state {transfer.worst_boundary_state}"
    )
    print(
        "rational upper bound <"
        f"{strict_decimal_upper_bound(certificate.contraction_upper_bound, 12)} "
        f"strict={certificate.is_strict}"
    )
    independent_coefficient = weighted_absolute_contraction(
        transfer.independent_chi_tensor_chi, (1.0,) * 4
    )
    print(
        "independent-environment counterfactual coefficient="
        f"{independent_coefficient:.12f}"
    )
    for depth in (1, 2, 3, 5, 10):
        print(
            f"depth={depth:2d} second_moment="
            f"{twisted_chain_second_moment(transfer.p, depth):.12g}"
        )
    polarization = boundary_polarization_diagnostic(transfer.p)
    print("polarized-boundary no-go")
    for field, value in zip(
        polarization.fields, polarization.second_moments, strict=True
    ):
        print(f"field={field:4.1f} second_moment={value:.12f}")
    print(
        "unbounded-field limit=1: retain polarization or use a centered "
        "annealed norm"
    )
    print(
        "not E2/T2: winning edge, Palm skeleton, pivotal state, and "
        "ancestral messages remain open"
    )


if __name__ == "__main__":
    main()
