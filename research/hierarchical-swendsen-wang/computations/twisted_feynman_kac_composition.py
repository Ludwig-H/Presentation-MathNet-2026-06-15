"""Finite twisted Feynman--Kac composition and backward Doob normalization.

For a positive lifted transfer ``T(z, z', epsilon)``, with
``epsilon in {-1, +1}``, set

    K(z, z') = sum_epsilon T(z, z', epsilon),
    U(z, z') = sum_epsilon epsilon T(z, z', epsilon).

This module assumes that every row of ``K`` already sums to one.  On the
support of ``K`` it then defines ``r = |U| / K``; at a zero of ``K``,
positivity forces ``U = 0`` and the convention is ``r = 0``.  Consequently
``0 <= r <= 1`` and, for a finite inhomogeneous path of kernels,

    |(U_1 ... U_N g)(z)|
        <= E_z[prod_j r_j(Z_{j-1}, Z_j) |g(Z_N)|].

The proof is the path expansion followed by the triangle inequality, since
``|U_j| = K_j r_j`` entry by entry.  If operators are numbered in their
order of application to ``g``, the left side is equivalently written
``|U_N ... U_1 g|``.  The public API numbers kernels in path order, from
``Z_0`` to ``Z_N``, to keep the Feynman--Kac expectation unambiguous.

The row-stochastic formulation is convenient but not necessary.  For a
finite inhomogeneous sequence of arbitrary nonnegative mass transfers
``K_j`` and dominated signed transfers ``|U_j| <= K_j``, backward partition
functions give an exact finite-horizon Doob transform.  The diagonal factors
telescope, so no common, time-homogeneous normalization has to be guessed.
Identifying those unnormalized transfers with the real Palm corridor and
controlling the resulting killed path measure remain open.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import isfinite
from numbers import Real

Scalar = int | float | Fraction
SignedPair = tuple[Scalar, Scalar]
LiftedKernel = tuple[tuple[SignedPair, ...], ...]
Matrix = tuple[tuple[Scalar, ...], ...]
Vector = tuple[Scalar, ...]

EPSILON_ORDER = (-1, 1)
DEFAULT_TOLERANCE = 1e-12
OPEN_STEPS = (
    "identification of the unnormalized transfers with the real Palm corridor",
    "thermodynamic control of the killed finite-horizon Doob path measure",
)


def _validate_tolerance(tolerance: float) -> None:
    if (
        isinstance(tolerance, bool)
        or not isinstance(tolerance, Real)
        or not isfinite(float(tolerance))
        or tolerance < 0
    ):
        raise ValueError("tolerance must be a finite nonnegative real")


def _validate_scalar(value: Scalar, label: str) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not isfinite(float(value))
    ):
        raise ValueError(f"{label} must be a finite real number")


def _is_exact(values: Sequence[Scalar]) -> bool:
    return all(isinstance(value, (int, Fraction)) for value in values)


def _divide(numerator: Scalar, denominator: Scalar) -> Scalar:
    """Divide while retaining rational arithmetic for integer inputs."""

    if _is_exact((numerator, denominator)):
        return Fraction(numerator) / Fraction(denominator)
    return numerator / denominator


def _equals_one(values: Sequence[Scalar], tolerance: float) -> bool:
    total = sum(values)
    if _is_exact(values):
        return total == 1
    return abs(float(total) - 1.0) <= tolerance


def _validate_matrix(
    matrix: Sequence[Sequence[Scalar]],
    label: str,
) -> tuple[int, int]:
    if not isinstance(matrix, Sequence) or not matrix:
        raise ValueError(f"{label} must have at least one row")
    if not isinstance(matrix[0], Sequence):
        raise ValueError(f"{label} must be a sequence of rows")
    width = len(matrix[0])
    if width == 0 or any(
        not isinstance(row, Sequence) or len(row) != width for row in matrix
    ):
        raise ValueError(f"{label} must be a nonempty rectangular matrix")
    for row in matrix:
        for value in row:
            _validate_scalar(value, f"an entry of {label}")
    return len(matrix), width


@dataclass(frozen=True)
class CommonNormalizedTransfer:
    """Mass, twisted sector, and local killing ratio of one finite block."""

    lifted: LiftedKernel
    mass: Matrix
    twisted: Matrix
    ratio: Matrix
    source_size: int
    target_size: int

    @property
    def is_exact(self) -> bool:
        return all(
            isinstance(value, (int, Fraction))
            for row in self.lifted
            for pair in row
            for value in pair
        )


@dataclass(frozen=True)
class DominatedTransfer:
    """One finite positive/signed pair without a row normalization.

    The sole structural condition is the entrywise domination
    ``|twisted| <= mass``.  It is exactly the condition needed to represent
    the signed entry as the difference of two nonnegative lifted masses.
    """

    mass: Matrix
    twisted: Matrix
    ratio: Matrix
    source_size: int
    target_size: int

    @property
    def is_exact(self) -> bool:
        return all(
            isinstance(value, (int, Fraction))
            for matrix in (self.mass, self.twisted, self.ratio)
            for row in matrix
            for value in row
        )


def dominated_transfer(
    mass: Sequence[Sequence[Scalar]],
    twisted: Sequence[Sequence[Scalar]],
) -> DominatedTransfer:
    """Validate an unnormalized pair ``(K,U)`` with ``|U| <= K``.

    No stochasticity or common eigenvector is assumed.  Integer and rational
    inputs retain exact arithmetic in the local ratios.
    """

    mass_shape = _validate_matrix(mass, "mass")
    twisted_shape = _validate_matrix(twisted, "twisted")
    if mass_shape != twisted_shape:
        raise ValueError("mass and twisted matrices must have the same shape")

    frozen_mass = []
    frozen_twisted = []
    ratios = []
    for mass_row, twisted_row in zip(mass, twisted, strict=True):
        if any(value < 0 for value in mass_row):
            raise ValueError("mass entries must be nonnegative")
        ratio_row = []
        for mass_value, twisted_value in zip(mass_row, twisted_row, strict=True):
            if abs(twisted_value) > mass_value:
                raise ValueError("the twist must satisfy |U| <= K entrywise")
            ratio_row.append(
                _divide(abs(twisted_value), mass_value)
                if mass_value != 0
                else type(mass_value)(0)
            )
        frozen_mass.append(tuple(mass_row))
        frozen_twisted.append(tuple(twisted_row))
        ratios.append(tuple(ratio_row))
    return DominatedTransfer(
        mass=tuple(frozen_mass),
        twisted=tuple(frozen_twisted),
        ratio=tuple(ratios),
        source_size=mass_shape[0],
        target_size=mass_shape[1],
    )


def normalize_lifted_transfer(
    lifted: Sequence[Sequence[Sequence[Scalar]]],
    *,
    tolerance: float = DEFAULT_TOLERANCE,
) -> CommonNormalizedTransfer:
    """Validate a lifted kernel whose common mass kernel is stochastic.

    Each innermost pair is ordered as ``(T_epsilon=-1, T_epsilon=+1)``.
    No normalization is silently performed: a non-stochastic mass kernel is
    rejected.  This makes the still-open common Doob normalization visible.
    """

    _validate_tolerance(tolerance)
    if not isinstance(lifted, Sequence) or not lifted:
        raise ValueError("a lifted transfer must have at least one row")
    if not isinstance(lifted[0], Sequence):
        raise ValueError("a lifted transfer must be a sequence of rows")
    width = len(lifted[0])
    if width == 0 or any(
        not isinstance(row, Sequence) or len(row) != width for row in lifted
    ):
        raise ValueError("a lifted transfer must be nonempty and rectangular")

    frozen_rows: list[tuple[SignedPair, ...]] = []
    mass_rows: list[tuple[Scalar, ...]] = []
    twisted_rows: list[tuple[Scalar, ...]] = []
    ratio_rows: list[tuple[Scalar, ...]] = []
    for row in lifted:
        frozen_row: list[SignedPair] = []
        mass_row: list[Scalar] = []
        twisted_row: list[Scalar] = []
        ratio_row: list[Scalar] = []
        for pair in row:
            if not isinstance(pair, Sequence) or len(pair) != 2:
                raise ValueError(
                    "each lifted entry must contain the masses for epsilon=-1,+1"
                )
            negative, positive = pair
            _validate_scalar(negative, "a lifted mass")
            _validate_scalar(positive, "a lifted mass")
            if negative < 0 or positive < 0:
                raise ValueError("lifted masses must be nonnegative")
            mass = negative + positive
            twisted = positive - negative
            ratio = _divide(abs(twisted), mass) if mass != 0 else type(mass)(0)
            if mass == 0 and twisted != 0:
                raise AssertionError("positivity should force U=0 wherever K=0")
            if ratio < 0 or ratio > 1:
                raise AssertionError("the local Radon--Nikodym ratio left [0,1]")
            frozen_row.append((negative, positive))
            mass_row.append(mass)
            twisted_row.append(twisted)
            ratio_row.append(ratio)
        if not _equals_one(mass_row, tolerance):
            raise ValueError("every row of the common mass kernel must sum to one")
        frozen_rows.append(tuple(frozen_row))
        mass_rows.append(tuple(mass_row))
        twisted_rows.append(tuple(twisted_row))
        ratio_rows.append(tuple(ratio_row))

    return CommonNormalizedTransfer(
        lifted=tuple(frozen_rows),
        mass=tuple(mass_rows),
        twisted=tuple(twisted_rows),
        ratio=tuple(ratio_rows),
        source_size=len(frozen_rows),
        target_size=width,
    )


def canonical_lift_from_mass_and_twist(
    mass: Sequence[Sequence[Scalar]],
    twisted: Sequence[Sequence[Scalar]],
    *,
    tolerance: float = DEFAULT_TOLERANCE,
) -> CommonNormalizedTransfer:
    """Build the unique binary-sign lift ``T_-/+ = (K -/+ U)/2``.

    Existence is equivalent to the entrywise domination ``|U| <= K``.  The
    check is strict, including for floating-point inputs; ``tolerance`` is
    used only for stochastic row sums.
    """

    _validate_tolerance(tolerance)
    mass_shape = _validate_matrix(mass, "mass")
    twisted_shape = _validate_matrix(twisted, "twisted")
    if mass_shape != twisted_shape:
        raise ValueError("mass and twisted matrices must have the same shape")

    lifted_rows = []
    for mass_row, twisted_row in zip(mass, twisted, strict=True):
        if any(value < 0 for value in mass_row):
            raise ValueError("mass entries must be nonnegative")
        if not _equals_one(mass_row, tolerance):
            raise ValueError("every row of the common mass kernel must sum to one")
        lifted_row = []
        for mass_value, twisted_value in zip(mass_row, twisted_row, strict=True):
            if abs(twisted_value) > mass_value:
                raise ValueError("the twist must satisfy |U| <= K entrywise")
            lifted_row.append(
                (
                    _divide(mass_value - twisted_value, 2),
                    _divide(mass_value + twisted_value, 2),
                )
            )
        lifted_rows.append(tuple(lifted_row))
    return normalize_lifted_transfer(tuple(lifted_rows), tolerance=tolerance)


def _validate_terminal(terminal: Sequence[Scalar]) -> Vector:
    if not terminal:
        raise ValueError("the terminal function must be nonempty")
    for value in terminal:
        _validate_scalar(value, "a terminal value")
    return tuple(terminal)


def _apply(matrix: Matrix, vector: Sequence[Scalar]) -> Vector:
    if not matrix or len(matrix[0]) != len(vector):
        raise ValueError("matrix and vector dimensions do not match")
    return tuple(
        sum(
            (entry * value for entry, value in zip(row, vector, strict=True)),
            0,
        )
        for row in matrix
    )


def _killed_mass_matrix(step: CommonNormalizedTransfer) -> Matrix:
    return tuple(
        tuple(mass * ratio for mass, ratio in zip(mass_row, ratio_row, strict=True))
        for mass_row, ratio_row in zip(step.mass, step.ratio, strict=True)
    )


def _validate_chain(
    steps: Sequence[CommonNormalizedTransfer], terminal_size: int
) -> None:
    if not steps:
        return
    for first, second in zip(steps[:-1], steps[1:], strict=True):
        if first.target_size != second.source_size:
            raise ValueError("successive transfer dimensions do not match")
    if steps[-1].target_size != terminal_size:
        raise ValueError("the last transfer does not match the terminal function")


def _validate_dominated_chain(
    steps: Sequence[DominatedTransfer], terminal_size: int
) -> None:
    if not steps:
        return
    for first, second in zip(steps[:-1], steps[1:], strict=True):
        if first.target_size != second.source_size:
            raise ValueError("successive transfer dimensions do not match")
    if steps[-1].target_size != terminal_size:
        raise ValueError("the last transfer does not match the terminal function")


def _pointwise_divide_on_support(
    numerator: Sequence[Scalar], denominator: Sequence[Scalar]
) -> Vector:
    if len(numerator) != len(denominator):
        raise ValueError("pointwise division requires vectors of equal length")
    result = []
    for top, bottom in zip(numerator, denominator, strict=True):
        if bottom < 0:
            raise ValueError("a reference potential cannot be negative")
        if bottom == 0:
            if top != 0:
                raise AssertionError("a dominated numerator survived off support")
            result.append(type(bottom)(0))
        else:
            result.append(_divide(top, bottom))
    return tuple(result)


@dataclass(frozen=True)
class FiniteHorizonDoobCertificate:
    """Exact telescoping normalization for unnormalized finite transfers.

    Rows attached to a zero backward potential are outside the finite-horizon
    support.  They are stored as zero rows; every active row of the Doob
    kernels is stochastic.
    """

    backward_potentials: tuple[Vector, ...]
    active_states: tuple[tuple[bool, ...], ...]
    doob_kernels: tuple[Matrix, ...]
    signed_values: Vector
    feynman_kac_envelope: Vector
    normalized_absolute_values: Vector
    normalized_feynman_kac_envelope: Vector
    doob_expectation: Vector
    normalized_slacks: Vector
    active_rows_are_stochastic: bool
    normalization_identity_holds: bool
    inequality_holds: bool
    is_exact: bool
    scope_label: str = (
        "finite inhomogeneous composition with an exact backward Doob transform"
    )
    open_steps: tuple[str, ...] = OPEN_STEPS


def finite_horizon_doob_certificate(
    steps: Sequence[DominatedTransfer],
    terminal: Sequence[Scalar],
    reference_terminal: Sequence[Scalar],
    *,
    tolerance: float = DEFAULT_TOLERANCE,
) -> FiniteHorizonDoobCertificate:
    """Normalize an arbitrary finite dominated chain by backward potentials.

    Let ``h_N`` be ``reference_terminal`` and recursively set
    ``h_{j-1}=K_j h_j``.  On ``{h_{j-1}>0}``,

    ``P_j(z,z') = K_j(z,z') h_j(z') / h_{j-1}(z)``

    is stochastic.  With ``r_j=|U_j|/K_j`` and ``|g|<=h_N``, diagonal
    factors telescope and give the normalized Feynman--Kac bound.  Zero
    backward potentials are retained as inactive states, so impossible rows
    need not be completed by an artificial stochastic kernel.
    """

    _validate_tolerance(tolerance)
    terminal_values = _validate_terminal(terminal)
    reference_values = _validate_terminal(reference_terminal)
    if len(terminal_values) != len(reference_values):
        raise ValueError("terminal and reference_terminal must have equal length")
    if any(value < 0 for value in reference_values):
        raise ValueError("reference_terminal must be nonnegative")
    if any(
        abs(value) > reference
        for value, reference in zip(terminal_values, reference_values, strict=True)
    ):
        raise ValueError("the terminal function must satisfy |g| <= h_N")

    frozen_steps = tuple(steps)
    if any(not isinstance(step, DominatedTransfer) for step in frozen_steps):
        raise ValueError("every step must be a DominatedTransfer")
    _validate_dominated_chain(frozen_steps, len(terminal_values))

    backward_reversed = [reference_values]
    for step in reversed(frozen_steps):
        backward_reversed.append(_apply(step.mass, backward_reversed[-1]))
    backward = tuple(reversed(backward_reversed))
    active_states = tuple(
        tuple(value > 0 for value in potential) for potential in backward
    )

    doob_kernels = []
    active_rows_are_stochastic = True
    for index, step in enumerate(frozen_steps):
        source_potential = backward[index]
        target_potential = backward[index + 1]
        rows = []
        for source, mass_row in enumerate(step.mass):
            if source_potential[source] == 0:
                rows.append(tuple(type(value)(0) for value in mass_row))
                continue
            row = tuple(
                _divide(
                    mass * target_potential[target],
                    source_potential[source],
                )
                for target, mass in enumerate(mass_row)
            )
            if not _equals_one(row, tolerance):
                active_rows_are_stochastic = False
            rows.append(row)
        doob_kernels.append(tuple(rows))

    signed = terminal_values
    envelope = tuple(abs(value) for value in terminal_values)
    doob_value = _pointwise_divide_on_support(envelope, reference_values)
    for index in reversed(range(len(frozen_steps))):
        step = frozen_steps[index]
        signed = _apply(step.twisted, signed)
        envelope = _apply(
            tuple(
                tuple(
                    mass * ratio
                    for mass, ratio in zip(mass_row, ratio_row, strict=True)
                )
                for mass_row, ratio_row in zip(step.mass, step.ratio, strict=True)
            ),
            envelope,
        )
        doob_value = _apply(
            tuple(
                tuple(
                    probability * ratio
                    for probability, ratio in zip(
                        probability_row, ratio_row, strict=True
                    )
                )
                for probability_row, ratio_row in zip(
                    doob_kernels[index], step.ratio, strict=True
                )
            ),
            doob_value,
        )

    absolute = tuple(abs(value) for value in signed)
    normalized_absolute = _pointwise_divide_on_support(absolute, backward[0])
    normalized_envelope = _pointwise_divide_on_support(envelope, backward[0])
    slacks = tuple(
        upper - lower
        for upper, lower in zip(normalized_envelope, normalized_absolute, strict=True)
    )
    exact = (
        all(step.is_exact for step in frozen_steps)
        and _is_exact(terminal_values)
        and _is_exact(reference_values)
    )
    if exact:
        normalization_holds = doob_value == normalized_envelope
        inequality_holds = all(slack >= 0 for slack in slacks)
    else:
        normalization_holds = all(
            abs(float(first) - float(second)) <= tolerance
            for first, second in zip(doob_value, normalized_envelope, strict=True)
        )
        inequality_holds = all(float(slack) >= -tolerance for slack in slacks)
    if not active_rows_are_stochastic:
        raise AssertionError("the backward Doob transform is not stochastic")
    if not normalization_holds:
        raise AssertionError("the Doob and unnormalized envelopes disagree")
    if not inequality_holds:
        raise AssertionError("the normalized Feynman--Kac domination failed")

    return FiniteHorizonDoobCertificate(
        backward_potentials=backward,
        active_states=active_states,
        doob_kernels=tuple(doob_kernels),
        signed_values=signed,
        feynman_kac_envelope=envelope,
        normalized_absolute_values=normalized_absolute,
        normalized_feynman_kac_envelope=normalized_envelope,
        doob_expectation=doob_value,
        normalized_slacks=slacks,
        active_rows_are_stochastic=active_rows_are_stochastic,
        normalization_identity_holds=normalization_holds,
        inequality_holds=inequality_holds,
        is_exact=exact,
    )


@dataclass(frozen=True)
class FiniteCompositionCertificate:
    """Pointwise certificate for a finite, commonly normalized sequence."""

    signed_values: Vector
    absolute_values: Vector
    feynman_kac_envelope: Vector
    slacks: Vector
    is_exact: bool
    inequality_holds: bool
    scope_label: str = "finite composition after common stochastic normalization"
    open_steps: tuple[str, ...] = OPEN_STEPS


def finite_composition_certificate(
    steps: Sequence[CommonNormalizedTransfer],
    terminal: Sequence[Scalar],
    *,
    tolerance: float = DEFAULT_TOLERANCE,
) -> FiniteCompositionCertificate:
    """Compute both sides of the finite twisted Feynman--Kac inequality."""

    _validate_tolerance(tolerance)
    terminal_values = _validate_terminal(terminal)
    frozen_steps = tuple(steps)
    if any(not isinstance(step, CommonNormalizedTransfer) for step in frozen_steps):
        raise ValueError("every step must be a CommonNormalizedTransfer")
    _validate_chain(frozen_steps, len(terminal_values))

    signed = terminal_values
    envelope = tuple(abs(value) for value in terminal_values)
    for step in reversed(frozen_steps):
        signed = _apply(step.twisted, signed)
        envelope = _apply(_killed_mass_matrix(step), envelope)
    absolute = tuple(abs(value) for value in signed)
    slacks = tuple(
        upper - lower for upper, lower in zip(envelope, absolute, strict=True)
    )
    exact = all(step.is_exact for step in frozen_steps) and _is_exact(terminal_values)
    if exact:
        holds = all(slack >= 0 for slack in slacks)
    else:
        holds = all(float(slack) >= -tolerance for slack in slacks)
    if not holds:
        raise AssertionError("the finite Feynman--Kac domination failed")
    return FiniteCompositionCertificate(
        signed_values=signed,
        absolute_values=absolute,
        feynman_kac_envelope=envelope,
        slacks=slacks,
        is_exact=exact,
        inequality_holds=holds,
    )


def brute_force_path_expansion(
    steps: Sequence[CommonNormalizedTransfer],
    terminal: Sequence[Scalar],
) -> tuple[Vector, Vector]:
    """Enumerate all state paths, independently of backward recursion.

    The first returned vector is the signed path sum.  The second is the
    killed Markov expectation with edge weight ``r``.  This routine is
    exponential and intended only as a finite counter-audit.
    """

    terminal_values = _validate_terminal(terminal)
    frozen_steps = tuple(steps)
    if any(not isinstance(step, CommonNormalizedTransfer) for step in frozen_steps):
        raise ValueError("every step must be a CommonNormalizedTransfer")
    _validate_chain(frozen_steps, len(terminal_values))
    if not frozen_steps:
        return terminal_values, tuple(abs(value) for value in terminal_values)

    signed_values = []
    envelope_values = []
    tail_spaces = tuple(range(step.target_size) for step in frozen_steps)
    for initial in range(frozen_steps[0].source_size):
        signed_total: Scalar = 0
        envelope_total: Scalar = 0
        for tail in product(*tail_spaces):
            signed_weight: Scalar = 1
            envelope_weight: Scalar = 1
            previous = initial
            for step, following in zip(frozen_steps, tail, strict=True):
                signed_weight *= step.twisted[previous][following]
                envelope_weight *= (
                    step.mass[previous][following] * step.ratio[previous][following]
                )
                previous = following
            signed_total += signed_weight * terminal_values[tail[-1]]
            envelope_total += envelope_weight * abs(terminal_values[tail[-1]])
        signed_values.append(signed_total)
        envelope_values.append(envelope_total)
    return tuple(signed_values), tuple(envelope_values)


@dataclass(frozen=True)
class FractionCompositionAudit:
    """Exact dynamic-versus-path audit for one small rational construction."""

    name: str
    certificate: FiniteCompositionCertificate
    brute_signed_values: Vector
    brute_envelope_values: Vector

    @property
    def matches_exactly(self) -> bool:
        return (
            self.certificate.signed_values == self.brute_signed_values
            and self.certificate.feynman_kac_envelope == self.brute_envelope_values
        )


def _fraction_audit(
    name: str,
    steps: Sequence[CommonNormalizedTransfer],
    terminal: Sequence[Fraction],
) -> FractionCompositionAudit:
    certificate = finite_composition_certificate(steps, terminal)
    brute_signed, brute_envelope = brute_force_path_expansion(steps, terminal)
    if not certificate.is_exact or not certificate.inequality_holds:
        raise AssertionError("the rational composition audit was not exact")
    if (
        certificate.signed_values != brute_signed
        or certificate.feynman_kac_envelope != brute_envelope
    ):
        raise AssertionError("dynamic and exhaustive rational audits disagree")
    return FractionCompositionAudit(
        name=name,
        certificate=certificate,
        brute_signed_values=brute_signed,
        brute_envelope_values=brute_envelope,
    )


def direct_fraction_lift_audit() -> FractionCompositionAudit:
    """Audit two inhomogeneous kernels specified directly at lifted level."""

    first = normalize_lifted_transfer(
        (
            (
                (Fraction(1, 12), Fraction(1, 4)),
                (Fraction(1, 6), Fraction(1, 6)),
                (Fraction(1, 4), Fraction(1, 12)),
            ),
            (
                (Fraction(0), Fraction(1, 2)),
                (Fraction(1, 8), Fraction(1, 8)),
                (Fraction(1, 4), Fraction(0)),
            ),
        )
    )
    second = normalize_lifted_transfer(
        (
            (
                (Fraction(1, 10), Fraction(2, 5)),
                (Fraction(1, 4), Fraction(1, 4)),
            ),
            (
                (Fraction(1, 3), Fraction(0)),
                (Fraction(1, 6), Fraction(1, 2)),
            ),
            (
                (Fraction(0), Fraction(1, 4)),
                (Fraction(3, 8), Fraction(3, 8)),
            ),
        )
    )
    return _fraction_audit(
        "direct positive lift on a 2-to-3-to-2 chain",
        (first, second),
        (Fraction(2, 3), Fraction(-3, 5)),
    )


def canonical_fraction_lift_audit() -> FractionCompositionAudit:
    """Audit an independent construction starting from dominated ``(K,U)``."""

    first = canonical_lift_from_mass_and_twist(
        (
            (Fraction(1, 3), Fraction(2, 3)),
            (Fraction(3, 4), Fraction(1, 4)),
        ),
        (
            (Fraction(1, 6), Fraction(-1, 3)),
            (Fraction(0), Fraction(1, 8)),
        ),
    )
    second = canonical_lift_from_mass_and_twist(
        (
            (Fraction(2, 5), Fraction(3, 5)),
            (Fraction(1, 2), Fraction(1, 2)),
        ),
        (
            (Fraction(-1, 5), Fraction(1, 10)),
            (Fraction(1, 4), Fraction(-1, 4)),
        ),
    )
    return _fraction_audit(
        "canonical binary-sign lift of two dominated (K,U) pairs",
        (first, second),
        (Fraction(1), Fraction(-2, 3)),
    )


def unnormalized_fraction_doob_audit() -> FiniteHorizonDoobCertificate:
    """Exact 2-to-3-to-2 audit with non-stochastic mass transfers."""

    first = dominated_transfer(
        (
            (Fraction(2), Fraction(1), Fraction(0)),
            (Fraction(1, 2), Fraction(3, 2), Fraction(1)),
        ),
        (
            (Fraction(1), Fraction(-1, 2), Fraction(0)),
            (Fraction(-1, 4), Fraction(1), Fraction(1, 2)),
        ),
    )
    second = dominated_transfer(
        (
            (Fraction(1), Fraction(2)),
            (Fraction(3), Fraction(1)),
            (Fraction(0), Fraction(5)),
        ),
        (
            (Fraction(1, 2), Fraction(-1)),
            (Fraction(-2), Fraction(1, 2)),
            (Fraction(0), Fraction(3)),
        ),
    )
    return finite_horizon_doob_certificate(
        (first, second),
        (Fraction(2), Fraction(1)),
        (Fraction(3), Fraction(3)),
    )


@dataclass(frozen=True)
class TriangularNeutralPowerDiagnostic:
    """Composition diagnostic for the existing neutral triangular cell."""

    p: float
    rank: float
    depths: tuple[int, ...]
    signed_values: tuple[float, ...]
    feynman_kac_envelopes: tuple[float, ...]
    uniform_coefficient: float
    uniform_power_bounds: tuple[float, ...]
    all_composition_bounds_hold: bool
    scope_label: str
    open_steps: tuple[str, ...] = OPEN_STEPS


def triangular_neutral_power_diagnostic(
    depths: Sequence[int] = (1, 2, 3, 5, 10),
    *,
    p: float = 0.805,
) -> TriangularNeutralPowerDiagnostic:
    """Apply the finite lemma to powers of the neutral all-closed cell.

    The import is local so this abstract module remains usable on its own.
    The canonical lift is valid because the mass and twisted Fourier blocks
    arise from one positive replicated transfer.
    """

    from triangular_band_collapsed_certificate import (  # local by design
        PARITY_PAIRS,
        triangular_cell_transfer,
    )

    requested_depths = tuple(depths)
    if not requested_depths or any(
        isinstance(depth, bool) or not isinstance(depth, int) or depth < 0
        for depth in requested_depths
    ):
        raise ValueError("depths must be a nonempty sequence of integers >= 0")

    transfer = triangular_cell_transfer(p=p)
    step = canonical_lift_from_mass_and_twist(transfer.mass, transfer.chi_tensor_chi)
    initial = PARITY_PAIRS.index((1, 1))
    signed_values = []
    envelopes = []
    for depth in requested_depths:
        certificate = finite_composition_certificate(
            (step,) * depth, (1.0,) * step.target_size
        )
        signed_values.append(float(certificate.signed_values[initial]))
        envelopes.append(float(certificate.feynman_kac_envelope[initial]))
    coefficient = transfer.uniform_weight_contraction
    power_bounds = tuple(coefficient**depth for depth in requested_depths)
    bounds_hold = all(
        abs(signed) <= envelope + DEFAULT_TOLERANCE
        and envelope <= power_bound + DEFAULT_TOLERANCE
        for signed, envelope, power_bound in zip(
            signed_values, envelopes, power_bounds, strict=True
        )
    )
    return TriangularNeutralPowerDiagnostic(
        p=transfer.p,
        rank=transfer.rank,
        depths=requested_depths,
        signed_values=tuple(signed_values),
        feynman_kac_envelopes=tuple(envelopes),
        uniform_coefficient=coefficient,
        uniform_power_bounds=power_bounds,
        all_composition_bounds_hold=bounds_hold,
        scope_label=(
            "finite powers of the E1+ neutral all-closed triangular cell; "
            "not a real Palm-corridor certificate"
        ),
    )


def main() -> None:
    for audit in (direct_fraction_lift_audit(), canonical_fraction_lift_audit()):
        print(
            f"{audit.name}: exact={audit.certificate.is_exact} "
            f"matches_paths={audit.matches_exactly} "
            f"slacks={audit.certificate.slacks}"
        )
    doob = unnormalized_fraction_doob_audit()
    print(
        "unnormalized finite-horizon Doob audit: "
        f"exact={doob.is_exact} "
        f"stochastic={doob.active_rows_are_stochastic} "
        f"telescopes={doob.normalization_identity_holds} "
        f"slacks={doob.normalized_slacks}"
    )
    diagnostic = triangular_neutral_power_diagnostic()
    print(diagnostic.scope_label)
    print(f"uniform coefficient={diagnostic.uniform_coefficient:.12f}")
    for depth, signed, envelope, power in zip(
        diagnostic.depths,
        diagnostic.signed_values,
        diagnostic.feynman_kac_envelopes,
        diagnostic.uniform_power_bounds,
        strict=True,
    ):
        print(
            f"depth={depth:2d} signed={signed:.12g} "
            f"FK={envelope:.12g} uniform={power:.12g}"
        )
    print("open: " + "; ".join(diagnostic.open_steps))


if __name__ == "__main__":
    main()
