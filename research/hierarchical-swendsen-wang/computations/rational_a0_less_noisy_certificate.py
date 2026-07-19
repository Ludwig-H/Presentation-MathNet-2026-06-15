"""Exact rational certificates for several triangle-channel milestones.

The first channel parameters are

    p = 161/200, q = 61/100, (a, s, e) = (41, 14, 42)/125,
    delta = 1/200.

The same proof also certifies the sharper rational candidate

    p = 809/1000, q = 309/500,
    (a, s, e) = (1660, 559, 1663)/5000, delta = 1/50000.

The third built-in candidate raises the bound further:

    p = 4047/5000, q = 1547/2500,
    (a, s, e) = (33317, 11118, 33329)/100000, delta = 1/1000000.

The sharpest built-in rational point is

    p = 809439/1000000, q = 309439/500000,
    (a, s, e) = (166642280, 55571811, 166642287)/500000000,
    delta = 1/50000000.

For each candidate the module proves the stronger fixed-prior quadratic
comparison

    Q_E(mu, f) - Q_Y(mu, f) >= delta Var_mu(f)

for every prior on the four relative-spin states.  It uses no floating-point
or sampled-simplex argument.  Four univariate constant-sign statements are
certified by exact Sturm sequences.  In the polarized sector they imply that
the exact three-by-three gap matrix is symmetric diagonally dominant.

The key cancellation is particularly useful.  Put

    d(t) = (a - delta - c_q(t)) / t.

After eliminating ``g[0]`` with ``sum(g)=0``, the gap matrix satisfies

    M[i,i] - M[i,j] - M[i,k] = d(mu[i]) - d(mu[0]).

Its off-diagonal entries have the form ``d(mu[0]) + s/D``.  Exact Sturm
certificates show that tail values of ``d`` dominate all values on
``[1/2,1]`` and that every off-diagonal entry is positive.  Gershgorin's
theorem then proves positive semidefiniteness.  Priors on proper faces follow
by support restriction (equivalently, by continuity of the finite-channel
quadratic forms).

Unlike a numerical stress test, :func:`a0_less_noisy_certificate` and
:func:`p809_less_noisy_certificate` and
:func:`p8094_less_noisy_certificate` and
:func:`p809439_less_noisy_certificate` return exhaustive certificates with
zero unresolved regions.
"""

from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import gcd, lcm
from typing import Sequence

Polynomial = tuple[Fraction, ...]  # coefficients in increasing degree
Matrix = tuple[tuple[Fraction, ...], ...]


@dataclass(frozen=True)
class RationalTriangleCandidate:
    """Rational triangle erasure channel and requested variance gap."""

    name: str
    p: Fraction
    q: Fraction
    full: Fraction
    single: Fraction
    empty: Fraction
    variance_gap: Fraction

    @property
    def reduced_full(self) -> Fraction:
        return self.full - self.variance_gap


A0_CANDIDATE = RationalTriangleCandidate(
    name="a0_p805",
    p=Fraction(161, 200),
    q=Fraction(61, 100),
    full=Fraction(41, 125),
    single=Fraction(14, 125),
    empty=Fraction(42, 125),
    variance_gap=Fraction(1, 200),
)

P809_CANDIDATE = RationalTriangleCandidate(
    name="p809",
    p=Fraction(809, 1000),
    q=Fraction(309, 500),
    full=Fraction(1660, 5000),
    single=Fraction(559, 5000),
    empty=Fraction(1663, 5000),
    variance_gap=Fraction(1, 50000),
)

P8094_CANDIDATE = RationalTriangleCandidate(
    name="p8094",
    p=Fraction(4047, 5000),
    q=Fraction(1547, 2500),
    full=Fraction(33317, 100000),
    single=Fraction(11118, 100000),
    empty=Fraction(33329, 100000),
    variance_gap=Fraction(1, 1000000),
)

P809439_CANDIDATE = RationalTriangleCandidate(
    name="p809439",
    p=Fraction(809439, 1000000),
    q=Fraction(309439, 500000),
    full=Fraction(166642280, 500000000),
    single=Fraction(55571811, 500000000),
    empty=Fraction(166642287, 500000000),
    variance_gap=Fraction(1, 50000000),
)

# Backward-compatible names for the original A0 certificate.
P_0 = A0_CANDIDATE.p
Q_0 = A0_CANDIDATE.q
R_0 = Q_0 * Q_0
A_0 = A0_CANDIDATE.full
S_0 = A0_CANDIDATE.single
E_0 = A0_CANDIDATE.empty
DELTA_0 = A0_CANDIDATE.variance_gap
REDUCED_FULL = A0_CANDIDATE.reduced_full

HALF = Fraction(1, 2)
ONE = Fraction(1)

STATES: tuple[tuple[int, int], ...] = tuple(product((1, -1), repeat=2))
CHARACTERS: tuple[tuple[int, ...], ...] = tuple(
    tuple((first, second, first * second)[index] for first, second in STATES)
    for index in range(3)
)


def _trim(polynomial: Sequence[Fraction | int]) -> Polynomial:
    values = [Fraction(value) for value in polynomial]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values) if values else (Fraction(0),)


def _add(first: Polynomial, second: Polynomial) -> Polynomial:
    size = max(len(first), len(second))
    return _trim(
        tuple(
            (first[index] if index < len(first) else 0)
            + (second[index] if index < len(second) else 0)
            for index in range(size)
        )
    )


def _subtract(first: Polynomial, second: Polynomial) -> Polynomial:
    return _add(first, tuple(-coefficient for coefficient in second))


def _multiply(first: Polynomial, second: Polynomial) -> Polynomial:
    answer = [Fraction(0)] * (len(first) + len(second) - 1)
    for first_degree, first_coefficient in enumerate(first):
        for second_degree, second_coefficient in enumerate(second):
            answer[first_degree + second_degree] += (
                first_coefficient * second_coefficient
            )
    return _trim(answer)


def _scale(polynomial: Polynomial, scalar: Fraction | int) -> Polynomial:
    scalar = Fraction(scalar)
    return _trim(tuple(scalar * coefficient for coefficient in polynomial))


def _derivative(polynomial: Polynomial) -> Polynomial:
    if len(polynomial) == 1:
        return (Fraction(0),)
    return _trim(
        tuple(degree * polynomial[degree] for degree in range(1, len(polynomial)))
    )


def _evaluate(polynomial: Polynomial, point: Fraction | int) -> Fraction:
    point = Fraction(point)
    answer = Fraction(0)
    for coefficient in reversed(polynomial):
        answer = answer * point + coefficient
    return answer


def _divide_polynomials(
    dividend: Polynomial, divisor: Polynomial
) -> tuple[Polynomial, Polynomial]:
    """Return exact quotient and remainder in ``Q[t]``."""

    dividend = _trim(dividend)
    divisor = _trim(divisor)
    if divisor == (0,):
        raise ZeroDivisionError("the polynomial divisor is zero")
    if len(dividend) < len(divisor):
        return (Fraction(0),), dividend

    remainder = list(dividend)
    quotient = [Fraction(0)] * (len(dividend) - len(divisor) + 1)
    while len(remainder) >= len(divisor) and any(remainder):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] = coefficient
        for index, divisor_coefficient in enumerate(divisor):
            remainder[index + shift] -= coefficient * divisor_coefficient
        remainder = list(_trim(remainder))
    return _trim(quotient), _trim(remainder)


def sturm_sequence(polynomial: Sequence[Fraction | int]) -> tuple[Polynomial, ...]:
    """Build the exact Sturm sequence of a nonconstant polynomial."""

    first = _trim(polynomial)
    second = _derivative(first)
    if len(first) <= 1 or second == (0,):
        raise ValueError("Sturm sequences require a nonconstant polynomial")
    sequence = [first, second]
    while sequence[-1] != (0,):
        _, remainder = _divide_polynomials(sequence[-2], sequence[-1])
        if remainder == (0,):
            break
        sequence.append(_scale(remainder, -1))
    return tuple(sequence)


def _sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def _sturm_variations(sequence: Sequence[Polynomial], point: Fraction) -> int:
    signs = tuple(
        sign
        for sign in (_sign(_evaluate(polynomial, point)) for polynomial in sequence)
        if sign
    )
    return sum(first != second for first, second in zip(signs, signs[1:]))


def sturm_root_count(
    polynomial: Sequence[Fraction | int],
    left: Fraction | int,
    right: Fraction | int,
) -> int:
    """Count distinct roots in ``(left,right)`` when endpoints are not roots."""

    left = Fraction(left)
    right = Fraction(right)
    if not left < right:
        raise ValueError("the Sturm interval must be nonempty")
    normalized = _trim(polynomial)
    if _evaluate(normalized, left) == 0 or _evaluate(normalized, right) == 0:
        raise ValueError("the current exact root counter excludes endpoint roots")
    sequence = sturm_sequence(normalized)
    return _sturm_variations(sequence, left) - _sturm_variations(sequence, right)


@dataclass(frozen=True)
class ConstantSignCertificate:
    """Exact constant-sign certificate for one integer polynomial."""

    name: str
    polynomial: tuple[int, ...]
    interval: tuple[Fraction, Fraction]
    witness: Fraction
    expected_sign: int
    root_count: int
    witness_value: Fraction

    @property
    def certified(self) -> bool:
        return (
            self.expected_sign in (-1, 1)
            and self.root_count == 0
            and _sign(self.witness_value) == self.expected_sign
        )


def _constant_sign_certificate(
    name: str,
    polynomial: Sequence[Fraction | int],
    interval: tuple[Fraction, Fraction],
    expected_sign: int,
    witness: Fraction | None = None,
) -> ConstantSignCertificate:
    normalized = _trim(polynomial)
    common_denominator = lcm(*(coefficient.denominator for coefficient in normalized))
    integer_coefficients = tuple(
        int(coefficient * common_denominator) for coefficient in normalized
    )
    content = 0
    for coefficient in integer_coefficients:
        content = gcd(content, abs(coefficient))
    if content:
        integer_coefficients = tuple(
            coefficient // content for coefficient in integer_coefficients
        )
    normalized = tuple(Fraction(coefficient) for coefficient in integer_coefficients)
    left, right = interval
    if witness is None:
        witness = (left + right) / 2
    root_count = sturm_root_count(normalized, left, right)
    result = ConstantSignCertificate(
        name=name,
        polynomial=integer_coefficients,
        interval=interval,
        witness=Fraction(witness),
        expected_sign=expected_sign,
        root_count=root_count,
        witness_value=_evaluate(normalized, Fraction(witness)),
    )
    if not result.certified:
        raise AssertionError(f"the constant-sign certificate failed for {name}")
    return result


def _profile_polynomials(q: Fraction) -> tuple[Polynomial, Polynomial]:
    """Return exact ``(C,H)`` such that ``c_q(t)=C(t)/H(t)``."""

    first_factor = ((1 - q) * (1 - q), 4 * q)
    second_factor = ((1 + q) * (1 + q), -4 * q)
    denominator = _multiply(first_factor, second_factor)
    bracket = _add(
        _scale(second_factor, 1 + q),
        _scale(first_factor, 1 - q),
    )
    numerator = _multiply((Fraction(0), 2 * q * q), bracket)
    return numerator, denominator


def _integer_profile_polynomials() -> tuple[Polynomial, Polynomial]:
    """Return the historical integer-scaled A0 profile polynomials."""

    first_factor = (Fraction(1521), Fraction(24400))
    second_factor = (Fraction(25921), Fraction(-24400))
    denominator = _multiply(first_factor, second_factor)
    numerator = _multiply(
        (Fraction(0), Fraction(14884)),
        (Fraction(21163), Fraction(-14884)),
    )
    return numerator, denominator


def prior_profile(mass: Fraction | int, q: Fraction = Q_0) -> Fraction:
    """Evaluate the exact triangle fixed-prior profile at channel quality ``q``."""

    mass = Fraction(mass)
    if not 0 <= mass <= 1:
        raise ValueError("mass must belong to [0,1]")
    if mass == 0:
        return Fraction(0)
    numerator, denominator = _profile_polynomials(Fraction(q))
    return _evaluate(numerator, mass) / _evaluate(denominator, mass)


def reduced_diagonal_coefficient(
    mass: Fraction | int,
    candidate: RationalTriangleCandidate = A0_CANDIDATE,
) -> Fraction:
    """Return ``d(t)=(a-delta-c_q(t))/t`` for positive ``t``."""

    mass = Fraction(mass)
    if not 0 < mass <= 1:
        raise ValueError("mass must belong to (0,1]")
    return (candidate.reduced_full - prior_profile(mass, candidate.q)) / mass


def _certificate_polynomials(
    candidate: RationalTriangleCandidate = A0_CANDIDATE,
) -> dict[str, Polynomial]:
    """Derive the four Sturm polynomials from the channel constants."""

    profile_numerator, profile_denominator = _profile_polynomials(candidate.q)

    # H(t) [a-delta+4st-c_q(t)].
    nonpolarized = _subtract(
        _multiply(
            (candidate.reduced_full, 4 * candidate.single),
            profile_denominator,
        ),
        profile_numerator,
    )

    # d(t)=R(t)/(t H(t)).
    d_numerator = _subtract(
        _scale(profile_denominator, candidate.reduced_full),
        profile_numerator,
    )
    d_denominator_core = _multiply((0, 1), profile_denominator)

    # Factor the endpoint zero from d(t)-d(1/2).  The quotient has a
    # constant positive sign on [0,1/2].
    half_core = _evaluate(d_denominator_core, HALF)
    half_numerator = _evaluate(d_numerator, HALF)
    tail_difference_numerator = _subtract(
        _scale(d_numerator, half_core),
        _scale(d_denominator_core, half_numerator),
    )
    tail_quotient, tail_remainder = _divide_polynomials(
        tail_difference_numerator, (1, -2)
    )
    if tail_remainder != (0,):
        raise AssertionError("the tail-separation endpoint factor was lost")

    # Numerator of -d'(t); its denominator is (tH(t))^2 > 0.
    decreasing_numerator = _subtract(
        _multiply(d_numerator, _derivative(d_denominator_core)),
        _multiply(_derivative(d_numerator), d_denominator_core),
    )

    # Numerator of d(t)+s/[t(1-t)] over t H(t)(1-t).
    off_diagonal_numerator = _add(
        _multiply(d_numerator, (1, -1)),
        _scale(profile_denominator, candidate.single),
    )

    return {
        "nonpolarized_affine_gap": nonpolarized,
        "tail_separation": tail_quotient,
        "dominant_decrease": decreasing_numerator,
        "off_diagonal_lower_bound": off_diagonal_numerator,
    }


def polarized_gap_matrix(
    prior: Sequence[Fraction | int],
    candidate: RationalTriangleCandidate = A0_CANDIDATE,
) -> Matrix:
    """Return the exact matrix of ``Q_E-Q_Y-delta*Var``.

    The first prior mass must be dominant and all four masses must be
    strictly positive.  Boundary priors are covered by support restriction in
    the theorem, rather than by inserting a zero into ``d(t)``.
    """

    masses = tuple(Fraction(mass) for mass in prior)
    if len(masses) != 4 or any(mass <= 0 for mass in masses):
        raise ValueError("the matrix requires four strictly positive masses")
    if sum(masses) != 1 or masses[0] <= HALF:
        raise ValueError("require a normalized prior with first mass > 1/2")

    diagonal = tuple(reduced_diagonal_coefficient(mass, candidate) for mass in masses)
    matrix = [[diagonal[0] for _ in range(3)] for _ in range(3)]
    for index in range(3):
        matrix[index][index] += diagonal[index + 1]

    vectors = ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    denominators = (
        (masses[0] + masses[1]) * (masses[2] + masses[3]),
        (masses[0] + masses[2]) * (masses[1] + masses[3]),
        (masses[0] + masses[3]) * (masses[1] + masses[2]),
    )
    for vector, denominator in zip(vectors, denominators, strict=True):
        for row in range(3):
            for column in range(3):
                matrix[row][column] += (
                    candidate.single * vector[row] * vector[column] / denominator
                )
    return tuple(tuple(row) for row in matrix)


def matrix_dominance_residuals(
    prior: Sequence[Fraction | int],
    matrix: Matrix | None = None,
    candidate: RationalTriangleCandidate = A0_CANDIDATE,
) -> tuple[Fraction, Fraction, Fraction]:
    """Return exact Gershgorin residuals ``M_ii-sum_{j!=i}|M_ij|``."""

    masses = tuple(Fraction(mass) for mass in prior)
    if matrix is None:
        matrix = polarized_gap_matrix(masses, candidate)
    return tuple(
        matrix[index][index]
        - sum(abs(matrix[index][other]) for other in range(3) if other != index)
        for index in range(3)
    )


def binary_face_margin(
    candidate: RationalTriangleCandidate = A0_CANDIDATE,
) -> Fraction:
    """Return the exact strengthened margin at the uniform binary prior."""

    squared_quality = candidate.q * candidate.q
    global_contraction = 2 * squared_quality / (1 + squared_quality)
    return candidate.reduced_full + 2 * candidate.single - global_contraction


@dataclass(frozen=True)
class A0LessNoisyCertificate:
    """Manifest of one exhaustive exact rational certificate."""

    candidate_name: str
    status: str
    proof_method: str
    unresolved_regions: int
    p: Fraction
    q: Fraction
    full: Fraction
    single: Fraction
    empty: Fraction
    variance_gap: Fraction
    binary_margin: Fraction
    order_slack: Fraction
    self_dual_slack: Fraction
    fkg_slack: Fraction
    density_squared_slack: Fraction
    nonpolarized: ConstantSignCertificate
    tail_separation: ConstantSignCertificate
    dominant_decrease: ConstantSignCertificate
    off_diagonal_lower_bound: ConstantSignCertificate
    boundary_method: str

    @property
    def exhaustive(self) -> bool:
        return self.status == "CERTIFIED_PSD" and self.unresolved_regions == 0


def rational_less_noisy_certificate(
    candidate: RationalTriangleCandidate,
) -> A0LessNoisyCertificate:
    """Return one exhaustive rational certificate, or raise on failure."""

    if candidate.q != 2 * candidate.p - 1:
        raise AssertionError("q must equal 2p-1")
    if candidate.full + 3 * candidate.single + candidate.empty != 1:
        raise AssertionError("the erasure parameters are not normalized")
    if not all(
        value > 0 for value in (candidate.full, candidate.single, candidate.empty)
    ):
        raise AssertionError("the erasure parameters must be strictly positive")
    if not candidate.full < candidate.empty:
        raise AssertionError("the strict Chayes-Lei order inequality failed")
    if not candidate.full * candidate.empty > 2 * candidate.single**2:
        raise AssertionError("the strict Chayes-Lei inequalities failed")
    if not 2 * candidate.full + 3 * candidate.single < 1:
        raise AssertionError(
            "the channel must lie strictly below the self-dual surface"
        )
    density_squared_slack = (candidate.full + candidate.empty + 8) ** 2 - 72
    if density_squared_slack <= 0:
        raise AssertionError("the Chayes-Lei density inequality failed")

    polynomials = _certificate_polynomials(candidate)
    nonpolarized = _constant_sign_certificate(
        "a-delta+4st-c_q(t)",
        polynomials["nonpolarized_affine_gap"],
        (Fraction(0), HALF),
        1,
    )
    tail_separation = _constant_sign_certificate(
        "d(t)-d(1/2) after removal of 1-2t",
        polynomials["tail_separation"],
        (Fraction(0), HALF),
        1,
    )
    dominant_decrease = _constant_sign_certificate(
        "-d'(t)",
        polynomials["dominant_decrease"],
        (HALF, ONE),
        1,
    )
    off_diagonal = _constant_sign_certificate(
        "d(t)+s/[t(1-t)] numerator",
        polynomials["off_diagonal_lower_bound"],
        (HALF, ONE),
        1,
    )
    binary = binary_face_margin(candidate)
    if binary <= 0:
        raise AssertionError("the binary face has no positive strengthened margin")

    return A0LessNoisyCertificate(
        candidate_name=candidate.name,
        status="CERTIFIED_PSD",
        proof_method="exact Sturm certificates plus diagonal dominance",
        unresolved_regions=0,
        p=candidate.p,
        q=candidate.q,
        full=candidate.full,
        single=candidate.single,
        empty=candidate.empty,
        variance_gap=candidate.variance_gap,
        binary_margin=binary,
        order_slack=candidate.empty - candidate.full,
        self_dual_slack=1 - (2 * candidate.full + 3 * candidate.single),
        fkg_slack=(candidate.full * candidate.empty - 2 * candidate.single**2),
        density_squared_slack=density_squared_slack,
        nonpolarized=nonpolarized,
        tail_separation=tail_separation,
        dominant_decrease=dominant_decrease,
        off_diagonal_lower_bound=off_diagonal,
        boundary_method="support restriction and continuity of finite-channel forms",
    )


def a0_less_noisy_certificate() -> A0LessNoisyCertificate:
    """Return the original exhaustive certificate at ``p=0.805``."""

    return rational_less_noisy_certificate(A0_CANDIDATE)


def p809_less_noisy_certificate() -> A0LessNoisyCertificate:
    """Return the sharper exhaustive certificate at ``p=0.809``."""

    return rational_less_noisy_certificate(P809_CANDIDATE)


def p8094_less_noisy_certificate() -> A0LessNoisyCertificate:
    """Return the exhaustive certificate at ``p=0.8094``."""

    return rational_less_noisy_certificate(P8094_CANDIDATE)


def p809439_less_noisy_certificate() -> A0LessNoisyCertificate:
    """Return the sharpest built-in certificate at ``p=0.809439``."""

    return rational_less_noisy_certificate(P809439_CANDIDATE)


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _print_certificate(certificate: A0LessNoisyCertificate) -> None:
    print(f"candidate: {certificate.candidate_name}")
    print(f"status: {certificate.status}")
    print(f"scope: {'exhaustive' if certificate.exhaustive else 'partial'}")
    print(f"unresolved_regions: {certificate.unresolved_regions}")
    print(f"method: {certificate.proof_method}")
    print(f"p: {_fraction_text(certificate.p)}")
    print(f"q: {_fraction_text(certificate.q)}")
    print(f"variance_gap: {_fraction_text(certificate.variance_gap)}")
    print(f"binary_margin: {_fraction_text(certificate.binary_margin)}")
    print(f"order_slack: {_fraction_text(certificate.order_slack)}")
    print(f"self_dual_slack: {_fraction_text(certificate.self_dual_slack)}")
    print(f"fkg_slack: {_fraction_text(certificate.fkg_slack)}")
    print(
        "density_squared_slack: " f"{_fraction_text(certificate.density_squared_slack)}"
    )
    for item in (
        certificate.nonpolarized,
        certificate.tail_separation,
        certificate.dominant_decrease,
        certificate.off_diagonal_lower_bound,
    ):
        print(
            f"sturm[{item.name}]: roots={item.root_count} "
            f"sign={item.expected_sign:+d} certified={str(item.certified).lower()}"
        )


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate",
        choices=("a0", "p809", "p8094", "p809439", "all"),
        default="a0",
        help="rational candidate to certify (default: a0)",
    )
    selection = parser.parse_args().candidate
    candidates = {
        "a0": (A0_CANDIDATE,),
        "p809": (P809_CANDIDATE,),
        "p8094": (P8094_CANDIDATE,),
        "p809439": (P809439_CANDIDATE,),
        "all": (
            A0_CANDIDATE,
            P809_CANDIDATE,
            P8094_CANDIDATE,
            P809439_CANDIDATE,
        ),
    }[selection]
    for index, candidate in enumerate(candidates):
        if index:
            print()
        _print_certificate(rational_less_noisy_certificate(candidate))


if __name__ == "__main__":
    main()
