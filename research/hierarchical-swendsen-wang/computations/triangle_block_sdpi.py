"""Exact local audits for the triangular three-terminal observation channel.

The module deliberately uses only the Python standard library.  It separates
proved channel identities from the conditional multi-state erasure candidate
described in ``11_TRIANGLE_BLOCK_SDPI.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import pi, sin, sqrt
from typing import Iterable, Sequence


Q_BOND_CRITICAL = 2.0 * sin(pi / 18.0)
P_INFORMATION = (1.0 + sqrt(Q_BOND_CRITICAL)) / 2.0
P_SCALAR_TRIANGLE = (1.0 + 1.0 / sqrt(3.0)) / 2.0
P_FIXED_PRIOR_TRIANGLE = (
    1.0 + sqrt((sqrt(13.0) - 1.0) / 6.0)
) / 2.0


def _check_q(q: float) -> None:
    if not 0.0 <= q < 1.0:
        raise ValueError("q must satisfy 0 <= q < 1")


def triangle_uniform_contraction(q: float) -> float:
    """Squared correlation for one isolated noisy triangle."""

    _check_q(q)
    square = q * q
    return square * (1.0 + 2.0 * square) / (
        1.0 + square + square * square
    )


def two_observation_contraction(q: float) -> float:
    """SDPI coefficient of two independent BSC observations of one bit."""

    _check_q(q)
    square = q * q
    return 2.0 * square / (1.0 + square)


def prior_profile(q: float, mass: float) -> float:
    """Diagonal fixed-prior chi-square profile of the triangle channel.

    If ``mu`` is a prior on the four relative-spin states and ``f`` has
    ``mu``-mean zero, the exact output variance is

        sum_x mu[x] * prior_profile(q, mu[x]) * f[x]**2.

    The formula remains useful at masses above one half, even though an
    individual diagonal coefficient can then exceed one: the mean-zero
    constraint prevents that coordinate from being used in isolation.
    """

    _check_q(q)
    if not 0.0 <= mass <= 1.0:
        raise ValueError("mass must satisfy 0 <= mass <= 1")
    if mass == 0.0 or q == 0.0:
        return 0.0
    first = (1.0 - q) ** 2 + 4.0 * q * mass
    second = (1.0 + q) ** 2 - 4.0 * q * mass
    return 2.0 * q * q * mass * (
        (1.0 + q) / first + (1.0 - q) / second
    )


def triangle_channel(q: float) -> tuple[tuple[float, ...], ...]:
    """Return the 4-by-8 transition matrix of the relative-spin channel."""

    _check_q(q)
    correct = (1.0 + q) / 2.0
    flipped = (1.0 - q) / 2.0
    inputs = tuple(product((1, -1), repeat=2))
    outputs = tuple(product((1, -1), repeat=3))
    rows: list[tuple[float, ...]] = []
    for first, second in inputs:
        codeword = (first, second, first * second)
        row = []
        for output in outputs:
            probability = 1.0
            for observed, truth in zip(output, codeword):
                probability *= correct if observed == truth else flipped
            row.append(probability)
        rows.append(tuple(row))
    return tuple(rows)


def direct_chi_square_ratio(
    q: float,
    prior: Sequence[float],
    function: Sequence[float],
) -> float:
    """Enumerate Var(E[f(X)|Y]) / Var(f(X)) for the triangle channel."""

    if len(prior) != 4 or len(function) != 4:
        raise ValueError("prior and function must have four entries")
    if any(mass < 0.0 for mass in prior):
        raise ValueError("prior masses must be nonnegative")
    if abs(sum(prior) - 1.0) > 1e-12:
        raise ValueError("prior must sum to one")
    mean = sum(mass * value for mass, value in zip(prior, function))
    centered = tuple(value - mean for value in function)
    variance = sum(
        mass * value * value for mass, value in zip(prior, centered)
    )
    if variance == 0.0:
        raise ValueError("function must have positive variance")

    channel = triangle_channel(q)
    output_variance = 0.0
    for output in range(8):
        output_mass = sum(prior[x] * channel[x][output] for x in range(4))
        numerator = sum(
            prior[x] * centered[x] * channel[x][output] for x in range(4)
        )
        output_variance += numerator * numerator / output_mass
    return output_variance / variance


def profile_chi_square_ratio(
    q: float,
    prior: Sequence[float],
    function: Sequence[float],
) -> float:
    """Evaluate the closed fixed-prior profile for comparison with enumeration."""

    mean = sum(mass * value for mass, value in zip(prior, function))
    centered = tuple(value - mean for value in function)
    variance = sum(
        mass * value * value for mass, value in zip(prior, centered)
    )
    if variance == 0.0:
        raise ValueError("function must have positive variance")
    numerator = sum(
        mass * prior_profile(q, mass) * value * value
        for mass, value in zip(prior, centered)
    )
    return numerator / variance


def projection_variance_sum(
    prior: Sequence[float], function: Sequence[float]
) -> float:
    """Sum variances after revealing each of the three nonzero characters."""

    if len(prior) != 4 or len(function) != 4:
        raise ValueError("prior and function must have four entries")
    states = tuple(product((1, -1), repeat=2))
    characters = tuple(
        tuple((first, second, first * second)[index] for first, second in states)
        for index in range(3)
    )
    mean = sum(mass * value for mass, value in zip(prior, function))
    centered = tuple(value - mean for value in function)
    total = 0.0
    for character in characters:
        positive_mass = sum(
            prior[index] for index, sign in enumerate(character) if sign == 1
        )
        negative_mass = 1.0 - positive_mass
        if positive_mass == 0.0 or negative_mass == 0.0:
            continue
        positive_sum = sum(
            prior[index] * centered[index]
            for index, sign in enumerate(character)
            if sign == 1
        )
        total += positive_sum * positive_sum / (
            positive_mass * negative_mass
        )
    return total


def projection_diagonal_lower_bound(
    prior: Sequence[float], function: Sequence[float]
) -> float:
    """The exact elementary lower bound used in the erasure comparison."""

    mean = sum(mass * value for mass, value in zip(prior, function))
    return 4.0 * sum(
        mass * mass * (value - mean) ** 2
        for mass, value in zip(prior, function)
    )


def less_noisy_gap_matrix(
    q: float,
    prior: Sequence[float],
    full: float,
    single: float,
) -> tuple[tuple[float, ...], ...]:
    """Return the exact 3-by-3 quadratic gap matrix in likelihood coordinates.

    Write ``g[x] = prior[x] * f[x]`` and eliminate ``g[0]`` using
    ``sum(g) = 0``.  The returned matrix represents ``Q_E - Q_Y`` in the
    remaining coordinates ``(g[1], g[2], g[3])``.  It is useful for finite
    counterexample searches; sampling priors does not turn the audit into a
    proof of the less-noisy comparison.
    """

    if len(prior) != 4 or any(mass <= 0.0 for mass in prior):
        raise ValueError("prior must have four strictly positive entries")
    if abs(sum(prior) - 1.0) > 1e-12:
        raise ValueError("prior must sum to one")
    diagonal = tuple(
        (full - prior_profile(q, mass)) / mass for mass in prior
    )
    matrix = [[diagonal[0] for _ in range(3)] for _ in range(3)]
    for index in range(3):
        matrix[index][index] += diagonal[index + 1]

    projection_vectors = ((0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 0.0))
    projection_denominators = (
        (prior[0] + prior[1]) * (prior[2] + prior[3]),
        (prior[0] + prior[2]) * (prior[1] + prior[3]),
        (prior[0] + prior[3]) * (prior[1] + prior[2]),
    )
    for vector, denominator in zip(
        projection_vectors, projection_denominators
    ):
        for row in range(3):
            for column in range(3):
                matrix[row][column] += (
                    single
                    * vector[row]
                    * vector[column]
                    / denominator
                )
    return tuple(tuple(row) for row in matrix)


def envelope_slope(q: float, mass: float) -> float:
    """Slope constraint for a line anchored at mass one half."""

    if not 0.0 <= mass < 0.5:
        raise ValueError("mass must satisfy 0 <= mass < 1/2")
    return (
        two_observation_contraction(q) - prior_profile(q, mass)
    ) / (2.0 - 4.0 * mass)


@dataclass(frozen=True)
class ErasureEnvelope:
    full: float
    single: float
    empty: float
    tangent_mass: float

    @property
    def self_dual_score(self) -> float:
        return 2.0 * self.full + 3.0 * self.single


def erasure_envelope(q: float) -> ErasureEnvelope:
    """Optimize the necessary affine profile of the three-state erasure model.

    This calculation does *not* prove the remaining polarized-prior less-noisy
    lemma.  It only returns the sharp affine candidate used in the note.
    """

    _check_q(q)
    if q == 0.0:
        return ErasureEnvelope(0.0, 0.0, 1.0, 0.0)
    square = q * q
    tangent = 0.5 - (1.0 + square) / (
        4.0 * (1.0 + sqrt(1.0 - square))
    )
    single = envelope_slope(q, tangent)
    global_coefficient = two_observation_contraction(q)
    full = global_coefficient - 2.0 * single
    empty = 1.0 - full - 3.0 * single
    return ErasureEnvelope(full, single, empty, tangent)


def conditional_candidate_q() -> float:
    """Solve the conditional self-dual equation 2*a(q)+3*s(q)=1."""

    left = 0.60
    right = 0.64
    for _ in range(70):
        middle = (left + right) / 2.0
        if erasure_envelope(middle).self_dual_score < 1.0:
            left = middle
        else:
            right = middle
    return (left + right) / 2.0


def conditional_candidate_polynomial(q: float) -> float:
    """Polynomial vanishing at the interior-tangency candidate."""

    square = q * q
    return (
        square**5
        + 46.0 * square**4
        + 45.0 * square**3
        - 20.0 * square**2
        - 12.0 * square
        + 4.0
    )


def sampled_line_defect(q: float, full: float, single: float) -> float:
    """Minimum of a+4*s*t-c_q(t) on a reproducible fine grid."""

    return min(
        full + 4.0 * single * mass - prior_profile(q, mass)
        for mass in (index / 20000.0 for index in range(10001))
    )


def _format_rows(rows: Iterable[tuple[str, float]]) -> str:
    return "\n".join(f"{label:<34} {value:.12f}" for label, value in rows)


def main() -> None:
    q_information = 2.0 * P_INFORMATION - 1.0
    candidate_q = conditional_candidate_q()
    candidate = erasure_envelope(candidate_q)
    naive_q = 2.0 * P_FIXED_PRIOR_TRIANGLE - 1.0

    print(
        _format_rows(
            (
                ("p_information", P_INFORMATION),
                ("p_scalar_triangle", P_SCALAR_TRIANGLE),
                ("p_fixed_prior_naive", P_FIXED_PRIOR_TRIANGLE),
                ("eta_triangle(q_information)", triangle_uniform_contraction(q_information)),
                ("gamma_2(q_information)", two_observation_contraction(q_information)),
                ("eta_triangle(naive_q)", triangle_uniform_contraction(naive_q)),
                ("conditional_candidate_q", candidate_q),
                ("conditional_candidate_p", (1.0 + candidate_q) / 2.0),
                ("candidate_a", candidate.full),
                ("candidate_s", candidate.single),
                ("candidate_e", candidate.empty),
                ("candidate_tangent_mass", candidate.tangent_mass),
                ("candidate_polynomial", conditional_candidate_polynomial(candidate_q)),
            )
        )
    )
    print("\nSTATUS: the 0.809909... value is conditional on the polarized-prior")
    print("less-noisy comparison; it is not a proved weak-recovery bound.")


if __name__ == "__main__":
    main()
