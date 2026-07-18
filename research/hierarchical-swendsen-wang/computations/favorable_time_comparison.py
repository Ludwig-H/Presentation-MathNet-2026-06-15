"""Exact local audit of the critical-versus-late comparison.

For a merger bucket of size ``m`` with ``k`` satisfied edges, the local
parity log-likelihood at time ``t`` is

    log(k / (m-k)) + (1-t) u_p (2k-m).

An outside ancestral message is then added.  Earlier times strengthen the
local majority, but they do not always increase the *absolute* total
message: an anti-aligned ancestor can cancel the critical local evidence.
This module gives the exact condition for that failure, a quantitative
one-node bound, and a valid counterexample at ``p=4/5``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, fsum, inf, log, pi, sin, tanh

from ancestral_lambda_chain import (
    AncestorBucket,
    closed_satisfaction_probability,
    coupling,
    local_result,
)


Q_TRIANGLE = 2.0 * sin(pi / 18.0)
RELIABILITY_LIPSCHITZ = 2.0 / (3.0 * (3.0**0.5))


def critical_beta(p: float) -> float:
    """Return the triangular bond-percolation clock level."""

    if Q_TRIANGLE >= 2.0 * p - 1.0:
        raise ValueError("the critical cut must occur strictly before time one")
    return -log(1.0 - Q_TRIANGLE / p) / coupling(p)


def local_parity_log_odds(
    p: float, size: int, satisfied: int, level: float
) -> float:
    """Return the exact local parity message of one unmarked bucket."""

    if size <= 0 or not 0 <= satisfied <= size:
        raise ValueError("require 0 <= satisfied <= size and size > 0")
    if not 0.0 <= level <= 1.0:
        raise ValueError("level must belong to [0, 1]")
    if satisfied == 0:
        return -inf
    if satisfied == size:
        return inf
    return log(satisfied / (size - satisfied)) + (
        (1.0 - level) * coupling(p) * (2 * satisfied - size)
    )


def local_persistence(
    outside_message: float,
    p: float,
    size: int,
    satisfied: int,
    level: float,
) -> float:
    """Squared parity persistence after adding an outside message."""

    local = local_parity_log_odds(p, size, satisfied, level)
    if local in (-inf, inf):
        return 1.0
    return tanh((outside_message + local) / 2.0) ** 2


def anti_alignment_violation(
    outside_message: float,
    p: float,
    size: int,
    satisfied: int,
    early_level: float,
    late_level: float,
) -> bool:
    """Whether the late node is more persistent than the early node.

    For a non-neutral local majority this is equivalent to

        sign(d) B < -(c_early + c_late) / 2,

    where ``d=2*k-m`` and ``c_t`` is the positive local message in the
    majority direction.  Thus every failure is caused by a sufficiently
    strong ancestral message pointing against the local majority.
    """

    if not 0.0 <= early_level <= late_level <= 1.0:
        raise ValueError("require 0 <= early_level <= late_level <= 1")
    difference = 2 * satisfied - size
    if difference == 0 or satisfied in (0, size):
        return False
    direction = 1.0 if difference > 0 else -1.0
    early = direction * local_parity_log_odds(
        p, size, satisfied, early_level
    )
    late = direction * local_parity_log_odds(
        p, size, satisfied, late_level
    )
    return direction * outside_message < -0.5 * (early + late)


def local_comparison_error_bound(
    outside_message: float,
    p: float,
    size: int,
    satisfied: int,
    early_level: float,
    late_level: float,
) -> float:
    """Bound the positive excess of late over early persistence."""

    if not anti_alignment_violation(
        outside_message,
        p,
        size,
        satisfied,
        early_level,
        late_level,
    ):
        return 0.0
    message_gap = (
        coupling(p)
        * abs(2 * satisfied - size)
        * (late_level - early_level)
    )
    return min(1.0, RELIABILITY_LIPSCHITZ * message_gap)


def screened_two_edge_contraction(
    p: float, level: float, message_bound: float
) -> float:
    """Worst mean persistence of a two-edge bucket with ``|B| <= b``."""

    if message_bound < 0.0:
        raise ValueError("message_bound must be nonnegative")
    residual = closed_satisfaction_probability(p, level)
    return residual + (1.0 - residual) * tanh(message_bound / 2.0) ** 2


def bucket_binary_experiment(
    p: float, size: int, level: float
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Return the two mirrored count laws of an unmarked merger bucket."""

    if size <= 0:
        raise ValueError("size must be positive")
    satisfaction = closed_satisfaction_probability(p, level)
    plus = [0.0] * (size + 1)
    minus = [0.0] * (size + 1)
    for count in range(1, size + 1):
        plus[count] = (
            comb(size - 1, count - 1)
            * satisfaction ** (count - 1)
            * (1.0 - satisfaction) ** (size - count)
        )
    for count in range(size):
        minus[count] = (
            comb(size - 1, count)
            * satisfaction ** (size - 1 - count)
            * (1.0 - satisfaction) ** count
        )
    return tuple(plus), tuple(minus)


def bucket_blackwell_minimum_call_gap(
    p: float, size: int, early_level: float, late_level: float
) -> float:
    """Return the finite convex-order certificate for critical dominance.

    For a binary state with uniform prior, one experiment Blackwell-dominates
    another iff its posterior belief is larger in convex order.  For finite
    support it suffices to compare the call functions at every posterior
    value in the union of the two supports.  A nonnegative return value is
    therefore an exact finite-support certificate, up to floating arithmetic.
    """

    if not 0.0 <= early_level <= late_level <= 1.0:
        raise ValueError("require 0 <= early_level <= late_level <= 1")
    early = bucket_binary_experiment(p, size, early_level)
    late = bucket_binary_experiment(p, size, late_level)
    return binary_experiment_blackwell_minimum_call_gap(early, late)


def binary_experiment_blackwell_minimum_call_gap(
    dominant: tuple[tuple[float, ...], tuple[float, ...]],
    degraded: tuple[tuple[float, ...], tuple[float, ...]],
) -> float:
    """Return the convex-order gap for two finite binary experiments."""

    def posterior_law(
        experiment: tuple[tuple[float, ...], tuple[float, ...]],
    ) -> tuple[tuple[float, float], ...]:
        plus, minus = experiment
        if len(plus) != len(minus) or not plus:
            raise ValueError("each experiment needs two rows on one support")
        if min((*plus, *minus)) < 0.0:
            raise ValueError("experiment masses must be nonnegative")
        if abs(fsum(plus) - 1.0) > 1e-12:
            raise ValueError("the plus row must sum to one")
        if abs(fsum(minus) - 1.0) > 1e-12:
            raise ValueError("the minus row must sum to one")
        law = []
        for plus_mass, minus_mass in zip(plus, minus, strict=True):
            total = plus_mass + minus_mass
            if total > 0.0:
                law.append((plus_mass / total, 0.5 * total))
        return tuple(law)

    early_law = posterior_law(dominant)
    late_law = posterior_law(degraded)
    strikes = {0.0, 1.0}
    strikes.update(belief for belief, _ in early_law)
    strikes.update(belief for belief, _ in late_law)
    return min(
        fsum(mass * max(belief - strike, 0.0) for belief, mass in early_law)
        - fsum(mass * max(belief - strike, 0.0) for belief, mass in late_law)
        for strike in strikes
    )


def cross_size_bucket_blackwell_minimum_call_gap(
    p: float,
    dominant_size: int,
    dominant_level: float,
    degraded_size: int,
    degraded_level: float,
) -> float:
    """Test Blackwell dominance when both bucket size and level may change."""

    dominant = bucket_binary_experiment(p, dominant_size, dominant_level)
    degraded = bucket_binary_experiment(p, degraded_size, degraded_level)
    return binary_experiment_blackwell_minimum_call_gap(dominant, degraded)


RationalInterval = tuple[Fraction, Fraction]


def _interval_constant(value: Fraction | int) -> RationalInterval:
    value = Fraction(value)
    return value, value


def _interval_add(
    first: RationalInterval, second: RationalInterval
) -> RationalInterval:
    return first[0] + second[0], first[1] + second[1]


def _interval_subtract(
    first: RationalInterval, second: RationalInterval
) -> RationalInterval:
    return first[0] - second[1], first[1] - second[0]


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
    if denominator[0] <= 0 <= denominator[1]:
        raise ZeroDivisionError("an interval denominator contains zero")
    reciprocal = Fraction(1, denominator[1]), Fraction(1, denominator[0])
    return _interval_multiply(numerator, reciprocal)


def _interval_power(value: RationalInterval, exponent: int) -> RationalInterval:
    if exponent < 0:
        raise ValueError("the interval exponent must be nonnegative")
    result = _interval_constant(1)
    for _ in range(exponent):
        result = _interval_multiply(result, value)
    return result


@dataclass(frozen=True)
class PEightCrossSizeCertificate:
    """Exact rational enclosures proving two experiments incomparable."""

    critical_satisfaction: RationalInterval
    late_satisfaction: RationalInterval
    forward_strike: RationalInterval
    critical_four_vs_late_two_gap: RationalInterval
    late_two_vs_critical_four_gap: RationalInterval


def p_eight_cross_size_incomparability_certificate(
) -> PEightCrossSizeCertificate:
    """Certify a cross-size Blackwell obstruction at ``p=t=4/5``.

    The critical experiment has size four.  The late experiment has size two
    and level ``t=4/5``.  A negative call-function gap in each direction
    proves that neither binary experiment Blackwell-dominates the other.

    Every operation below is on ``Fraction`` intervals.  The critical bond
    threshold ``q`` is enclosed via ``q^3-3q+1=0``; ``y=4^(-1/5)`` is
    enclosed via ``y^5=1/4``.  The displayed decimal endpoints are therefore
    input rational numbers, not floating-point approximations used in the
    proof.
    """

    scale = 10**15
    q_lower = Fraction(347296355333860, scale)
    q_upper = Fraction(347296355333861, scale)
    q_polynomial = lambda value: value**3 - 3 * value + 1
    if not 0 < q_lower < q_upper < 1:
        raise AssertionError("the q_triangle interval must lie in (0, 1)")
    if not q_polynomial(q_lower) > 0 > q_polynomial(q_upper):
        raise AssertionError("the rational interval does not enclose q_triangle")

    y_lower = Fraction(757858283255199, scale)
    y_upper = Fraction(757858283255200, scale)
    if not y_lower**5 < Fraction(1, 4) < y_upper**5:
        raise AssertionError("the rational interval does not enclose 4^(-1/5)")

    one = _interval_constant(1)
    half = _interval_constant(Fraction(1, 2))
    q_interval = q_lower, q_upper
    y_interval = y_lower, y_upper

    critical_satisfaction = _interval_divide(
        _interval_subtract(_interval_constant(Fraction(4, 5)), q_interval),
        _interval_subtract(one, q_interval),
    )
    y_fourth = _interval_power(y_interval, 4)
    late_satisfaction = _interval_divide(
        _interval_multiply(_interval_constant(4), y_fourth),
        _interval_add(
            one, _interval_multiply(_interval_constant(4), y_fourth)
        ),
    )

    # Posterior law of the critical size-four experiment.  Its beliefs are
    # 0, z, 1/2, 1-z, 1 with symmetric masses a0, a1, a2, a1, a0.
    residual = _interval_subtract(one, critical_satisfaction)
    residual_square = _interval_power(residual, 2)
    satisfaction_square = _interval_power(critical_satisfaction, 2)
    posterior_denominator = _interval_add(
        residual_square,
        _interval_multiply(_interval_constant(3), satisfaction_square),
    )
    strike = _interval_divide(residual_square, posterior_denominator)
    if not 0 < strike[0] <= strike[1] < Fraction(1, 2):
        raise AssertionError("the call-function branch assumptions failed")
    mass_a0 = _interval_multiply(
        half, _interval_power(critical_satisfaction, 3)
    )
    mass_a1 = _interval_multiply(
        half, _interval_multiply(residual, posterior_denominator)
    )
    mass_a2 = _interval_multiply(
        _interval_constant(3),
        _interval_multiply(critical_satisfaction, residual_square),
    )

    # The late size-two posterior law has beliefs 0, 1/2, 1 with masses
    # b0, b1, b0.
    mass_b0 = _interval_multiply(half, late_satisfaction)
    mass_b1 = _interval_subtract(one, late_satisfaction)

    critical_call_at_strike = _interval_add(
        _interval_add(
            _interval_multiply(
                mass_a2, _interval_subtract(half, strike)
            ),
            _interval_multiply(
                mass_a1,
                _interval_subtract(
                    one, _interval_multiply(_interval_constant(2), strike)
                ),
            ),
        ),
        _interval_multiply(mass_a0, _interval_subtract(one, strike)),
    )
    late_call_at_strike = _interval_add(
        _interval_multiply(mass_b1, _interval_subtract(half, strike)),
        _interval_multiply(mass_b0, _interval_subtract(one, strike)),
    )
    forward_gap = _interval_subtract(
        critical_call_at_strike, late_call_at_strike
    )

    late_call_at_half = _interval_multiply(mass_b0, half)
    critical_call_at_half = _interval_add(
        _interval_multiply(mass_a1, _interval_subtract(half, strike)),
        _interval_multiply(mass_a0, half),
    )
    reverse_gap = _interval_subtract(
        late_call_at_half, critical_call_at_half
    )

    if not (forward_gap[1] < 0 and reverse_gap[1] < 0):
        raise AssertionError("the rational intervals do not prove incomparability")
    return PEightCrossSizeCertificate(
        critical_satisfaction=critical_satisfaction,
        late_satisfaction=late_satisfaction,
        forward_strike=strike,
        critical_four_vs_late_two_gap=forward_gap,
        late_two_vs_critical_four_gap=reverse_gap,
    )


def two_edge_blackwell_erasure_probability(
    p: float, early_level: float, late_level: float
) -> float:
    """Extra erasure that degrades an early two-edge count experiment."""

    if not 0.0 <= early_level <= late_level <= 1.0:
        raise ValueError("require 0 <= early_level <= late_level <= 1")
    early = closed_satisfaction_probability(p, early_level)
    late = closed_satisfaction_probability(p, late_level)
    return 1.0 - late / early


@dataclass(frozen=True)
class PEightCounterexample:
    """A valid one-ancestor failure of pointwise time monotonicity."""

    critical_level: float
    late_level: float
    critical_log_odds: float
    late_log_odds: float
    critical_persistence: float
    late_persistence: float


def p_eight_counterexample() -> PEightCounterexample:
    """Return an explicit admissible counterexample at ``p=4/5``."""

    p = 4.0 / 5.0
    rate = coupling(p)
    ancestor = AncestorBucket(
        beta=0.81,
        total=(rate, 6.0 * rate, 6.0 * rate),
        satisfied=(rate, 0.0, 6.0 * rate),
    )
    critical_level = critical_beta(p)
    late_level = 0.8
    critical = local_result(
        2.0 * rate,
        3.0 * rate,
        critical_level,
        (ancestor,),
    )
    late = local_result(
        2.0 * rate,
        3.0 * rate,
        late_level,
        (ancestor,),
    )
    return PEightCounterexample(
        critical_level=critical_level,
        late_level=late_level,
        critical_log_odds=critical.log_odds,
        late_log_odds=late.log_odds,
        critical_persistence=critical.reliability,
        late_persistence=late.reliability,
    )


def main() -> None:
    p = 4.0 / 5.0
    example = p_eight_counterexample()
    print(example)
    for bound in (0.0, 0.5, 1.0, 1.5, 2.0):
        contraction = screened_two_edge_contraction(
            p, example.critical_level, bound
        )
        print(f"message_bound={bound:.1f} contraction={contraction:.12f}")
    certificate = p_eight_cross_size_incomparability_certificate()
    print(
        "critical m=4 vs late m=2 call gap in "
        f"[{float(certificate.critical_four_vs_late_two_gap[0]):.12g}, "
        f"{float(certificate.critical_four_vs_late_two_gap[1]):.12g}]"
    )
    print(
        "late m=2 vs critical m=4 call gap in "
        f"[{float(certificate.late_two_vs_critical_four_gap[0]):.12g}, "
        f"{float(certificate.late_two_vs_critical_four_gap[1]):.12g}]"
    )


if __name__ == "__main__":
    main()
