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

    def posterior_law(
        experiment: tuple[tuple[float, ...], tuple[float, ...]],
    ) -> tuple[tuple[float, float], ...]:
        plus, minus = experiment
        law = []
        for plus_mass, minus_mass in zip(plus, minus, strict=True):
            total = plus_mass + minus_mass
            if total > 0.0:
                law.append((plus_mass / total, 0.5 * total))
        return tuple(law)

    early_law = posterior_law(early)
    late_law = posterior_law(late)
    strikes = {0.0, 1.0}
    strikes.update(belief for belief, _ in early_law)
    strikes.update(belief for belief, _ in late_law)
    return min(
        fsum(mass * max(belief - strike, 0.0) for belief, mass in early_law)
        - fsum(mass * max(belief - strike, 0.0) for belief, mass in late_law)
        for strike in strikes
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


if __name__ == "__main__":
    main()
