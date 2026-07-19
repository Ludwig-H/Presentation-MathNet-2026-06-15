"""Exact two-port obstruction to naive corridor criticalization.

Consider a merger bucket containing two labelled physical edges.  Conditional
on the unmarked merger time, one winning edge is uniform and certainly
satisfied.  Its identity is *not* observed.  If ``x in {+1,-1}^2`` is the
vector of true edge relations and ``y`` the labelled residual signs, the
channel is

    choose G uniformly in {0,1};
    set y_G = x_G;
    send the other coordinate through a BSC with success probability s.

After marginalizing ``G`` this is an additive channel on the group
``{+1,-1}^2``.  Its noise law is

    P_s(++ ) = s,
    P_s(+-) = P_s(-+) = (1-s)/2,
    P_s(-- ) = 0.

For a single binary bucket state, every edge relation is complemented
together and the scalar count experiment is Blackwell-monotone in ``s``.
The present four-state experiment is different: descendant orientations can
change the two labelled relations separately.  It has no degradation from an
earlier ``s_c`` to any strictly later ``s_l`` in ``(1/2, s_c)``.

Indeed, any putative post-processing can be averaged over simultaneous group
translations of its input and output.  The averaged kernel remains stochastic
and remains a degradation because both channels are translation covariant.
It is therefore convolution by a probability ``R``.  Fourier division makes
``R`` unique, and Fourier inversion gives

    R(-- ) = (1 - 2 s_l/s_c + h_l/h_c)/4
            = (1-h_c)(h_l-h_c)/(4 h_c (1+h_c)) < 0,

where ``h=2s-1``.  This module certifies that sign at ``p=161/200`` and also
certifies a posterior-variance reversal for an exact rational prior and the
character ``f(x)=x_1*x_2``.  The critical rank is enclosed using rational
endpoints and ``q^3-3q+1=0``; all certificate operations use ``Fraction``
intervals.

This obstruction concerns a marginalized winner.  Revealing the winner would
permit coordinatewise degradation, but it would also reveal a satisfied edge
at every Kruskal merger and is therefore a much stronger, generally useless
oracle for a decorrelation proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

State = tuple[int, int]
Matrix = tuple[tuple[Fraction, ...], ...]
Interval = tuple[Fraction, Fraction]
IntervalMatrix = tuple[tuple[Interval, ...], ...]

STATES: tuple[State, ...] = (
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1),
)

P_805 = Fraction(161, 200)
LATE_RANK = Fraction(11, 20)
LATE_SATISFACTION = Fraction(17, 30)
PRIOR: tuple[Fraction, ...] = (
    Fraction(19, 20),
    Fraction(1, 50),
    Fraction(1, 40),
    Fraction(1, 200),
)
TARGET: tuple[int, ...] = tuple(first * second for first, second in STATES)


def _multiply_states(first: State, second: State) -> State:
    return first[0] * second[0], first[1] * second[1]


def _interval(value: Fraction | int) -> Interval:
    scalar = Fraction(value)
    return scalar, scalar


def _add(first: Interval, second: Interval) -> Interval:
    return first[0] + second[0], first[1] + second[1]


def _subtract(first: Interval, second: Interval) -> Interval:
    return first[0] - second[1], first[1] - second[0]


def _multiply(first: Interval, second: Interval) -> Interval:
    products = tuple(
        first[index] * second[other] for index in (0, 1) for other in (0, 1)
    )
    return min(products), max(products)


def _divide(numerator: Interval, denominator: Interval) -> Interval:
    if denominator[0] <= 0 <= denominator[1]:
        raise ZeroDivisionError("an interval denominator contains zero")
    reciprocal = Fraction(1, denominator[1]), Fraction(1, denominator[0])
    return _multiply(numerator, reciprocal)


def _square(value: Interval) -> Interval:
    upper = max(value[0] * value[0], value[1] * value[1])
    lower = (
        Fraction(0)
        if value[0] <= 0 <= value[1]
        else min(value[0] * value[0], value[1] * value[1])
    )
    return lower, upper


def _sum_intervals(values: Sequence[Interval]) -> Interval:
    answer = _interval(0)
    for value in values:
        answer = _add(answer, value)
    return answer


def marginalized_winner_noise(satisfaction: Fraction) -> tuple[Fraction, ...]:
    """Return ``P_s`` in ``STATES`` order, with the winner marginalized."""

    satisfaction = Fraction(satisfaction)
    if not Fraction(1, 2) <= satisfaction <= 1:
        raise ValueError("satisfaction must belong to [1/2, 1]")
    half_failure = (1 - satisfaction) / 2
    return satisfaction, half_failure, half_failure, Fraction(0)


def marginalized_winner_channel(satisfaction: Fraction) -> Matrix:
    """Return the exact four-state channel ``W_s(y|x)=P_s(y*x)``."""

    noise = marginalized_winner_noise(satisfaction)
    noise_index = {state: index for index, state in enumerate(STATES)}
    return tuple(
        tuple(noise[noise_index[_multiply_states(output, latent)]] for output in STATES)
        for latent in STATES
    )


def winner_noise_fourier(
    satisfaction: Fraction,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Return Fourier coefficients ``(1, chi_1, chi_2, chi_1*chi_2)``."""

    satisfaction = Fraction(satisfaction)
    marginalized_winner_noise(satisfaction)
    return (
        Fraction(1),
        satisfaction,
        satisfaction,
        2 * satisfaction - 1,
    )


def symmetrize_postprocessing(kernel: Matrix) -> Matrix:
    """Average a four-state post-processing over group conjugations.

    If ``G`` degrades one translation-covariant channel into another, this
    averaged kernel does too.  Its rows are translates of a single noise law,
    so it is a convolution kernel.
    """

    if len(kernel) != 4 or any(len(row) != 4 for row in kernel):
        raise ValueError("kernel must be a 4 by 4 matrix")
    if any(entry < 0 for row in kernel for entry in row):
        raise ValueError("kernel entries must be nonnegative")
    if any(sum(row) != 1 for row in kernel):
        raise ValueError("kernel rows must sum exactly to one")

    index = {state: position for position, state in enumerate(STATES)}
    rows = []
    for source in STATES:
        row = []
        for target in STATES:
            row.append(
                sum(
                    kernel[index[_multiply_states(shift, source)]][
                        index[_multiply_states(shift, target)]
                    ]
                    for shift in STATES
                )
                / 4
            )
        rows.append(tuple(row))
    return tuple(rows)


def candidate_covariant_degradation_noise(
    dominant_satisfaction: Fraction,
    degraded_satisfaction: Fraction,
) -> tuple[Fraction, ...]:
    """Fourier-invert the unique candidate convolution degradation.

    A negative entry proves that no stochastic degradation exists, including
    a non-covariant one, because any such degradation could be symmetrized.
    """

    dominant = Fraction(dominant_satisfaction)
    degraded = Fraction(degraded_satisfaction)
    if not Fraction(1, 2) < degraded < dominant < 1:
        raise ValueError("require 1/2 < degraded < dominant < 1")
    dominant_margin = 2 * dominant - 1
    degraded_margin = 2 * degraded - 1
    first_ratio = degraded / dominant
    product_ratio = degraded_margin / dominant_margin
    return tuple(
        (
            1
            + first_ratio * state[0]
            + first_ratio * state[1]
            + product_ratio * state[0] * state[1]
        )
        / 4
        for state in STATES
    )


def posterior_variance(
    channel: Matrix,
    prior: Sequence[Fraction],
    target: Sequence[int | Fraction],
) -> Fraction:
    """Return ``Var(E[f(X)|Y])`` under an exact finite channel."""

    if len(channel) != 4 or any(len(row) != 4 for row in channel):
        raise ValueError("channel must be a 4 by 4 matrix")
    prior_values = tuple(Fraction(value) for value in prior)
    target_values = tuple(Fraction(value) for value in target)
    if len(prior_values) != 4 or len(target_values) != 4:
        raise ValueError("prior and target must have four entries")
    if any(mass < 0 for mass in prior_values) or sum(prior_values) != 1:
        raise ValueError("prior must be a probability vector")
    if any(entry < 0 for row in channel for entry in row):
        raise ValueError("channel entries must be nonnegative")
    if any(sum(row) != 1 for row in channel):
        raise ValueError("channel rows must sum exactly to one")

    mean = sum(
        mass * value for mass, value in zip(prior_values, target_values, strict=True)
    )
    centered_masses = tuple(
        mass * (value - mean)
        for mass, value in zip(prior_values, target_values, strict=True)
    )
    answer = Fraction(0)
    for output in range(4):
        output_mass = sum(
            prior_values[latent] * channel[latent][output] for latent in range(4)
        )
        numerator = sum(
            centered_masses[latent] * channel[latent][output] for latent in range(4)
        )
        if output_mass <= 0:
            raise AssertionError("the selected prior gives a zero output mass")
        answer += numerator * numerator / output_mass
    return answer


def _interval_channel(satisfaction: Interval) -> IntervalMatrix:
    failure = _subtract(_interval(1), satisfaction)
    half_failure = _divide(failure, _interval(2))
    noise = satisfaction, half_failure, half_failure, _interval(0)
    noise_index = {state: index for index, state in enumerate(STATES)}
    return tuple(
        tuple(noise[noise_index[_multiply_states(output, latent)]] for output in STATES)
        for latent in STATES
    )


def _posterior_variance_interval(
    channel: IntervalMatrix,
    prior: Sequence[Fraction],
    target: Sequence[int | Fraction],
) -> Interval:
    prior_values = tuple(Fraction(value) for value in prior)
    target_values = tuple(Fraction(value) for value in target)
    mean = sum(
        mass * value for mass, value in zip(prior_values, target_values, strict=True)
    )
    centered_masses = tuple(
        mass * (value - mean)
        for mass, value in zip(prior_values, target_values, strict=True)
    )
    terms = []
    for output in range(4):
        output_mass = _sum_intervals(
            tuple(
                _multiply(_interval(prior_values[latent]), channel[latent][output])
                for latent in range(4)
            )
        )
        numerator = _sum_intervals(
            tuple(
                _multiply(
                    _interval(centered_masses[latent]),
                    channel[latent][output],
                )
                for latent in range(4)
            )
        )
        terms.append(_divide(_square(numerator), output_mass))
    return _sum_intervals(tuple(terms))


@dataclass(frozen=True)
class MultiportBlackwellCounterexample:
    """Rational-interval certificate at ``p=161/200``."""

    p: Fraction
    critical_rank: Interval
    late_rank: Fraction
    critical_satisfaction: Interval
    late_satisfaction: Fraction
    critical_margin: Interval
    late_margin: Fraction
    candidate_noise_minus_minus: Interval
    prior: tuple[Fraction, ...]
    target: tuple[int, ...]
    target_mean: Fraction
    target_variance: Fraction
    critical_posterior_variance: Interval
    late_posterior_variance: Fraction
    variance_reversal_gap: Interval
    critical_posterior_second_moment: Interval
    late_posterior_second_moment: Fraction

    @property
    def certifies_no_blackwell_degradation(self) -> bool:
        return self.candidate_noise_minus_minus[1] < 0

    @property
    def certifies_posterior_variance_reversal(self) -> bool:
        return self.variance_reversal_gap[0] > 0


def p_805_multiport_blackwell_counterexample() -> MultiportBlackwellCounterexample:
    """Build the exact interval certificate for the documented example."""

    scale = 10**15
    critical_rank = (
        Fraction(347296355333860, scale),
        Fraction(347296355333861, scale),
    )
    polynomial = lambda value: value**3 - 3 * value + 1
    if not polynomial(critical_rank[0]) > 0 > polynomial(critical_rank[1]):
        raise AssertionError("the rational interval does not enclose q_c")

    critical_satisfaction = (
        (P_805 - critical_rank[1]) / (1 - critical_rank[1]),
        (P_805 - critical_rank[0]) / (1 - critical_rank[0]),
    )
    critical_margin = _subtract(
        _multiply(_interval(2), critical_satisfaction), _interval(1)
    )
    late_margin = 2 * LATE_SATISFACTION - 1

    # Fourier inversion of the unique covariant degradation candidate.
    first_ratio = _divide(_interval(LATE_SATISFACTION), critical_satisfaction)
    product_ratio = _divide(_interval(late_margin), critical_margin)
    candidate_minus_minus = _divide(
        _add(
            _subtract(_interval(1), _multiply(_interval(2), first_ratio)),
            product_ratio,
        ),
        _interval(4),
    )

    target_mean = sum(mass * value for mass, value in zip(PRIOR, TARGET, strict=True))
    target_variance = sum(
        mass * (value - target_mean) ** 2
        for mass, value in zip(PRIOR, TARGET, strict=True)
    )
    late_variance = posterior_variance(
        marginalized_winner_channel(LATE_SATISFACTION), PRIOR, TARGET
    )
    critical_variance = _posterior_variance_interval(
        _interval_channel(critical_satisfaction), PRIOR, TARGET
    )
    reversal_gap = _subtract(_interval(late_variance), critical_variance)
    critical_second_moment = _add(
        _interval(target_mean * target_mean), critical_variance
    )
    late_second_moment = target_mean * target_mean + late_variance

    certificate = MultiportBlackwellCounterexample(
        p=P_805,
        critical_rank=critical_rank,
        late_rank=LATE_RANK,
        critical_satisfaction=critical_satisfaction,
        late_satisfaction=LATE_SATISFACTION,
        critical_margin=critical_margin,
        late_margin=late_margin,
        candidate_noise_minus_minus=candidate_minus_minus,
        prior=PRIOR,
        target=TARGET,
        target_mean=target_mean,
        target_variance=target_variance,
        critical_posterior_variance=critical_variance,
        late_posterior_variance=late_variance,
        variance_reversal_gap=reversal_gap,
        critical_posterior_second_moment=critical_second_moment,
        late_posterior_second_moment=late_second_moment,
    )
    if not certificate.certifies_no_blackwell_degradation:
        raise AssertionError("the candidate degradation noise was not negative")
    if not certificate.certifies_posterior_variance_reversal:
        raise AssertionError("the posterior-variance reversal was not certified")
    return certificate


def _decimal(value: Fraction, places: int = 15) -> str:
    return f"{float(value):.{places}f}"


def _decimal_interval(value: Interval, places: int = 15) -> str:
    return f"[{_decimal(value[0], places)}, {_decimal(value[1], places)}]"


def main() -> None:
    certificate = p_805_multiport_blackwell_counterexample()
    print("two-port marginalized-winner Blackwell counterexample")
    print(f"states={STATES}")
    print(f"p={certificate.p} late rank={certificate.late_rank}")
    print(
        "critical satisfaction="
        f"{_decimal_interval(certificate.critical_satisfaction)}"
    )
    print(f"late satisfaction={certificate.late_satisfaction}")
    print(
        "unique covariant degradation candidate R(-- )="
        f"{_decimal_interval(certificate.candidate_noise_minus_minus)}"
    )
    print(
        "negative upper endpoint proves that no Blackwell degradation exists: "
        f"{certificate.certifies_no_blackwell_degradation}"
    )
    print(f"prior={certificate.prior}")
    print(f"target x1*x2={certificate.target}")
    print(
        f"target mean={certificate.target_mean} "
        f"variance={certificate.target_variance}"
    )
    print(
        "critical Var(E[f|Y])="
        f"{_decimal_interval(certificate.critical_posterior_variance)}"
    )
    print("late Var(E[f|Y])=" f"{_decimal(certificate.late_posterior_variance)}")
    print(
        "late-minus-critical gap="
        f"{_decimal_interval(certificate.variance_reversal_gap)}"
    )
    print(
        "critical E[E[f|Y]^2]="
        f"{_decimal_interval(certificate.critical_posterior_second_moment)}"
    )
    print("late E[E[f|Y]^2]=" f"{_decimal(certificate.late_posterior_second_moment)}")
    print(
        "posterior-variance reversal certified: "
        f"{certificate.certifies_posterior_variance_reversal}"
    )
    print(
        "scope: the scalar global-complement bucket lemma remains valid; "
        "the obstruction starts when two edge relations vary independently"
    )


if __name__ == "__main__":
    main()
