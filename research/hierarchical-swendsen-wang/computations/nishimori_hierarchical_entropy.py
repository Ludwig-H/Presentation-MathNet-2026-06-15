"""Exact entropy identities behind the triangular Nishimori conjecture.

The module is deliberately dependency-free.  It checks that equation (28) of
Nishimori--Ohzeki (2006) is exactly the conditional entropy balance

    H(Z_1, Z_2, Z_3 | Z_1 Z_2 Z_3) = 1 bit,

and exposes the equivalent four-state exponential-race law.  This recovers
the conjectured number, not a proof that it is the weak-recovery threshold.
"""

from __future__ import annotations

from itertools import product
from math import atanh, log, log2
from typing import Iterable


NoiseWord = tuple[int, int, int]


def _check_probability(p: float, *, strict: bool = False) -> None:
    lower_ok = p > 0.0 if strict else p >= 0.0
    upper_ok = p < 1.0 if strict else p <= 1.0
    if not (lower_ok and upper_ok):
        interval = "0 < p < 1" if strict else "0 <= p <= 1"
        raise ValueError(f"p must satisfy {interval}")


def binary_entropy_bits(p: float) -> float:
    """Return the binary entropy h_2(p), with the endpoint convention 0 log 0=0."""

    _check_probability(p)
    if p in (0.0, 1.0):
        return 0.0
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p)


def face_syndrome_probability(p: float) -> float:
    """Return P(Z_1 Z_2 Z_3=+1) for iid signs with P(Z_e=+1)=p."""

    _check_probability(p)
    return p**3 + 3.0 * p * (1.0 - p) ** 2


def face_noise_probability(word: NoiseWord, p: float) -> float:
    """Return the probability of one three-edge noise word."""

    _check_probability(p)
    if len(word) != 3 or any(sign not in (-1, 1) for sign in word):
        raise ValueError("word must contain exactly three signs in {-1,+1}")
    plus = sum(sign == 1 for sign in word)
    return p**plus * (1.0 - p) ** (3 - plus)


def conditional_face_noise_laws(
    p: float,
) -> dict[int, dict[NoiseWord, float]]:
    """Return the two four-state laws of Z conditioned on its face syndrome."""

    _check_probability(p, strict=True)
    joint: dict[int, dict[NoiseWord, float]] = {1: {}, -1: {}}
    for word in product((-1, 1), repeat=3):
        typed_word: NoiseWord = word
        syndrome = word[0] * word[1] * word[2]
        joint[syndrome][typed_word] = face_noise_probability(typed_word, p)
    for syndrome in (-1, 1):
        mass = sum(joint[syndrome].values())
        joint[syndrome] = {
            word: probability / mass
            for word, probability in joint[syndrome].items()
        }
    return joint


def face_residual_entropy_bits(p: float) -> float:
    """Return H(Z_1,Z_2,Z_3 | Z_1 Z_2 Z_3), in bits."""

    _check_probability(p)
    syndrome_plus = face_syndrome_probability(p)
    return 3.0 * binary_entropy_bits(p) - binary_entropy_bits(syndrome_plus)


def nishimori_entropy_balance_bits(p: float) -> float:
    """Return H(Z_1,Z_2,Z_3 | product)-1; its upper-half zero is p_N."""

    return face_residual_entropy_bits(p) - 1.0


def nishimori_ohzeki_equation_nats(p: float) -> float:
    """Return left minus right in equation (28) of Nishimori--Ohzeki (2006)."""

    _check_probability(p, strict=True)
    q = 1.0 - p
    a = 4.0 * p * p - 6.0 * p + 3.0
    b = 4.0 * p * p - 2.0 * p + 1.0
    left = (
        2.0 * p * p * (3.0 - 2.0 * p) * log(p)
        + 2.0 * q * q * (1.0 + 2.0 * p) * log(q)
        + log(2.0)
    )
    right = p * a * log(a) + q * b * log(b)
    return left - right


def entropy_balance_derivative_nats(p: float) -> float:
    """Derivative of 3 H(p)-H((1+(2p-1)^3)/2), using natural logs."""

    if not 0.5 < p < 1.0:
        raise ValueError("p must satisfy 0.5 < p < 1")
    x = 2.0 * p - 1.0
    return -6.0 * (atanh(x) - x * x * atanh(x**3))


def conjectured_nishimori_root(tolerance: float = 5e-16) -> float:
    """Return the unique zero of the face-entropy balance in (1/2,1)."""

    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    lower, upper = 0.5, 1.0
    while upper - lower > tolerance:
        midpoint = (lower + upper) / 2.0
        if nishimori_entropy_balance_bits(midpoint) > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def race_winner_probabilities(rates: Iterable[float]) -> tuple[float, ...]:
    """Winner law for independent exponential clocks with the given rates."""

    positive_rates = tuple(float(rate) for rate in rates)
    if not positive_rates or any(rate < 0.0 for rate in positive_rates):
        raise ValueError("rates must be a non-empty family of non-negative values")
    total = sum(positive_rates)
    if total <= 0.0:
        raise ValueError("at least one rate must be positive")
    return tuple(rate / total for rate in positive_rates)


def exponential_race_entropy_bits(rates: Iterable[float]) -> float:
    """Return E[-log_2 pi_W] for the winner W of an exponential race."""

    probabilities = race_winner_probabilities(rates)
    return -sum(
        probability * log2(probability)
        for probability in probabilities
        if probability > 0.0
    )


def main() -> None:
    root = conjectured_nishimori_root()
    syndrome_plus = face_syndrome_probability(root)
    laws = conditional_face_noise_laws(root)
    conditional_entropies = {
        syndrome: exponential_race_entropy_bits(law.values())
        for syndrome, law in laws.items()
    }
    weighted_entropy = (
        syndrome_plus * conditional_entropies[1]
        + (1.0 - syndrome_plus) * conditional_entropies[-1]
    )

    print(f"p_N^(0)                  = {root:.12f}")
    print(f"P(face syndrome = +1)   = {syndrome_plus:.12f}")
    print(f"H(Z | syndrome=+1)      = {conditional_entropies[1]:.12f} bits")
    print(f"H(Z | syndrome=-1)      = {conditional_entropies[-1]:.12f} bits")
    print(f"H(Z | syndrome)         = {weighted_entropy:.12f} bits")
    print(
        "Nishimori--Ohzeki eq. 28 = "
        f"{nishimori_ohzeki_equation_nats(root):.3e} nats"
    )


if __name__ == "__main__":
    main()
