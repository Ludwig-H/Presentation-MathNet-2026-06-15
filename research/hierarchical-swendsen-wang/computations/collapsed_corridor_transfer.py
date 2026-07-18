"""Exact finite transfer for a collapsed hierarchical corridor.

The latent variables are binary parities carried by a fixed corridor.  Given
those parities, merger-bucket counts are independent and follow the mirrored
experiments of ``favorable_time_comparison``.  The module enumerates the
posterior second moment of any bounded target, including the product parity
between the two endpoints.

This is exact for the stated fixed-skeleton experiment.  It is not a claim
that the triangular-grid Palm corridor has independent latent parities or a
bounded boundary state.
"""

from __future__ import annotations

from itertools import product
from math import exp, fsum
from typing import Mapping, Sequence

from favorable_time_comparison import (
    bucket_binary_experiment,
    critical_beta,
    screened_two_edge_contraction,
)


SpinState = tuple[int, ...]


def binary_experiment_reliability(
    plus: Sequence[float], minus: Sequence[float]
) -> float:
    """Return E[E[X | Y]^2] for a uniform binary input X."""

    if len(plus) != len(minus) or not plus:
        raise ValueError("the two channel rows must have the same support")
    if min((*plus, *minus)) < 0.0:
        raise ValueError("channel masses must be nonnegative")
    if abs(fsum(plus) - 1.0) > 1e-12 or abs(fsum(minus) - 1.0) > 1e-12:
        raise ValueError("each channel row must sum to one")
    return 0.5 * fsum(
        (plus_mass - minus_mass) ** 2 / (plus_mass + minus_mass)
        for plus_mass, minus_mass in zip(plus, minus, strict=True)
        if plus_mass + minus_mass > 0.0
    )


def uniform_spin_prior(length: int) -> dict[SpinState, float]:
    """Return the uniform law on ``{-1,+1}^length``."""

    if length <= 0:
        raise ValueError("length must be positive")
    states = tuple(product((-1, 1), repeat=length))
    mass = 1.0 / len(states)
    return {state: mass for state in states}


def ising_chain_prior(
    length: int, interaction: float, field: float = 0.0
) -> dict[SpinState, float]:
    """Return a finite correlated prior used to counter-audit tensorization."""

    if length <= 0:
        raise ValueError("length must be positive")
    unnormalized: dict[SpinState, float] = {}
    for state in product((-1, 1), repeat=length):
        energy = interaction * sum(
            state[index] * state[index + 1]
            for index in range(length - 1)
        ) + field * sum(state)
        unnormalized[state] = exp(energy)
    normalizer = fsum(unnormalized.values())
    return {
        state: weight / normalizer for state, weight in unnormalized.items()
    }


def _normalized_prior(
    length: int, prior: Mapping[SpinState, float]
) -> tuple[tuple[SpinState, float], ...]:
    expected = set(product((-1, 1), repeat=length))
    if set(prior) != expected:
        raise ValueError("the prior must contain every binary corridor state")
    if min(prior.values()) < 0.0:
        raise ValueError("prior masses must be nonnegative")
    normalizer = fsum(prior.values())
    if normalizer <= 0.0:
        raise ValueError("the prior must have positive mass")
    return tuple((state, prior[state] / normalizer) for state in sorted(prior))


def collapsed_corridor_second_moment(
    p: float,
    sizes: Sequence[int],
    levels: Sequence[float],
    prior: Mapping[SpinState, float] | None = None,
    target: Mapping[SpinState, float] | None = None,
) -> float:
    """Enumerate E[E[F(X) | K_1,...,K_h]^2] exactly in finite volume."""

    if not sizes or len(sizes) != len(levels):
        raise ValueError("sizes and levels must have the same positive length")
    if min(sizes) <= 0:
        raise ValueError("bucket sizes must be positive")
    length = len(sizes)
    if prior is None:
        prior = uniform_spin_prior(length)
    normalized = _normalized_prior(length, prior)
    if target is None:
        target = {
            state: float(_spin_product(state)) for state, _ in normalized
        }
    if set(target) != {state for state, _ in normalized}:
        raise ValueError("the target must be defined on every corridor state")
    if max(abs(value) for value in target.values()) > 1.0 + 1e-12:
        raise ValueError("the target must take values in [-1, 1]")

    channels = tuple(
        bucket_binary_experiment(p, size, level)
        for size, level in zip(sizes, levels, strict=True)
    )
    result = 0.0
    for counts in product(*(range(size + 1) for size in sizes)):
        mixture_terms = []
        numerator_terms = []
        for state, prior_mass in normalized:
            likelihood = 1.0
            for parity, channel, count in zip(
                state, channels, counts, strict=True
            ):
                row = channel[0] if parity == 1 else channel[1]
                likelihood *= row[count]
            joint_mass = prior_mass * likelihood
            mixture_terms.append(joint_mass)
            numerator_terms.append(joint_mass * target[state])
        mixture = fsum(mixture_terms)
        if mixture > 0.0:
            numerator = fsum(numerator_terms)
            result += numerator * numerator / mixture
    return result


def _spin_product(state: SpinState) -> int:
    value = 1
    for spin in state:
        value *= spin
    return value


def factorized_corridor_reliability(
    p: float, sizes: Sequence[int], levels: Sequence[float]
) -> float:
    """Return the exact product for independent uniform corridor parities."""

    if not sizes or len(sizes) != len(levels):
        raise ValueError("sizes and levels must have the same positive length")
    reliability = 1.0
    for size, level in zip(sizes, levels, strict=True):
        reliability *= binary_experiment_reliability(
            *bucket_binary_experiment(p, size, level)
        )
    return reliability


def screened_two_edge_corridor_bound(
    p: float, level: float, message_bound: float, block_count: int
) -> float:
    """Return kappa_2(b)^N for an explicitly screened two-edge corridor."""

    if block_count < 0:
        raise ValueError("block_count must be nonnegative")
    coefficient = screened_two_edge_contraction(p, level, message_bound)
    return coefficient**block_count


def main() -> None:
    p = 4.0 / 5.0
    critical = critical_beta(p)
    sizes = (2, 3, 2, 4)
    critical_levels = (critical,) * len(sizes)
    late_levels = (0.55, 0.7, 0.85, 1.0)
    for name, prior in (
        ("uniform", uniform_spin_prior(len(sizes))),
        ("correlated", ising_chain_prior(len(sizes), interaction=0.6)),
    ):
        critical_value = collapsed_corridor_second_moment(
            p, sizes, critical_levels, prior
        )
        late_value = collapsed_corridor_second_moment(
            p, sizes, late_levels, prior
        )
        print(
            f"{name}: critical={critical_value:.12f} "
            f"late={late_value:.12f} gap={critical_value-late_value:.12f}"
        )
    for count in (5, 10, 20, 40):
        value = screened_two_edge_corridor_bound(
            p, critical, message_bound=0.0, block_count=count
        )
        print(f"neutral m=2 blocks={count:2d} bound={value:.12g}")


if __name__ == "__main__":
    main()
