"""Calibrate hierarchical recovery on the symmetric two-community SBM.

The local weak limit is a Poisson Galton--Watson broadcast tree.  Conditional
on a ``+1`` root spin, its exact posterior-magnetization recursion is

    M_{t+1} = tanh(sum_i atanh(theta S_i M_{t,i})),

where ``K ~ Poisson(d)``, ``E[S_i] = theta``, and the ``M_{t,i}`` are
independent copies of ``M_t`` conditional on a ``+1`` child spin.  Revealed
leaves give ``M_0 = 1``.

The Monte Carlo population dynamics is only a diagnostic.  The exact
Kesten--Stigum transition is certified by the deterministic sandwich

    ell_t(lambda) <= E[M_t^2] <= r_t(lambda),

where ``ell_t`` is the linear boundary-score bound and ``r_t`` is the
information-percolation survival probability:

    ell_t = 1 / sum_{s=0}^t lambda^{-s},
    r_0 = 1,   r_{t+1} = 1 - exp(-lambda r_t).

Thus the local tree is nonreconstructible for ``lambda <= 1`` and
reconstructible for ``lambda > 1``.  This module does not by itself prove the
transfer from the tree to graph weak recovery.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from statistics import fmean
from typing import Sequence


PASS_STATUS = "PASS_EMPIRICAL_CONSISTENCY"
WARNING_STATUS = "WARNING_EMPIRICAL_INCONSISTENCY"
DETERMINISTIC_CERTIFICATE = (
    "exact PGW broadcast threshold certified by the linear-score lower bound "
    "and chi-square information-percolation upper bound"
)


@dataclass(frozen=True)
class DepthDiagnostic:
    """Population estimate and exact deterministic bounds at one depth."""

    depth: int
    conditional_mean_hat: float
    conditional_mean_standard_error: float | None
    q_hat: float
    q_standard_error: float | None
    nishimori_gap: float
    ell_t: float
    r_t: float
    nishimori_consistent_up_to_tolerance: bool
    sandwich_consistent_up_to_tolerance: bool


@dataclass(frozen=True)
class LambdaDiagnostic:
    """One signal-to-noise trajectory."""

    signal_to_noise: float
    theta: float
    deterministic_regime: str
    asymptotic_linear_score_lower_bound: float
    empirical_consistency_status: str
    estimates: tuple[DepthDiagnostic, ...]


@dataclass(frozen=True)
class BroadcastDiagnosticSummary:
    """Reproducible collection of PGW broadcast trajectories."""

    degree: float
    max_depth: int
    particles: int
    batches: int
    seed: int
    consistency_tolerance: float
    deterministic_certificate: str
    diagnostic_status: str
    lambda_diagnostics: tuple[LambdaDiagnostic, ...]
    graph_weak_recovery_claimed: bool
    interpretation: str


def linear_score_lower_bound(signal_to_noise: float, depth: int) -> float:
    """Return ``ell_t = 1 / sum_{s=0}^t lambda^{-s}`` stably."""

    if signal_to_noise < 0.0:
        raise ValueError("signal_to_noise must be non-negative")
    if depth < 0:
        raise ValueError("depth must be non-negative")
    if depth == 0:
        return 1.0
    if signal_to_noise == 0.0:
        return 0.0
    if signal_to_noise == 1.0:
        return 1.0 / (depth + 1)
    if signal_to_noise < 1.0:
        power = signal_to_noise**depth
        return (
            (1.0 - signal_to_noise)
            * power
            / (1.0 - signal_to_noise * power)
        )
    inverse = 1.0 / signal_to_noise
    return (1.0 - inverse) / (1.0 - inverse ** (depth + 1))


def information_percolation_trajectory(
    signal_to_noise: float,
    max_depth: int,
) -> tuple[float, ...]:
    """Return ``r_0,...,r_T`` for PGW information percolation."""

    if signal_to_noise < 0.0:
        raise ValueError("signal_to_noise must be non-negative")
    if max_depth < 0:
        raise ValueError("max_depth must be non-negative")
    values = [1.0]
    for _ in range(max_depth):
        values.append(-math.expm1(-signal_to_noise * values[-1]))
    return tuple(values)


def deterministic_regime(signal_to_noise: float) -> str:
    """State the rigorous PGW reconstruction conclusion at one lambda."""

    if signal_to_noise < 0.0:
        raise ValueError("signal_to_noise must be non-negative")
    if signal_to_noise <= 1.0:
        return (
            "NONRECONSTRUCTION_CERTIFIED: r_t tends to zero for lambda <= 1"
        )
    return (
        "RECONSTRUCTION_CERTIFIED: liminf q_t >= "
        f"{(signal_to_noise - 1.0) / signal_to_noise:.12g}"
    )


def posterior_update(
    signed_child_magnetizations: Sequence[float],
    theta: float,
) -> float:
    """Combine already signed child messages by the exact tree recursion."""

    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must belong to [0, 1]")
    field = 0.0
    positive_infinity = False
    negative_infinity = False
    for magnetization in signed_child_magnetizations:
        if not -1.0 <= magnetization <= 1.0:
            raise ValueError("child magnetizations must belong to [-1, 1]")
        weighted = theta * magnetization
        if weighted >= 1.0:
            positive_infinity = True
        elif weighted <= -1.0:
            negative_infinity = True
        else:
            field += math.atanh(weighted)
    if positive_infinity and negative_infinity:
        raise ValueError("incompatible perfectly informative child messages")
    if positive_infinity:
        return 1.0
    if negative_infinity:
        return -1.0
    return math.tanh(field)


def _sample_poisson(mean: float, rng: random.Random) -> int:
    """Sample Poisson(mean) exactly by summing stable Knuth chunks."""

    if mean < 0.0:
        raise ValueError("a Poisson mean must be non-negative")
    if mean == 0.0:
        return 0
    chunk_count = max(1, math.ceil(mean / 20.0))
    chunk_mean = mean / chunk_count
    cutoff = math.exp(-chunk_mean)
    result = 0
    for _ in range(chunk_count):
        product = 1.0
        count = -1
        while product > cutoff:
            product *= rng.random()
            count += 1
        result += count
    return result


def _population_step(
    population: Sequence[float],
    *,
    degree: float,
    theta: float,
    rng: random.Random,
) -> list[float]:
    """Advance one conditional-root population-dynamics step."""

    size = len(population)
    sign_plus_probability = 0.5 * (1.0 + theta)
    result: list[float] = []
    append = result.append
    for _ in range(size):
        field = 0.0
        positive_infinity = False
        negative_infinity = False
        for _ in range(_sample_poisson(degree, rng)):
            child = population[rng.randrange(size)]
            if rng.random() >= sign_plus_probability:
                child = -child
            weighted = theta * child
            if weighted >= 1.0:
                positive_infinity = True
            elif weighted <= -1.0:
                negative_infinity = True
            else:
                field += math.atanh(weighted)
        if positive_infinity and negative_infinity:
            raise AssertionError(
                "the broadcast recursion produced incompatible exact messages"
            )
        if positive_infinity:
            append(1.0)
        elif negative_infinity:
            append(-1.0)
        else:
            append(math.tanh(field))
    return result


def _simulate_batch(
    *,
    degree: float,
    theta: float,
    particle_count: int,
    max_depth: int,
    seed: int,
) -> tuple[tuple[float, float], ...]:
    """Return ``(conditional mean, second moment)`` at every depth."""

    rng = random.Random(seed)
    population = [1.0] * particle_count
    trajectory = []
    for depth in range(max_depth + 1):
        conditional_mean = math.fsum(population) / particle_count
        q_hat = math.fsum(value * value for value in population) / particle_count
        trajectory.append((conditional_mean, q_hat))
        if depth < max_depth:
            population = _population_step(
                population,
                degree=degree,
                theta=theta,
                rng=rng,
            )
    return tuple(trajectory)


def _mean_and_cluster_standard_error(
    values: Sequence[float],
) -> tuple[float, float | None]:
    """Aggregate independent batch estimates."""

    if not values:
        raise ValueError("at least one batch estimate is required")
    mean = fmean(values)
    if len(values) == 1:
        return mean, None
    standard_error = math.sqrt(
        math.fsum((value - mean) ** 2 for value in values)
        / (len(values) * (len(values) - 1))
    )
    return mean, standard_error


def _particle_allocation(particles: int, batches: int) -> tuple[int, ...]:
    if particles <= 0:
        raise ValueError("particles must be positive")
    if batches <= 0:
        raise ValueError("batches must be positive")
    if particles < batches:
        raise ValueError("particles must be at least batches")
    quotient, remainder = divmod(particles, batches)
    return tuple(
        quotient + (batch < remainder)
        for batch in range(batches)
    )


def _validate_parameters(
    *,
    degree: float,
    signal_to_noise_values: Sequence[float],
    max_depth: int,
) -> tuple[float, ...]:
    if degree < 0.0:
        raise ValueError("degree must be non-negative")
    if max_depth < 0:
        raise ValueError("max_depth must be non-negative")
    values = tuple(float(value) for value in signal_to_noise_values)
    if not values:
        raise ValueError("at least one signal-to-noise value is required")
    if any(value < 0.0 for value in values):
        raise ValueError("signal-to-noise values must be non-negative")
    if degree == 0.0:
        if any(value != 0.0 for value in values):
            raise ValueError("degree zero only permits lambda zero")
    elif any(value > degree for value in values):
        raise ValueError("lambda must not exceed d because theta <= 1")
    return values


def run_diagnostic(
    *,
    degree: float,
    signal_to_noise_values: Sequence[float],
    max_depth: int,
    particles: int,
    batches: int,
    seed: int,
) -> BroadcastDiagnosticSummary:
    """Run reproducible conditional-root population dynamics."""

    values = _validate_parameters(
        degree=degree,
        signal_to_noise_values=signal_to_noise_values,
        max_depth=max_depth,
    )
    allocations = _particle_allocation(particles, batches)
    tolerance = min(1.0, 6.0 / math.sqrt(particles))
    master_rng = random.Random(seed)
    lambda_diagnostics = []

    for signal_to_noise in values:
        theta = (
            0.0
            if degree == 0.0
            else math.sqrt(signal_to_noise / degree)
        )
        trajectories = tuple(
            _simulate_batch(
                degree=degree,
                theta=theta,
                particle_count=particle_count,
                max_depth=max_depth,
                seed=master_rng.getrandbits(64),
            )
            for particle_count in allocations
        )
        upper_bounds = information_percolation_trajectory(
            signal_to_noise, max_depth
        )
        estimates = []
        for depth in range(max_depth + 1):
            conditional_mean, conditional_mean_se = (
                _mean_and_cluster_standard_error(
                    tuple(trajectory[depth][0] for trajectory in trajectories)
                )
            )
            q_hat, q_standard_error = _mean_and_cluster_standard_error(
                tuple(trajectory[depth][1] for trajectory in trajectories)
            )
            lower_bound = linear_score_lower_bound(signal_to_noise, depth)
            upper_bound = upper_bounds[depth]
            gap = conditional_mean - q_hat
            estimates.append(
                DepthDiagnostic(
                    depth=depth,
                    conditional_mean_hat=conditional_mean,
                    conditional_mean_standard_error=conditional_mean_se,
                    q_hat=q_hat,
                    q_standard_error=q_standard_error,
                    nishimori_gap=gap,
                    ell_t=lower_bound,
                    r_t=upper_bound,
                    nishimori_consistent_up_to_tolerance=(
                        abs(gap) <= tolerance
                    ),
                    sandwich_consistent_up_to_tolerance=(
                        lower_bound - tolerance
                        <= q_hat
                        <= upper_bound + tolerance
                    ),
                )
            )
        empirical_pass = all(
            estimate.nishimori_consistent_up_to_tolerance
            and estimate.sandwich_consistent_up_to_tolerance
            for estimate in estimates
        )
        lambda_diagnostics.append(
            LambdaDiagnostic(
                signal_to_noise=signal_to_noise,
                theta=theta,
                deterministic_regime=deterministic_regime(signal_to_noise),
                asymptotic_linear_score_lower_bound=(
                    0.0
                    if signal_to_noise <= 1.0
                    else (signal_to_noise - 1.0) / signal_to_noise
                ),
                empirical_consistency_status=(
                    PASS_STATUS if empirical_pass else WARNING_STATUS
                ),
                estimates=tuple(estimates),
            )
        )

    global_pass = all(
        diagnostic.empirical_consistency_status == PASS_STATUS
        for diagnostic in lambda_diagnostics
    )
    return BroadcastDiagnosticSummary(
        degree=degree,
        max_depth=max_depth,
        particles=particles,
        batches=batches,
        seed=seed,
        consistency_tolerance=tolerance,
        deterministic_certificate=DETERMINISTIC_CERTIFICATE,
        diagnostic_status=PASS_STATUS if global_pass else WARNING_STATUS,
        lambda_diagnostics=tuple(lambda_diagnostics),
        graph_weak_recovery_claimed=False,
        interpretation=(
            "Monte Carlo checks the exact conditional-root recursion; the "
            "ell_t/r_t sandwich, not the samples, certifies the PGW threshold "
            "lambda=1. A separate tree-to-graph argument is required."
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=float, default=3.0)
    parser.add_argument(
        "--lambdas",
        type=float,
        nargs="+",
        default=(0.8, 0.95, 1.0, 1.05, 1.2),
    )
    parser.add_argument("--depth", type=int, default=30)
    parser.add_argument(
        "--particles",
        type=int,
        default=20_000,
        help="total population, split across independent batches",
    )
    parser.add_argument("--batches", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260726)
    arguments = parser.parse_args()

    summary = run_diagnostic(
        degree=arguments.degree,
        signal_to_noise_values=arguments.lambdas,
        max_depth=arguments.depth,
        particles=arguments.particles,
        batches=arguments.batches,
        seed=arguments.seed,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
