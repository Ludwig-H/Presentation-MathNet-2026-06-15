"""Exact full-dendrogram global-port convolution for the ferromagnetic SBM.

For the ordinary assortative SBM, every edge used by the full dendrogram is
satisfied exactly when its endpoint spins are equal.  Consequently every
final dendrogram root ``R`` is monochromatic.  If its size is ``s_R``, its
only remaining state is an orientation ``z_R in {-1,+1}``, contributing
``z_R s_R`` to the total magnetization

    M = sum_R z_R s_R.

The convolution in this module computes the exact multiplicity of every
possible ``M``.  Common internal root factors, which are invariant under a
global flip of one root, are divided out.  They would multiply every
partition function below by the same product and do not change the port law.

For the i.i.d.-label SBM, the non-edge port has field

    h0 = 0.5 * log((1-a_n/n)/(1-b_n/n)) < 0

and the orientation partition function is

    Z_iid = sum_M multiplicity(M) exp(h0 M^2 / 2).

For the exactly balanced planted bisection,

    Z_bal = multiplicity(0).

An exhaustive ``2^r`` enumeration is included only as a small-instance
audit.  The module neither samples the law of the dendrogram nor proves a
weak-, almost-exact-, or exact-recovery threshold.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from dataclasses import asdict, dataclass
from typing import Mapping, Sequence


EXACT_STATUS = "EXACT_FULL_D_GLOBAL_PORT_CONVOLUTION"
AUDIT_PASS_STATUS = "PASS_EXACT_ENUMERATION_IDENTITY"
AUDIT_SKIPPED_STATUS = "SKIPPED_ENUMERATION_ROOT_LIMIT"
_LOG_MAX_FLOAT = math.log(sys.float_info.max)


@dataclass(frozen=True)
class MagnetizationPortTerm:
    """One coefficient of the convolved root-orientation measure."""

    magnetization: int
    orientation_multiplicity: int
    iid_port_factor: float
    log_iid_weighted_multiplicity: float
    iid_weighted_multiplicity: float | None
    balanced_weighted_multiplicity: int


@dataclass(frozen=True)
class PortPartitionFunctions:
    """Partition functions after common internal root factors are removed."""

    iid_partition_function: float | None
    iid_log_partition_function: float
    iid_partition_normalized_by_uniform_orientations: float
    balanced_partition_function: int
    balanced_log_partition_function: float | None
    balanced_partition_normalized_by_uniform_orientations: float
    uniform_orientation_normalizer: int
    balanced_model_feasible: bool


@dataclass(frozen=True)
class DirectEnumerationAudit:
    """Comparison of convolution with exhaustive root-sign enumeration."""

    enumeration_performed: bool
    enumerated_orientation_count: int | None
    multiplicities_match_exactly: bool | None
    iid_partition_absolute_error: float | None
    balanced_partition_matches_exactly: bool | None
    audit_status: str


@dataclass(frozen=True)
class GlobalPortConvolutionDiagnostic:
    """Complete JSON-serializable audit of one collection of roots."""

    n: int
    within_parameter_a_n: float
    between_parameter_b_n: float
    nonedge_field_h0: float
    root_sizes: tuple[int, ...]
    root_count: int
    roots_are_monochromatic: bool
    common_internal_root_factors_omitted: bool
    magnetization_terms: tuple[MagnetizationPortTerm, ...]
    partition_functions: PortPartitionFunctions
    direct_enumeration_audit: DirectEnumerationAudit
    diagnostic_status: str
    dendrogram_law_sampled: bool
    recovery_threshold_claimed: bool
    interpretation: str


def _validated_root_sizes(root_sizes: Sequence[int]) -> tuple[int, ...]:
    sizes = tuple(root_sizes)
    if not sizes:
        raise ValueError("at least one root size is required")
    for size in sizes:
        if isinstance(size, bool) or not isinstance(size, int):
            raise ValueError("root sizes must be integers")
        if size <= 0:
            raise ValueError("root sizes must be positive")
    return sizes


def finite_nonedge_field(n: int, a_n: float, b_n: float) -> float:
    """Return the exact finite-``n`` field carried by all non-edges."""

    if n < 2:
        raise ValueError("n must be at least 2")
    if not math.isfinite(a_n) or not math.isfinite(b_n):
        raise ValueError("a_n and b_n must be finite")
    if not 0.0 < b_n < a_n < n:
        raise ValueError("ferromagnetic parameters must satisfy 0<b_n<a_n<n")
    return 0.5 * math.log((1.0 - a_n / n) / (1.0 - b_n / n))


def magnetization_multiplicities(
    root_sizes: Sequence[int],
) -> dict[int, int]:
    """Convolve ``delta_{-s_R}+delta_{+s_R}`` over all roots."""

    sizes = _validated_root_sizes(root_sizes)
    multiplicities = {0: 1}
    for size in sizes:
        following: dict[int, int] = {}
        for magnetization, multiplicity in multiplicities.items():
            following[magnetization - size] = (
                following.get(magnetization - size, 0) + multiplicity
            )
            following[magnetization + size] = (
                following.get(magnetization + size, 0) + multiplicity
            )
        multiplicities = following
    return dict(sorted(multiplicities.items()))


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("at least one logarithmic term is required")
    maximum = max(values)
    return maximum + math.log(
        math.fsum(math.exp(value - maximum) for value in values)
    )


def _safe_exp(value: float) -> float | None:
    if value > _LOG_MAX_FLOAT:
        return None
    return math.exp(value)


def iid_log_partition_function(
    multiplicities: Mapping[int, int],
    h0: float,
) -> float:
    """Return ``log sum_M N(M) exp(h0 M^2/2)`` stably."""

    if not multiplicities:
        raise ValueError("multiplicities must not be empty")
    if not math.isfinite(h0) or h0 > 0.0:
        raise ValueError("the non-edge field h0 must be finite and non-positive")
    logarithmic_terms = []
    for magnetization, multiplicity in multiplicities.items():
        if multiplicity <= 0:
            raise ValueError("multiplicities must be positive")
        logarithmic_terms.append(
            math.log(multiplicity)
            + 0.5 * h0 * magnetization * magnetization
        )
    return _logsumexp(logarithmic_terms)


def iid_partition_function(
    multiplicities: Mapping[int, int],
    h0: float,
) -> float:
    """Return the unnormalized i.i.d. orientation partition when finite."""

    log_partition = iid_log_partition_function(multiplicities, h0)
    partition = _safe_exp(log_partition)
    if partition is None:
        raise OverflowError(
            "the partition exceeds the floating-point range; use its log"
        )
    return partition


def balanced_partition_function(
    multiplicities: Mapping[int, int],
) -> int:
    """Return the exact number of sign assignments with total spin zero."""

    return multiplicities.get(0, 0)


def port_partition_functions(
    root_sizes: Sequence[int],
    h0: float,
) -> PortPartitionFunctions:
    """Compute both port partition functions from one convolution."""

    sizes = _validated_root_sizes(root_sizes)
    multiplicities = magnetization_multiplicities(sizes)
    log_iid_partition = iid_log_partition_function(multiplicities, h0)
    iid_partition = _safe_exp(log_iid_partition)
    orientation_normalizer = 1 << len(sizes)
    log_orientation_normalizer = len(sizes) * math.log(2.0)
    normalized_iid = math.exp(
        log_iid_partition - log_orientation_normalizer
    )
    balanced_partition = balanced_partition_function(multiplicities)
    if balanced_partition:
        log_balanced_partition = math.log(balanced_partition)
        normalized_balanced = math.exp(
            log_balanced_partition - log_orientation_normalizer
        )
    else:
        log_balanced_partition = None
        normalized_balanced = 0.0
    return PortPartitionFunctions(
        iid_partition_function=iid_partition,
        iid_log_partition_function=log_iid_partition,
        iid_partition_normalized_by_uniform_orientations=normalized_iid,
        balanced_partition_function=balanced_partition,
        balanced_log_partition_function=log_balanced_partition,
        balanced_partition_normalized_by_uniform_orientations=(
            normalized_balanced
        ),
        uniform_orientation_normalizer=orientation_normalizer,
        balanced_model_feasible=(balanced_partition > 0),
    )


def direct_orientation_enumeration(
    root_sizes: Sequence[int],
    h0: float,
) -> tuple[dict[int, int], float, int]:
    """Enumerate all signs; intended only for small audit instances."""

    sizes = _validated_root_sizes(root_sizes)
    if not math.isfinite(h0) or h0 > 0.0:
        raise ValueError("the non-edge field h0 must be finite and non-positive")
    multiplicities: dict[int, int] = {}
    iid_terms = []
    balanced_partition = 0
    for signs in itertools.product((-1, 1), repeat=len(sizes)):
        magnetization = sum(
            sign * size for sign, size in zip(signs, sizes)
        )
        multiplicities[magnetization] = (
            multiplicities.get(magnetization, 0) + 1
        )
        iid_terms.append(
            math.exp(0.5 * h0 * magnetization * magnetization)
        )
        if magnetization == 0:
            balanced_partition += 1
    return (
        dict(sorted(multiplicities.items())),
        math.fsum(iid_terms),
        balanced_partition,
    )


def audit_global_port_convolution(
    root_sizes: Sequence[int],
    h0: float,
    *,
    max_enumerated_roots: int = 22,
) -> DirectEnumerationAudit:
    """Compare convolution and direct enumeration whenever affordable."""

    sizes = _validated_root_sizes(root_sizes)
    if max_enumerated_roots < 0:
        raise ValueError("max_enumerated_roots must be non-negative")
    if len(sizes) > max_enumerated_roots:
        return DirectEnumerationAudit(
            enumeration_performed=False,
            enumerated_orientation_count=None,
            multiplicities_match_exactly=None,
            iid_partition_absolute_error=None,
            balanced_partition_matches_exactly=None,
            audit_status=AUDIT_SKIPPED_STATUS,
        )

    convolved = magnetization_multiplicities(sizes)
    direct_counts, direct_iid, direct_balanced = (
        direct_orientation_enumeration(sizes, h0)
    )
    convolved_iid = iid_partition_function(convolved, h0)
    convolved_balanced = balanced_partition_function(convolved)
    counts_match = convolved == direct_counts
    iid_error = abs(convolved_iid - direct_iid)
    balanced_match = convolved_balanced == direct_balanced
    tolerance = 1e-12 * max(1.0, convolved_iid, direct_iid)
    passed = counts_match and iid_error <= tolerance and balanced_match
    return DirectEnumerationAudit(
        enumeration_performed=True,
        enumerated_orientation_count=1 << len(sizes),
        multiplicities_match_exactly=counts_match,
        iid_partition_absolute_error=iid_error,
        balanced_partition_matches_exactly=balanced_match,
        audit_status=(
            AUDIT_PASS_STATUS
            if passed
            else "FAILED_CONVOLUTION_ENUMERATION_IDENTITY"
        ),
    )


def run_diagnostic(
    *,
    root_sizes: Sequence[int],
    a_n: float,
    b_n: float,
    max_enumerated_roots: int = 22,
) -> GlobalPortConvolutionDiagnostic:
    """Run the exact convolution and its optional exhaustive audit."""

    sizes = _validated_root_sizes(root_sizes)
    n = sum(sizes)
    h0 = finite_nonedge_field(n, a_n, b_n)
    multiplicities = magnetization_multiplicities(sizes)
    terms = []
    for magnetization, multiplicity in multiplicities.items():
        log_weighted_multiplicity = (
            math.log(multiplicity)
            + 0.5 * h0 * magnetization * magnetization
        )
        terms.append(
            MagnetizationPortTerm(
                magnetization=magnetization,
                orientation_multiplicity=multiplicity,
                iid_port_factor=math.exp(
                    0.5 * h0 * magnetization * magnetization
                ),
                log_iid_weighted_multiplicity=(
                    log_weighted_multiplicity
                ),
                iid_weighted_multiplicity=_safe_exp(
                    log_weighted_multiplicity
                ),
                balanced_weighted_multiplicity=(
                    multiplicity if magnetization == 0 else 0
                ),
            )
        )
    audit = audit_global_port_convolution(
        sizes,
        h0,
        max_enumerated_roots=max_enumerated_roots,
    )
    return GlobalPortConvolutionDiagnostic(
        n=n,
        within_parameter_a_n=a_n,
        between_parameter_b_n=b_n,
        nonedge_field_h0=h0,
        root_sizes=sizes,
        root_count=len(sizes),
        roots_are_monochromatic=True,
        common_internal_root_factors_omitted=True,
        magnetization_terms=tuple(terms),
        partition_functions=port_partition_functions(sizes, h0),
        direct_enumeration_audit=audit,
        diagnostic_status=EXACT_STATUS,
        dendrogram_law_sampled=False,
        recovery_threshold_claimed=False,
        interpretation=(
            "The convolution is exact conditional on the supplied full-D "
            "ferromagnetic root sizes. It does not generate those sizes, "
            "sample the dendrogram law, analyze mixing, or prove recovery."
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root-sizes",
        type=int,
        nargs="+",
        default=(3, 2, 1),
    )
    parser.add_argument("--a", dest="a_n", type=float, default=4.0)
    parser.add_argument("--b", dest="b_n", type=float, default=1.0)
    parser.add_argument("--max-enumerated-roots", type=int, default=22)
    arguments = parser.parse_args()
    diagnostic = run_diagnostic(
        root_sizes=arguments.root_sizes,
        a_n=arguments.a_n,
        b_n=arguments.b_n,
        max_enumerated_roots=arguments.max_enumerated_roots,
    )
    print(
        json.dumps(
            asdict(diagnostic),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
