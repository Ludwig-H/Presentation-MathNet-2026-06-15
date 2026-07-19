"""Exact critical-boundary calculus for hierarchical Swendsen--Wang.

Only edges crossing the current merge cut enter a Lambda rate.  The module
therefore separates boundary geometry, residual boundary marks, pair-Palm
selection, and the four grouped ancestor rates.  Internal edges are never
silently assigned the residual boundary law.

The formulas do not assert a global weak-recovery threshold.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from math import comb, exp, fsum, log, log1p, prod, sqrt, tanh

from critical_band_thresholds import Q_CRITICAL, coupling


P_SW = (1.0 + Q_CRITICAL) / 2.0
P_BOUNDARY_LATE = (2.0 + Q_CRITICAL) / 3.0
P_INFO = (1.0 + sqrt(Q_CRITICAL)) / 2.0

Edge = tuple[int, int]
Partition = tuple[frozenset[int], ...]


@dataclass(frozen=True)
class ClosedCategories:
    """Categories for a boundary edge conditional on no opening by ``time``."""

    p: float
    time: float
    late_true: float
    censored_true: float
    false: float
    true_probability: float
    signed_margin: float


@dataclass(frozen=True)
class CriticalMasses:
    """Unconditional four-category masses at the percolation time."""

    early_true: float
    late_true: float
    censored_true: float
    false: float


@dataclass(frozen=True)
class InternalProportions:
    """Targets among all potential internal edges at fixed open density."""

    open_true: float
    late_true: float
    censored_true: float
    false: float
    unactivated_true: float


def _validate_p(p: float) -> None:
    if not 0.5 < p < 1.0:
        raise ValueError("p must satisfy 0.5 < p < 1")


def critical_time_untruncated(p: float) -> float:
    """Solve q_p(t)=q_c without restricting t to the censoring interval."""

    _validate_p(p)
    return -log1p(-Q_CRITICAL / p) / coupling(p)


def future_activation_crossover_time(p: float) -> float:
    """Solve late-true mass = false mass without truncating to [0, 1].

    Here ``late true`` means a true edge whose clock lies in ``(time, 1]``.
    This diagnostic deliberately excludes true edges censored after time 1.
    The returned value belongs to ``[0, 1]`` exactly when ``p >= 2/3``.
    """

    _validate_p(p)
    return 1.0 - log(2.0) / coupling(p)


def unactivated_true_minus_false_mass(p: float, time: float) -> float:
    """Unconditional mass gap: all unactivated true edges minus false edges."""

    _validate_p(p)
    if not 0.0 <= time <= 1.0:
        raise ValueError("time must belong to [0, 1]")
    return p * exp(-coupling(p) * time) - (1.0 - p)


def future_true_minus_false_mass(p: float, time: float) -> float:
    """Mass gap: true clocks in ``(time, 1]`` minus false edges."""

    _validate_p(p)
    if not 0.0 <= time <= 1.0:
        raise ValueError("time must belong to [0, 1]")
    return p * exp(-coupling(p) * time) - 2.0 * (1.0 - p)


def closed_categories(p: float, time: float) -> ClosedCategories:
    """Return the exact categories conditional on no opening before ``time``."""

    _validate_p(p)
    if not -1e-14 <= time <= 1.0 + 1e-14:
        raise ValueError("time must belong to [0, 1]")
    time = min(1.0, max(0.0, time))

    q = 1.0 - p
    surviving_true = p * exp(-coupling(p) * time)
    normalizer = q + surviving_true
    late_true = (surviving_true - q) / normalizer
    censored_true = q / normalizer
    false = q / normalizer
    true_probability = late_true + censored_true
    signed_margin = 2.0 * true_probability - 1.0
    return ClosedCategories(
        p=p,
        time=time,
        late_true=late_true,
        censored_true=censored_true,
        false=false,
        true_probability=true_probability,
        signed_margin=signed_margin,
    )


def internal_proportions(
    p: float, time: float, open_density: float
) -> InternalProportions:
    """Return conditional targets among all potential edges inside a cluster.

    The cluster selection may be any event measurable by the open graph at
    ``time``.  Its only input here is the realized or limiting fraction of
    already-open internal edges.  The residual split is the universal closed
    edge law returned by :func:`closed_categories`.
    """

    if not 0.0 <= open_density <= 1.0:
        raise ValueError("open_density must belong to [0, 1]")
    residual = closed_categories(p, time)
    closed_density = 1.0 - open_density
    late_true = closed_density * residual.late_true
    censored_true = closed_density * residual.censored_true
    false = closed_density * residual.false
    return InternalProportions(
        open_true=open_density,
        late_true=late_true,
        censored_true=censored_true,
        false=false,
        unactivated_true=late_true + censored_true,
    )


def critical_masses(p: float) -> CriticalMasses:
    """Return the four unconditional masses at beta_c, for beta_c <= 1."""

    _validate_p(p)
    if p < P_SW:
        raise ValueError("the critical time lies after the censoring time 1")
    q = 1.0 - p
    return CriticalMasses(
        early_true=Q_CRITICAL,
        late_true=2.0 * p - 1.0 - Q_CRITICAL,
        censored_true=q,
        false=q,
    )


def critical_closed_categories(p: float) -> ClosedCategories:
    """Return conditional quality of a critical boundary edge."""

    return closed_categories(p, critical_time_untruncated(p))


def boundary_late_true_fraction(p: float) -> float:
    """Boundary diagnostic: late-true mass among late-true and false marks."""

    masses = critical_masses(p)
    return masses.late_true / (masses.late_true + masses.false)


def connected_components(
    vertices: Iterable[int], open_edges: Iterable[Edge]
) -> Partition:
    """Return the canonical component partition of a finite open graph."""

    vertex_set = set(vertices)
    parent = {vertex: vertex for vertex in vertex_set}

    def find(vertex: int) -> int:
        root = vertex
        while parent[root] != root:
            root = parent[root]
        while parent[vertex] != vertex:
            next_vertex = parent[vertex]
            parent[vertex] = root
            vertex = next_vertex
        return root

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left, right in open_edges:
        if left not in vertex_set or right not in vertex_set:
            raise ValueError("every edge endpoint must belong to vertices")
        union(left, right)

    blocks: dict[int, set[int]] = {}
    for vertex in vertex_set:
        blocks.setdefault(find(vertex), set()).add(vertex)
    return tuple(
        sorted(
            (frozenset(block) for block in blocks.values()),
            key=lambda block: (min(block), len(block)),
        )
    )


def partition_boundary(
    edges: Iterable[Edge], partition: Partition
) -> tuple[Edge, ...]:
    """Return exactly the edges whose endpoints lie in distinct blocks."""

    membership: dict[int, int] = {}
    for block_index, block in enumerate(partition):
        if not block:
            raise ValueError("partition blocks must be nonempty")
        for vertex in block:
            if vertex in membership:
                raise ValueError("partition blocks must be disjoint")
            membership[vertex] = block_index

    boundary: list[Edge] = []
    for left, right in edges:
        if left not in membership or right not in membership:
            raise ValueError("every edge endpoint must belong to the partition")
        if membership[left] != membership[right]:
            boundary.append((left, right))
    return tuple(boundary)


def pair_palm_weights(
    partition: Partition,
    distance: Callable[[int, int], float],
    minimum_distance: float,
) -> dict[frozenset[int], float]:
    """Return cluster weights induced by an ordered distant vertex pair."""

    if minimum_distance < 0.0:
        raise ValueError("minimum_distance must be nonnegative")
    counts = {
        block: sum(
            left != right and distance(left, right) >= minimum_distance
            for left in block
            for right in block
        )
        for block in partition
    }
    total = sum(counts.values())
    if total == 0:
        raise ValueError("the partition contains no eligible ordered pair")
    return {block: count / total for block, count in counts.items()}


def strict_majority_probability(size: int, p: float, time: float) -> float:
    """P(1 + Bin(size-1, s_p(time)) > size/2)."""

    if size <= 0:
        raise ValueError("size must be positive")
    s = closed_categories(p, time).true_probability
    start = size // 2
    return fsum(
        comb(size - 1, successes)
        * s**successes
        * (1.0 - s) ** (size - 1 - successes)
        for successes in range(start, size)
    )


def bucket_vote_moments(size: int, p: float, time: float) -> tuple[float, float]:
    """Return mean and variance of V=2K-size for K=1+Bin(size-1,s)."""

    if size <= 0:
        raise ValueError("size must be positive")
    h = closed_categories(p, time).signed_margin
    mean = 1.0 + (size - 1) * h
    variance = (size - 1) * (1.0 - h * h)
    return mean, variance


def margin_failure_bound(
    residual_count: int,
    winner_indicator: int,
    signed_margin: float,
    require_strict: bool = False,
) -> float:
    """Hoeffding bound for failure of a nonnegative or strict margin."""

    if residual_count < 0:
        raise ValueError("residual_count must be nonnegative")
    if winner_indicator not in (0, 1):
        raise ValueError("winner_indicator must be 0 or 1")
    if not 0.0 <= signed_margin <= 1.0:
        raise ValueError("signed_margin must belong to [0, 1]")
    if residual_count == 0:
        fails = winner_indicator <= 0 if require_strict else winner_indicator < 0
        return float(fails)
    deviation = residual_count * signed_margin + winner_indicator
    return exp(-deviation * deviation / (2.0 * residual_count))


def margin_success_probability(
    residual_count: int,
    winner_indicator: int,
    true_probability: float,
    require_strict: bool = False,
) -> float:
    """Exact P(g + 2 Bin(n,s) - n >= 0), or >0 when requested."""

    if residual_count < 0:
        raise ValueError("residual_count must be nonnegative")
    if winner_indicator not in (0, 1):
        raise ValueError("winner_indicator must be 0 or 1")
    if not 0.0 <= true_probability <= 1.0:
        raise ValueError("true_probability must belong to [0, 1]")
    difference = residual_count - winner_indicator
    if require_strict:
        start = difference // 2 + 1
    else:
        start = (difference + 1) // 2
    start = max(0, start)
    return fsum(
        comb(residual_count, successes)
        * true_probability**successes
        * (1.0 - true_probability) ** (residual_count - successes)
        for successes in range(start, residual_count + 1)
    )


def grouped_rates(
    sizes: tuple[int, int, int],
    true_counts: tuple[int, int, int],
    weight: float,
) -> dict[tuple[int, int], float]:
    """Return Lambda_v^{ab} for the three-group ancestor decomposition."""

    if weight <= 0.0:
        raise ValueError("weight must be positive")
    if any(size < 0 for size in sizes):
        raise ValueError("group sizes must be nonnegative")
    if any(not 0 <= count <= size for count, size in zip(true_counts, sizes)):
        raise ValueError("each true count must belong to its group")

    _, m1, m2 = sizes
    k0, k1, k2 = true_counts
    return {
        (a, b): weight
        * (
            k0
            + (k1 if a == 0 else m1 - k1)
            + (k2 if b == 0 else m2 - k2)
        )
        for a in (0, 1)
        for b in (0, 1)
    }


def clock_factor_weights(
    rates: dict[tuple[int, int], float], beta: float
) -> dict[tuple[int, int], float]:
    """Apply F_beta(x)=x exp((1-beta)x) to four nonnegative rates."""

    expected = {(0, 0), (0, 1), (1, 0), (1, 1)}
    if set(rates) != expected:
        raise ValueError("rates must contain the four flip states")
    if any(value < 0.0 for value in rates.values()):
        raise ValueError("rates must be nonnegative")
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must belong to [0, 1]")
    return {
        state: value * exp((1.0 - beta) * value)
        for state, value in rates.items()
    }


def multiply_four_weights(
    *tables: dict[tuple[int, int], float],
) -> dict[tuple[int, int], float]:
    """Pointwise product of four-state factors."""

    expected = {(0, 0), (0, 1), (1, 0), (1, 1)}
    if not tables:
        return {state: 1.0 for state in expected}
    if any(set(table) != expected for table in tables):
        raise ValueError("each table must contain the four flip states")
    return {
        state: prod(table[state] for table in tables)
        for state in expected
    }


def walsh_coefficients(
    weights: dict[tuple[int, int], float],
) -> dict[tuple[int, ...], float]:
    """Return the four normalized Walsh coefficients on {0,1}^2."""

    expected = {(0, 0), (0, 1), (1, 0), (1, 1)}
    if set(weights) != expected:
        raise ValueError("weights must contain the four flip states")
    f00 = weights[(0, 0)]
    f10 = weights[(1, 0)]
    f01 = weights[(0, 1)]
    f11 = weights[(1, 1)]
    return {
        (): (f00 + f10 + f01 + f11) / 4.0,
        (1,): (f00 + f01 - f10 - f11) / 4.0,
        (2,): (f00 + f10 - f01 - f11) / 4.0,
        (1, 2): (f00 + f11 - f10 - f01) / 4.0,
    }


def parity_keep_probability(weights: dict[tuple[int, int], float]) -> float:
    """Return (q00+q11)/sum(qab) for four nonnegative heat-bath weights."""

    expected = {(0, 0), (0, 1), (1, 0), (1, 1)}
    if set(weights) != expected:
        raise ValueError("weights must contain the four flip states")
    if any(value < 0.0 for value in weights.values()):
        raise ValueError("weights must be nonnegative")
    total = fsum(weights.values())
    if total <= 0.0:
        raise ValueError("at least one weight must be positive")
    return (weights[(0, 0)] + weights[(1, 1)]) / total


def parity_log_odds(weights: dict[tuple[int, int], float]) -> float:
    """Return the exact even-versus-odd log odds for positive weights."""

    if any(value <= 0.0 for value in weights.values()):
        raise ValueError("log odds require four positive weights")
    return log(
        (weights[(0, 0)] + weights[(1, 1)])
        / (weights[(1, 0)] + weights[(0, 1)])
    )


def main() -> None:
    print("critical boundary thresholds")
    print(f"q_c       = {Q_CRITICAL:.12f}")
    print(f"p_SW      = {P_SW:.12f}")
    print(f"p_dlate   = {P_BOUNDARY_LATE:.12f}")
    print(f"p_info    = {P_INFO:.12f}")
    p_target = 4.0 / 5.0
    beta_critical = critical_time_untruncated(p_target)
    beta_active = future_activation_crossover_time(p_target)
    print("\nsimple residual-balance audit at p=0.8")
    print(f"beta_c    = {beta_critical:.12f}")
    print(f"beta_act  = {beta_active:.12f}")
    print(
        "all-unactivated true minus false at beta_c = "
        f"{unactivated_true_minus_false_mass(p_target, beta_critical):.12f}"
    )
    print(
        "future-ringing true minus false at beta_act = "
        f"{future_true_minus_false_mass(p_target, beta_active):.12f}"
    )
    print("\ncritical category audit")
    print("p              beta_c        late          censored      false         s_c")
    for p in (P_SW, P_BOUNDARY_LATE, P_INFO, 0.835805792367):
        beta = critical_time_untruncated(p)
        categories = critical_closed_categories(p)
        print(
            f"{p:.10f}   {beta:.10f}   {categories.late_true:.10f}   "
            f"{categories.censored_true:.10f}   {categories.false:.10f}   "
            f"{categories.true_probability:.10f}"
        )


if __name__ == "__main__":
    main()
