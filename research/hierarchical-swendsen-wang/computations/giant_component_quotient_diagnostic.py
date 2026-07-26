"""Final-giant quotient-path diagnostic above the critical Kruskal cut.

For one signal parameter ``p``, every physical edge receives an independent
uniform percolation rank.  Kruskal is run up to the final rank

    q_1 = 2p - 1.

The largest final tree ``R_*`` is then cut at the triangular-lattice bond
threshold ``q_c``.  This gives the critical blocks ``Pi_{q_c}`` contained in
``R_*``.  Endpoint pairs are sampled uniformly from all ordered pairs in
``R_*``, conditional on being macroscopically distant and lying in distinct
critical blocks.

Only merger nodes with realized rank strictly above ``q_c`` are retained on
the two-arm path.  Thus the reported path is exactly the path in the tree
obtained by contracting every critical block.  Bucket sizes remain the
physical Kruskal cut sizes from :mod:`critical_pair_path_geometry`.

The final attenuation is deliberately labelled ``PATH-FAC local oracle
(non-proof)``.  It multiplies ancestor-neutral local reliabilities along the
quotient path.  It does not establish independence of the node updates,
control exterior messages, or prove weak-recovery impossibility.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from math import fsum, sqrt
from statistics import fmean
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import (
    CriticalKruskalForest,
    RankedEdge,
    RankedPairPathGeometry,
    clock_time_from_rank,
    evaluate_ranked_geometry,
    sample_ranked_edges,
    triangular_torus_distance,
)


BASELINE_P = 0.809439
P_SW = (1.0 + Q_CRITICAL) / 2.0
PATH_FAC_STATUS = (
    "PATH-FAC local factorized oracle (non-proof): ancestor-neutral bucket "
    "reliabilities are multiplied; dependence, exterior messages, and the "
    "joint hierarchical Gibbs sweep are not controlled"
)


class NoEligibleCrossBlockPairError(RuntimeError):
    """The sampled final tree has no pair satisfying the diagnostic filter."""


@dataclass(frozen=True)
class CriticalPartition:
    """Critical blocks inside one final Kruskal tree."""

    final_root: int
    component: tuple[int, ...]
    block_roots: tuple[int, ...]
    blocks: tuple[tuple[int, ...], ...]
    block_index_by_vertex: tuple[int, ...]

    def block_index(self, vertex: int) -> int:
        """Return the block index of a vertex in the final component."""

        if not 0 <= vertex < len(self.block_index_by_vertex):
            raise ValueError("vertex is outside the torus")
        result = self.block_index_by_vertex[vertex]
        if result < 0:
            raise ValueError("vertex is outside the selected final component")
        return result

    def block_size(self, vertex: int) -> int:
        """Return the size of the critical block containing ``vertex``."""

        return len(self.blocks[self.block_index(vertex)])


@dataclass(frozen=True)
class UpperQuotientPathGeometry(RankedPairPathGeometry):
    """Ranked path after contracting both endpoint critical blocks."""

    first_critical_block_size: int
    second_critical_block_size: int

    def unit_bucket_count(self) -> int:
        """Count buckets consisting only of the winning Kruskal edge."""

        return sum(size == 1 for size in self.bucket_sizes)

    def exact_size_two_bucket_count(self) -> int:
        """Count quotient-path buckets of physical size two."""

        return sum(size == 2 for size in self.bucket_sizes)


@dataclass(frozen=True)
class EnvironmentDiagnostic:
    """Pair means within one independent rank environment."""

    repetition: int
    final_root: int
    pair_count: int
    final_giant_fraction: float
    second_largest_final_component_fraction: float
    critical_block_count: float
    mean_critical_block_size: float
    same_critical_block_pair_probability_in_giant: float
    largest_critical_block_fraction_of_giant: float
    mean_endpoint_critical_block_size: float
    mean_pair_distance: float
    mean_upper_quotient_node_count: float
    mean_upper_quotient_bucket_size: float
    mean_upper_quotient_maximum_bucket_size: float
    mean_upper_quotient_lca_bucket_size: float
    mean_upper_quotient_lca_rank: float
    mean_upper_quotient_unit_bucket_count: float
    mean_upper_quotient_exact_size_two_bucket_count: float
    mean_upper_quotient_nontrivial_bucket_count_at_most_four: float
    mean_upper_quotient_nontrivial_bucket_count_at_most_eight: float
    mean_path_fac_local_oracle_attenuation: float
    mean_path_fac_local_oracle_correlation: float
    mean_path_fac_local_oracle_same_relation_probability: float


@dataclass(frozen=True)
class ScalarEstimate:
    """Mean and standard error across independent rank environments."""

    mean: float
    standard_error: float | None


@dataclass(frozen=True)
class GiantQuotientDiagnosticSummary:
    """Reproducible environment-clustered summary at one torus size."""

    side_length: int
    vertex_count: int
    p: float
    critical_rank: float
    critical_beta: float
    final_rank: float
    final_beta: float
    repetitions: int
    eligible_environment_count: int
    ineligible_environment_count: int
    pairs_per_environment: int
    pair_count: int
    distance_fraction: float
    seed: int
    geometry_estimates: dict[str, ScalarEstimate]
    path_fac_local_oracle_nonproof_estimates: dict[str, ScalarEstimate]
    uncertainty_note: str
    estimates_conditioned_on_eligible_environments: bool
    path_fac_status: str
    weak_recovery_claimed: bool


def _validate_p(p: float) -> float:
    """Return ``q_1`` when the final forest is strictly supercritical."""

    if not P_SW < p < 1.0:
        raise ValueError(f"p must satisfy {P_SW} < p < 1")
    final_rank = 2.0 * p - 1.0
    if final_rank <= Q_CRITICAL:
        raise ValueError("the final rank must lie strictly above q_c")
    return final_rank


def _tree_root(forest: CriticalKruskalForest, vertex: int) -> int:
    """Return the dendrogram root containing one leaf."""

    if not 0 <= vertex < forest.vertex_count:
        raise ValueError("vertex is outside the torus")
    node = vertex
    while forest.tree_parent[node] != -1:
        node = forest.tree_parent[node]
    return node


def largest_final_component(
    forest: CriticalKruskalForest,
) -> tuple[int, ...]:
    """Choose the largest final component, breaking ties by its first leaf."""

    if not forest.components:
        raise AssertionError("a finite forest must contain a component")
    return max(
        forest.components,
        key=lambda component: (len(component), -min(component)),
    )


def critical_partition_of_largest_tree(
    forest: CriticalKruskalForest,
    critical_rank: float = Q_CRITICAL,
) -> CriticalPartition:
    """Cut the largest final tree at ``critical_rank``."""

    if not 0.0 <= critical_rank < 1.0:
        raise ValueError("critical_rank must belong to [0, 1)")
    component = largest_final_component(forest)
    final_roots = {_tree_root(forest, vertex) for vertex in component}
    if len(final_roots) != 1:
        raise AssertionError("one final component has several tree roots")
    final_root = next(iter(final_roots))

    vertices_by_block_root: dict[int, list[int]] = {}
    for vertex in component:
        node = vertex
        parent = forest.tree_parent[node]
        while parent != -1 and forest.merge_rank[parent] <= critical_rank:
            node = parent
            parent = forest.tree_parent[node]
        vertices_by_block_root.setdefault(node, []).append(vertex)

    ordered = sorted(
        vertices_by_block_root.items(),
        key=lambda item: (min(item[1]), len(item[1])),
    )
    block_roots = tuple(root for root, _ in ordered)
    blocks = tuple(tuple(vertices) for _, vertices in ordered)
    membership = [-1] * forest.vertex_count
    for block_index, block in enumerate(blocks):
        for vertex in block:
            membership[vertex] = block_index
    if sum(len(block) for block in blocks) != len(component):
        raise AssertionError("critical blocks do not partition the final tree")
    return CriticalPartition(
        final_root=final_root,
        component=component,
        block_roots=block_roots,
        blocks=blocks,
        block_index_by_vertex=tuple(membership),
    )


def sample_uniform_far_cross_block_pair(
    partition: CriticalPartition,
    side_length: int,
    distance_fraction: float,
    rng: random.Random,
    maximum_attempts: int = 100_000,
) -> tuple[int, int]:
    """Sample uniformly conditional on distance and distinct critical blocks.

    The proposal is uniform on ordered distinct pairs in ``R_*``.  Rejection
    therefore gives the exact conditional uniform law requested here.
    """

    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0, 1)")
    if maximum_attempts <= 0:
        raise ValueError("maximum_attempts must be positive")
    vertices = partition.component
    if len(vertices) < 2 or len(partition.blocks) < 2:
        raise NoEligibleCrossBlockPairError(
            "the largest final tree has no cross-block pair"
        )

    minimum_distance = distance_fraction * side_length
    for _ in range(maximum_attempts):
        first_index = rng.randrange(len(vertices))
        second_index = rng.randrange(len(vertices) - 1)
        if second_index >= first_index:
            second_index += 1
        first = vertices[first_index]
        second = vertices[second_index]
        if partition.block_index(first) == partition.block_index(second):
            continue
        if (
            triangular_torus_distance(first, second, side_length)
            >= minimum_distance
        ):
            return first, second
    raise NoEligibleCrossBlockPairError(
        "failed to sample a distant pair from distinct critical blocks"
    )


def upper_quotient_path_geometry(
    forest: CriticalKruskalForest,
    partition: CriticalPartition,
    first: int,
    second: int,
    critical_rank: float = Q_CRITICAL,
) -> UpperQuotientPathGeometry:
    """Return only the strict postcritical nodes on the endpoint path."""

    if partition.block_index(first) == partition.block_index(second):
        raise ValueError("endpoints must lie in distinct critical blocks")
    if not forest.connected(first, second):
        raise ValueError("endpoints must lie in one final tree")
    path_nodes = forest.pair_path_nodes(first, second)
    upper_nodes = tuple(
        node
        for node in path_nodes
        if forest.merge_rank[node] > critical_rank
    )
    if not upper_nodes:
        raise AssertionError("a cross-block path must contain an upper node")
    lca = forest.tree_lca(first, second)
    if upper_nodes[-1] != lca:
        raise AssertionError("the upper quotient path must end at the LCA")
    return UpperQuotientPathGeometry(
        distance=triangular_torus_distance(
            first, second, forest.side_length
        ),
        component_size=len(partition.component),
        lca_rank=forest.merge_rank[lca],
        path_length=len(upper_nodes),
        bucket_sizes=tuple(forest.bucket_size[node] for node in upper_nodes),
        bucket_ranks=tuple(forest.merge_rank[node] for node in upper_nodes),
        first_critical_block_size=partition.block_size(first),
        second_critical_block_size=partition.block_size(second),
    )


def diagnose_environment(
    *,
    forest: CriticalKruskalForest,
    repetition: int,
    p: float,
    pairs_per_environment: int,
    distance_fraction: float,
    rng: random.Random,
) -> EnvironmentDiagnostic:
    """Compute pair means in one final-forest environment."""

    _validate_p(p)
    if pairs_per_environment <= 0:
        raise ValueError("pairs_per_environment must be positive")
    partition = critical_partition_of_largest_tree(forest)
    paths: list[UpperQuotientPathGeometry] = []
    attenuations = []
    correlations = []
    for _ in range(pairs_per_environment):
        first, second = sample_uniform_far_cross_block_pair(
            partition,
            forest.side_length,
            distance_fraction,
            rng,
        )
        path = upper_quotient_path_geometry(
            forest, partition, first, second
        )
        evaluation = evaluate_ranked_geometry(path, p)
        paths.append(path)
        attenuations.append(evaluation.attenuation)
        correlations.append(evaluation.correlation)

    critical_block_sizes = tuple(len(block) for block in partition.blocks)
    ordered_final_sizes = sorted(
        (len(component) for component in forest.components),
        reverse=True,
    )
    second_largest_final_size = (
        ordered_final_sizes[1] if len(ordered_final_sizes) >= 2 else 0
    )
    endpoint_block_sizes = tuple(
        0.5
        * (
            path.first_critical_block_size
            + path.second_critical_block_size
        )
        for path in paths
    )
    return EnvironmentDiagnostic(
        repetition=repetition,
        final_root=partition.final_root,
        pair_count=len(paths),
        final_giant_fraction=len(partition.component) / forest.vertex_count,
        second_largest_final_component_fraction=(
            second_largest_final_size / forest.vertex_count
        ),
        critical_block_count=float(len(partition.blocks)),
        mean_critical_block_size=fmean(critical_block_sizes),
        same_critical_block_pair_probability_in_giant=(
            fsum(
                size * (size - 1)
                for size in critical_block_sizes
            )
            / (
                len(partition.component)
                * (len(partition.component) - 1)
            )
        ),
        largest_critical_block_fraction_of_giant=(
            max(critical_block_sizes) / len(partition.component)
        ),
        mean_endpoint_critical_block_size=fmean(endpoint_block_sizes),
        mean_pair_distance=fmean(path.distance for path in paths),
        mean_upper_quotient_node_count=fmean(
            path.path_length for path in paths
        ),
        mean_upper_quotient_bucket_size=fmean(
            fmean(path.bucket_sizes) for path in paths
        ),
        mean_upper_quotient_maximum_bucket_size=fmean(
            max(path.bucket_sizes) for path in paths
        ),
        mean_upper_quotient_lca_bucket_size=fmean(
            path.bucket_sizes[-1] for path in paths
        ),
        mean_upper_quotient_lca_rank=fmean(
            path.lca_rank for path in paths
        ),
        mean_upper_quotient_unit_bucket_count=fmean(
            path.unit_bucket_count() for path in paths
        ),
        mean_upper_quotient_exact_size_two_bucket_count=fmean(
            path.exact_size_two_bucket_count() for path in paths
        ),
        mean_upper_quotient_nontrivial_bucket_count_at_most_four=fmean(
            path.bounded_bucket_count(4) for path in paths
        ),
        mean_upper_quotient_nontrivial_bucket_count_at_most_eight=fmean(
            path.bounded_bucket_count(8) for path in paths
        ),
        mean_path_fac_local_oracle_attenuation=fmean(attenuations),
        mean_path_fac_local_oracle_correlation=fmean(correlations),
        mean_path_fac_local_oracle_same_relation_probability=0.5
        * (1.0 + fmean(correlations)),
    )


def mean_with_environment_standard_error(
    values: Sequence[float],
) -> ScalarEstimate:
    """Aggregate independent environment values without pair pseudo-replication."""

    if not values:
        raise ValueError("at least one environment value is required")
    mean = fmean(values)
    if len(values) < 2:
        standard_error = None
    else:
        standard_error = sqrt(
            fsum((value - mean) ** 2 for value in values)
            / (len(values) * (len(values) - 1))
        )
    return ScalarEstimate(mean=mean, standard_error=standard_error)


GEOMETRY_METRICS = (
    "final_giant_fraction",
    "second_largest_final_component_fraction",
    "critical_block_count",
    "mean_critical_block_size",
    "same_critical_block_pair_probability_in_giant",
    "largest_critical_block_fraction_of_giant",
    "mean_endpoint_critical_block_size",
    "mean_pair_distance",
    "mean_upper_quotient_node_count",
    "mean_upper_quotient_bucket_size",
    "mean_upper_quotient_maximum_bucket_size",
    "mean_upper_quotient_lca_bucket_size",
    "mean_upper_quotient_lca_rank",
    "mean_upper_quotient_unit_bucket_count",
    "mean_upper_quotient_exact_size_two_bucket_count",
    "mean_upper_quotient_nontrivial_bucket_count_at_most_four",
    "mean_upper_quotient_nontrivial_bucket_count_at_most_eight",
)

PATH_FAC_METRICS = (
    "mean_path_fac_local_oracle_attenuation",
    "mean_path_fac_local_oracle_correlation",
    "mean_path_fac_local_oracle_same_relation_probability",
)


def _aggregate_metrics(
    environments: Sequence[EnvironmentDiagnostic],
    metric_names: Sequence[str],
) -> dict[str, ScalarEstimate]:
    return {
        name: mean_with_environment_standard_error(
            tuple(float(getattr(environment, name)) for environment in environments)
        )
        for name in metric_names
    }


def run_diagnostic(
    *,
    side_length: int,
    repetitions: int,
    pairs_per_environment: int,
    p: float,
    distance_fraction: float,
    seed: int,
) -> GiantQuotientDiagnosticSummary:
    """Run the final-giant quotient diagnostic at one torus size."""

    final_rank = _validate_p(p)
    if side_length < 4:
        raise ValueError("side_length must be at least 4")
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if pairs_per_environment <= 0:
        raise ValueError("pairs_per_environment must be positive")
    if not 0.0 < distance_fraction < 1.0:
        raise ValueError("distance_fraction must belong to (0, 1)")

    seed_rng = random.Random(seed)
    environments = []
    for repetition in range(repetitions):
        rank_rng = random.Random(seed_rng.getrandbits(64))
        pair_rng = random.Random(seed_rng.getrandbits(64))
        ranked_edges: tuple[RankedEdge, ...] = sample_ranked_edges(
            side_length, rank_rng
        )
        forest = CriticalKruskalForest(
            side_length,
            ranked_edges,
            critical_rank=final_rank,
        )
        try:
            environments.append(
                diagnose_environment(
                    forest=forest,
                    repetition=repetition,
                    p=p,
                    pairs_per_environment=pairs_per_environment,
                    distance_fraction=distance_fraction,
                    rng=pair_rng,
                )
            )
        except NoEligibleCrossBlockPairError:
            continue

    geometry_estimates = (
        _aggregate_metrics(environments, GEOMETRY_METRICS)
        if environments
        else {}
    )
    path_fac_estimates = (
        _aggregate_metrics(environments, PATH_FAC_METRICS)
        if environments
        else {}
    )

    return GiantQuotientDiagnosticSummary(
        side_length=side_length,
        vertex_count=side_length * side_length,
        p=p,
        critical_rank=Q_CRITICAL,
        critical_beta=clock_time_from_rank(Q_CRITICAL, p),
        final_rank=final_rank,
        final_beta=1.0,
        repetitions=repetitions,
        eligible_environment_count=len(environments),
        ineligible_environment_count=repetitions - len(environments),
        pairs_per_environment=pairs_per_environment,
        pair_count=len(environments) * pairs_per_environment,
        distance_fraction=distance_fraction,
        seed=seed,
        geometry_estimates=geometry_estimates,
        path_fac_local_oracle_nonproof_estimates=path_fac_estimates,
        uncertainty_note=(
            "standard errors use independent rank-environment means; sampled "
            "pairs within one environment are not treated as iid clusters; "
            "environments with no eligible distant cross-block pair are "
            "counted explicitly and excluded from the conditional estimates"
        ),
        estimates_conditioned_on_eligible_environments=True,
        path_fac_status=PATH_FAC_STATUS,
        weak_recovery_claimed=False,
    )


def _parse_side_lengths(raw: str) -> tuple[int, ...]:
    result = tuple(int(value) for value in raw.split(","))
    if not result:
        raise ValueError("at least one side length is required")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sides", default="16,32,64")
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--pairs", type=int, default=100)
    parser.add_argument("--p", type=float, default=BASELINE_P)
    parser.add_argument("--distance-fraction", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260726)
    arguments = parser.parse_args()

    summaries = []
    for side_length in _parse_side_lengths(arguments.sides):
        side_seed = (
            arguments.seed
            ^ (side_length * 0x9E3779B97F4A7C15)
        ) & ((1 << 64) - 1)
        summaries.append(
            run_diagnostic(
                side_length=side_length,
                repetitions=arguments.repetitions,
                pairs_per_environment=arguments.pairs,
                p=arguments.p,
                distance_fraction=arguments.distance_fraction,
                seed=side_seed,
            )
        )
    payload = {
        "diagnostic": "final-giant critical-quotient PATH-FAC oracle/non-proof",
        "summaries": [asdict(summary) for summary in summaries],
    }
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
