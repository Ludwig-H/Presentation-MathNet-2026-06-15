"""Unselected population audit of exact two-step hierarchical L2 cells.

For each realized side-four hierarchy, this module samples one pair at
maximal torus distance independently of the dendrogram.  If the pair lies in
one final root, it enumerates every pair of consecutive strict-arm merger
updates.  Each cell is then passed to :mod:`two_step_projective_l2_cell`, so
the full posterior and every attained exterior coset are treated exactly.

The sampling convention removes the post-hoc choice used to find the D1
witness.  It is still only a tiny-volume diagnostic: cells are pooled by
occurrence along the sampled pair, overlap one another, and are neither
independent nor an asymptotic Palm sample.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from math import fsum
from statistics import fmean, median
from typing import Sequence

from critical_band_thresholds import Q_CRITICAL
from critical_pair_path_geometry import sample_ranked_edges
from joint_hierarchical_sweep import HierarchicalSweepEnvironment
from nested_projection_l2_diagnostic import _sample_uniform_distant_pair
from two_step_projective_l2_cell import analyze_two_step_cell


@dataclass(frozen=True)
class PopulationCellDiagnostic:
    """One unselected consecutive strict-arm cell."""

    repetition: int
    pair_first: int
    pair_second: int
    pair_distance: int
    pair_lca: int
    lower_node: int
    upper_node: int
    lower_rank: float
    upper_rank: float
    upper_beta: float
    spine_size: int
    attachment_size: int
    ancestral_sibling_size: int
    lower_bucket_size: int
    upper_bucket_size: int
    incoming_energy: float
    second_absolute_loss: float
    second_relative_loss: float
    strictly_positive_exterior_posterior_mass: float
    strictly_positive_exterior_energy_share: float
    strictly_positive_exterior_energy_ratio: float | None
    minimum_strictly_positive_exterior_second_relative_loss: float | None
    median_strictly_positive_exterior_second_relative_loss: float | None
    minimum_all_attainable_second_relative_loss: float
    maximum_pythagorean_error: float
    maximum_projective_factorization_error: float


@dataclass(frozen=True)
class TwoStepPopulationSummary:
    """Finite population statistics under the declared pair sampling law."""

    side_length: int
    repetitions: int
    p: float
    distance_fraction: float
    critical_rank: float
    critical_rank_window: float
    seed: int
    connected_pair_count: int
    pair_with_eligible_cell_count: int
    eligible_cell_count: int
    mean_eligible_cells_per_connected_pair: float | None
    energy_weighted_second_relative_loss: float | None
    mean_cell_second_relative_loss: float | None
    median_cell_second_relative_loss: float | None
    ninetieth_percentile_cell_second_relative_loss: float | None
    positive_second_loss_fraction: float | None
    second_loss_above_one_percent_fraction: float | None
    top_ten_percent_second_absolute_loss_share: float | None
    top_twenty_percent_second_absolute_loss_share: float | None
    near_critical_cell_count: int
    near_critical_cell_fraction: float | None
    near_critical_incoming_energy_share: float | None
    near_critical_second_absolute_loss_share: float | None
    near_critical_energy_weighted_second_relative_loss: float | None
    median_near_critical_second_relative_loss: float | None
    mean_strictly_positive_exterior_posterior_mass: float | None
    energy_weighted_strictly_positive_exterior_share: float | None
    boundary_zero_margin_cell_fraction: float | None
    cell_with_strictly_positive_exterior_potential_count: int
    interior_margin_above_one_per_mille_fraction_among_interior_cells: float | None
    interior_margin_above_one_percent_fraction_among_interior_cells: float | None
    near_critical_boundary_zero_margin_cell_fraction: float | None
    near_critical_energy_weighted_strictly_positive_exterior_share: float | None
    maximum_pythagorean_error: float
    maximum_projective_factorization_error: float
    pair_sample_is_independent_of_dendrogram: bool
    witness_selected_after_scan: bool
    exact_full_posterior_enumeration: bool
    weak_recovery_claimed: bool
    interpretation: str


def eligible_two_step_nodes(
    environment: HierarchicalSweepEnvironment,
    first: int,
    second: int,
) -> tuple[int, ...]:
    """List lower nodes whose parent is also strictly before the pair LCA."""

    forest = environment.forest
    if first == second or not forest.connected(first, second):
        return ()
    lca = forest.tree_lca(first, second)
    result = []
    for endpoint in (first, second):
        node = endpoint
        while True:
            node = forest.tree_parent[node]
            if node == -1 or node == lca:
                break
            if node in forest.internal_nodes and forest.tree_parent[node] != lca:
                result.append(node)
    if len(result) != len(set(result)):
        raise AssertionError("the two strict arms must be disjoint below the LCA")
    return tuple(result)


def _linear_quantile(values: Sequence[float], probability: float) -> float:
    if not values:
        raise ValueError("a quantile requires at least one value")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must belong to [0,1]")
    ordered = sorted(values)
    location = probability * (len(ordered) - 1)
    lower = int(location)
    upper = math.ceil(location)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (location - lower) * (ordered[upper] - ordered[lower])


def _upper_tail_share(values: Sequence[float], fraction: float) -> float | None:
    total = fsum(values)
    if not values or total == 0.0:
        return None
    count = max(1, math.ceil(fraction * len(values)))
    return fsum(sorted(values, reverse=True)[:count]) / total


def _summarize(
    cells: Sequence[PopulationCellDiagnostic],
    *,
    side_length: int,
    repetitions: int,
    p: float,
    distance_fraction: float,
    critical_rank_window: float,
    seed: int,
    connected_pair_count: int,
    pair_with_eligible_cell_count: int,
) -> TwoStepPopulationSummary:
    incoming = tuple(cell.incoming_energy for cell in cells)
    losses = tuple(cell.second_absolute_loss for cell in cells)
    relative = tuple(cell.second_relative_loss for cell in cells)
    total_incoming = fsum(incoming)
    total_loss = fsum(losses)
    near_critical = tuple(
        cell
        for cell in cells
        if abs(cell.upper_rank - Q_CRITICAL) <= critical_rank_window
    )
    near_critical_incoming = fsum(cell.incoming_energy for cell in near_critical)
    near_critical_loss = fsum(cell.second_absolute_loss for cell in near_critical)
    near_critical_interior_energy = fsum(
        cell.incoming_energy * cell.strictly_positive_exterior_energy_share
        for cell in near_critical
    )
    interior_energy = fsum(
        cell.incoming_energy * cell.strictly_positive_exterior_energy_share
        for cell in cells
    )
    interior_minima = tuple(
        cell.minimum_strictly_positive_exterior_second_relative_loss
        for cell in cells
        if cell.minimum_strictly_positive_exterior_second_relative_loss is not None
    )
    return TwoStepPopulationSummary(
        side_length=side_length,
        repetitions=repetitions,
        p=p,
        distance_fraction=distance_fraction,
        critical_rank=Q_CRITICAL,
        critical_rank_window=critical_rank_window,
        seed=seed,
        connected_pair_count=connected_pair_count,
        pair_with_eligible_cell_count=pair_with_eligible_cell_count,
        eligible_cell_count=len(cells),
        mean_eligible_cells_per_connected_pair=(
            None if connected_pair_count == 0 else len(cells) / connected_pair_count
        ),
        energy_weighted_second_relative_loss=(
            None if total_incoming == 0.0 else fsum(losses) / total_incoming
        ),
        mean_cell_second_relative_loss=(None if not relative else fmean(relative)),
        median_cell_second_relative_loss=(None if not relative else median(relative)),
        ninetieth_percentile_cell_second_relative_loss=(
            None if not relative else _linear_quantile(relative, 0.9)
        ),
        positive_second_loss_fraction=(
            None if not losses else sum(loss > 1e-14 for loss in losses) / len(losses)
        ),
        second_loss_above_one_percent_fraction=(
            None
            if not relative
            else sum(loss > 0.01 for loss in relative) / len(relative)
        ),
        top_ten_percent_second_absolute_loss_share=_upper_tail_share(losses, 0.1),
        top_twenty_percent_second_absolute_loss_share=_upper_tail_share(losses, 0.2),
        near_critical_cell_count=len(near_critical),
        near_critical_cell_fraction=(
            None if not cells else len(near_critical) / len(cells)
        ),
        near_critical_incoming_energy_share=(
            None if total_incoming == 0.0 else near_critical_incoming / total_incoming
        ),
        near_critical_second_absolute_loss_share=(
            None if total_loss == 0.0 else near_critical_loss / total_loss
        ),
        near_critical_energy_weighted_second_relative_loss=(
            None
            if near_critical_incoming == 0.0
            else near_critical_loss / near_critical_incoming
        ),
        median_near_critical_second_relative_loss=(
            None
            if not near_critical
            else median(cell.second_relative_loss for cell in near_critical)
        ),
        mean_strictly_positive_exterior_posterior_mass=(
            None
            if not cells
            else fmean(cell.strictly_positive_exterior_posterior_mass for cell in cells)
        ),
        energy_weighted_strictly_positive_exterior_share=(
            None if total_incoming == 0.0 else interior_energy / total_incoming
        ),
        boundary_zero_margin_cell_fraction=(
            None
            if not cells
            else sum(
                cell.minimum_all_attainable_second_relative_loss <= 1e-14
                for cell in cells
            )
            / len(cells)
        ),
        cell_with_strictly_positive_exterior_potential_count=len(interior_minima),
        interior_margin_above_one_per_mille_fraction_among_interior_cells=(
            None
            if not interior_minima
            else sum(value > 0.001 for value in interior_minima) / len(interior_minima)
        ),
        interior_margin_above_one_percent_fraction_among_interior_cells=(
            None
            if not interior_minima
            else sum(value > 0.01 for value in interior_minima) / len(interior_minima)
        ),
        near_critical_boundary_zero_margin_cell_fraction=(
            None
            if not near_critical
            else sum(
                cell.minimum_all_attainable_second_relative_loss <= 1e-14
                for cell in near_critical
            )
            / len(near_critical)
        ),
        near_critical_energy_weighted_strictly_positive_exterior_share=(
            None
            if near_critical_incoming == 0.0
            else near_critical_interior_energy / near_critical_incoming
        ),
        maximum_pythagorean_error=max(
            (cell.maximum_pythagorean_error for cell in cells), default=0.0
        ),
        maximum_projective_factorization_error=max(
            (cell.maximum_projective_factorization_error for cell in cells),
            default=0.0,
        ),
        pair_sample_is_independent_of_dendrogram=True,
        witness_selected_after_scan=False,
        exact_full_posterior_enumeration=True,
        weak_recovery_claimed=False,
        interpretation=(
            "finite unselected maximal-distance pair audit; cells are pooled by "
            "strict-arm occurrence, overlap, and yield no thermodynamic or Palm "
            "conclusion"
        ),
    )


def run_population_diagnostic(
    *,
    side_length: int = 4,
    repetitions: int = 12,
    p: float = 0.805,
    distance_fraction: float = 0.5,
    critical_rank_window: float = 0.02,
    seed: int = 20260729,
) -> tuple[TwoStepPopulationSummary, tuple[PopulationCellDiagnostic, ...]]:
    if side_length != 4:
        raise ValueError("exact posterior enumeration requires side_length=4")
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if critical_rank_window <= 0.0:
        raise ValueError("critical_rank_window must be positive")
    rng = random.Random(seed)
    cells = []
    connected_pair_count = 0
    pair_with_eligible_cell_count = 0
    for repetition in range(repetitions):
        environment = HierarchicalSweepEnvironment(
            side_length, sample_ranked_edges(side_length, rng), p
        )
        first, second = _sample_uniform_distant_pair(
            side_length, distance_fraction, rng
        )
        if not environment.forest.connected(first, second):
            continue
        connected_pair_count += 1
        nodes = eligible_two_step_nodes(environment, first, second)
        if nodes:
            pair_with_eligible_cell_count += 1
        for lower_node in nodes:
            result = analyze_two_step_cell(
                environment,
                first,
                second,
                lower_node,
                witness_selected_after_exploratory_scan=False,
            )
            if result.witness_selected_after_exploratory_scan:
                raise AssertionError("population cells cannot be post-hoc witnesses")
            cells.append(
                PopulationCellDiagnostic(
                    repetition=repetition,
                    pair_first=first,
                    pair_second=second,
                    pair_distance=result.geometry.pair_distance,
                    pair_lca=result.geometry.pair_lca,
                    lower_node=result.geometry.lower_node,
                    upper_node=result.geometry.upper_node,
                    lower_rank=result.geometry.lower_rank,
                    upper_rank=result.geometry.upper_rank,
                    upper_beta=result.geometry.upper_beta,
                    spine_size=result.geometry.spine_size,
                    attachment_size=result.geometry.attachment_size,
                    ancestral_sibling_size=(result.geometry.ancestral_sibling_size),
                    lower_bucket_size=result.geometry.lower_bucket_size,
                    upper_bucket_size=result.geometry.upper_bucket_size,
                    incoming_energy=result.after_first_norm_square,
                    second_absolute_loss=result.second_absolute_loss,
                    second_relative_loss=result.second_relative_loss,
                    strictly_positive_exterior_posterior_mass=(
                        result.strictly_positive_exterior_posterior_mass
                    ),
                    strictly_positive_exterior_energy_share=(
                        0.0
                        if result.after_first_norm_square == 0.0
                        else result.strictly_positive_exterior_incoming_energy
                        / result.after_first_norm_square
                    ),
                    strictly_positive_exterior_energy_ratio=(
                        result.strictly_positive_exterior_energy_ratio
                    ),
                    minimum_strictly_positive_exterior_second_relative_loss=(
                        result.minimum_strictly_positive_exterior_second_relative_loss
                    ),
                    median_strictly_positive_exterior_second_relative_loss=(
                        result.median_strictly_positive_exterior_second_relative_loss
                    ),
                    minimum_all_attainable_second_relative_loss=(
                        result.minimum_all_attainable_second_relative_loss
                    ),
                    maximum_pythagorean_error=result.maximum_pythagorean_error,
                    maximum_projective_factorization_error=(
                        result.maximum_projective_factorization_error
                    ),
                )
            )
    result_cells = tuple(cells)
    return (
        _summarize(
            result_cells,
            side_length=side_length,
            repetitions=repetitions,
            p=p,
            distance_fraction=distance_fraction,
            critical_rank_window=critical_rank_window,
            seed=seed,
            connected_pair_count=connected_pair_count,
            pair_with_eligible_cell_count=pair_with_eligible_cell_count,
        ),
        result_cells,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=4)
    parser.add_argument("--repetitions", type=int, default=12)
    parser.add_argument("--p", type=float, default=0.805)
    parser.add_argument("--distance-fraction", type=float, default=0.5)
    parser.add_argument("--critical-rank-window", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--include-cells", action="store_true")
    arguments = parser.parse_args()
    summary, cells = run_population_diagnostic(
        side_length=arguments.side,
        repetitions=arguments.repetitions,
        p=arguments.p,
        distance_fraction=arguments.distance_fraction,
        critical_rank_window=arguments.critical_rank_window,
        seed=arguments.seed,
    )
    payload: dict[str, object] = {"summary": asdict(summary)}
    if arguments.include_cells:
        payload["cells"] = [asdict(cell) for cell in cells]
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
