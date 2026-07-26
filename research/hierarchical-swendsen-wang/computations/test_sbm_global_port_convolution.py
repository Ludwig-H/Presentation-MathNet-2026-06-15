"""Tests for the exact full-D ferromagnetic-SBM global port."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

from sbm_global_port_convolution import (
    AUDIT_PASS_STATUS,
    AUDIT_SKIPPED_STATUS,
    EXACT_STATUS,
    audit_global_port_convolution,
    balanced_partition_function,
    direct_orientation_enumeration,
    finite_nonedge_field,
    iid_partition_function,
    magnetization_multiplicities,
    port_partition_functions,
    run_diagnostic,
)


class MagnetizationConvolutionTests(unittest.TestCase):
    def test_two_singletons_have_the_expected_coefficients(self) -> None:
        self.assertEqual(
            magnetization_multiplicities((1, 1)),
            {-2: 1, 0: 2, 2: 1},
        )

    def test_coefficients_are_symmetric_and_sum_to_two_to_r(self) -> None:
        root_sizes = (4, 3, 2, 1)
        multiplicities = magnetization_multiplicities(root_sizes)
        self.assertEqual(sum(multiplicities.values()), 2 ** len(root_sizes))
        for magnetization, multiplicity in multiplicities.items():
            self.assertEqual(
                multiplicity,
                multiplicities[-magnetization],
            )

    def test_convolution_matches_direct_counts_on_small_cases(self) -> None:
        for root_sizes in (
            (1,),
            (1, 1),
            (3, 2, 1),
            (4, 3, 2, 1),
            (2, 2, 2, 2),
        ):
            with self.subTest(root_sizes=root_sizes):
                convolved = magnetization_multiplicities(root_sizes)
                direct, _, _ = direct_orientation_enumeration(
                    root_sizes,
                    -0.1,
                )
                self.assertEqual(convolved, direct)


class PartitionFunctionTests(unittest.TestCase):
    def test_finite_nonedge_field_has_the_exact_negative_value(self) -> None:
        h0 = finite_nonedge_field(10, 4.0, 1.0)
        self.assertAlmostEqual(
            h0,
            0.5 * math.log((1.0 - 4.0 / 10.0) / (1.0 - 1.0 / 10.0)),
        )
        self.assertLess(h0, 0.0)

    def test_two_singleton_partition_functions_have_closed_forms(self) -> None:
        h0 = -0.2
        multiplicities = magnetization_multiplicities((1, 1))
        self.assertAlmostEqual(
            iid_partition_function(multiplicities, h0),
            2.0 + 2.0 * math.exp(2.0 * h0),
        )
        self.assertEqual(
            balanced_partition_function(multiplicities),
            2,
        )
        partitions = port_partition_functions((1, 1), h0)
        self.assertAlmostEqual(
            partitions.iid_partition_normalized_by_uniform_orientations,
            0.5 + 0.5 * math.exp(2.0 * h0),
        )
        self.assertEqual(partitions.balanced_partition_function, 2)
        self.assertAlmostEqual(
            partitions.balanced_partition_normalized_by_uniform_orientations,
            0.5,
        )

    def test_odd_total_size_has_no_balanced_assignment(self) -> None:
        partitions = port_partition_functions((2, 1), -0.1)
        self.assertFalse(partitions.balanced_model_feasible)
        self.assertEqual(partitions.balanced_partition_function, 0)
        self.assertIsNone(partitions.balanced_log_partition_function)
        self.assertEqual(
            partitions.balanced_partition_normalized_by_uniform_orientations,
            0.0,
        )

    def test_even_total_size_need_not_admit_a_balanced_root_signing(self) -> None:
        partitions = port_partition_functions((4, 2), -0.1)
        self.assertFalse(partitions.balanced_model_feasible)
        self.assertEqual(partitions.balanced_partition_function, 0)


class EnumerationAuditTests(unittest.TestCase):
    def test_exhaustive_audit_passes_exactly(self) -> None:
        for root_sizes in ((3, 2, 1), (4, 2, 2, 1, 1)):
            with self.subTest(root_sizes=root_sizes):
                audit = audit_global_port_convolution(
                    root_sizes,
                    -0.13,
                    max_enumerated_roots=10,
                )
                self.assertTrue(audit.enumeration_performed)
                self.assertTrue(audit.multiplicities_match_exactly)
                self.assertTrue(audit.balanced_partition_matches_exactly)
                self.assertLessEqual(
                    audit.iid_partition_absolute_error,
                    1e-12 * (2 ** len(root_sizes)),
                )
                self.assertEqual(audit.audit_status, AUDIT_PASS_STATUS)

    def test_enumeration_is_skipped_above_the_explicit_limit(self) -> None:
        audit = audit_global_port_convolution(
            (1, 1, 1),
            -0.1,
            max_enumerated_roots=2,
        )
        self.assertFalse(audit.enumeration_performed)
        self.assertIsNone(audit.multiplicities_match_exactly)
        self.assertEqual(audit.audit_status, AUDIT_SKIPPED_STATUS)


class CombinedDiagnosticTests(unittest.TestCase):
    def test_diagnostic_states_assumptions_and_nonclaims(self) -> None:
        diagnostic = run_diagnostic(
            root_sizes=(3, 2, 1),
            a_n=4.0,
            b_n=1.0,
        )
        self.assertEqual(diagnostic.diagnostic_status, EXACT_STATUS)
        self.assertEqual(diagnostic.n, 6)
        self.assertTrue(diagnostic.roots_are_monochromatic)
        self.assertTrue(diagnostic.common_internal_root_factors_omitted)
        self.assertFalse(diagnostic.dendrogram_law_sampled)
        self.assertFalse(diagnostic.recovery_threshold_claimed)
        self.assertEqual(
            sum(
                term.orientation_multiplicity
                for term in diagnostic.magnetization_terms
            ),
            2 ** diagnostic.root_count,
        )
        self.assertEqual(
            diagnostic.direct_enumeration_audit.audit_status,
            AUDIT_PASS_STATUS,
        )

    def test_cli_emits_json_with_a_successful_identity_audit(self) -> None:
        module_path = Path(__file__).with_name(
            "sbm_global_port_convolution.py"
        )
        completed = subprocess.run(
            (
                sys.executable,
                str(module_path),
                "--root-sizes",
                "3",
                "2",
                "1",
                "--a",
                "4",
                "--b",
                "1",
            ),
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["diagnostic_status"], EXACT_STATUS)
        self.assertEqual(
            payload["direct_enumeration_audit"]["audit_status"],
            AUDIT_PASS_STATUS,
        )
        self.assertFalse(payload["dendrogram_law_sampled"])
        self.assertFalse(payload["recovery_threshold_claimed"])


class ValidationTests(unittest.TestCase):
    def test_invalid_root_sizes_are_rejected(self) -> None:
        for root_sizes in ((), (0, 1), (-1, 2), (True, 1), (1.5, 2)):
            with self.subTest(root_sizes=root_sizes):
                with self.assertRaises(ValueError):
                    magnetization_multiplicities(root_sizes)

    def test_invalid_sbm_parameters_are_rejected(self) -> None:
        for n, a_n, b_n in (
            (1, 0.5, 0.1),
            (10, 1.0, 1.0),
            (10, 0.5, 1.0),
            (10, 10.0, 1.0),
            (10, 4.0, 0.0),
        ):
            with self.subTest(n=n, a_n=a_n, b_n=b_n):
                with self.assertRaises(ValueError):
                    finite_nonedge_field(n, a_n, b_n)

    def test_positive_port_field_is_rejected(self) -> None:
        multiplicities = magnetization_multiplicities((1, 1))
        with self.assertRaises(ValueError):
            iid_partition_function(multiplicities, 0.1)
        with self.assertRaises(ValueError):
            direct_orientation_enumeration((1, 1), 0.1)


if __name__ == "__main__":
    unittest.main()
