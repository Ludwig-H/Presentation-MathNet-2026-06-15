from __future__ import annotations

import unittest
from fractions import Fraction
from math import comb

from rational_a0_less_noisy_certificate import (
    A_0,
    A0_CANDIDATE,
    CHARACTERS,
    DELTA_0,
    E_0,
    HALF,
    Q_0,
    P809_CANDIDATE,
    P8094_CANDIDATE,
    P809439_CANDIDATE,
    RationalTriangleCandidate,
    REDUCED_FULL,
    S_0,
    STATES,
    a0_less_noisy_certificate,
    binary_face_margin,
    matrix_dominance_residuals,
    p809_less_noisy_certificate,
    p8094_less_noisy_certificate,
    p809439_less_noisy_certificate,
    polarized_gap_matrix,
    prior_profile,
    reduced_diagonal_coefficient,
    sturm_root_count,
)


def _centered(
    prior: tuple[Fraction, ...], function: tuple[Fraction, ...]
) -> tuple[Fraction, ...]:
    mean = sum(mass * value for mass, value in zip(prior, function, strict=True))
    return tuple(value - mean for value in function)


def _variance(prior: tuple[Fraction, ...], function: tuple[Fraction, ...]) -> Fraction:
    centered = _centered(prior, function)
    return sum(
        mass * value * value for mass, value in zip(prior, centered, strict=True)
    )


def _projection_variance_sum(
    prior: tuple[Fraction, ...], function: tuple[Fraction, ...]
) -> Fraction:
    centered = _centered(prior, function)
    total = Fraction(0)
    for character in CHARACTERS:
        positive_mass = sum(
            prior[index] for index, sign in enumerate(character) if sign == 1
        )
        negative_mass = 1 - positive_mass
        if not positive_mass or not negative_mass:
            continue
        numerator = sum(
            prior[index] * centered[index]
            for index, sign in enumerate(character)
            if sign == 1
        )
        total += numerator * numerator / (positive_mass * negative_mass)
    return total


def _direct_strengthened_gap(
    prior: tuple[Fraction, ...],
    function: tuple[Fraction, ...],
    candidate: RationalTriangleCandidate = A0_CANDIDATE,
) -> Fraction:
    centered = _centered(prior, function)
    variance = _variance(prior, function)
    physical = sum(
        mass * prior_profile(mass, candidate.q) * value * value
        for mass, value in zip(prior, centered, strict=True)
    )
    erasure = candidate.full * variance + candidate.single * (
        _projection_variance_sum(prior, function)
    )
    return erasure - physical - candidate.variance_gap * variance


def _matrix_quadratic(matrix, coordinates) -> Fraction:
    return sum(
        coordinates[row] * matrix[row][column] * coordinates[column]
        for row in range(len(coordinates))
        for column in range(len(coordinates))
    )


def _bernstein_coefficients(
    polynomial: tuple[int, ...], left: Fraction, right: Fraction
) -> tuple[Fraction, ...]:
    """Return exact Bernstein coefficients after mapping ``[left,right]`` to ``[0,1]``."""

    degree = len(polynomial) - 1
    width = right - left
    power = [Fraction(0)] * (degree + 1)
    for source_degree, coefficient in enumerate(polynomial):
        for target_degree in range(source_degree + 1):
            power[target_degree] += (
                coefficient
                * comb(source_degree, target_degree)
                * left ** (source_degree - target_degree)
                * width**target_degree
            )
    return tuple(
        sum(
            power[source_degree]
            * Fraction(comb(index, source_degree), comb(degree, source_degree))
            for source_degree in range(index + 1)
        )
        for index in range(degree + 1)
    )


def _positive_bernstein_cover(
    polynomial: tuple[int, ...],
    left: Fraction,
    right: Fraction,
    depth: int = 0,
    maximum_depth: int = 12,
) -> tuple[int, int]:
    coefficients = _bernstein_coefficients(polynomial, left, right)
    if all(coefficient > 0 for coefficient in coefficients):
        return 1, depth
    if depth >= maximum_depth:
        raise AssertionError("the exact Bernstein positivity cover did not close")
    middle = (left + right) / 2
    left_count, left_depth = _positive_bernstein_cover(
        polynomial, left, middle, depth + 1, maximum_depth
    )
    right_count, right_depth = _positive_bernstein_cover(
        polynomial, middle, right, depth + 1, maximum_depth
    )
    return left_count + right_count, max(left_depth, right_depth)


class RationalA0LessNoisyCertificateTests(unittest.TestCase):
    def test_constants_and_exact_profile_formula(self) -> None:
        self.assertEqual(Q_0, Fraction(61, 100))
        self.assertEqual(A_0 + 3 * S_0 + E_0, 1)
        self.assertEqual(REDUCED_FULL, Fraction(323, 1000))
        q = Q_0
        for mass in (
            Fraction(1, 100),
            Fraction(1, 4),
            Fraction(1, 2),
            Fraction(3, 4),
            Fraction(99, 100),
        ):
            direct = (
                2
                * q
                * q
                * mass
                * (
                    (1 + q) / ((1 - q) ** 2 + 4 * q * mass)
                    + (1 - q) / ((1 + q) ** 2 - 4 * q * mass)
                )
            )
            self.assertEqual(prior_profile(mass), direct)

    def test_sturm_engine_counts_known_roots(self) -> None:
        # (t-1/3)(t-2/3) = t^2-t+2/9.
        polynomial = (Fraction(2, 9), Fraction(-1), Fraction(1))
        self.assertEqual(sturm_root_count(polynomial, 0, 1), 2)
        self.assertEqual(sturm_root_count(polynomial, 0, Fraction(1, 2)), 1)
        with self.assertRaises(ValueError):
            sturm_root_count((Fraction(-1), Fraction(1)), 0, 1)

    def test_exact_diagonal_dominance_ceiling_above_p809439(self) -> None:
        q = Fraction(309439019, 500000000)
        tail_point = Fraction(6127079, 20000000)
        off_diagonal_point = Fraction(56557217, 100000000)
        tail_lower_bound = (
            prior_profile(tail_point, q) - 2 * tail_point * prior_profile(HALF, q)
        ) / (1 - 2 * tail_point)
        off_diagonal_upper_bound = (
            1 - 3 * (1 - off_diagonal_point) * prior_profile(off_diagonal_point, q)
        ) / (3 * off_diagonal_point - 1)
        self.assertEqual(
            tail_lower_bound,
            Fraction(
                8022626566837929926456876745349269849352750411481961087,
                24071400765421871584513133172002858211131157280382500000,
            ),
        )
        self.assertEqual(
            off_diagonal_upper_bound,
            Fraction(
                90054211950004859732870915240433502464745259388712786399,
                270202164286984350616107038750357334771130625000000000000,
            ),
        )
        self.assertEqual(
            tail_lower_bound - off_diagonal_upper_bound,
            Fraction(
                14065222857405417517792295205369626373869821016600100484402322184760321915204457808629397159792953,
                2601657833694544556860086082446992935025579384744178481259486514689111374965381788339862985625000000000000,
            ),
        )

        # Primitive numerator of the same bound difference as a function of z=q^2.
        numerator = (
            -1562500000000000000000000000000000000,
            -3723156998744256875000000000000000000,
            4474518809248178047275159638775000000,
            25105318581993286191938002251843285299,
            34353196848684139192702888251025000000,
            20387125116385832068338520637956714701,
            4504071041702544145298475000000000000,
        )
        squared_quality = q * q
        self.assertEqual(sturm_root_count(numerator, squared_quality, 1), 0)
        value = sum(
            Fraction(coefficient) * squared_quality**degree
            for degree, coefficient in enumerate(numerator)
        )
        self.assertGreater(value, 0)

    def test_exhaustive_certificate_has_no_unresolved_region(self) -> None:
        certificate = a0_less_noisy_certificate()
        self.assertTrue(certificate.exhaustive)
        self.assertEqual(certificate.status, "CERTIFIED_PSD")
        self.assertEqual(certificate.unresolved_regions, 0)
        self.assertIn("diagonal dominance", certificate.proof_method)
        self.assertTrue(certificate.nonpolarized.certified)
        self.assertTrue(certificate.tail_separation.certified)
        self.assertTrue(certificate.dominant_decrease.certified)
        self.assertTrue(certificate.off_diagonal_lower_bound.certified)
        self.assertTrue(
            all(
                item.root_count == 0
                for item in (
                    certificate.nonpolarized,
                    certificate.tail_separation,
                    certificate.dominant_decrease,
                    certificate.off_diagonal_lower_bound,
                )
            )
        )

    def test_p809_exhaustive_certificate_has_no_unresolved_region(self) -> None:
        certificate = p809_less_noisy_certificate()
        self.assertTrue(certificate.exhaustive)
        self.assertEqual(certificate.candidate_name, "p809")
        self.assertEqual(certificate.p, Fraction(809, 1000))
        self.assertEqual(certificate.q, Fraction(309, 500))
        self.assertEqual(certificate.variance_gap, Fraction(1, 50000))
        self.assertTrue(
            all(
                item.certified and item.root_count == 0
                for item in (
                    certificate.nonpolarized,
                    certificate.tail_separation,
                    certificate.dominant_decrease,
                    certificate.off_diagonal_lower_bound,
                )
            )
        )

    def test_p8094_exhaustive_certificate_has_no_unresolved_region(self) -> None:
        certificate = p8094_less_noisy_certificate()
        self.assertTrue(certificate.exhaustive)
        self.assertEqual(certificate.candidate_name, "p8094")
        self.assertEqual(certificate.p, Fraction(4047, 5000))
        self.assertEqual(certificate.q, Fraction(1547, 2500))
        self.assertEqual(certificate.variance_gap, Fraction(1, 1000000))
        self.assertTrue(
            all(
                item.certified and item.root_count == 0
                for item in (
                    certificate.nonpolarized,
                    certificate.tail_separation,
                    certificate.dominant_decrease,
                    certificate.off_diagonal_lower_bound,
                )
            )
        )

    def test_p809439_exhaustive_certificate_has_no_unresolved_region(self) -> None:
        certificate = p809439_less_noisy_certificate()
        self.assertTrue(certificate.exhaustive)
        self.assertEqual(certificate.candidate_name, "p809439")
        self.assertEqual(certificate.p, Fraction(809439, 1000000))
        self.assertEqual(certificate.q, Fraction(309439, 500000))
        self.assertEqual(certificate.variance_gap, Fraction(1, 50000000))
        self.assertTrue(
            all(
                item.certified and item.root_count == 0
                for item in (
                    certificate.nonpolarized,
                    certificate.tail_separation,
                    certificate.dominant_decrease,
                    certificate.off_diagonal_lower_bound,
                )
            )
        )

    def test_derived_sturm_polynomials_are_stable(self) -> None:
        certificate = a0_less_noisy_certificate()
        self.assertEqual(
            certificate.nonpolarized.polynomial,
            (
                12734546643,
                -105026035232,
                295953456000,
                -266721280000,
            ),
        )
        self.assertEqual(
            certificate.tail_separation.polynomial,
            (
                174730714488603,
                -1096598972008000,
                1792103257120000,
            ),
        )
        self.assertEqual(
            certificate.off_diagonal_lower_bound.polynomial,
            (
                17150240835,
                -68743038643,
                85240668000,
                -29232176000,
            ),
        )

    def test_p809_derived_sturm_polynomials_are_stable(self) -> None:
        certificate = p809_less_noisy_certificate()
        self.assertEqual(
            certificate.nonpolarized.polynomial,
            (
                396319738471239,
                -3370596266968040,
                9493561252800000,
                -8539820640000000,
            ),
        )
        self.assertEqual(
            certificate.tail_separation.polynomial,
            (
                45640313188927373653,
                -297658582360966800000,
                485484077838348000000,
            ),
        )
        self.assertEqual(
            certificate.dominant_decrease.polynomial,
            (
                3154192724499694324012093,
                100909346530592989224000000,
                -656024362891036720509600000,
                994139601778209600000000000,
                -121418809934342400000000000,
            ),
        )
        self.assertEqual(
            certificate.off_diagonal_lower_bound.polynomial,
            (
                529787256879229,
                -2165830919071239,
                2723251793400000,
                -953740612800000,
            ),
        )

    def test_p8094_derived_sturm_polynomials_are_stable(self) -> None:
        certificate = p8094_less_noisy_certificate()
        self.assertEqual(
            certificate.nonpolarized.polynomial,
            (
                1651944835903320363,
                -14070058150170968560,
                39445224932932000000,
                -35476930216000000000,
            ),
        )
        self.assertEqual(
            certificate.tail_separation.polynomial,
            (
                14278104473183101691364867,
                -93204333895741039244000000,
                152109867179442963700000000,
            ),
        )
        self.assertEqual(
            certificate.dominant_decrease.polynomial,
            (
                24572409865872707621348365738203,
                790689849757469884522973400000000,
                -5140035038301541643972134792000000,
                7789945054613015403200000000000000,
                -949695863121411478800000000000000,
            ),
        )
        self.assertEqual(
            certificate.off_diagonal_lower_bound.polynomial,
            (
                734402096122192741,
                -3019272080642440121,
                3791388707652000000,
                -1322764905644000000,
            ),
        )

    def test_p8094_has_an_independent_exact_bernstein_cover(self) -> None:
        certificate = p8094_less_noisy_certificate()
        expected_covers = ((3, 2), (9, 8), (1, 0), (4, 3))
        for item, expected in zip(
            (
                certificate.nonpolarized,
                certificate.tail_separation,
                certificate.dominant_decrease,
                certificate.off_diagonal_lower_bound,
            ),
            expected_covers,
            strict=True,
        ):
            self.assertEqual(
                _positive_bernstein_cover(
                    item.polynomial, item.interval[0], item.interval[1]
                ),
                expected,
            )

    def test_p809439_polynomials_and_independent_bernstein_cover(self) -> None:
        certificate = p809439_less_noisy_certificate()
        expected_polynomials = (
            (
                1982400847413554713024317119535,
                -16886556652293525468619903939898,
                47330451139550895979364000000000,
                -42569116315311277848000000000000,
            ),
            (
                137084007706052200483756974963331002694947,
                -894938462691318390938821881638944400000000,
                1460628246646072375716946304372000000000000,
            ),
            (
                9433172315345871934969498145675028982913299390874787,
                303711722682995718660328951892648021559576000000000000,
                -1974331516564512256114987698060786064097177024800000000,
                2992214636015683935852477362179332800000000000000000000,
                -364727750098333998369989835623555200000000000000000000,
            ),
            (
                5286982498277589323500182269121,
                -21742080055560441029094634239070,
                27299948009212567865774000000000,
                -9522669648479236262728000000000,
            ),
        )
        items = (
            certificate.nonpolarized,
            certificate.tail_separation,
            certificate.dominant_decrease,
            certificate.off_diagonal_lower_bound,
        )
        expected_covers = ((3, 2), (12, 11), (1, 0), (11, 10))
        for item, polynomial, cover in zip(
            items, expected_polynomials, expected_covers, strict=True
        ):
            self.assertEqual(item.polynomial, polynomial)
            self.assertEqual(
                _positive_bernstein_cover(
                    item.polynomial,
                    item.interval[0],
                    item.interval[1],
                    maximum_depth=14,
                ),
                cover,
            )

    def test_binary_face_has_an_exact_strengthened_margin(self) -> None:
        expected = Fraction(547, 1000) - Fraction(7442, 13721)
        self.assertEqual(binary_face_margin(), expected)
        self.assertEqual(expected, Fraction(63387, 13721000))
        self.assertGreater(expected, 0)

    def test_p809_binary_face_has_an_exact_strengthened_margin(self) -> None:
        expected = Fraction(49016699, 17274050000)
        self.assertEqual(binary_face_margin(P809_CANDIDATE), expected)
        self.assertGreater(expected, 0)

    def test_p8094_binary_face_has_an_exact_strengthened_margin(self) -> None:
        expected = Fraction(15135252561, 8643209000000)
        self.assertEqual(binary_face_margin(P8094_CANDIDATE), expected)
        self.assertGreater(expected, 0)

    def test_p809439_binary_face_has_an_exact_strengthened_margin(self) -> None:
        expected = Fraction(
            73167609074569033,
            43219061840125000000,
        )
        self.assertEqual(binary_face_margin(P809439_CANDIDATE), expected)
        self.assertGreater(expected, 0)

    def test_matrix_has_exact_diagonal_dominance_decomposition(self) -> None:
        priors = (
            (
                Fraction(3, 5),
                Fraction(1, 5),
                Fraction(1, 10),
                Fraction(1, 10),
            ),
            (
                Fraction(51, 100),
                Fraction(47, 100),
                Fraction(1, 100),
                Fraction(1, 100),
            ),
            (
                Fraction(9, 10),
                Fraction(1, 20),
                Fraction(3, 100),
                Fraction(1, 50),
            ),
        )
        for prior in priors:
            matrix = polarized_gap_matrix(prior)
            off_diagonal = (
                matrix[0][1],
                matrix[0][2],
                matrix[1][2],
            )
            self.assertTrue(all(value > 0 for value in off_diagonal))
            residuals = matrix_dominance_residuals(prior, matrix)
            expected = tuple(
                reduced_diagonal_coefficient(prior[index])
                - reduced_diagonal_coefficient(prior[0])
                for index in range(1, 4)
            )
            self.assertEqual(residuals, expected)
            self.assertTrue(all(value >= 0 for value in residuals))

    def test_matrix_matches_both_quadratic_forms_exactly(self) -> None:
        prior = (
            Fraction(13, 20),
            Fraction(1, 5),
            Fraction(1, 10),
            Fraction(1, 20),
        )
        functions = (
            (Fraction(0), Fraction(1), Fraction(-2), Fraction(3)),
            (Fraction(5), Fraction(-1), Fraction(4), Fraction(2)),
            tuple(Fraction(first * second) for first, second in STATES),
        )
        matrix = polarized_gap_matrix(prior)
        for function in functions:
            centered = _centered(prior, function)
            coordinates = tuple(prior[index] * centered[index] for index in range(1, 4))
            self.assertEqual(
                _direct_strengthened_gap(prior, function),
                _matrix_quadratic(matrix, coordinates),
            )

    def test_p809_matrix_matches_forms_and_dominance_exactly(self) -> None:
        priors = (
            (
                Fraction(13, 20),
                Fraction(1, 5),
                Fraction(1, 10),
                Fraction(1, 20),
            ),
            (
                Fraction(51, 100),
                Fraction(47, 100),
                Fraction(1, 100),
                Fraction(1, 100),
            ),
            (
                Fraction(9, 10),
                Fraction(1, 20),
                Fraction(3, 100),
                Fraction(1, 50),
            ),
        )
        functions = (
            (Fraction(0), Fraction(1), Fraction(-2), Fraction(3)),
            (Fraction(5), Fraction(-1), Fraction(4), Fraction(2)),
        )
        for prior in priors:
            matrix = polarized_gap_matrix(prior, P809_CANDIDATE)
            self.assertTrue(
                all(value > 0 for value in (matrix[0][1], matrix[0][2], matrix[1][2]))
            )
            residuals = matrix_dominance_residuals(prior, matrix, P809_CANDIDATE)
            expected = tuple(
                reduced_diagonal_coefficient(prior[index], P809_CANDIDATE)
                - reduced_diagonal_coefficient(prior[0], P809_CANDIDATE)
                for index in range(1, 4)
            )
            self.assertEqual(residuals, expected)
            self.assertTrue(all(value >= 0 for value in residuals))
            for function in functions:
                centered = _centered(prior, function)
                coordinates = tuple(
                    prior[index] * centered[index] for index in range(1, 4)
                )
                self.assertEqual(
                    _direct_strengthened_gap(prior, function, P809_CANDIDATE),
                    _matrix_quadratic(matrix, coordinates),
                )

    def test_p8094_matrix_matches_forms_and_dominance_exactly(self) -> None:
        prior = (
            Fraction(13, 20),
            Fraction(1, 5),
            Fraction(1, 10),
            Fraction(1, 20),
        )
        functions = (
            (Fraction(0), Fraction(1), Fraction(-2), Fraction(3)),
            (Fraction(5), Fraction(-1), Fraction(4), Fraction(2)),
        )
        matrix = polarized_gap_matrix(prior, P8094_CANDIDATE)
        self.assertTrue(
            all(value > 0 for value in (matrix[0][1], matrix[0][2], matrix[1][2]))
        )
        residuals = matrix_dominance_residuals(prior, matrix, P8094_CANDIDATE)
        expected = tuple(
            reduced_diagonal_coefficient(prior[index], P8094_CANDIDATE)
            - reduced_diagonal_coefficient(prior[0], P8094_CANDIDATE)
            for index in range(1, 4)
        )
        self.assertEqual(residuals, expected)
        self.assertTrue(all(value >= 0 for value in residuals))
        for function in functions:
            centered = _centered(prior, function)
            coordinates = tuple(prior[index] * centered[index] for index in range(1, 4))
            self.assertEqual(
                _direct_strengthened_gap(prior, function, P8094_CANDIDATE),
                _matrix_quadratic(matrix, coordinates),
            )

    def test_p809439_matrix_matches_forms_and_dominance_exactly(self) -> None:
        priors = (
            (
                Fraction(13, 20),
                Fraction(1, 5),
                Fraction(1, 10),
                Fraction(1, 20),
            ),
            (
                Fraction(51, 100),
                Fraction(47, 100),
                Fraction(1, 100),
                Fraction(1, 100),
            ),
        )
        functions = (
            (Fraction(0), Fraction(1), Fraction(-2), Fraction(3)),
            (Fraction(5), Fraction(-1), Fraction(4), Fraction(2)),
        )
        for prior in priors:
            matrix = polarized_gap_matrix(prior, P809439_CANDIDATE)
            self.assertTrue(
                all(value > 0 for value in (matrix[0][1], matrix[0][2], matrix[1][2]))
            )
            residuals = matrix_dominance_residuals(prior, matrix, P809439_CANDIDATE)
            expected = tuple(
                reduced_diagonal_coefficient(prior[index], P809439_CANDIDATE)
                - reduced_diagonal_coefficient(prior[0], P809439_CANDIDATE)
                for index in range(1, 4)
            )
            self.assertEqual(residuals, expected)
            self.assertTrue(all(value >= 0 for value in residuals))
            for function in functions:
                centered = _centered(prior, function)
                coordinates = tuple(
                    prior[index] * centered[index] for index in range(1, 4)
                )
                self.assertEqual(
                    _direct_strengthened_gap(prior, function, P809439_CANDIDATE),
                    _matrix_quadratic(matrix, coordinates),
                )

    def test_boundary_faces_are_consistent_with_the_certificate(self) -> None:
        examples = (
            (
                (HALF, HALF, Fraction(0), Fraction(0)),
                (Fraction(1), Fraction(-1), Fraction(0), Fraction(0)),
            ),
            (
                (
                    Fraction(3, 5),
                    Fraction(1, 4),
                    Fraction(3, 20),
                    Fraction(0),
                ),
                (Fraction(0), Fraction(2), Fraction(-3), Fraction(7)),
            ),
            (
                (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
                (Fraction(2), Fraction(-5), Fraction(3), Fraction(9)),
            ),
        )
        for candidate in (
            A0_CANDIDATE,
            P809_CANDIDATE,
            P8094_CANDIDATE,
            P809439_CANDIDATE,
        ):
            for prior, function in examples:
                self.assertGreaterEqual(
                    _direct_strengthened_gap(prior, function, candidate), 0
                )

    def test_chayes_lei_conditions_have_exact_slack(self) -> None:
        self.assertLess(A_0, E_0)
        self.assertEqual(2 * A_0 + 3 * S_0, Fraction(124, 125))
        self.assertGreater(A_0 * E_0, 2 * S_0 * S_0)
        self.assertEqual(A_0 * E_0, Fraction(1722, 15625))
        self.assertEqual(2 * S_0 * S_0, Fraction(392, 15625))

    def test_p809_chayes_lei_conditions_have_exact_slack(self) -> None:
        candidate = P809_CANDIDATE
        certificate = p809_less_noisy_certificate()
        self.assertEqual(candidate.full + 3 * candidate.single + candidate.empty, 1)
        self.assertEqual(certificate.order_slack, Fraction(3, 5000))
        self.assertEqual(certificate.self_dual_slack, Fraction(3, 5000))
        self.assertEqual(certificate.fkg_slack, Fraction(1067809, 12500000))
        self.assertEqual(
            certificate.density_squared_slack,
            Fraction(76882329, 25000000),
        )
        self.assertEqual(
            candidate.full * candidate.empty,
            Fraction(138029, 1250000),
        )
        self.assertEqual(
            2 * candidate.single * candidate.single,
            Fraction(312481, 12500000),
        )

    def test_p8094_chayes_lei_conditions_have_exact_slack(self) -> None:
        candidate = P8094_CANDIDATE
        certificate = p8094_less_noisy_certificate()
        self.assertEqual(candidate.full + 3 * candidate.single + candidate.empty, 1)
        self.assertEqual(certificate.order_slack, Fraction(3, 25000))
        self.assertEqual(certificate.self_dual_slack, Fraction(3, 25000))
        self.assertEqual(
            certificate.fkg_slack,
            Fraction(172640489, 2000000000),
        )
        self.assertEqual(
            certificate.density_squared_slack,
            Fraction(7768822329, 2500000000),
        )

    def test_p809439_chayes_lei_conditions_have_exact_slack(self) -> None:
        candidate = P809439_CANDIDATE
        certificate = p809439_less_noisy_certificate()
        self.assertEqual(candidate.full + 3 * candidate.single + candidate.empty, 1)
        self.assertEqual(certificate.order_slack, Fraction(7, 500000000))
        self.assertEqual(certificate.self_dual_slack, Fraction(7, 500000000))
        self.assertEqual(
            certificate.fkg_slack,
            Fraction(10796599147227459, 125000000000000000),
        )
        self.assertEqual(
            certificate.density_squared_slack,
            Fraction(777355138600377489, 250000000000000000),
        )


if __name__ == "__main__":
    unittest.main()
