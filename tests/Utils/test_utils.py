import numpy as np

from caat.utils import (
    MIN_PHOTOMETRY_ERROR,
    ab_zero_point,
    filter_zero_point_function,
    floor_photometry_errors,
)

# Zero points that are deliberately not exactly AB, so the nodes matter
ZPS = {"U": 8.0e-9, "B": 6.0e-9, "V": 3.0e-9, "r": 2.5e-9}
FILTERS = np.asarray(["U", "U", "B", "B", "V", "V"])
WAVELENGTHS = np.asarray([3290.0, 3310.0, 3990.0, 4010.0, 5600.0, 5640.0])


class TestAbZeroPoint:
    def test_matches_the_ab_formula(self):
        assert np.isclose(ab_zero_point(5000.0), 3631e-23 * 3e18 / 5000.0**2)

    def test_is_an_array_for_array_input(self):
        assert ab_zero_point(np.asarray([3000.0, 6000.0])).shape == (2,)

    def test_falls_by_a_factor_of_four_when_wavelength_doubles(self):
        assert np.isclose(ab_zero_point(4000.0) / ab_zero_point(8000.0), 4.0)


class TestFilterZeroPointFunction:
    def test_reproduces_the_filter_zero_points_at_their_median_wavelength(self):
        zero_point = filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)
        for filt, wl in (("U", 3300.0), ("B", 4000.0), ("V", 5620.0)):
            # median in log space of the two points is within a fraction of an Angstrom of wl
            assert np.isclose(zero_point(wl), ZPS[filt], rtol=1e-2)

    def test_is_exact_at_the_median_wavelength_in_log_space(self):
        zero_point = filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)
        median_wl = 10 ** np.median(np.log10(WAVELENGTHS[FILTERS == "B"]))
        assert np.isclose(zero_point(median_wl), ZPS["B"])

    def test_is_log_linear_between_filters(self):
        zero_point = filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)
        wl_u, wl_b = (10 ** np.median(np.log10(WAVELENGTHS[FILTERS == f])) for f in ("U", "B"))
        mid = 10 ** np.mean(np.log10([wl_u, wl_b]))
        assert np.isclose(np.log10(zero_point(mid)), np.mean(np.log10([ZPS["U"], ZPS["B"]])))

    def test_continues_with_the_ab_slope_beyond_the_ends(self):
        zero_point = filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)
        # Beyond the blue end, a factor of 2 in wavelength changes the zero point by a factor of 4
        assert np.isclose(zero_point(1500.0) / zero_point(3000.0), 4.0)
        # and the same beyond the red end
        assert np.isclose(zero_point(9000.0) / zero_point(18000.0), 4.0)

    def test_is_continuous_at_the_end_nodes(self):
        zero_point = filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)
        wl_u = 10 ** np.median(np.log10(WAVELENGTHS[FILTERS == "U"]))
        assert np.isclose(zero_point(wl_u * (1 - 1e-9)), zero_point(wl_u * (1 + 1e-9)), rtol=1e-6)

    def test_does_not_depend_on_the_order_of_the_input(self):
        order = np.random.default_rng(0).permutation(len(FILTERS))
        wl = np.linspace(2000.0, 8000.0, 25)
        a = filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)(wl)
        b = filter_zero_point_function(ZPS, FILTERS[order], WAVELENGTHS[order])(wl)
        assert np.allclose(a, b)

    def test_ignores_filters_without_a_zero_point(self):
        filters = np.concatenate((FILTERS, ["zz", "zz"]))
        wavelengths = np.concatenate((WAVELENGTHS, [9000.0, 9100.0]))
        wl = np.linspace(2000.0, 8000.0, 25)
        assert np.allclose(
            filter_zero_point_function(ZPS, filters, wavelengths)(wl),
            filter_zero_point_function(ZPS, FILTERS, WAVELENGTHS)(wl),
        )

    def test_falls_back_to_the_ab_zero_point_with_fewer_than_two_filters(self):
        wl = np.linspace(2000.0, 8000.0, 25)
        one_filter = filter_zero_point_function(ZPS, np.asarray(["B", "B"]), np.asarray([3990.0, 4010.0]))
        assert np.allclose(one_filter(wl), ab_zero_point(wl))
        no_known_filter = filter_zero_point_function(ZPS, np.asarray(["zz", "yy"]), np.asarray([3990.0, 4010.0]))
        assert np.allclose(no_known_filter(wl), ab_zero_point(wl))


class TestFloorPhotometryErrors:
    def test_raises_small_and_zero_errors_to_the_floor(self):
        floored = floor_photometry_errors(np.asarray([0.0, 0.005, MIN_PHOTOMETRY_ERROR / 2]))
        assert np.allclose(floored, MIN_PHOTOMETRY_ERROR)

    def test_keeps_errors_above_the_floor(self):
        errors = np.asarray([0.05, 0.1, 0.3])
        assert np.allclose(floor_photometry_errors(errors), errors)

    def test_sets_missing_errors_to_the_floor(self):
        assert np.allclose(floor_photometry_errors(np.asarray([np.nan, 0.1])), [MIN_PHOTOMETRY_ERROR, 0.1])

    def test_does_not_modify_the_input(self):
        errors = np.asarray([0.0, 0.1])
        floor_photometry_errors(errors)
        assert np.allclose(errors, [0.0, 0.1])

    def test_uses_a_given_floor(self):
        assert np.allclose(floor_photometry_errors(np.asarray([0.0, 0.2]), floor=0.1), [0.1, 0.2])

    def test_floor_is_small_compared_with_typical_uncertainties(self):
        assert 0.0 < MIN_PHOTOMETRY_ERROR <= 0.05

