from types import SimpleNamespace

import matplotlib
import numpy as np
import pandas as pd
import pytest
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

from caat import SNModel
from caat.SNModel import SurfaceArray
from caat.utils import (
    WLE,
    ab_zero_point,
    convert_shifted_fluxes_to_shifted_mags,
    MIN_PHOTOMETRY_ERROR,
    filter_zero_point_function,
)

matplotlib.use("Agg")

LOG_TRANSFORM = 22
PHASES = np.linspace(-10.0, 50.0, 61)
WAVELENGTHS = np.linspace(3000.0, 9000.0, 31)


def _surface_grid():
    """A smooth, fully known surface in magnitudes from peak, shape (n_wl, n_phase)"""
    phase_term = -0.03 * (PHASES[None, :] + 10.0)
    wl_term = -0.0004 * (WAVELENGTHS[:, None] - 3000.0)
    return phase_term + wl_term


def _template_grid(scale=1.0):
    """A template in log flux, shape (n_phase, n_wl), deliberately far from zero"""
    return scale * (-0.5 - 0.02 * np.arange(len(PHASES))[:, None] + 0.03 * np.arange(len(WAVELENGTHS))[None, :])


def _mock_sn():
    # Like a real SN, with a zero point for every filter (AB f_lambda at the nominal wavelength)
    return SimpleNamespace(
        name="mock_sn",
        zps={filt: 3631e-23 * 3e18 / wl**2 for filt, wl in WLE.items()},
        info={"peak_filt": "V", "peak_mag": 15.0},
    )


def _surface_array_model(template_scale=1.0, iqr=0.0):
    surface = SurfaceArray(
        surface=np.asarray([_surface_grid(), np.full((len(WAVELENGTHS), len(PHASES)), iqr)]),
        phase_grid=np.log(PHASES + LOG_TRANSFORM),
        wl_grid=np.log10(WAVELENGTHS),
        kernel=RBF(),
    )
    return SNModel(
        surface=surface,
        template_mags=_template_grid(template_scale),
        phase_grid=PHASES,
        wl_grid=WAVELENGTHS,
        filters_fit=["B", "V"],
        sncollection=SimpleNamespace(type="Other", subtype="TDE", sne=[]),
        log_transform=LOG_TRANSFORM,
    )


FIT_PHASE_MIN = -5.0


def _residual_gp_model():
    """A per-SN style model: a GP on residuals in log flux, fit in log(phase - first phase + 0.1), with an SN to convert to magnitudes"""
    rng = np.random.default_rng(0)
    phases = np.append(rng.uniform(FIT_PHASE_MIN, 45, 19), FIT_PHASE_MIN)
    x = np.column_stack(
        (np.log(phases - FIT_PHASE_MIN + 0.1), np.log10(rng.uniform(3500, 8500, 20)))
    )
    y = 0.05 * rng.normal(size=20)
    gp = GaussianProcessRegressor(kernel=RBF([1.0, 0.3]), alpha=0.01, optimizer=None).fit(x, y)
    model = SNModel(
        surface=gp,
        template_mags=_template_grid(),
        phase_grid=PHASES,
        wl_grid=WAVELENGTHS,
        filters_fit=["B", "V"],
        sn=_mock_sn(),
        log_transform=LOG_TRANSFORM,
    )
    model._fit_phase_min = FIT_PHASE_MIN
    return model


def _expected_surface(phases, wavelengths):
    model = _surface_array_model()
    X = np.column_stack((np.log(np.asarray(phases) + LOG_TRANSFORM), np.log10(wavelengths)))
    return model.surface.predict(X)[0]


class TestSurfaceArrayModelReturnsMagnitudes:
    """The stored SurfaceArray surface already contains the template, in magnitudes from peak"""

    def test_predict_lightcurve_returns_the_surface(self):
        phases, mags, _ = _surface_array_model().predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        assert np.allclose(mags, _expected_surface(phases, np.full(len(phases), 5000.0)))

    def test_predict_lightcurve_does_not_depend_on_the_template(self):
        _, mags, _ = _surface_array_model(template_scale=1.0).predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        _, mags_other, _ = _surface_array_model(template_scale=10.0).predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        assert np.allclose(mags, mags_other)

    def test_predict_sed_returns_the_surface(self):
        waves, mags, _ = _surface_array_model().predict_sed(3500.0, 8500.0, 10.0, show=False)
        assert np.allclose(mags, _expected_surface(np.full(len(waves), 10.0), waves))

    def test_predict_sed_does_not_depend_on_the_template(self):
        _, mags, _ = _surface_array_model(template_scale=1.0).predict_sed(3500.0, 8500.0, 10.0, show=False)
        _, mags_other, _ = _surface_array_model(template_scale=10.0).predict_sed(3500.0, 8500.0, 10.0, show=False)
        assert np.allclose(mags, mags_other)

    def test_predict_photometry_points_returns_the_surface(self):
        wavelengths = np.asarray([4500.0, 4500.0, 7000.0, 7000.0])
        phases = np.asarray([0.0, 10.0, 0.0, 10.0])
        _, mags, _ = _surface_array_model().predict_photometry_points(wavelengths, phases)
        assert np.allclose(mags, _expected_surface(phases, wavelengths))

    def test_methods_agree_with_each_other(self):
        model = _surface_array_model()
        _, lc_mags, _ = model.predict_lightcurve(10.0, 10.0 + 1e-9, 5000.0, show=False)
        _, sed_mags, _ = model.predict_sed(5000.0, 5000.0 + 1e-9, 10.0, show=False)
        _, point_mags, _ = model.predict_photometry_points(np.asarray([5000.0]), np.asarray([10.0]))
        assert np.allclose(lc_mags[0], sed_mags[0], atol=1e-3)
        assert np.allclose(lc_mags[0], point_mags[0], atol=1e-3)

    def test_uncertainties_are_the_surface_spread_in_magnitudes(self):
        _, _, dev = _surface_array_model(iqr=0.2).predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        assert np.allclose(dev, 0.2)


class TestResidualGPModelReturnsMagnitudes:
    """A GP on residuals needs the template added back and a conversion from log flux to magnitudes"""

    @staticmethod
    def _expected(model, phases, wavelengths):
        X = np.column_stack(
            (np.log(np.maximum(np.asarray(phases) - model._fit_phase_min + 0.1, 0.1)), np.log10(wavelengths))
        )
        residual, std = model.surface.predict(X, return_std=True)
        template = np.asarray(
            [
                model.template[np.argmin(abs(PHASES - p)), np.argmin(abs(WAVELENGTHS - w))]
                for p, w in zip(phases, wavelengths)
            ]
        )
        zps = ab_zero_point(np.asarray(wavelengths))
        mags = 2.5 * (residual + template + np.log10(model.sn.zps["V"] / zps))
        return mags, 2.5 * std

    def test_conversion_formula(self):
        """Magnitudes from peak are 2.5 * (log flux + log10(zp of peak filter / zp))"""
        sn = _mock_sn()
        f = np.asarray([-0.3, 0.0, 0.4])
        zp = sn.zps["B"]
        assert np.allclose(
            convert_shifted_fluxes_to_shifted_mags(f, sn, zp),
            2.5 * (f + np.log10(sn.zps["V"] / zp)),
        )

    def test_predict_lightcurve_returns_magnitudes(self):
        model = _residual_gp_model()
        phases, mags, dev = model.predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        expected_mags, expected_dev = self._expected(model, phases, np.full(len(phases), 5000.0))
        assert np.allclose(mags, expected_mags)
        assert np.allclose(dev, expected_dev)

    def test_predict_sed_returns_magnitudes(self):
        model = _residual_gp_model()
        waves, mags, dev = model.predict_sed(3500.0, 8500.0, 10.0, show=False)
        expected_mags, expected_dev = self._expected(model, np.full(len(waves), 10.0), waves)
        assert np.allclose(mags, expected_mags)
        assert np.allclose(dev, expected_dev)

    def test_predict_photometry_points_returns_magnitudes(self):
        model = _residual_gp_model()
        wavelengths = np.asarray([4500.0, 4500.0, 7000.0, 7000.0])
        phases = np.asarray([0.0, 10.0, 0.0, 10.0])
        np.random.seed(0)
        draws = np.asarray([model.predict_photometry_points(wavelengths, phases)[1] for _ in range(300)])
        expected_mags, expected_dev = self._expected(model, phases, wavelengths)
        assert np.allclose(draws.mean(axis=0), expected_mags, atol=4 * expected_dev.max() / np.sqrt(300) + 1e-3)

    def test_raises_without_an_sn_to_convert_with(self):
        model = _residual_gp_model()
        model.sn = None
        for call in (
            lambda: model.predict_lightcurve(0.0, 40.0, 5000.0, show=False),
            lambda: model.predict_sed(3500.0, 8500.0, 10.0, show=False),
            lambda: model.predict_photometry_points(np.asarray([5000.0]), np.asarray([10.0])),
        ):
            with pytest.raises(ValueError, match="convert"):
                call()


    def test_raises_without_a_phase_reference(self):
        model = _residual_gp_model()
        model._fit_phase_min = None
        for call in (
            lambda: model.predict_lightcurve(0.0, 40.0, 5000.0, show=False),
            lambda: model.predict_sed(3500.0, 8500.0, 10.0, show=False),
            lambda: model.predict_photometry_points(np.asarray([5000.0]), np.asarray([10.0])),
        ):
            with pytest.raises(ValueError, match="phase"):
                call()


class TestGPPhaseCoordinate:
    """A GP on residuals is fit in log(phase - first phase + 0.1); a SurfaceArray is on a log(phase + log_transform) grid"""

    def test_values(self):
        model = _residual_gp_model()
        phases = np.asarray([-5.0, 0.0, 20.0])
        assert np.allclose(model._gp_phase_coordinate(phases), np.log(phases - FIT_PHASE_MIN + 0.1))

    def test_phases_at_or_before_the_first_epoch_map_to_the_first_epoch(self):
        model = _residual_gp_model()
        coordinate = model._gp_phase_coordinate(np.asarray([FIT_PHASE_MIN, FIT_PHASE_MIN - 3.0, FIT_PHASE_MIN - 30.0]))
        assert np.all(np.isfinite(coordinate))
        assert np.allclose(coordinate, np.log(0.1))

    def test_raises_without_a_phase_reference(self):
        model = _residual_gp_model()
        model._fit_phase_min = None
        with pytest.raises(ValueError, match="phase"):
            model._gp_phase_coordinate(np.asarray([0.0]))

    def test_surface_array_models_do_not_depend_on_it(self):
        model = _surface_array_model()
        _, mags, _ = model.predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        model._fit_phase_min = 123.0
        _, mags_other, _ = model.predict_lightcurve(0.0, 40.0, 5000.0, show=False)
        assert np.allclose(mags, mags_other)


class TestFitPhotometryKeepsTheFitSN:
    PHOTOMETRY = {
        "Filter": ["B", "B", "B", "V", "V", "V"],
        "Phase": [-5.0, 5.0, 20.0, -5.0, 5.0, 20.0],
        "Mag": [-1.0, 0.0, -0.8, -0.6, 0.0, -0.5],
        "MagErr": [0.1] * 6,
    }

    def test_keep_new_fit_stores_the_sn_used_for_the_fit(self):
        model = _surface_array_model()
        model.fit_photometry(photometry=dict(self.PHOTOMETRY), keep_new_fit=True, show=False)
        assert isinstance(model.surface, GaussianProcessRegressor)
        assert model._fit_sn is not None

    def test_predictions_after_keep_new_fit_are_in_magnitudes(self):
        model = _surface_array_model()
        model.fit_photometry(photometry=dict(self.PHOTOMETRY), keep_new_fit=True, show=False)
        phases, mags, dev = model.predict_lightcurve(-5.0, 20.0, WLE["V"], show=False)

        X = np.column_stack(
            (np.log(phases - model._fit_phase_min + 0.1), np.full(len(phases), np.log10(WLE["V"])))
        )
        residual, std = model.surface.predict(X, return_std=True)
        template = np.asarray(
            [model.template[np.argmin(abs(PHASES - p)), np.argmin(abs(WAVELENGTHS - WLE["V"]))] for p in phases]
        )
        sn = model._fit_sn
        zp = sn.zps["V"]  # WLE["V"] is the V filter, so that is the closest filter
        expected = 2.5 * (residual + template + np.log10(sn.zps[sn.info["peak_filt"]] / zp))
        assert np.allclose(mags, expected)
        assert np.allclose(dev, 2.5 * std)

    def test_keep_new_fit_stores_the_first_phase_and_fits_in_the_shifted_coordinate(self):
        model = _surface_array_model()
        model.fit_photometry(photometry=dict(self.PHOTOMETRY), keep_new_fit=True, show=False)
        first = min(self.PHOTOMETRY["Phase"])
        assert model._fit_phase_min == first
        assert np.allclose(
            np.sort(model.surface.X_train_[:, 0]),
            np.sort(np.log(np.asarray(self.PHOTOMETRY["Phase"]) - first + 0.1)),
        )

    def test_noise_of_the_fit_is_the_floored_variance(self):
        photometry = dict(self.PHOTOMETRY)
        photometry["MagErr"] = [0.0, 0.005, 0.1, 0.2, 0.0, 0.1]
        model = _surface_array_model()
        model.fit_photometry(photometry=photometry, keep_new_fit=True, show=False)
        # The residuals are in the order of the filters in `filters_to_fit`, so compare as sets
        expected = np.maximum(np.asarray(photometry["MagErr"]), MIN_PHOTOMETRY_ERROR) ** 2
        assert np.allclose(np.sort(model.surface.alpha), np.sort(expected))

    def test_duplicate_points_without_a_reported_error_do_not_make_the_fit_singular(self):
        photometry = {
            "Filter": ["B", "B", "B", "V", "V", "V"],
            "Phase": [0.0, 0.0, 10.0, 0.0, 0.0, 10.0],
            "Mag": [-0.5, -0.5, -0.2, -0.4, -0.4, -0.1],
            "MagErr": [0.0] * 6,
        }
        model = _surface_array_model()
        model.fit_photometry(photometry=photometry, keep_new_fit=True, show=False)
        _, mags, dev = model.predict_lightcurve(-5.0, 20.0, WLE["V"], show=False)
        assert np.all(np.isfinite(mags)) and np.all(np.isfinite(dev))

    def test_no_fit_sn_is_stored_without_keep_new_fit(self):
        model = _surface_array_model()
        model.fit_photometry(photometry=dict(self.PHOTOMETRY), keep_new_fit=False, show=False)
        assert isinstance(model.surface, SurfaceArray)
        assert getattr(model, "_fit_sn", None) is None


class TestZeroPointOfTheConversion:
    """The conversion to magnitudes uses a zero point interpolated between the fitted filters"""

    @staticmethod
    def _sn_with_non_ab_zero_points():
        sn = _mock_sn()
        sn.zps = {**sn.zps, "B": 6.5e-9, "V": 3.0e-9}  # not AB, so the filters' own zero points matter
        return sn

    def test_conversion_reproduces_the_filter_zero_point_at_the_wavelength_of_its_photometry(self):
        model = _residual_gp_model()
        model.sn = self._sn_with_non_ab_zero_points()
        model._zero_point = filter_zero_point_function(
            model.sn.zps, np.asarray(["B", "B", "V", "V"]), np.asarray([4000.0, 4000.0, 5600.0, 5600.0])
        )
        f = np.asarray([-0.3])
        for filt, wl in (("B", 4000.0), ("V", 5600.0)):
            expected = 2.5 * (f + np.log10(model.sn.zps["V"] / model.sn.zps[filt]))
            assert np.allclose(model._log_flux_to_mags(f, np.asarray([wl])), expected)

    def test_conversion_does_not_step_as_the_wavelength_crosses_between_filters(self):
        model = _residual_gp_model()
        model.sn = self._sn_with_non_ab_zero_points()
        model._zero_point = filter_zero_point_function(
            model.sn.zps, np.asarray(["B", "B", "V", "V"]), np.asarray([4000.0, 4000.0, 5600.0, 5600.0])
        )
        wavelengths = np.arange(2000.0, 9000.0, 10.0)
        mags = model._log_flux_to_mags(np.zeros(len(wavelengths)), wavelengths)
        assert np.abs(np.diff(mags)).max() < 0.05

    def test_conversion_uses_the_ab_zero_point_without_a_zero_point_function(self):
        model = _residual_gp_model()
        assert model._zero_point is None
        f, wl = np.asarray([-0.3, 0.1]), np.asarray([4000.0, 7000.0])
        assert np.allclose(
            model._log_flux_to_mags(f, wl), 2.5 * (f + np.log10(model.sn.zps["V"] / ab_zero_point(wl)))
        )

    def test_the_sn_of_the_last_fit_is_used_instead_of_the_sn_of_the_model(self):
        model = _residual_gp_model()
        fit_sn = _mock_sn()
        fit_sn.info = {"peak_filt": "r", "peak_mag": 16.0}
        model._fit_sn = fit_sn
        f, wl = np.asarray([-0.3]), np.asarray([5000.0])
        assert np.allclose(
            model._log_flux_to_mags(f, wl), 2.5 * (f + np.log10(fit_sn.zps["r"] / ab_zero_point(wl)))
        )

    def test_initializing_the_surface_fit_builds_the_zero_point_from_the_photometry(self, tmp_path):
        sn = self._sn_with_non_ab_zero_points()
        sn.base_path, sn.classification, sn.subtype = str(tmp_path), "Other", "TDE"
        (tmp_path / "Other" / "TDE" / sn.name).mkdir(parents=True)
        rows = [
            dict(Phase=p, ShiftedFilter=f, ShiftedWavelength=wl, ShiftedFlux=-0.1 * p / 10, ShiftedFluxerr=0.1, Nondetection=False)
            for f, wl in (("B", 4000.0), ("V", 5600.0))
            for p in (0.0, 10.0, 20.0, 30.0)
        ]
        rows[0]["ShiftedFluxerr"] = 0.0  # a detection with no reported uncertainty
        # A non-detection of another filter must not create a node
        rows.append(dict(Phase=15.0, ShiftedFilter="r", ShiftedWavelength=7000.0, ShiftedFlux=-0.5, ShiftedFluxerr=0.1, Nondetection=True))
        pd.DataFrame(rows).to_csv(tmp_path / "Other" / "TDE" / sn.name / f"{sn.name}_datacube_mangled.csv")

        model = _residual_gp_model()
        model.sn = sn
        model._initialize_surface_fit()

        # The noise is the variance (0.1**2), and the point without an uncertainty has the floor
        assert np.allclose(np.sort(np.unique(model.surface.alpha)), np.sort(np.unique([0.1**2, MIN_PHOTOMETRY_ERROR**2])))
        assert model._fit_phase_min == 0.0
        assert np.isclose(model.surface.X_train_[:, 0].min(), np.log(0.1))
        assert np.isclose(model._zero_point(4000.0), sn.zps["B"])
        assert np.isclose(model._zero_point(5600.0), sn.zps["V"])
        # With only B and V as nodes, the red end continues with the AB slope from V (not from r)
        assert np.isclose(model._zero_point(7000.0) / model._zero_point(14000.0), 4.0)
