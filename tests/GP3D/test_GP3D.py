import inspect

import numpy as np
import pandas as pd
import pytest
from unittest.mock import Mock, patch

from sklearn.gaussian_process.kernels import RBF

from caat import GP3D, Kernel
from caat.utils import MIN_PHOTOMETRY_ERROR, ab_zero_point, filter_zero_point_function


class TestGP3D:
    
    @pytest.fixture(autouse=True)
    def setup(self, mock_gp3d):
        self.phasemin = -20
        self.phasemax = 50
        self.gp = mock_gp3d
        self.phase_grid = np.arange(0.5, 5.0, 0.1)
        self.wl_grid = np.log10(np.arange(2000.0, 8000.0, 100.0))
    
    def test_build_samples(self):
        """Test build_samples method"""
        ### Test all boolean combinations of log_transform and use_flux
        with patch("caat.GP.GP._process_dataset", Mock(
            return_value=(
                np.asarray([0.5, 3.5]), 
                np.asarray([-1.0, 3.0]),
                np.asarray([0.1, 0.1]), 
                np.asarray([3.5, 3.5])
            )
        )):

            phases, wls, mags, err_grid = (
                self.gp._build_samples('B')
            )

            assert all([isinstance(var, np.ndarray) for var in [phases, mags, wls, err_grid]])
            assert len(phases) == len(mags) == len(wls) == len(err_grid)

    def test_process_dataset(self):
        """Test process_dataset method"""
        template_df = self.gp._process_dataset()
        assert isinstance(template_df, pd.DataFrame)

    def test_median_grid(self, mock_datacube):
        """Test construct_median_grid method"""
        phase_grid, wl_grid, mag_grid, err_grid = (
            self.gp._construct_median_grid(
                self.phasemin,
                self.phasemax,
                ['B'],
                mock_datacube,
                log_transform=22,
                plot=False,
            )
        )

        assert mag_grid.shape == (len(phase_grid), len(wl_grid))
        assert err_grid.shape == (len(phase_grid), len(wl_grid))

    def test_polynomial_grid(self, mock_datacube):
        """Test construct_polynomial_grid method"""
        phase_grid, wl_grid, mag_grid, err_grid = (
            self.gp._construct_polynomial_grid(
                self.phasemin,
                self.phasemax,
                ['B'],
                mock_datacube,
                log_transform=22,
                plot=False,
            )
        )

        assert mag_grid.shape == (len(phase_grid), len(wl_grid))
        assert err_grid.shape == (len(phase_grid), len(wl_grid))

    def test_subtract_from_grid(self, mock_sn):
        """Test subtract_from_grid method"""
        residuals = self.gp._subtract_data_from_grid(
            mock_sn,
            ['B'],
            self.phase_grid,
            self.wl_grid,
            np.random.random((len(self.phase_grid), len(self.wl_grid))),
            np.ones((len(self.phase_grid), len(self.wl_grid))) * 0.01,
        )

        assert isinstance(residuals, pd.DataFrame)

    def test_run_gp_full_sample_without_specifying_subtract(self):
        """Should raise an Exception when running GP without specifying a subtract method"""
        with pytest.raises(Exception, match=r'Must toggle either .*'):
            self.gp.run_gp_on_full_sample(
                plot=False,
            )

    def test_run_gp_individually_without_specifying_subtract(self):
        """Should raise an Exception when running GP without specifying a subtract method"""
        with pytest.raises(Exception, match=r'Must toggle either .*'):
            self.gp.run_gp_individually(
                plot=False,
            )

    def test_build_test_wavelength_phase_grid_from_photometry(self, mock_datacube):
        """Test build_test_wavelength_phase_grid_from_photometry method"""
        (
            x, 
            y, 
            wl_inds_fitted, 
            phase_inds_fitted, 
            min_phase
        ) = self.gp._build_test_wavelength_phase_grid_from_photometry(
            mock_datacube["LogShiftedWavelength"].values,
            mock_datacube["LogPhase"].values,
            self.wl_grid,
            self.phase_grid,
        )

        assert len(x) == len(y)
        assert isinstance(wl_inds_fitted, list | np.ndarray)
        assert isinstance(phase_inds_fitted, list | np.ndarray)
        assert isinstance(min_phase, float | None)

    def _optimize_hyperparams_with_residuals(self, residuals):
        """Run optimize_hyperparams with the grid construction and subtraction mocked"""
        with patch(
            "caat.GP3D.GP3D._construct_polynomial_grid", Mock(
                return_value=(
                    self.phase_grid,
                    self.wl_grid,
                    np.random.random((len(self.phase_grid), len(self.wl_grid))),
                    np.ones((len(self.phase_grid), len(self.wl_grid))) * 0.01,
                )
            )
        ), patch(
            "caat.GP3D.GP3D._subtract_data_from_grid", Mock(return_value=residuals)
        ):
            return self.gp.optimize_hyperparams(subtract_polynomial=True)

    @staticmethod
    def _mock_residuals():
        return pd.DataFrame(
            [
                {
                    "Filter": "B",
                    "Phase": phase,
                    "Wavelength": wl,
                    "MagResidual": residual,
                    "MagErr": 0.1,
                    "Mag": -0.5,
                    "Nondetection": False,
                }
                for phase, wl, residual in [
                    (np.log(30.0), 3.70, 0.1),
                    (np.log(35.0), 3.72, 0.2),
                    (np.log(42.0), 3.70, -0.1),
                    (np.log(55.0), 3.72, 0.0),
                ]
            ]
        )

    def test_optimize_hyperparameters(self):
        """optimize_hyperparams should fix the kernel at the exponentiated median theta"""
        kernel_params = self._optimize_hyperparams_with_residuals(self._mock_residuals())

        assert isinstance(kernel_params, list)
        # kernel_params are the log-transformed thetas, not the hyperparameters themselves
        expected_hyperparams = np.exp(np.median(np.asarray(kernel_params), axis=0))
        assert np.allclose(self.gp.kernel.length_scale, expected_hyperparams)
        assert self.gp.kernel.length_scale_bounds == "fixed"
        assert len(self.gp.kernel.theta) == 0

    def test_optimize_hyperparameters_modifies_kernel_object_in_place(self):
        """A Kernel object should be modified in place, so other references see the change"""
        kernel = Kernel(RBF())
        self.gp.kernel = kernel
        kernel_params = self._optimize_hyperparams_with_residuals(self._mock_residuals())

        assert self.gp.kernel is kernel
        expected_hyperparams = np.exp(np.median(np.asarray(kernel_params), axis=0))
        assert np.allclose(kernel.kernel.length_scale, expected_hyperparams)
        assert len(kernel.kernel.theta) == 0

    def test_optimize_hyperparameters_raises_for_fixed_kernel(self):
        """optimize_hyperparams should raise if the kernel has no free hyperparameters"""
        self.gp.kernel = Kernel(RBF(1.0, "fixed"))
        with pytest.raises(ValueError, match=r"no free hyperparameters"):
            self._optimize_hyperparams_with_residuals(self._mock_residuals())

    def test_optimize_hyperparameters_raises_when_called_twice(self):
        """The first call fixes the hyperparameters, so a second call should raise"""
        self._optimize_hyperparams_with_residuals(self._mock_residuals())
        with pytest.raises(ValueError, match=r"no free hyperparameters"):
            self._optimize_hyperparams_with_residuals(self._mock_residuals())

    def test_optimize_hyperparameters_raises_without_fittable_data(self):
        """optimize_hyperparams should raise if no transient had data to fit"""
        with pytest.raises(ValueError, match=r"enough data"):
            self._optimize_hyperparams_with_residuals(pd.DataFrame())

    def test_build_gaussian_process_uses_the_noise_and_kernel(self):
        gp = self.gp._build_gaussian_process(np.asarray([0.01, 0.04]))
        assert np.allclose(gp.alpha, [0.01, 0.04])
        assert gp.n_restarts_optimizer == 10
        assert gp.kernel is self.gp.kernel

    def test_build_gaussian_process_uses_the_inner_kernel_of_a_kernel_object(self):
        kernel = Kernel(RBF())
        self.gp.kernel = kernel
        assert self.gp._build_gaussian_process(np.asarray([0.01])).kernel is kernel.kernel

    def test_optimize_hyperparameters_floors_the_noise(self):
        """A point without a reported uncertainty must still have noise, as a variance"""
        residuals = self._mock_residuals()
        residuals.loc[0, "MagErr"] = 0.0
        with patch.object(
            self.gp, "_build_gaussian_process", wraps=self.gp._build_gaussian_process
        ) as build:
            self._optimize_hyperparams_with_residuals(residuals)

        alphas = [call.args[0] for call in build.call_args_list]
        assert len(alphas) > 0
        for alpha in alphas:
            assert np.isclose(alpha[0], MIN_PHOTOMETRY_ERROR**2)
            assert np.allclose(alpha[1:], 0.1**2)

    def _run_individually(self, mock_sn, residuals, **kwargs):
        """
        Run run_gp_individually on mocked grids and residuals. The warp, smoothing and random
        sampling are replaced so that the returned grid is the converted surface itself
        (the mean prediction in magnitudes from peak), at the fitted wavelengths and phases.
        """
        mock_sn.zps = {"B": 6.0e-9, "V": 3.0e-9}
        mock_sn.info = {**mock_sn.info, "peak_filt": "B", "peak_mag": 16.0}
        np.random.seed(0)
        with patch(
            "caat.GP3D.GP3D._construct_polynomial_grid", Mock(
                return_value=(
                    self.phase_grid,
                    self.wl_grid,
                    np.zeros((len(self.phase_grid), len(self.wl_grid))),
                    np.ones((len(self.phase_grid), len(self.wl_grid))) * 0.01,
                )
            )
        ), patch(
            "caat.GP3D.GP3D._subtract_data_from_grid", Mock(return_value=residuals)
        ), patch(
            "caat.GP3D.GP3D._iteratively_warp_sed", lambda self, _res, pred, *a, **k: pred
        ), patch(
            "caat.GP3D.GP3D._smooth_predicted_model", lambda self, pred, *a, **k: pred
        ), patch(
            "caat.GP3D.GP3D._sample_predicted_sed", lambda self, grid, std: grid.copy()
        ):
            return self.gp.run_gp_individually(subtract_polynomial=True, **kwargs)

    def _two_filter_residuals(self):
        """Residuals in two filters, at wavelengths that sit exactly on rows of the wavelength grid"""
        return pd.DataFrame(
            [
                {
                    "Filter": filt,
                    "Phase": np.log(phase + 22),
                    "Wavelength": np.log10(wl),
                    "MagResidual": 0.1 * i,
                    "MagErr": 0.1,
                    "Mag": -0.5,
                    "Nondetection": False,
                }
                for filt, wl in (("B", 5000.0), ("V", 6000.0))
                for i, phase in enumerate((5.0, 10.0, 20.0, 30.0))
            ]
        )

    def test_filter_zero_points_are_the_default_for_fitting_individually(self):
        for method in (GP3D.run_gp_individually, GP3D.predict):
            assert inspect.signature(method).parameters["use_filter_zero_points"].default is True

    def test_run_gp_individually_converts_with_the_filter_zero_points_by_default(self, mock_sn):
        """
        By default, the surface should be converted with a zero point that reproduces each
        filter's own zero point at the wavelength of its photometry, not the AB zero point.
        The two conversions differ by -2.5 log10(zp_filter / zp_AB) at each wavelength.
        """
        residuals = self._two_filter_residuals()
        default, _, _, wl_grid = self._run_individually(mock_sn, residuals)
        filter_zps, _, _, _ = self._run_individually(mock_sn, residuals, use_filter_zero_points=True)
        ab_zps, _, _, _ = self._run_individually(mock_sn, residuals, use_filter_zero_points=False)

        assert len(default) > 0
        assert np.allclose(default[0], filter_zps[0], equal_nan=True)

        for filt, wl in (("B", 5000.0), ("V", 6000.0)):
            row = int(np.argmin(abs(10**wl_grid - wl)))
            assert np.isclose(10**wl_grid[row], wl)
            shift = default[0][row] - ab_zps[0][row]
            expected = -2.5 * np.log10(mock_sn.zps[filt] / ab_zero_point(wl))
            # The grid is only filled at the fitted phases, and it is filled in the same places
            assert np.array_equal(np.isfinite(default[0][row]), np.isfinite(ab_zps[0][row]))
            filled = np.isfinite(shift)
            assert filled.sum() > 0
            assert np.allclose(shift[filled], expected)

    def test_run_gp_individually_uses_the_ab_zero_point_when_asked(self, mock_sn):
        residuals = self._two_filter_residuals()
        with patch("caat.GP3D.ab_zero_point", wraps=ab_zero_point) as ab, patch(
            "caat.GP3D.filter_zero_point_function", wraps=filter_zero_point_function
        ) as filter_zp:
            self._run_individually(mock_sn, residuals, use_filter_zero_points=False)
            assert ab.called
            assert not filter_zp.called

            ab.reset_mock()
            self._run_individually(mock_sn, residuals)
            assert filter_zp.called
            assert not ab.called

    # def test_run_gp3d(self):

    #     gp = GP3D(sncollection, kernel)
    #     gaussian_processes, phase_grid, kernel_params, wl_grid = (
    #         gp.run_gp(
    #             ['B', 'g', 'V'],
    #             -20,
    #             50,
    #             log_transform=30,
    #             fit_residuals=True,
    #             set_to_normalize=sncollection.sne,
    #             subtract_polynomial=True,
    #             use_fluxes=use_flux
    #         )
    #     )

    #     assert len(gaussian_processes) > 0 and len(phase_grid) > 0 and len(kernel_params) > 0 and len(wl_grid) > 0
    #     assert gaussian_processes[0].shape <= (len(phase_grid), len(wl_grid))

    def test_interpolate_grid_fills_nans(self):
        """interpolate_grid should replace NaNs between valid values (requires >4 valid points)"""
        # filter_window must be > polyorder (3), so use 5; need >4 non-NaN values per row
        row = np.array([1.0, 1.5, 2.0, np.nan, 3.0, 3.5, 4.0, 4.5, 5.0], dtype=float)
        grid = np.vstack([row, row + 1.0])
        interp_array = np.linspace(0, 1, grid.shape[1])
        result = self.gp.interpolate_grid(grid.copy(), interp_array, filter_window=5)
        # The NaN at index 3 is between valid values and should be filled
        assert not np.isnan(result[0, 3])
        assert not np.isnan(result[1, 3])

    def test_interpolate_grid_preserves_endpoints(self):
        """interpolate_grid should leave NaNs outside the range of valid values"""
        grid = np.array([[np.nan, 1.0, 2.0, np.nan]], dtype=float)
        interp_array = np.linspace(0, 1, 4)
        result = self.gp.interpolate_grid(grid.copy(), interp_array, filter_window=3)
        # Endpoints outside valid data range stay NaN
        assert np.isnan(result[0, 0])
        assert np.isnan(result[0, 3])

    def test_interpolate_grid_all_nan_row_unchanged(self):
        """Rows with all NaNs should remain all NaN"""
        grid = np.array([[np.nan, np.nan, np.nan]], dtype=float)
        interp_array = np.linspace(0, 1, 3)
        result = self.gp.interpolate_grid(grid.copy(), interp_array, filter_window=3)
        assert np.all(np.isnan(result))

    def test_process_dataset_with_set_to_normalize(self):
        """_process_dataset with set_to_normalize returns a DataFrame"""
        with patch("caat.GP.GP._process_dataset", Mock(
            return_value=(
                np.asarray([0.5, 3.5]),
                np.asarray([-1.0, 3.0]),
                np.asarray([0.1, 0.1]),
                np.asarray([3.5, 3.5])
            )
        )):
            template_df = self.gp._process_dataset(
                set_to_normalize=self.gp.set_to_normalize
            )
        assert isinstance(template_df, pd.DataFrame)
        assert set(template_df.columns) >= {"Phase", "Wavelength", "Mag", "MagErr"}

    def test_build_samples_returns_empty_arrays_when_no_data(self):
        """_build_samples should return empty arrays when _process_dataset returns nothing"""
        with patch("caat.GP.GP._process_dataset", Mock(
            return_value=(
                np.asarray([]),
                np.asarray([]),
                np.asarray([]),
                np.asarray([])
            )
        )):
            phases, wls, mags, err_grid = self.gp._build_samples('B')

        assert len(phases) == 0
        assert len(wls) == 0
        assert len(mags) == 0
        assert len(err_grid) == 0

    def test_build_test_wavelength_phase_grid_returns_empty_when_no_phases_in_range(
        self, mock_datacube
    ):
        """Returns empty lists when no phases in the grid match the measured phases"""
        # Put measured phases far outside the wl_grid range
        measured_phases = np.asarray([100.0, 101.0])  # >> phase_grid max
        result = self.gp._build_test_wavelength_phase_grid_from_photometry(
            mock_datacube["LogShiftedWavelength"].values,
            measured_phases,
            self.wl_grid,
            self.phase_grid,
        )
        x, y, wl_inds_fitted, phase_inds_fitted, min_phase = result
        assert len(x) == 0
        assert len(y) == 0
        assert min_phase is None

    def test_sample_predicted_sed_shape(self):
        """_sample_predicted_sed returns array with same shape as input"""
        mean = np.random.random((10, 5))
        std = np.random.random((10, 5)) * 0.1
        result = self.gp._sample_predicted_sed(mean, std)
        assert result.shape == mean.shape

    def test_sample_predicted_sed_within_bounds(self):
        """Sampled SED should be within ±1 sigma of the mean"""
        mean = np.zeros((20, 20))
        std = np.ones((20, 20))
        result = self.gp._sample_predicted_sed(mean, std)
        assert np.all(result >= mean - std)
        assert np.all(result <= mean + std)

    def test_smooth_predicted_model_shape_preserved(self):
        """_smooth_predicted_model should preserve the shape of the input array"""
        model = np.random.random((5, 10))
        result = self.gp._smooth_predicted_model(model, window_size=3)
        assert result.shape == model.shape

    def test_smooth_predicted_model_transpose(self):
        """_smooth_predicted_model with transpose=True should still return original shape"""
        model = np.random.random((5, 10))
        result = self.gp._smooth_predicted_model(model, window_size=3, transpose=True)
        assert result.shape == model.shape

    def test_smooth_predicted_model_even_window_becomes_odd(self):
        """Even window_size should be incremented to odd without raising an error"""
        model = np.random.random((4, 8))
        result = self.gp._smooth_predicted_model(model, window_size=4)
        assert result.shape == model.shape

    def test_optimize_hyperparams_raises_without_subtract_flag(self):
        """optimize_hyperparams should raise when neither subtract flag is set"""
        with patch("caat.GP3D.GP3D._process_dataset", Mock(return_value=pd.DataFrame(
            {"Phase": [], "Wavelength": [], "Mag": [], "MagErr": []}
        ))):
            with pytest.raises(Exception, match=r'Must toggle either .*'):
                self.gp.optimize_hyperparams()

    def test_subtract_data_from_grid_skips_nan_grid_values(self, mock_sn):
        """Residuals for data points where mag_grid is NaN should be skipped"""
        nan_mag_grid = np.full(
            (len(self.phase_grid), len(self.wl_grid)), np.nan
        )
        residuals = self.gp._subtract_data_from_grid(
            mock_sn,
            ['B'],
            self.phase_grid,
            self.wl_grid,
            nan_mag_grid,
            np.ones((len(self.phase_grid), len(self.wl_grid))) * 0.01,
        )
        assert isinstance(residuals, pd.DataFrame)
        assert len(residuals) == 0

    def test_subtract_data_from_grid_columns(self, mock_sn):
        """Residuals DataFrame should contain the expected columns"""
        mag_grid = np.random.random((len(self.phase_grid), len(self.wl_grid)))
        residuals = self.gp._subtract_data_from_grid(
            mock_sn,
            ['B'],
            self.phase_grid,
            self.wl_grid,
            mag_grid,
            np.ones((len(self.phase_grid), len(self.wl_grid))) * 0.01,
        )
        assert isinstance(residuals, pd.DataFrame)
        if len(residuals) > 0:
            expected_cols = {"Filter", "Phase", "Wavelength", "MagResidual", "MagErr", "Mag", "Nondetection"}
            assert expected_cols <= set(residuals.columns)