import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

from caat import SN
from caat.utils import ROOT_DIR, WLE
from caat.utils import colors as FILTER_COLORS



CUBE_ROOT = os.path.join(ROOT_DIR, "data/")


def _nightly_bin(df: pd.DataFrame) -> pd.DataFrame:
    """
    Bin intra-night repeat exposures down to one row per filter per night,
    using an inverse-variance weighted mean for magnitude and combined error.
    Rows with Magerr == 0 are excluded from weighting (treated as missing error).
    """
    df = df.copy()
    if df.empty:
        return df
    
    df["night"] = df["MJD"].astype(int)

    rows = []
    for (filt, night), grp in df.groupby(["ShiftedFilter", "night"]):
        if len(grp) == 1:
            rows.append(grp.iloc[0].to_dict())
            continue

        good = grp[grp["Magerr"] > 0]
        if len(good) == 0:
            # No usable errors — fall back to simple mean
            row = grp.iloc[0].to_dict()
            row["ShiftedMag"] = grp["ShiftedMag"].mean()
            row["MJD"] = grp["MJD"].mean()
            row["Phase"] = grp["Phase"].mean()
        else:
            w = 1.0 / good["Magerr"].values ** 2
            row = good.iloc[0].to_dict()
            row["ShiftedMag"] = float(np.average(good["ShiftedMag"].values, weights=w))
            row["Magerr"] = float(1.0 / np.sqrt(w.sum()))
            row["MJD"] = float(np.average(good["MJD"].values, weights=w))
            row["Phase"] = float(np.average(good["Phase"].values, weights=w))
        rows.append(row)

    return pd.DataFrame(rows).drop(columns=["night"]).reset_index(drop=True)


def load_datacube(
    sn_name: str,
    sn_type: str,
    sn_subtype: str,
    min_phase: float,
    max_phase: float,
) -> pd.DataFrame:
    """Load the mangled datacube for a single SN, keeping only real detections
    within [min_phase, max_phase], and binning any intra-night repeat exposures
    down to one point per night."""
    path = os.path.join(CUBE_ROOT, sn_type, sn_subtype, sn_name, f"{sn_name}_datacube_mangled.csv")
    cube = pd.read_csv(path, index_col=0)
    # Cast explicitly: pd.read_csv reads booleans as strings or ints from CSV
    cube = cube[~cube["Nondetection"].astype(bool)].reset_index(drop=True)
    cube = cube[(cube["Phase"] >= min_phase) & (cube["Phase"] <= max_phase)].reset_index(drop=True)
    cube = _nightly_bin(cube)
    cube["MagFromPeak"] = SN(name=sn_name).info["peak_mag"] - cube["Mag"]
    return cube


def get_observed_points(cube: pd.DataFrame, band_map: dict) -> pd.DataFrame:
    """
    Return a tidy DataFrame of observed detections for filters in band_map.

    ShiftedMag = peak_mag - Mag  (magnitude units, 0 at peak, negative for fading).
    """
    rows = []
    for gopreaux_filt, comp_filt in band_map.items():
        subset = cube[cube["ShiftedFilter"] == gopreaux_filt]
        for _, row in subset.iterrows():
            rows.append({
                "Filter": gopreaux_filt,
                "CompFilter": comp_filt,
                "Phase": row["Phase"],
                "Wavelength": row["ShiftedWavelength"],
                "ShiftedMag": row["MagFromPeak"],
                "Magerr": row["Magerr"],
            })
    return pd.DataFrame(rows)


def _gopreaux_filter_peak(model, wavelength: float) -> float:
    """
    Evaluate the GOPREAUX model on a dense phase grid at the given effective
    wavelength and return the peak log10(flux) value (maximum = peak brightness).

    Used to compute per-filter offsets relative to any reference filter.
    For the model's native reference (typically g-band), returns ~0.
    Bluer/brighter filters return positive values; redder/dimmer return negative.
    """
    phases_c = np.clip(
        np.linspace(model.min_phase + 0.5, model.max_phase - 0.5, 500),
        model.min_phase + 0.01,
        model.max_phase - 0.01,
    )
    wl_c = np.clip(wavelength, model.min_wl + 1, model.max_wl - 1)

    log_phases = np.log(phases_c + model.log_transform)
    log_waves  = np.full_like(log_phases, np.log10(wl_c))

    mean_log10_flux, _ = model.surface.predict(
        np.vstack((log_phases, log_waves)).T, return_std=True
    )

    wl_idx = np.argmin(np.abs(model.wl_grid - wl_c))
    template_vals = np.array([
        model.template[np.argmin(np.abs(model.phase_grid - p)), wl_idx]
        for p in phases_c
    ])

    return float(np.nanmax(mean_log10_flux + template_vals))


def _gopreaux_dense_prediction(model, wavelength: float, phase_min: float, phase_max: float, n: int = 200):
    """
    Evaluate the GOPREAUX model on a dense phase grid at a single wavelength.

    Returns:
        phases      : 1-D array of linear phases (length n)
        log10_flux  : mean log10(flux) prediction at each phase
        dev         : 1-sigma uncertainty at each phase
    """
    phases_c = np.clip(
        np.linspace(phase_min, phase_max, n),
        model.min_phase + 0.01,
        model.max_phase - 0.01,
    )
    wl_c = np.clip(wavelength, model.min_wl + 1, model.max_wl - 1)

    log_phases = np.log(phases_c + model.log_transform)
    log_waves  = np.full_like(log_phases, np.log10(wl_c))

    mean_log10_flux, dev = model.surface.predict(
        np.vstack((log_phases, log_waves)).T, return_std=True
    )

    wl_idx = np.argmin(np.abs(model.wl_grid - wl_c))
    template_vals = np.array([
        model.template[np.argmin(np.abs(model.phase_grid - p)), wl_idx]
        for p in phases_c
    ])

    return phases_c, mean_log10_flux + template_vals, dev


def predict_gopreaux(
        model,
        obs: pd.DataFrame,
        peak_filt: str | None = None,
        z: float = 0.0,
        plot: bool = True
    ) -> pd.DataFrame:
    """
    Evaluate the GOPREAUX model against observed photometry, normalised to
    the per-SN peak filter rather than a fixed g-band reference.

    For each filter the model is evaluated on a dense phase grid (smooth curve
    for plotting), then interpolated to the observed phases for residuals.

    The reference offset  gop_ref_peak = _gopreaux_filter_peak(model, WLE[peak_filt])
    is subtracted from every filter's prediction so that the model zero-point
    matches the observed MagFromPeak zero-point (0 at peak_filt peak).

    Returns obs with added columns (log10(flux) units, relative to peak_filt):
        GOPREAUX_mag  – model prediction interpolated to each observed phase
        GOPREAUX_err  – 1-sigma uncertainty interpolated to each observed phase
    """
    if obs.empty:
        return obs.assign(GOPREAUX_mag=np.nan, GOPREAUX_err=np.nan)

    # Compute once: the model's predicted peak in the SN's reference filter
    gop_ref_peak = _gopreaux_filter_peak(model, WLE[peak_filt]) if peak_filt is not None else 0

    result = obs.copy()
    result["GOPREAUX_mag"] = np.nan
    result["GOPREAUX_err"] = np.nan

    filters = obs["Filter"].unique()

    if plot:
        fig, axes = plt.subplots(
            1, len(filters),
            figsize=(4 * len(filters), 4),
            sharey=True,
            squeeze=False,
        )

    for idx, filt in enumerate(filters):
        mask = obs["Filter"] == filt
        grp  = obs[mask].sort_values("Phase")

        wavelength = float(grp["Wavelength"].median()) / (1 + z)
        phase_min  = float(grp["Phase"].min())
        phase_max  = float(grp["Phase"].max())

        dense_phases, dense_log10_flux, dense_dev = _gopreaux_dense_prediction(
            model, wavelength, phase_min, phase_max
        )

        # Shift so model zero-point aligns with the SN's peak_filt reference
        dense_log10_flux_shifted = dense_log10_flux - gop_ref_peak

        interp_mag = interp1d(dense_phases, dense_log10_flux_shifted, bounds_error=False, fill_value=np.nan)
        interp_err = interp1d(dense_phases, dense_dev,                bounds_error=False, fill_value=np.nan)

        result.loc[mask, "GOPREAUX_mag"] = interp_mag(obs.loc[mask, "Phase"].values)
        result.loc[mask, "GOPREAUX_err"] = interp_err(obs.loc[mask, "Phase"].values)

        if plot:
            ax    = axes[0][idx]
            color = FILTER_COLORS.get(filt, "steelblue")

            ax.errorbar(
                grp["Phase"], grp["ShiftedMag"], yerr=grp["Magerr"],
                fmt="o", color=color, label="Observed", zorder=3,
            )
            ax.plot(dense_phases, dense_log10_flux_shifted, color="k", label="GOPREAUX", zorder=2)
            ax.fill_between(
                dense_phases,
                dense_log10_flux_shifted - dense_dev,
                dense_log10_flux_shifted + dense_dev,
                color="k", alpha=0.15, zorder=1,
            )

            ax.set_xlabel("Phase (days)")
            ax.set_title(filt)
            ax.legend(fontsize=8)

    if plot:
        axes[0][0].set_ylabel(f"log₁₀(flux) relative to {peak_filt}-peak")
        fig.suptitle("GOPREAUX predictions vs. observed", y=1.02)
        plt.tight_layout()
        plt.show()

    return result


def compute_metrics(residuals: np.ndarray, errors: np.ndarray) -> dict:
    """
    Compute comparison metrics for a set of residuals (model - observed).

    Args:
        residuals: model_mag - observed_mag at each point
        errors:    photometric uncertainty at each point

    Returns:
        dict with keys: n, bias, mad, iqr, chi2_reduced
    """
    # Use a joint mask so residuals and errors stay the same length
    mask = ~np.isnan(residuals) & ~np.isnan(errors) & (errors > 0)
    residuals = residuals[mask]
    errors    = errors[mask]

    # Use a floor of 0.05 mag for errors, to avoid inflated chi sq. values
    small_error_mask = errors < 0.05
    errors[small_error_mask] = 0.05

    if len(residuals) == 0:
        return {"n": 0, "bias": np.nan, "mad": np.nan, "iqr": np.nan, "chi2_reduced": np.nan}

    bias = np.median(residuals)
    mad  = np.median(np.abs(residuals - bias))
    iqr  = np.subtract(*np.percentile(residuals, [75, 25]))
    chi2_reduced = np.mean((residuals / errors) ** 2)

    return {"n": len(residuals), "bias": bias, "mad": mad, "iqr": iqr, "chi2_reduced": chi2_reduced}