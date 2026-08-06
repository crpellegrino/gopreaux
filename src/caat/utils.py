import os
from io import BytesIO
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
from astropy.io.votable import parse

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


FILT_TEL_CONVERSION = {
    "UVW2": "Swift",
    "UVM2": "Swift",
    "UVW1": "Swift",
    "U": "Swift",
    "c": "Atlas",
    "o": "Atlas",
    "B": "Swift",
    "V": "Swift",
    "g": "ZTF",
    "r": "ZTF",
    "i": "ZTF",
    "R": "CTIO",
    "I": "CTIO",
    "J": "CTIO",
    "H": "CTIO",
    "K": "CTIO",
    "Y": "DECam",
    "u": "DECam",
    "G": "GAIA",
    "y": "PAN-STARRS",
    "z": "PAN-STARRS",
    "w": "PAN-STARRS",
}

colors = {
    "U": "purple",
    "B": "blue",
    "V": "lime",
    "g": "cyan",
    "r": "orange",
    "i": "red",
    "UVW2": "#FE0683",
    "UVM2": "#BF01BC",
    "UVW1": "#8B06FF",
    "c": "turquoise",
    "o": "salmon",
}

WLE = {
    "u": 3560,
    "g": 4830,
    "r": 6260,
    "i": 7670,
    "z": 8890,
    "y": 9600,
    "w": 5985,
    "u'": 3560,
    "g'": 4830,
    "r'": 6260,
    "i'": 7670,
    "z'": 8890,
    "y'": 9600,
    "w'": 5985,
    "Y": 9600,
    "U": 3600,
    "B": 4380,
    "V": 5450,
    "R": 6410,
    "G": 6730,
    "E": 6730,
    "I": 7980,
    "J": 12200,
    "H": 16300,
    "K": 21900,
    "Ks": 21900,
    "UVW2": 2030,
    "UVM2": 2231,
    "UVW1": 2634,
    "F": 1516,
    "N": 2267,
    "o": 6790,
    "c": 5330,
    "W": 33526,
    "Q": 46028,
}


def query_svo_service(instrument: str, filter: str):
    """
    Query the SVO Filter Service to retrieve filter curves
    for a given instrument and filter.

    NOTE: Many of these are approximations for the true data in
    our sample. 

    Args:
        instrument (str): The name of the instrument that took the data.
        filter (str): The name of the filter of the data.

    Returns:
        astropy.Table: An astropy Table object containing the
            filter curve information.
    """
    base_url = "http://svo2.cab.inta-csic.es/theory/fps/fps.php?"
    if instrument.lower() == "swift":
        url = base_url + f"ID={instrument}/UVOT.{filter}"
    elif instrument.lower() == "atlas":
        filter_dict = {"o": "orange", "c": "cyan"}
        url = base_url + f"ID=Misc/ATLAS.{filter_dict[filter]}"
    elif instrument.lower() == "ztf":
        url = base_url + f"ID=Palomar/ZTF.{filter}"
    elif instrument.lower() == "ctio" and filter not in ["J", "H", "K"]:
        url = base_url + f"ID=CTIO/ANDICAM.{filter}_KPNO"
    elif instrument.lower() == "ctio":
        url = base_url + f"ID=CTIO/ANDICAM.{filter}"
    elif instrument.lower() == "decam":
        url = base_url + f"ID=CTIO/DECam.{filter}"
    elif instrument.lower() == "gaia":
        url = base_url + "ID=GAIA/GAIA0.G"
    elif instrument.lower() == "pan-starrs":
        url = base_url + f"ID=PAN-STARRS/PS1.{filter}"
    else:
        url = base_url + f"ID={instrument}/{instrument}.{filter}"
    s = BytesIO(urlopen(url).read())
    table = parse(s).get_first_table().to_table(use_names_over_ids=True)
    return table["Wavelength"], table["Transmission"]


def bin_spec(
    wl: np.ndarray,
    flux: np.ndarray,
    wl2: list | np.ndarray,
    plot: bool = False
):
    """
    Bin a spectrum to a certain resolution
    
    Args:
        wl (np.ndarray): The input wavelength array
        flux (np.ndarray): The input flux array
        wl2 (list | np.ndarray): The wavelength array to bin to
        plt (bool, optional): Plot the binned and unbinned arrays to compare.
            Defaults to False.

    Returns:
        tuple(np.ndarray): The binned wavelength and flux arrays.
    """

    binned_wl = []
    binned_flux = []
    for i in range(len(wl2)):
        ind = np.argmin(abs(wl - wl2[i]))
        if wl[ind] not in binned_wl:
            binned_wl.append(wl[ind])
            binned_flux.append(flux[ind])

    if plot:
        plt.plot(binned_wl, binned_flux, color="blue", label="Binned", alpha=0.5)
        plt.plot(wl, flux, color="orange", label="Original", alpha=0.5)
        plt.legend()
        plt.show()
    return np.asarray(binned_wl), np.asarray(binned_flux)


def convert_shifted_fluxes_to_shifted_mags(
    fluxes: np.ndarray,
    sn,
    zp_at_wl: float
):
    """
    Convert between log flux relative to peak to mags relative to peak.

    Args:
        fluxes (np.ndarray): An array of shifted (relative to peak)
            flux measurements.
        sn (SN): The SN object the fluxes belong to.
        zp_at_wl (float): The zeropoint of the filter at the given effective wavelength.
    """
    shifted_peak_mag = np.log10(
        sn.zps[sn.info["peak_filt"]] * 10 ** (-0.4 * sn.info["peak_mag"])
    )
    shifted_mags = -1 * (
        (np.log10(10 ** (fluxes + shifted_peak_mag) / zp_at_wl) / -0.4)
        - sn.info["peak_mag"]
    )
    return shifted_mags
