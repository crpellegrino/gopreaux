from caat import SN, CAAT
from caat.utils import WLE
from extinction import ccm89 as ccm
import numpy as np
import pytest


def test_bad_sn_name():
    try:
        sn = SN(name='bogus_name')
        assert False
    except Exception as e:
        assert True

def test_good_sn_name():
    caat = CAAT().caat
    sn = SN(name=caat.Name.values[0])
    assert sn.name != ""

def test_create_sn_from_data_dict():
    sn = SN(data={'B': [{'mjd': 60000.0, 'mag': 15.0, 'err': 0.1}]})
    assert len(sn.data) > 0


def test_read_sn_info_from_caat():
    caat = CAAT().caat
    row = caat.sample()
    sn = SN(name=row["Name"].values[0])

    assert (sn.info.get("peak_mjd", np.nan) == row["Tmax"].values[0]) or (np.isnan(sn.info.get("peak_mjd", np.nan)) and np.isnan(row["Tmax"].values[0]))
    assert (sn.info.get("peak_mag", np.nan) == row["Magmax"].values[0]) or (np.isnan(sn.info.get("peak_mag", np.nan)) and np.isnan(row["Magmax"].values[0]))
    assert (sn.info.get("peak_filt", np.nan) == row["Filtmax"].values[0]) or (np.isnan(sn.info.get("peak_filt", np.nan)) and np.isnan(row["Filtmax"].values[0]))

def test_load_swift_file():
    sn = SN(name='SN2022acko')
    sn.data = {}
    sn.load_swift_data()
    assert len(sn.data) > 0

def test_load_json_data():
    sn = SN(name='SN2022acko')
    sn.data = {}
    sn.load_json_data()
    assert len(sn.data) > 0

# def test_convert_to_fluxes():
#     sn = SN(name='SN2022acko')
#     sn.load_json_data()
#     sn.load_swift_data()
#     sn.convert_to_fluxes()
#     assert all([d.get('flux', False) for f in sn.data.keys() for d in sn.data[f]])

def test_extinction_correction():
    sn = SN(name='SN2022acko')
    sn.load_json_data()
    sn.load_swift_data()
    sn.correct_for_galactic_extinction()
    assert all([d.get('ext_corrected', False) for f in sn.data.keys() for d in sn.data[f]])

def test_shift_to_max():
    sn = SN(name='SN2022acko')
    sn.load_json_data()
    sn.load_swift_data()
    mjds, _, _, _ = sn.shift_to_max(list(sn.data.keys())[0])
    assert len(mjds) > 0 and all([mjd < 50000.0 for mjd in mjds])


# ---------------------------------------------------------------------------
# correct_for_host_extinction tests
# ---------------------------------------------------------------------------

def test_host_ext_reduces_mag(fake_sn_with_data):
    """Extinction correction should brighten (reduce) magnitudes."""
    mag_before = fake_sn_with_data.data["B"][0]["mag"]
    fake_sn_with_data.correct_for_host_extinction(a_v=1.0)
    assert fake_sn_with_data.data["B"][0]["mag"] < mag_before


def test_host_ext_magnitude_matches_ccm(fake_sn_with_data):
    """Corrected mag should equal original minus CCM89 extinction at that wavelength."""
    filt = "B"
    a_v, r_v = 1.0, 3.1
    mag_before = fake_sn_with_data.data[filt][0]["mag"]
    expected_ext = ccm(np.asarray([float(WLE[filt])]), a_v, r_v)[0]
    fake_sn_with_data.correct_for_host_extinction(a_v=a_v, r_v=r_v)
    assert fake_sn_with_data.data[filt][0]["mag"] == pytest.approx(mag_before - expected_ext)


def test_host_ext_sets_flag(fake_sn_with_data):
    """All photometry dicts should carry host_ext_corrected=True after the call."""
    fake_sn_with_data.correct_for_host_extinction(a_v=1.0)
    assert all(
        d.get("host_ext_corrected", False)
        for f in fake_sn_with_data.data
        for d in fake_sn_with_data.data[f]
    )


def test_host_ext_should_not_double_correct(fake_sn_with_data):
    """Calling correct_for_host_extinction twice must not double-correct."""
    fake_sn_with_data.correct_for_host_extinction(a_v=1.0)
    mag_after_first = fake_sn_with_data.data["B"][0]["mag"]
    fake_sn_with_data.correct_for_host_extinction(a_v=1.0)
    assert fake_sn_with_data.data["B"][0]["mag"] == pytest.approx(mag_after_first)


def test_host_ext_zero_av(fake_sn_with_data):
    """a_v=0 should leave all magnitudes unchanged."""
    mags_before = {f: [d["mag"] for d in fake_sn_with_data.data[f]] for f in fake_sn_with_data.data}
    fake_sn_with_data.correct_for_host_extinction(a_v=0.0)
    for f in fake_sn_with_data.data:
        for i, d in enumerate(fake_sn_with_data.data[f]):
            assert d["mag"] == pytest.approx(mags_before[f][i])


def test_host_ext_unknown_filter_cleared(fake_sn_with_data):
    """Filters absent from self.wle should be cleared to an empty list."""
    fake_sn_with_data.data["BOGUS_FILT"] = [{"mjd": 60000.0, "mag": 15.0, "err": 0.1}]
    fake_sn_with_data.correct_for_host_extinction(a_v=1.0)
    assert fake_sn_with_data.data["BOGUS_FILT"] == []


def test_host_ext_corrects_shifted_data(fake_sn_with_data):
    """shifted_data should be corrected alongside data when present."""
    fake_sn_with_data.shifted_data = {
        "B": [{"mjd": 60000.0 + i, "mag": 15.0, "err": 0.1} for i in range(3)]
    }
    fake_sn_with_data.correct_for_host_extinction(a_v=1.0)
    assert all(d.get("host_ext_corrected", False) for d in fake_sn_with_data.shifted_data["B"])
    assert all(d["mag"] < 15.0 for d in fake_sn_with_data.shifted_data["B"])