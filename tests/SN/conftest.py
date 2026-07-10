import pytest

from caat import SN


@pytest.fixture
def fake_sn_with_data(filts=("B", "V"), n_obs=3, mag=15.0):
    """Return an SN built from a minimal in-memory data dict (no file I/O).
    type/subtype are passed to bypass the archive name lookup."""
    data = {
        f: [{"mjd": 60000.0 + i, "mag": mag, "err": 0.1} for i in range(n_obs)]
        for f in filts
    }
    return SN(name="test_sn", data=data, type="test", subtype="test", info={"z": 0.0})