"""Smoke test for the thermorawfile Python bindings.

Build first, then run:
    cd thermorawfile-connector && maturin develop --release
    pytest tests/
"""
import os

import numpy as np
import pytest
import thermorawfile

FIXTURE = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "data", "small2.RAW")


def test_open_and_metadata():
    rf = thermorawfile.RawFile(FIXTURE)
    assert rf.version == 66
    assert rf.n_scans == 95 == len(rf)
    assert rf.instrument_model == "LTQ Orbitrap Velos"
    assert rf.acquired == "2018-04-03 19:15:49"
    assert rf.checksum_valid()
    assert "LTQ Orbitrap Velos" in repr(rf)


def test_peaks_are_numpy():
    rf = thermorawfile.RawFile(FIXTURE)
    mz, intensity = rf.peaks(2)
    assert mz.dtype == np.float64 and intensity.dtype == np.float32
    assert mz.shape == (196,) and intensity.shape == (196,)
    assert abs(mz[0] - 116.0264) < 1e-3
    assert abs(intensity[0] - 12.1) < 0.5


def test_filter_event_and_params():
    rf = thermorawfile.RawFile(FIXTURE)
    assert rf.scan_filter(1) == "FTMS + p NSI Full ms [350.00-1200.00]"
    assert rf.scan_filter(2) == "ITMS + c NSI d w Full ms2 398.54@cid35.00 [95.00-1210.00]"

    ev = rf.scan_event(2)
    assert ev["ms_order"] == 2
    assert abs(ev["isolation_center"] - 398.5411) < 1e-3

    sp = rf.scan_params(2)
    assert sp["Charge State:"] == 3
    assert sp["Ion Injection Time (ms):"] == 50.0
    assert abs(sp["Monoisotopic M/Z:"] - 398.5411) < 1e-3


def test_author_centroids_roundtrip(tmp_path):
    rf = thermorawfile.RawFile(FIXTURE)
    mz = np.array([200.12, 400.34, 600.56], dtype=np.float64)  # fits scan-2's packet budget
    intensity = np.array([5000.0, 2500.0, 1250.0], dtype=np.float32)
    rf.author_centroids(2, mz, intensity)
    out = str(tmp_path / "authored.raw")
    rf.save(out)

    rf2 = thermorawfile.RawFile(out)
    m2, i2 = rf2.peaks(2)
    assert np.allclose(m2, mz)
    assert np.allclose(i2, intensity)
    assert rf2.checksum_valid()  # the keyless checksum recomputes on save — hence the need for mzprov


def test_write_input_validation():
    rf = thermorawfile.RawFile(FIXTURE)
    with pytest.raises(ValueError):  # length mismatch
        rf.author_centroids(2, np.array([1.0]), np.array([1.0, 2.0], dtype=np.float32))
    with pytest.raises(ValueError):  # non-positive m/z
        rf.author_centroids(2, np.array([-5.0]), np.array([1.0], dtype=np.float32))
