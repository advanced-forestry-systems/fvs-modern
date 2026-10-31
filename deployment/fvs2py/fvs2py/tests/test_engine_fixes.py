"""Tests for the fvs2py engine-interface fixes of October 30, 2026 (PR 143).

1. fvsTreeAttr is called with its six arguments (ntrees included).
2. The summary buffer holds the engine's 22 columns per row.
3. The calibrated overlay is applied once per keyword file (at stop point 7 of
   the first stand in the file), not once per object.
4. A second live FVS object on the same library is detected; closing or
   collecting an object closes its library handle, so the next object gets a
   fresh engine.

They need a built variant library: set FVS_LIB_DIR (default <repo>/lib) so that
FVSne.so is found; otherwise they are skipped. The keyword file is the package's
own tests/keyfiles/NE.key (inline TREEDATA, no external data).
Run: python3 -m pytest deployment/fvs2py/fvs2py/tests/test_engine_fixes.py -q
"""
from __future__ import annotations

import gc
import importlib.resources
import os
import sys
import threading
import warnings
from pathlib import Path

import numpy as np
import pytest

from fvs2py import _core, common
from fvs2py._base import FVS
from fvs2py._core import SharedEngineWarning
from fvs2py.constants import STR_NTREES, SUMMARY_COLS

REPO = Path(__file__).resolve().parents[4]
# config.config_loader (the calibrated overlay) is imported from the repository root.
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
LIB = Path(os.environ.get("FVS_LIB_DIR", REPO / "lib")) / "FVSne.so"
KEY_SRC = importlib.resources.files("fvs2py.tests.keyfiles").joinpath("NE.key")

needs_lib = pytest.mark.skipif(not LIB.exists(), reason=f"variant library not found: {LIB}")
needs_maps = pytest.mark.skipif(not Path("/proc/self/maps").exists(), reason="needs /proc/self/maps")


@pytest.fixture(autouse=True)
def _no_live_engines():
    gc.collect()
    _core._LIVE_ENGINES.clear()
    yield
    gc.collect()
    _core._LIVE_ENGINES.clear()


@pytest.fixture
def keyfile(tmp_path):
    p = tmp_path / "ne_test.key"
    p.write_text(KEY_SRC.read_text())
    return p


def _finish(fvs):
    for _ in range(5):
        if fvs.itrncd == 2:
            break
        fvs.run()


def _lib_mapped() -> bool:
    return os.path.realpath(LIB) in Path("/proc/self/maps").read_text()


def _project(keyfile):
    """One complete projection of the test stand on a new object; the object is dropped on return."""
    fvs = FVS(LIB)
    fvs.load_keyfile(keyfile)
    fvs.run()
    _finish(fvs)
    return fvs.summary.copy()


# --------------------------------------------------------------------------- 2
def test_summary_cols_match_engine_width():
    assert len(SUMMARY_COLS) == 22
    assert SUMMARY_COLS[-2:] == ("scuft", "rscuft")


@needs_lib
def test_summary_returns_22_columns(keyfile):
    fvs = FVS(LIB)
    try:
        fvs.load_keyfile(keyfile)
        fvs.run()
        s = fvs.summary
        assert list(s.columns) == list(SUMMARY_COLS)
        assert s.shape[1] == 22 and len(s) >= 2
        assert (s["year"].diff().dropna() > 0).all()
        assert (s[["scuft", "rscuft"]] >= 0).all().all()
    finally:
        fvs.close()


# --------------------------------------------------------------------------- 1
@needs_lib
def test_tree_attr_get_and_set_with_six_arguments(keyfile):
    fvs = FVS(LIB)
    try:
        fvs.load_keyfile(keyfile)
        fvs.run(stop_point_code=1, stop_point_year=-1)  # before the first Event Monitor call, cycle 1
        n = fvs.dims[STR_NTREES]
        assert n > 0
        dbh = fvs.get_tree_attr("dbh").copy()
        assert dbh.shape == (n,) and np.all(np.isfinite(dbh)) and np.all(dbh > 0)
        tab = fvs.tree_table(("id", "species", "dbh", "ht"))
        assert len(tab) == n
        np.testing.assert_allclose(tab["dbh"].to_numpy(), dbh)
        # a real write: every diameter changes, is read back, then is restored
        fvs.set_tree_attr("dbh", dbh + 0.1)
        np.testing.assert_allclose(fvs.get_tree_attr("dbh"), dbh + 0.1, rtol=0, atol=1e-5)
        fvs.set_tree_attr("dbh", dbh)
        np.testing.assert_allclose(fvs.get_tree_attr("dbh"), dbh, rtol=0, atol=1e-5)
        with pytest.raises(NameError):
            fvs.get_tree_attr("nosuch")
        fvs.run()
        _finish(fvs)
    finally:
        fvs.close()


# --------------------------------------------------------------------------- 3
@needs_lib
def test_overlay_applied_once_per_keyfile(keyfile, monkeypatch):
    calls = []
    orig = FVS.apply_calibrated_config

    def spy(self):
        calls.append(self.keyfile_path)
        return orig(self)

    monkeypatch.setattr(FVS, "apply_calibrated_config", spy)
    fvs = FVS(LIB, config_version="calibrated", config_dir=str(REPO / "config"))
    try:
        fvs.load_keyfile(keyfile)
        assert fvs._config_applied is False
        fvs.run()
        _finish(fvs)
        first = fvs.summary
        assert len(calls) == 1 and fvs._config_applied is True
        fvs.load_keyfile(keyfile)
        assert fvs._config_applied is False
        fvs.run()
        _finish(fvs)
        second = fvs.summary
        assert len(calls) == 2 and fvs._config_applied is True
        assert first.equals(second)
    finally:
        fvs.close()


# --------------------------------------------------------------------------- 4
@needs_lib
def test_second_live_object_on_same_library_warns(monkeypatch):
    monkeypatch.delenv("FVS2PY_SHARED_ENGINE", raising=False)
    a = FVS(LIB)
    with pytest.warns(SharedEngineWarning) as rec:
        b = FVS(LIB)
    # the warning points at the line that built the second object, not inside fvs2py
    assert rec[0].filename == __file__
    b.close()
    a.close()
    with warnings.catch_warnings():
        warnings.simplefilter("error", SharedEngineWarning)
        c = FVS(LIB)  # nothing live any more: no warning
        c.close()


@needs_lib
def test_second_live_object_raises_when_strict(monkeypatch):
    a = FVS(LIB)
    try:
        monkeypatch.setenv("FVS2PY_SHARED_ENGINE", "error")
        with pytest.raises(RuntimeError, match="share one engine"):
            FVS(LIB)
    finally:
        a.close()


@needs_lib
def test_concurrent_construction_registers_one_owner(monkeypatch):
    monkeypatch.setenv("FVS2PY_SHARED_ENGINE", "error")
    barrier = threading.Barrier(4)
    built, refused = [], []

    def build():
        barrier.wait()
        try:
            built.append(FVS(LIB))
        except RuntimeError:
            refused.append(1)

    threads = [threading.Thread(target=build) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    try:
        assert len(built) == 1 and len(refused) == 3
    finally:
        for f in built:
            f.close()


@needs_lib
def test_failed_load_does_not_leave_a_registration(monkeypatch):
    monkeypatch.setenv("FVS2PY_SHARED_ENGINE", "error")

    def broken(path):
        raise OSError("simulated load failure")

    with monkeypatch.context() as m:
        m.setattr(_core, "load_dll", broken)
        with pytest.raises(OSError, match="simulated"):
            FVS(LIB)
    assert os.path.realpath(LIB) not in _core._LIVE_ENGINES
    a = FVS(LIB)  # strict mode: would raise if the failed object were still registered
    a.close()


@needs_lib
def test_reload_and_garbage_collected_object_do_not_warn(keyfile, monkeypatch):
    monkeypatch.delenv("FVS2PY_SHARED_ENGINE", raising=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error", SharedEngineWarning)
        a = FVS(LIB)
        a.load_keyfile(keyfile)
        a.run()
        a.load_keyfile(keyfile)  # reload of the same object
        a.run()
        a.close()
        b = FVS(LIB)
        del b
        gc.collect()
        c = FVS(LIB)  # the previous object was collected
        c.close()


def test_dlclose_comes_from_the_process_not_the_library():
    f = common._dlclose()
    assert callable(f) and f.restype is not None


@needs_lib
@needs_maps
def test_close_and_collection_unload_the_library():
    a = FVS(LIB)
    assert _lib_mapped()
    a.close()
    assert not _lib_mapped()
    b = FVS(LIB)
    assert _lib_mapped()
    del b
    gc.collect()
    assert not _lib_mapped()


@needs_lib
def test_fresh_engine_after_collection_reproduces_run(keyfile):
    # one object per projection, never closed: the pattern of run_scenario() in
    # calibration/python/bakuzis_*_comparison.py. Each projection must start from
    # a fresh engine (fresh COMMON blocks) and give the same result.
    first = _project(keyfile)
    gc.collect()
    if Path("/proc/self/maps").exists():
        assert not _lib_mapped()
    second = _project(keyfile)
    gc.collect()
    third = _project(keyfile)
    assert len(first) >= 2
    assert first.equals(second) and first.equals(third)
