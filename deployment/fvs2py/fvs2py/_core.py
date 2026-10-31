from __future__ import annotations

import ctypes as ct
import logging
import os
import sys
import threading
import warnings
import weakref
from pathlib import Path

from fvs2py.common import load_dll, unload_dll
from fvs2py.constants import NEEDED_ROUTINES


class SharedEngineWarning(RuntimeWarning):
    """A second live FVS object was created on a library another live object holds.

    ctypes opens a library path once per process, so the two objects share one
    engine (one set of Fortran COMMON blocks). Their simulations, species
    attributes and stop points interfere, and load_keyfile() cannot reset an
    engine that another handle still holds open. Set the environment variable
    FVS2PY_SHARED_ENGINE=error to raise instead of warning, or =ignore to silence.
    """


# One live owner per resolved library path in this process (October 30, 2026).
# A reentrant lock: a finalizer (_release) can run during garbage collection
# triggered inside a locked section of the same thread.
_LIVE_ENGINES: dict[str, weakref.ref] = {}
_LIVE_ENGINES_LOCK = threading.RLock()


def _is_internal_frame(frame) -> bool:
    name = frame.f_globals.get("__name__", "") or ""
    if name == "fvs2py.tests" or name.startswith("fvs2py.tests."):
        return False
    return name == "fvs2py" or name.startswith("fvs2py.")


def _external_stacklevel() -> int:
    """stacklevel for a warnings.warn() issued in the caller: the first frame outside fvs2py."""
    frame = sys._getframe(1)  # the function that will call warnings.warn
    level = 1
    while frame is not None and _is_internal_frame(frame):
        frame = frame.f_back
        level += 1
    return level if frame is not None else 1


def _check_and_register(owner, lib_path: Path) -> None:
    key = os.path.realpath(lib_path)
    msg = None
    with _LIVE_ENGINES_LOCK:
        ref = _LIVE_ENGINES.get(key)
        other = ref() if ref is not None else None
        if other is not None and other is not owner and getattr(other, "_lib", None) is not None:
            mode = os.environ.get("FVS2PY_SHARED_ENGINE", "warn").strip().lower()
            msg = (
                f"Another live FVS object already holds {key}. Both objects share one engine, "
                "so their runs interfere and load_keyfile() cannot reset it. Call close() on "
                "the first object (or let it be garbage collected) before creating a second, "
                "reuse one object for successive keyword files, or run stands in separate processes."
            )
            if mode == "error":
                raise RuntimeError(msg)
            if mode == "ignore":
                msg = None
        _LIVE_ENGINES[key] = weakref.ref(owner)
    if msg is not None:
        warnings.warn(msg, SharedEngineWarning, stacklevel=_external_stacklevel())


def _unregister(owner, lib_path: Path) -> None:
    key = os.path.realpath(lib_path)
    with _LIVE_ENGINES_LOCK:
        ref = _LIVE_ENGINES.get(key)
        obj = ref() if ref is not None else None
        if ref is not None and (obj is owner or obj is None):
            del _LIVE_ENGINES[key]


def _release(lib_path: str, lib: ct.CDLL) -> None:
    """Finalizer of a garbage-collected FvsCore: close its library handle.

    When no other handle holds the library, it leaves the process, so the next
    object on that path loads a fresh engine (fresh COMMON blocks). Never raises.
    """
    try:
        unload_dll(lib)
    except Exception as e:  # a finalizer must not raise
        logging.error(f"fvs2py: could not unload {lib_path} for a collected FVS object: {e}")
    key = os.path.realpath(lib_path)
    with _LIVE_ENGINES_LOCK:
        ref = _LIVE_ENGINES.get(key)
        if ref is not None and ref() is None:
            del _LIVE_ENGINES[key]


class FvsCore:
    """Base class for FVS API wrapper."""

    def __init__(self, lib_path: str | os.PathLike):
        """Loads FVS shared library and checks to ensure needed routines exist.

        Args:
          lib_path : path to FVS library
        """
        self.lib_path: Path = Path(os.path.abspath(lib_path))
        self._lib: ct.CDLL | None = None
        self._finalizer: weakref.finalize | None = None
        self.variant: str = (
            os.path.basename(self.lib_path)
            .split(".")[0]
            .split("FVS")[-1]
            .upper()
        )

        # Declare function attributes with type annotations for mypy
        # These are ctypes foreign function pointers loaded from the shared library
        self._fvs: ct._FuncPointer
        self._fvsAddActivity: ct._FuncPointer
        self._fvsAddTrees: ct._FuncPointer
        self._fvsDimSizes: ct._FuncPointer
        self._fvsEvmonAttr: ct._FuncPointer
        self._fvsFFEAttrs: ct._FuncPointer
        self._fvsGetRestartCode: ct._FuncPointer
        self._fvsGetRtnCode: ct._FuncPointer
        self._fvsGetICCode: ct._FuncPointer
        self._fvsSVSDimSizes: ct._FuncPointer
        self._fvsSetStoppointCodes: ct._FuncPointer
        self._fvsSetCmdLine: ct._FuncPointer
        self._fvsSVSObjData: ct._FuncPointer
        self._fvsSpeciesAttr: ct._FuncPointer
        self._fvsSpeciesCode: ct._FuncPointer
        self._fvsStandID: ct._FuncPointer
        self._fvsSummary: ct._FuncPointer
        self._fvsTreeAttr: ct._FuncPointer
        self._fvsUnitConversion: ct._FuncPointer

        self._load_fvs()
        assert self._lib is not None  # to satisfy type checker

        # check for needed routines that are missing
        missing = []

        for routine in NEEDED_ROUTINES:
            if hasattr(self._lib, routine) and callable(
                getattr(self._lib, routine)
            ):
                logging.debug(f"Found {routine} as expected.")
                # anticipate subroutine name changes depending upon compiler and OS
                # unix pattern on fortran functions
            elif hasattr(self._lib, f"{routine.lower()}_") and callable(
                getattr(self._lib, f"{routine.lower()}_")
            ):
                logging.debug(f"Found {routine} renamed as {routine.lower()}_.")
            else:
                missing.append(routine)

        if len(missing) > 0:
            msg = " ".join(
                [
                    ", ".join(missing),
                    "are needed routines that are not available in library, "
                    "(maybe they weren't exported when library was built)",
                ]
            )
            raise ImportError(msg)
        try:
            self._fvs = self._lib.fvs_
            self._fvsAddActivity = self._lib.fvsaddactivity_
            self._fvsAddTrees = self._lib.fvsaddtrees_
            self._fvsDimSizes = self._lib.fvsdimsizes_
            self._fvsEvmonAttr = self._lib.fvsevmonattr_
            self._fvsFFEAttrs = self._lib.fvsffeattrs_
            self._fvsGetRestartCode = self._lib.fvsgetrestartcode_
            self._fvsGetRtnCode = self._lib.fvsgetrtncode_
            self._fvsGetICCode = self._lib.fvsgeticcode_
            self._fvsSVSDimSizes = self._lib.fvssvsdimsizes_
            self._fvsSetStoppointCodes = self._lib.fvssetstoppointcodes_
            self._fvsSetCmdLine = self._lib.fvssetcmdline_
            self._fvsSVSObjData = self._lib.fvssvsobjdata_
            self._fvsSpeciesAttr = self._lib.fvsspeciesattr_
            self._fvsSpeciesCode = self._lib.fvsspeciescode_
            self._fvsStandID = self._lib.fvsstandid_
            self._fvsSummary = self._lib.fvssummary_
            self._fvsTreeAttr = self._lib.fvstreeattr_
            self._fvsUnitConversion = self._lib.fvsunitconversion_
        except AttributeError:
            self._fvs = self._lib.fvs
            self._fvsAddActivity = self._lib.fvsAddActivity
            self._fvsAddTrees = self._lib.fvsAddTrees
            self._fvsDimSizes = self._lib.fvsDimSizes
            self._fvsEvmonAttr = self._lib.fvsEvmonAttr
            self._fvsFFEAttrs = self._lib.fvsFFEAttrs
            self._fvsGetRestartCode = self._lib.fvsGetRestartCode
            self._fvsGetRtnCode = self._lib.fvsGetRtnCode
            self._fvsGetICCode = self._lib.fvsGetICCode
            self._fvsSVSDimSizes = self._lib.fvsSVSDimSizes
            self._fvsSetStoppointCodes = self._lib.fvsSetStoppointCodes
            self._fvsSetCmdLine = self._lib.fvsSetCmdLine
            self._fvsSVSObjData = self._lib.fvsSVSObjData
            self._fvsSpeciesAttr = self._lib.fvsSpeciesAttr
            self._fvsSpeciesCode = self._lib.fvsSpeciesCode
            self._fvsStandID = self._lib.fvsStandID
            self._fvsSummary = self._lib.fvsSummary
            self._fvsTreeAttr = self._lib.fvsTreeAttr
            self._fvsUnitConversion = self._lib.fvsUnitConversion

        return

    def _load_fvs(self) -> None:
        if self._lib:
            logging.debug("Unloading existing library.")
            self._unload_fvs()
        _check_and_register(self, self.lib_path)
        lib = load_dll(self.lib_path)
        self._lib = lib
        # ctypes never closes a library, so before October 30, 2026 a collected
        # object left its engine loaded and the next object silently inherited its
        # state. The finalizer closes the handle when the object is collected. It is
        # skipped at interpreter exit (atexit False): the process is ending and
        # unloading a Fortran library during shutdown buys nothing.
        self._finalizer = weakref.finalize(self, _release, str(self.lib_path), lib)
        self._finalizer.atexit = False
        return

    def _unload_fvs(self) -> None:
        fin = getattr(self, "_finalizer", None)
        if fin is not None:
            fin.detach()
            self._finalizer = None
        if self._lib:
            unload_dll(self._lib)
            self._lib = None
            _unregister(self, self.lib_path)
        else:
            logging.debug("No library to unload.")
        return

    def close(self) -> None:
        """Release the engine library held by this object.

        After close() another FVS object on the same library can be created
        without a SharedEngineWarning. The object cannot run again until a new
        object is constructed. Garbage collection of the object does the same.
        """
        self._unload_fvs()
        return
