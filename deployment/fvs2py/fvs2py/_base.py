from __future__ import annotations
import ctypes as ct
import logging
import os
import warnings
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd

from fvs2py._core import FvsCore
from fvs2py.common import class_requires_fvs_library, fvs_property
from fvs2py.constants import (
    FVS_ITRNCD_FINISHED_ALL_STANDS,
    FVS_ITRNCD_NOT_STARTED,
    MGMT_ID_COLUMN_NAME,
    SPECIES_ATTRS,
    SPECIES_COLUMN_NAMES,
    STAND_CN_COLUMN_NAME,
    STAND_ID_COLUMN_NAME,
    STR_C_CONTIGUOUS,
    STR_MAXCYCLES,
    STR_MAXPLOTS,
    STR_MAXSPECIES,
    STR_MAXTREES,
    STR_NCYCLES,
    STR_NPLOTS,
    STR_NTREES,
    SUMMARY_COLS,
)
from fvs2py.enums import FvsAttributeAccessor


@class_requires_fvs_library
class FVS(FvsCore):
    """Main class for interacting with FVS at runtime.

    Args:
        lib_path: Path to the FVS shared library (.so or .dll)
        config_version: Which parameter set to use:
            'default'    = original FVS parameters (compiled in)
            'calibrated' = Bayesian posterior estimates from FIA data
            'custom'     = user supplied JSON calibrated to independent data
            None         = do not apply any config overlay (same as 'default')
        config_dir: Override path to the config directory containing
            variant JSON files and calibrated/ subdirectory
        custom_config: Path to a user supplied JSON config file.
            Required when config_version='custom'. Allows calibration
            to independent datasets (cooperative plots, regional
            inventories, silvicultural trials, etc.)
        uncertainty: Enable Monte Carlo uncertainty estimation using
            posterior draws from the Bayesian calibration. When True,
            run_ensemble() will execute FVS n_draws times with different
            parameter sets sampled from the posterior distribution.
        n_draws: Number of posterior draws to use for uncertainty
            estimation. Ignored when uncertainty=False. Default is 100.
        seed: Random seed for reproducible uncertainty draws.
    """

    def __init__(
        self,
        lib_path: str | os.PathLike,
        config_version: str | None = None,
        config_dir: str | os.PathLike | None = None,
        custom_config: str | os.PathLike | None = None,
        uncertainty: bool = False,
        n_draws: int = 100,
        seed: int | None = None,
    ):
        super().__init__(lib_path=lib_path)
        self._initialize_attributes()
        self._config_version = config_version
        self._config_dir = config_dir
        self._custom_config = custom_config
        self._config_applied = False
        self._uncertainty = uncertainty
        self._n_draws = n_draws
        self._seed = seed
        self._uncertainty_engine = None
        self._ensemble_results: list[pd.DataFrame] = []
        return

    def _initialize_attributes(self) -> None:
        self.keyfile_path: Path | None = None
        self.keyfile: str | None = None
        self._exit_code = ct.c_int(0)
        self._itrncd = ct.c_int(-1)
        self._maxcycles = ct.c_int(0)
        self._maxplots = ct.c_int(0)
        self._maxspecies = ct.c_int(0)
        self._maxtrees = ct.c_int(0)
        self._mgmt_id = ct.create_string_buffer(4)
        self._ncycles = ct.c_int(0)
        self._nplots = ct.c_int(0)
        self._ntrees = ct.c_int(0)
        self._restart_code = ct.c_int(0)
        self._species_attrs = dict.fromkeys(SPECIES_ATTRS)
        self._stand_cn = ct.create_string_buffer(40)
        self._stand_id = ct.create_string_buffer(26)
        self._stop_point_code = None
        self._stop_point_year = None

    @fvs_property
    def dims(self) -> dict:
        """Return the max dimensions of important FVS data storage."""
        self._fvsDimSizes.argtypes = [
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
        ]
        self._fvsDimSizes.restype = None

        self._dims = {
            STR_NTREES: self._ntrees,
            STR_NCYCLES: self._ncycles,
            STR_NPLOTS: self._nplots,
            STR_MAXTREES: self._maxtrees,
            STR_MAXSPECIES: self._maxspecies,
            STR_MAXPLOTS: self._maxplots,
            STR_MAXCYCLES: self._maxcycles,
        }

        self._fvsDimSizes(
            self._ntrees,
            self._ncycles,
            self._nplots,
            self._maxtrees,
            self._maxspecies,
            self._maxplots,
            self._maxcycles,
        )
        return {key: val.value for key, val in self._dims.items()}

    @fvs_property
    def exit_code(self) -> int:
        """Gets the integer code returned when FVS exits.

        Possible values are:
          0 - No serious errors occurred.
          1 - Input data error.
          2 - Keyword or expression error.
          3 - Extension or group activities error.
          4 - Scratch file error.
        """
        self._fvsGetICCode.argtypes = [ct.POINTER(ct.c_int)]
        self._fvsGetICCode.restype = None
        self._fvsGetICCode(self._exit_code)

        return self._exit_code.value

    @fvs_property
    def itrncd(self) -> int:
        """Returns with the current return code value in FVS.

        -1: indicates that FVS has not been started.
         0: indicates that FVS is in good running state.
         1: indicates that FVS has detected an error of some kind and should not
                be used until reset by specifying new input.
         2: indicates that FVS has finished processing all the stands; new input
                can be specified.
        """
        self._fvsGetRtnCode.argtypes = [ct.POINTER(ct.c_int)]
        self._fvsGetRtnCode.restype = None

        self._fvsGetRtnCode(self._itrncd)

        return self._itrncd.value

    @fvs_property
    def restart_code(self) -> int:
        """A code indicating when FVS stopped.

          1: Stop was done just before the first call to the Event Monitor.
          2: Stop was done just after the first call to the Event Monitor.
          3: Stop was done just before the second call to the Event Monitor.
          4: Stop was done just after the second call to the Event Monitor.
          5: Stop was done after growth and mortality has been computed, but
                prior to applying them.
          6: Stop was done just before the ESTAB routines are called.
        100: Stop was done after a stand has been simulated but prior to
                starting a subsequent stand.
        """
        self._fvsGetRestartCode.argtypes = [ct.POINTER(ct.c_int)]
        self._fvsGetRestartCode.restype = None
        self._fvsGetRestartCode(self._restart_code)

        return self._restart_code.value

    @fvs_property
    def species(self) -> pd.DataFrame:
        """Returns species codes and attributes for all species."""
        codes = self.species_codes
        attrs = self.species_attrs
        return codes.merge(attrs, left_index=True, right_index=True, copy=False)

    @fvs_property
    def species_codes(self) -> pd.DataFrame:
        """Fetch the various codes used to refer to different tree species."""
        self._fvsSpeciesCode.argtypes = [
            ct.c_char_p,  # FVS species code
            ct.c_char_p,  # FIA species code
            ct.c_char_p,  # PLANTS code
            ct.POINTER(ct.c_int),  # species index for this variant
            ct.POINTER(ct.c_int),  # num char for fvs code
            ct.POINTER(ct.c_int),  # num char for fia code
            ct.POINTER(ct.c_int),  # num char for plants code
            ct.POINTER(ct.c_int),  # return code
        ]
        self._fvsSpeciesCode.restype = None
        dims = self.dims

        _fvs_spp = ct.create_string_buffer(4)
        _fia_spp = ct.create_string_buffer(4)
        _plants_spp = ct.create_string_buffer(6)
        returncd = ct.c_int(0)

        spp_codes = pd.DataFrame(
            index=range(dims[STR_MAXSPECIES]),
            columns=SPECIES_COLUMN_NAMES,
        )

        for i in range(dims[STR_MAXSPECIES]):
            self._fvsSpeciesCode(
                _fvs_spp,
                _fia_spp,
                _plants_spp,
                ct.c_int(i + 1),
                ct.c_int(0),
                ct.c_int(0),
                ct.c_int(0),
                ct.c_int(0),
            )
            if returncd.value != 0:
                msg = f"Index {i + 1} out of range"
                raise IndexError(msg)
            spp_codes.iloc[i] = (
                i + 1,
                _fvs_spp.value.decode().strip(),
                _fia_spp.value.decode().strip(),
                _plants_spp.value.decode().strip(),
            )

        return spp_codes

    @fvs_property
    def species_attrs(self) -> pd.DataFrame:
        """Returns a dataframe of species attributes.

        Fields returned are:
            spccf: CCF for each species, recomputed in FVS so setting will
                likely have no effect
            spsdi: SDI maximums for each species
            spsiteindx: Species site indices
            bfmind: Min diameter related to BFVOLUME keyword
            bftopd: Top diameter related to BFVOLUME keyword
            bfstmp: Stump height related to BFVOLUME keyword
            frmcls: Form class related to BFVOLUME keyword
            bfmeth: Volume calculation code related to BFVOLUME keyword
                (internal FVS variable methb)
            mcmind: Min diameter related to VOLUME keyword (internal FVS
                variable dbhmin)
            mctopd: Top diameter related to VOLUME keyword (internal FVS
                variable topd)
            mcstmp: Stump height related to VOLUME keyword (internal FVS
                variable stmp)
            mcmeth: Volume calculation code related to VOLUME keyword (internal
                FVS variable methc)
            baimult: Basal area increment multiplier for large trees (internal
                FVS variable xdmult)
            htgmult: Height growth multiplier for large trees (internal FVS
                variable xhmult)
            mortmult: Mortality rate multiplier (internal FVS variable xmmult)
            mortdia1: Lower diameter limit for mortality multiplier (internal
                FVS variable xmdia1)
            mortdia2: Upper diameter limit for mortality multiplier (internal
                FVS variable xmdia2)
            regdmult: Diameter growth mulitplier for regeneration (internal FVS
                variable xrdmlt)
            reghmult: Height growth multiplier for regeneration (internal FVS
                variable xrhmlt)
        """
        for attr in self._species_attrs:
            _ = self.get_species_attr(attr)

        attrs = pd.DataFrame(self._species_attrs, copy=False)
        if (attrs == 0).all().all():
            warnings.warn("No species attributes initialized yet.")
            return attrs.replace(0, None)

        return attrs

    @fvs_property
    def stand_ids(self) -> dict:
        """Return stand identification codes."""
        self._fvsStandID.argtypes = [
            ct.c_char_p,  # stand id
            ct.c_char_p,  # database control number
            ct.c_char_p,  # management id
            ct.POINTER(ct.c_int),  # length of stand id
            ct.POINTER(ct.c_int),  # length of control number
            ct.POINTER(ct.c_int),  # length of management id
        ]
        self._fvsStandID.restype = None

        if self.keyfile is None:
            msg = "Keyfile not loaded yet."
            raise AttributeError(msg)
        if self.stop_point_code is None:
            msg = "No inventory data loaded yet. Call `run` method."
            raise RuntimeError(msg)

        self._fvsStandID(
            self._stand_id,
            self._stand_cn,
            self._mgmt_id,
            ct.c_int(0),
            ct.c_int(0),
            ct.c_int(0),
        )

        return {
            STAND_ID_COLUMN_NAME: self._stand_id.value.decode().strip(),
            STAND_CN_COLUMN_NAME: self._stand_cn.value.decode().strip(),
            MGMT_ID_COLUMN_NAME: self._mgmt_id.value.decode().strip(),
        }

    @property
    def stop_point_code(self) -> int | None:
        """A code used to instruct FVS when to stop during a cycle.

        -1 : Stop at every stop location.
         0 : Never stop.
         1 : Stop just before the first call to the Event Monitor.
         2 : Stop just after the first call to the Event Monitor.
         3 : Stop just before the second call to the Event Monitor.
         4 : Stop just after the second call to the Event Monitor.
         5 : Stop after growth and mortality has been computed, but prior to
                applying them.
         6 : Stop just before the ESTAB routines are called.
         7 : Stop just after input is read but before missing values are imputed
                (tree heights and crown ratios, for example) and model
                calibration (argument stptyr is ignored).
        """
        if self._stop_point_code is not None:
            return self._stop_point_code.value
        return None

    @property
    def stop_point_year(self) -> int | None:
        """A code indicating which cycles FVS should stop at.

        0 : Never stop.
        1 : Stop at every cycle.
        YYYY : A specific year during the simulation period.
        """
        if self._stop_point_year is not None:
            return self._stop_point_year.value
        return None

    @fvs_property
    def summary(self) -> pd.DataFrame:
        """Return a dataframe with FVS Summary Statistics for all initiated cycles.

        The returned dataframe omits cycles that have not yet been initiated, which are
        identifiable where all values in that row are zero.
        """
        self._fvsSummary.argtypes = [
            np.ctypeslib.ndpointer(np.intc, flags=STR_C_CONTIGUOUS),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
        ]
        self._fvsSummary.restype = None

        dims = self.dims
        if dims[STR_NCYCLES] == 0:
            return None
        summary = np.zeros(
            dtype=np.intc,
            shape=(dims[STR_NCYCLES] + 1, len(SUMMARY_COLS)),
        )
        for i in range(dims[STR_NCYCLES] + 1):
            self._fvsSummary(
                summary[i],
                ct.c_int(i + 1),  # icycle
                ct.c_int(dims[STR_NCYCLES]),  # ncycles
                ct.c_int(0),  # maxrow
                ct.c_int(0),  # maxcol
                ct.c_int(0),  # rtncode
            )

        empty_years = (summary == 0).all(axis=1)
        return pd.DataFrame(
            summary[~empty_years, :], columns=SUMMARY_COLS
        ).copy()

    def load_keyfile(self, keywordfile: str | os.PathLike) -> None:
        """Sets the keywordfile as a command line argument to FVS.

        Args:
          keywordfile (str | os.PathLike): path to the FVS keyword file
        """
        if self.itrncd != FVS_ITRNCD_NOT_STARTED:
            if self.itrncd != FVS_ITRNCD_FINISHED_ALL_STANDS:
                msg = (
                    "FVS had not completed the previous simulation. "
                    "Outputs from that simulation may be incomplete."
                )
                warnings.warn(msg)
            logging.debug("FVS was already started. Resetting.")
            self._reload_fvs()

        self._fvsSetCmdLine.argtypes = [
            ct.c_char_p,
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
        ]
        self._fvsSetCmdLine.restype = None

        self.keyfile_path = Path(os.path.abspath(keywordfile))
        with open(self.keyfile_path) as f:
            self.keyfile = f.read()

        cmdline = f"--keywordfile={self.keyfile_path}"
        nch = len(cmdline)

        self._fvsSetCmdLine(cmdline.encode(), ct.c_int(nch), self._itrncd)
        logging.debug(f"Return code updated to {self.itrncd}")

        return

    def set_stop_point_codes(
        self,
        stop_point_code: int | None = None,
        stop_point_year: int | None = None,
    ) -> None:
        """Sets FVS stop point codes.

        Args:
            stop_point_code (int): Optional code for when FVS should stop during
                a cycle:
               -1 : Stop at every stop location
                0 : Never stop
                1 : Stop just before the first call to the Event Monitor
                2 : Stop just after the first call to the Event Monitor
                3 : Stop just before the second call to the Event Monitor
                4 : Stop just after the second call to the Event Monitor
                5 : Stop after growth and mortality has been computed, but
                        prior to applying them
                6 : Stop just before the ESTAB routines are called
                7 : Stop just after input is read but before missing values
                        are imputed
            stop_point_year (int): Optional, years FVS should stop, options are:
                0 : Never stop
               -1 : Stop at every cycle
               YYYY : A specific year during the simulation period
        """
        self._fvsSetStoppointCodes.argtypes = [
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int),
        ]
        self._fvsSetStoppointCodes.restype = None

        if stop_point_code is not None:
            if stop_point_code in range(-1, 8):
                self._stop_point_code = ct.c_int(stop_point_code)  # type: ignore[assignment]
            else:
                msg = "Invalid value for stop_point_code"
                raise ValueError(msg)
        elif self._stop_point_code is None:
            self._stop_point_code = ct.c_int(0)  # type: ignore[assignment]

        if stop_point_year is not None:
            if stop_point_code is not None:
                self._stop_point_year = ct.c_int(stop_point_year)  # type: ignore[assignment]
            else:
                msg = (
                    "Must specify stop_point_year if also specifying "
                    "stop_point_code"
                )
                raise ValueError(msg)
        elif self._stop_point_year is None:
            self._stop_point_year = ct.c_int(0)  # type: ignore[assignment]

        self._fvsSetStoppointCodes(self._stop_point_code, self._stop_point_year)

        return

    def get_species_attr(self, attr: str) -> npt.NDArray[np.float64]:
        """Gets a single attribute for all existing species.

        Args:
            attr (str): name of species attribute to fetch

        Returns:
            array with values of requested attribute for all trees

        """
        self._species_attr(attr, FvsAttributeAccessor.GET)
        return self._species_attrs[attr]

    def set_species_attr(self, attr: str, arr: npt.NDArray[np.float64]) -> None:
        """Sets a single attribute for all existing species.

        Args:
            attr (str): name of species attribute to fetch
            arr (npt.NDArray[np.float64]): array of values to set

        """
        return self._species_attr(attr, FvsAttributeAccessor.SET, arr)

    def get_tree_attr(self, attr: str) -> npt.NDArray[np.float64]:
        """Get a single attribute for all live tree records, in process.

        Reads per-tree state directly from the engine's memory. Valid names
        include 'id', 'species', 'dbh', 'dg', 'ht', 'htg', 'crwdth', 'cratio',
        'tpa', 'mort', 'plot', 'age'. The returned length follows the current
        live tree count (dims['ntrees']), which changes cycle to cycle.

        Args:
            attr (str): name of the tree attribute to fetch.

        Returns:
            array of length ntrees with the requested attribute.
        """
        return self._tree_attr(attr, FvsAttributeAccessor.GET)

    def set_tree_attr(self, attr: str, arr: npt.NDArray[np.float64]) -> None:
        """Set a single attribute for all live tree records, in process.

        This is the write side that enables injecting external growth
        predictions. Stop FVS at restart code 5 (after diameter growth and
        mortality are computed but before they are applied), overwrite 'dg'
        (and 'htg' for height growth) with predictions from another model such
        as the fvs-conus equations, then resume so the engine applies them.

        Args:
            attr (str): name of the tree attribute to set.
            arr (npt.NDArray[np.float64]): values to set, length ntrees.
        """
        self._tree_attr(attr, FvsAttributeAccessor.SET, arr)
        return

    def _tree_attr(
        self,
        attr: str,
        action: FvsAttributeAccessor,
        arr: npt.NDArray[np.float64] | None = None,
    ) -> npt.NDArray[np.float64] | None:
        """Get or set a single per-tree attribute via the FVS fvsTreeAttr API.

        Mirrors `_species_attr` but the array is sized to the live tree count,
        which is dynamic, so dims are read fresh on every call.
        """
        dims = self.dims
        ntrees = dims[STR_NTREES]
        maxtrees = dims[STR_MAXTREES]
        # Allocate to the maximum tree capacity (mirrors the species accessor, which
        # allocates to maxspecies): the engine may index the full array, so a buffer
        # sized only to the current live count can be written past, segfaulting.
        buf = np.zeros(dtype=np.float64, shape=(maxtrees,))
        if action == FvsAttributeAccessor.SET:
            if arr is None:
                msg = "Must provide `arr` if `action` is 'set'"
                raise TypeError(msg)
            src = np.ascontiguousarray(arr, dtype=np.float64)
            if src.shape != (ntrees,):
                msg = f"`arr` must have shape (ntrees,)=({ntrees},), got {src.shape}"
                raise ValueError(msg)
            buf[:ntrees] = src

        self._fvsTreeAttr.argtypes = [
            ct.POINTER(ct.c_char),  # attr name
            ct.POINTER(ct.c_int),  # number of characters in attr name
            ct.POINTER(ct.c_char),  # action ('get' or 'set')
            np.ctypeslib.ndpointer(
                np.float64, ndim=1, flags=STR_C_CONTIGUOUS
            ),  # per-tree array
            ct.POINTER(ct.c_int),  # return code
        ]
        self._fvsTreeAttr.restype = None

        rtncode = ct.c_int(0)
        self._fvsTreeAttr(
            ct.c_char_p(attr.encode()),
            ct.c_int(len(attr)),
            ct.c_char_p(action.encode()),
            buf,
            rtncode,
        )
        if rtncode.value != 0:
            if rtncode.value == 1:
                msg = f"{attr} not found among tree attributes"
                raise NameError(msg)
            raise RuntimeError(f"fvsTreeAttr returned code {rtncode.value} for '{attr}'")

        return buf[:ntrees].copy() if action == FvsAttributeAccessor.GET else None

    def tree_table(
        self, attrs: tuple[str, ...] = ("id", "species", "dbh", "ht", "dg", "htg", "crwdth")
    ) -> pd.DataFrame:
        """Read several per-tree attributes into a DataFrame (one row per live tree)."""
        return pd.DataFrame({a: self.get_tree_attr(a) for a in attrs}, copy=False)

    def add_trees(self, dbh, species, ht, cratio, plot, tpa) -> int:
        """Add tree records to the engine in memory, via the fvsAddTrees API.

        Loads a stand without a keyword-file database, which avoids the in-process DBS SQLite read
        that otherwise leaves the example-tree arrays unpopulated and segfaults extree. Call after the
        engine has read its (database-free) keyword file, then project. All six inputs are equal-length;
        species and plot are stored as doubles per the API but hold integer indices.

        Args:
            dbh: DBH per tree (inches). species: FVS species index. ht: total height (ft).
            cratio: crown ratio (0-100). plot: plot/point number. tpa: trees per acre (expansion).
        Returns:
            number of trees added.
        """
        n = len(dbh)
        a_dbh, a_sp, a_ht, a_cr, a_plot, a_tpa = (
            np.ascontiguousarray(x, dtype=np.float64) for x in (dbh, species, ht, cratio, plot, tpa)
        )
        nd = np.ctypeslib.ndpointer(np.float64, ndim=1, flags=STR_C_CONTIGUOUS)
        self._fvsAddTrees.argtypes = [nd, nd, nd, nd, nd, nd, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
        self._fvsAddTrees.restype = None
        ntrees = ct.c_int(n)
        rtncode = ct.c_int(0)
        self._fvsAddTrees(a_dbh, a_sp, a_ht, a_cr, a_plot, a_tpa, ct.byref(ntrees), ct.byref(rtncode))
        if rtncode.value != 0:
            raise RuntimeError(
                f"fvsAddTrees returned {rtncode.value} (1 = no room for the trees, or ntrees is zero)"
            )
        return n

    def apply_calibrated_config(self) -> dict[str, bool]:
        """Apply calibrated parameters from a JSON config file.

        Uses the config_version and config_dir set at initialization.
        This is called automatically by run() when config_version='calibrated',
        but can also be called manually after load_keyfile() for fine grained
        control.

        Returns:
            Dictionary indicating which parameter groups were applied
        """
        if self._config_version is None or self._config_version == "default":
            return {}

        try:
            # Import here to avoid circular dependency
            from config.config_loader import FvsConfigLoader

            loader = FvsConfigLoader(
                self.variant.lower(),
                version=self._config_version,
                config_dir=self._config_dir,
                custom_config=self._custom_config,
            )
            result = loader.apply_to_fvs(self)
            self._config_applied = True
            logging.info(
                f"Applied {self._config_version} config for variant {self.variant}: {result}"
            )
            return result

        except ImportError:
            logging.warning(
                "config_loader not found. Install fvs-modern or add config/ to PYTHONPATH."
            )
            return {}
        except FileNotFoundError as e:
            logging.warning(f"Config file not found: {e}")
            return {}
        except Exception as e:
            logging.warning(f"Could not apply calibrated config: {e}")
            return {}

    def run(
        self,
        stop_point_code: int = 0,
        stop_point_year: int = 0,
    ) -> None:
        """Runs FVS.

        If config_version was set to 'calibrated' at initialization, the
        calibrated parameters will be applied automatically after reading
        inventory data (stop point 7) on the first call to run().

        Note that stopping after the simulation of each stand in a simulation is
        done even when no stop request has been scheduled (that is, FVS will
        return at the end of each stand in a simulation even if there are no
        stop codes specified). Once a stand has been fully processed by FVS, the
        FVS `restart_code` is set to 100 and the call to run() returns.

        If there are multiple stands in a single keyfile, the simulation of the
        next stand can be triggered by calling run() again.

        The main output text file may be truncated even after the last stand has
        been simulated. To conclude FVS writing to the main output file, call
        run() one last time. The `itrncd` attribute should then change to a
        value of 2, indicating all stands have been processed.

        Args:
            stop_point_code (optional, int): when FVS should stop during a cycle:
               -1 : Stop at every stop location
                0 : Never stop
                1 : Stop just before the first call to the Event Monitor
                2 : Stop just after the first call to the Event Monitor
                3 : Stop just before the second call to the Event Monitor
                4 : Stop just after the second call to the Event Monitor
                5 : Stop after growth and mortality has been computed, but
                        prior to applying them
                6 : Stop just before the ESTAB routines are called
                7 : Stop just after input is read but before missing values
                        are imputed
            stop_point_year (optional, int): years FVS should stop, options are:
                0 : Never stop
               -1 : Stop at every cycle
               YYYY : A specific year during the simulation period
        """
        self._fvs.argtypes = [ct.POINTER(ct.c_int)]
        self._fvs.restype = None

        if self.keyfile is None:
            msg = "No keyfile loaded yet."
            raise AttributeError(msg)
        logging.debug("Found keyfile.")

        # Apply calibrated or custom config on first run if requested
        if self._config_version in ("calibrated", "custom") and not self._config_applied:
            # First pass with stop point 7 to let FVS read input,
            # then apply calibrated parameters before imputation
            self.set_stop_point_codes(7, 0)
            while self.itrncd == 0:
                self._fvs(self._itrncd)
                if self.restart_code != 0:
                    break
            self.apply_calibrated_config()
            # Now continue with the user's requested stop points
            self.set_stop_point_codes(stop_point_code, stop_point_year)
        else:
            self.set_stop_point_codes(stop_point_code, stop_point_year)

        logging.debug(
            f"Set stop point codes, {stop_point_code}:{self.stop_point_code}, {stop_point_year}:{self.stop_point_year}"
        )
        while self.itrncd == 0:
            logging.debug("itrncd still zero.")
            self._fvs(self._itrncd)
            logging.debug(f"Ran _fvs routine, itrncd is {self.itrncd}")
            if self.restart_code != 0:
                logging.debug("restart code not zero... halting run.")
                break

        return

    def run_ensemble(
        self,
        stop_point_code: int = 0,
        stop_point_year: int = 0,
    ) -> list[pd.DataFrame]:
        """Run FVS multiple times with different posterior parameter draws.

        This method implements Monte Carlo uncertainty propagation by:
          1. Loading posterior draws from the calibration pipeline
          2. For each draw, reloading the keyfile, applying that draw's
             parameters, running FVS, and collecting the summary table
          3. Returning the list of summary DataFrames

        The ensemble can then be summarized into credible intervals using
        UncertaintyEngine.summarize_ensemble().

        Requires uncertainty=True and config_version='calibrated' or 'custom'
        to have been set at initialization.

        Args:
            stop_point_code: FVS stop point code (see run() docstring)
            stop_point_year: FVS stop point year

        Returns:
            List of summary DataFrames, one per posterior draw
        """
        if not self._uncertainty:
            raise RuntimeError(
                "Uncertainty mode not enabled. Set uncertainty=True at initialization."
            )

        if self._config_version not in ("calibrated", "custom"):
            raise RuntimeError(
                "Uncertainty estimation requires config_version='calibrated' or 'custom'."
            )

        if self.keyfile_path is None:
            raise AttributeError("No keyfile loaded yet. Call load_keyfile() first.")

        # Lazy import to avoid circular dependency at module load
        try:
            from config.uncertainty import UncertaintyEngine
            from config.config_loader import FvsConfigLoader
        except ImportError:
            raise ImportError(
                "config.uncertainty and config.config_loader must be importable. "
                "Add the fvs-modern root to PYTHONPATH."
            )

        # Initialize the uncertainty engine and config loader
        if self._uncertainty_engine is None:
            self._uncertainty_engine = UncertaintyEngine(
                self.variant.lower(),
                config_dir=self._config_dir,
                seed=self._seed,
            )

        loader = FvsConfigLoader(
            self.variant.lower(),
            version="default",
            config_dir=self._config_dir,
        )
        default_config = loader.config

        keyfile_path = self.keyfile_path
        self._ensemble_results = []
        n = min(self._n_draws, self._uncertainty_engine.n_draws)

        logging.info(f"Running uncertainty ensemble: {n} draws for variant {self.variant}")

        for i in range(n):
            draw_idx = self._uncertainty_engine.sample_draw_index()
            draw = self._uncertainty_engine.get_draw(draw_idx)

            # Reload FVS for a fresh simulation
            self._reload_fvs()
            self.load_keyfile(keyfile_path)

            # Advance to stop point 7 (after input read, before imputation)
            self.set_stop_point_codes(7, 0)
            while self.itrncd == 0:
                self._fvs(self._itrncd)
                if self.restart_code != 0:
                    break

            # Apply this draw's parameters
            self._uncertainty_engine.apply_draw_to_fvs(
                self, draw, default_config
            )

            # Run the full simulation
            self.set_stop_point_codes(stop_point_code, stop_point_year)
            while self.itrncd == 0:
                self._fvs(self._itrncd)
                if self.restart_code != 0:
                    break

            # Collect summary
            summ = self.summary
            if summ is not None:
                summ = summ.copy()
                summ["_draw"] = i
                summ["_draw_idx"] = draw_idx
                self._ensemble_results.append(summ)

            if (i + 1) % 10 == 0 or i == n - 1:
                logging.info(f"  Completed draw {i + 1}/{n}")

        logging.info(f"Ensemble complete: {len(self._ensemble_results)} successful runs")
        return self._ensemble_results

    @property
    def uncertainty_results(self) -> list[pd.DataFrame]:
        """Results from the most recent ensemble run.

        Returns:
            List of FVS summary DataFrames, one per posterior draw
        """
        return self._ensemble_results

    @property
    def uncertainty_summary(self) -> pd.DataFrame | None:
        """Summarized ensemble results with credible intervals.

        Convenience property that calls UncertaintyEngine.summarize_ensemble()
        on the stored ensemble results.

        Returns:
            DataFrame with mean, std, and quantile columns per variable,
            or None if no ensemble has been run yet
        """
        if not self._ensemble_results:
            return None

        try:
            from config.uncertainty import UncertaintyEngine
            return UncertaintyEngine.summarize_ensemble(self._ensemble_results)
        except ImportError:
            logging.warning("config.uncertainty not importable for summary")
            return None

    def _species_attr(
        self,
        attr: str,
        action: FvsAttributeAccessor,
        arr: npt.NDArray[np.float64] | None = None,
    ) -> None:
        """Gets or sets a single attribute for all existing species.

        Args:
            attr (str): name of species attribute to get or set
            action (FvsAttributeAccessor): 'get' or 'set'
            arr (optional, npt.NDArray[np.float64]): array of values to set
        """
        if attr not in self._species_attrs:
            msg = "Invalid variable requested. Valid options are"
            raise NameError(msg, self._species_attrs)

        dims = self.dims
        if (
            action == FvsAttributeAccessor.GET
            and self._species_attrs[attr] is None
        ):
            self._species_attrs[attr] = np.empty(
                dtype=np.float64, shape=(dims[STR_MAXSPECIES])
            )
        elif action == FvsAttributeAccessor.SET:
            if arr is None:
                msg = "Must provide `arr` if `action` is 'set'"
                raise TypeError(msg)
            if arr.shape != (dims[STR_MAXSPECIES],):
                msg = (
                    "`arr` must be same shape as `maxspecies` "
                    f"({dims[STR_MAXSPECIES]},)"
                )
                raise ValueError(msg)
            self._species_attrs[attr] = arr

        self._fvsSpeciesAttr.argtypes = [
            ct.POINTER(ct.c_char),  # attr name
            ct.POINTER(ct.c_int),  # number of characters in attr name
            ct.POINTER(ct.c_char),  # action (set or get)
            np.ctypeslib.ndpointer(
                np.float64, shape=(dims[STR_MAXSPECIES]), flags=STR_C_CONTIGUOUS
            ),  # array to fill
            ct.POINTER(ct.c_int),  # return code
        ]
        self._fvsSpeciesAttr.restype = None

        rtncode = ct.c_int(0)
        self._fvsSpeciesAttr(
            ct.c_char_p(attr.encode()),  # attribute requested
            ct.c_int(len(attr)),  # number of characters in attr name
            ct.c_char_p(action.encode()),  # "set" or "get"
            self._species_attrs[attr],  # array to fill
            rtncode,  # return code of setting/getting operation
        )
        ERRS = {
            0: "OK",
            1: "name not found",
            4: "length of name string was too large or small",
        }

        if rtncode.value != 0:
            if rtncode.value == 1:
                msg = f"{attr} not found among species attributes"
                raise NameError(msg)
            raise RuntimeError(ERRS[rtncode.value])

    def _reload_fvs(self) -> None:
        self._unload_fvs()
        self._load_fvs()
        self._initialize_attributes()
        return
