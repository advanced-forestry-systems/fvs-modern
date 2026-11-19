STR_C_CONTIGUOUS = "C_CONTIGUOUS"
STR_MAXCYCLES = "maxcycles"
STR_MAXPLOTS = "maxplots"
STR_MAXSPECIES = "maxspecies"
STR_MAXTREES = "maxtrees"
STR_NCYCLES = "ncycles"
STR_NPLOTS = "nplots"
STR_NTREES = "ntrees"

FVS_ITRNCD_NOT_STARTED = -1
FVS_ITRNCD_GOOD_RUNNING_STATE = 0
FVS_ITRNCD_FINISHED_ALL_STANDS = 2

STAND_CN_COLUMN_NAME = "stand_cn"
STAND_ID_COLUMN_NAME = "stand_id"
MGMT_ID_COLUMN_NAME = "mgmt_id"

NEEDED_ROUTINES = (
    "fvs",
    "fvsAddActivity",
    "fvsAddTrees",
    "fvsDimSizes",
    "fvsEvmonAttr",
    "fvsFFEAttrs",
    "fvsGetRestartCode",
    "fvsGetRtnCode",
    "fvsGetICCode",
    "fvsSVSDimSizes",
    "fvsSetStoppointCodes",
    "fvsSetCmdLine",
    "fvsSVSObjData",
    "fvsSpeciesAttr",
    "fvsSpeciesCode",
    "fvsStandID",
    "fvsSummary",
    "fvsTreeAttr",
    "fvsUnitConversion",
)

SPECIES_INDEX_COLUMN_NAME = "fvs_index"
SPECIES_ALPHA_COLUMN_NAME = "fvs_alpha"
SPECIES_FIA_COLUMN_NAME = "fia"
SPECIES_PLANTS_COLUMN_NAME = "plants"
SPECIES_COLUMN_NAMES = [
    SPECIES_INDEX_COLUMN_NAME,
    SPECIES_ALPHA_COLUMN_NAME,
    SPECIES_FIA_COLUMN_NAME,
    SPECIES_PLANTS_COLUMN_NAME,
]
SPECIES_ATTRS = (
    "spccf",
    "spsdi",
    "spsiteindx",
    "bfmind",
    "bftopd",
    "bfstmp",
    "frmcls",
    "bfmeth",
    "mcmind",
    "mctopd",
    "mcstmp",
    "mcmeth",
    "baimult",
    "htgmult",
    "mortmult",
    "mortdia1",
    "mortdia2",
    "regdmult",
    "reghmult",
)

SUMMARY_COLS = (
    "year",
    "age",
    "tpa",
    "tcuft",
    "mcuft",
    "bdft",
    "rtpa",
    "rtcuft",
    "rmcuft",
    "rbdft",
    "atba",
    "atccf",
    "attopht",
    "prdlen",
    "acc",
    "mort",
    "sampwt",
    "fortyp",
    "sizecls",
    "stkcls",
    # Columns 21 and 22 of the engine's summary(22) (apisubs.f90 fvsSummary; IOSUM(21) and IOSUM(22)
    # in vbase/disply.f90): sawlog cubic volume and sawlog cubic volume removed. Before October 30, 2026
    # fvs2py allocated 20 integers per row, so these were dropped and the last row's two extra integers
    # were written past the end of the buffer.
    "scuft",
    "rscuft",
)
