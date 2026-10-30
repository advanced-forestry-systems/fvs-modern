# Quickstart

This page takes the shortest path from a fresh checkout to a 100 year projection of one stand with the Northeast (NE) variant and prints the stand summary table. It assumes `gfortran`, `gcc`, Python 3.11 or later, `numpy` and `pandas` are installed ([Installation](installation.md)).

## 1. Build the NE library

From the repository root:

```bash
bash deployment/scripts/build_fvs_libraries.sh src-converted lib ne
ls lib/FVSne.so
```

## 2. Put fvs2py and config on the Python path

```bash
export PYTHONPATH="$PWD/deployment/fvs2py:$PWD"
mkdir -p runs
```

## 3. Write a keyword file

The keyword file below describes one stand with six tree records read inline (`TREEDATA` with 15 in field 1). It follows the inventory-mode layout used by `calibration/python/inventory_keyfile.py`: a `DESIGN` record that makes each record's expansion factor equal to trees per acre on a single 1 acre fixed plot, `STDINFO` with the NE forest code 919, stand age 45, aspect 0, slope 10 percent and elevation 800 ft (entered in hundreds of feet), and a `SITECODE` record giving balsam fir a site index of 60 ft. `TIMEINT` sets every cycle to 10 years and `NUMCYCLE` asks for 10 cycles, so the run covers 100 years.

Keyword records are fixed width: columns 1 to 10 hold the keyword and each following 10 column block holds one field. The tree records use the engine's default tree format: plot (columns 1 to 4), tree number (5 to 7), trees per acre (8 to 13), history code 1 for live (14), species code (15 to 17), DBH in inches (18 to 21), total height in feet (25 to 27) and a crown ratio class from 1 to 9 (35). Save the text exactly as shown, spaces included, as `runs/quickstart.key`:

```text
STDIDENT
QS001     quickstart stand
DESIGN          -1.0       1.0     999.0       1.0       0.0
STDINFO          919                  45       0.0      10.0       8.0
SITECODE          BF      60.0       1.0
INVYEAR         2025
TIMEINT            0        10
NUMCYCLE          10
TREEDATA          15
   1  1 120.01BF  5.2    38       6
   1  2  80.01RS  7.8    52       5
   1  3  60.01RM  6.5    45       4
   1  4  40.01YB  9.1    58       5
   1  5  30.01WP 12.4    70       5
   1  6  20.01BF  3.1    24       6
-999
PROCESS
STOP
```

The species codes are NE alpha codes (`BF` balsam fir, *Abies balsamea* (L.) Mill.; `RS` red spruce, *Picea rubens* Sarg.; `RM` red maple, *Acer rubrum* L.; `YB` yellow birch, *Betula alleghaniensis* Britton; `WP` eastern white pine, *Pinus strobus* L.).

## 4. Run it

Save this as `runs/quickstart.py`:

```python
from fvs2py import FVS

fvs = FVS("lib/FVSne.so")
fvs.load_keyfile("runs/quickstart.key")
fvs.run()

print("exit code:", fvs.exit_code)
cols = ["year", "age", "tpa", "atba", "atccf", "attopht", "prdlen"]
print(fvs.summary[cols].to_string(index=False))
```

and run it from the repository root:

```bash
python3 runs/quickstart.py
```

FVS writes its main output listing next to the keyword file (`runs/quickstart.out`). If the process ends with a segmentation fault, see [Troubleshooting](../user-guide/troubleshooting.md#crashes-segmentation-faults).

## 5. Read the output

The run used to write this page printed:

```text
exit code: 1
 year  age  tpa  atba  atccf  attopht  prdlen
 2025   45  350   102    191       67      10
 2035   55  341   127    217       74      10
 2045   65  336   154    243       80      10
 2055   75  313   172    253       84      10
 2065   85  293   189    262       87      10
 2075   95  276   205    270       89      10
 2085  105  262   221    277       91      10
 2095  115  249   237    286       92      10
 2105  125  229   242    283       93      10
 2115  135  208   242    276       95      10
 2125  145  191   242    269       97       0
```

There is one row per cycle boundary, 11 rows for 10 cycles. The values are integers in the engine's imperial units: `tpa` is trees per acre, `atba` is basal area after treatment in ft² ac⁻¹, `atccf` is crown competition factor, `attopht` is top height in feet and `prdlen` is the length of the period that starts in that year. [Outputs](../user-guide/outputs.md) lists every column.

Exit code 1 means the run finished with informational notes. In this stand the note is `FVS43 WARNING: FIAVBC REQUESTED BUT ECODIVISION CODE NOT RECOGNIZED OR MISSING`, which NE raises because it switches FIA volume, biomass and carbon on by default and the keyword file carries no ecodivision; see [Troubleshooting](../user-guide/troubleshooting.md).

Basal area stops at 242 ft² ac⁻¹ from 2105 onward. That plateau is the native NE mortality clamp described in `KNOWN_ISSUES.md` (mortality is inflated until residual basal area is within 1 ft² ac⁻¹ of BAMAX), not a property of this stand.

!!! note "How these numbers were produced"
    The table was produced from a library built from a text-only snapshot of commit 9b5aa90 that did not include the NVEL `.for` volume sources or the C sources (SQLite, FOFEM, `apisubsc.c`). Volume columns came out as zero in that build and are omitted above. Growth and mortality do not call the volume library, so the columns shown should match a full build, but that was not checked against one.

## 6. Next steps

To run the configuration the catalog recommends for NE (default parameters with GOMPIT survival), select it through the catalog before constructing `FVS`:

```python
import os
from fvs2py import FVS
from fvs2py.catalog import select_config

choice = select_config("ne")
os.environ.update(choice.env)   # sets FVS_GOMPIT and FVS_GOMPIT_COEF
fvs = FVS("lib/FVSne.so", config_version=choice.config_version)
fvs.load_keyfile("runs/quickstart.key")
fvs.run()
print(fvs.summary[["year", "tpa", "atba"]].to_string(index=False))
```

In the documentation build this lowered trees per acre in 2125 from 191 to 162 and basal area from 242 to 226 ft² ac⁻¹. Run it in a fresh Python process: the engine reads the GOMPIT switches once per process. [Choosing a configuration](../user-guide/choosing-a-configuration.md) explains the catalog, and [Running projections](../user-guide/running-projections.md) covers the rest of the API.
