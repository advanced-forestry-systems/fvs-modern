# FVS diagnostic guards

These standalone Python modules preserve the August 2026 Northeast driver audit. They are optional checks; importing them does not change FVS projections.

- `fvs_input_guards.py` classifies FIA species codes against the NE species table and builds collision-free stand and tree identifiers. An unknown code such as 9999 is rejected; a real code outside NE is reported before FVS silently maps it to other hardwood. Its source tables are checked against `config/ne.json` by the tests.
- `fvs_run_guards.py` checks `NUMCYCLE` against the compiled `MAXCYC`, scans the listing for FVS errors, and verifies the process return code and expected summary row count. FVS normally returns code 10. A rejected `NUMCYCLE` can otherwise leave a valid-looking two-row database after only one cycle.

Run `python3 -m unittest discover -s tools/tests` from the repository root. The optional collision reproduction on a private plot list is available with `FVS_GUARD_PLOT_LIST=/private/path/list.csv python3 tools/fvs_input_guards.py`; plot identifiers and lookup outputs stay outside this repository.

Example before and after one projection:

```sh
python3 tools/fvs_run_guards.py --key stand.key
python3 tools/fvs_run_guards.py --out stand.out --rc 10 --expected-cycles 10 --db fvsout.db
```

A calling driver must use these checks explicitly. The modules do not replace the driver's own input validation or prove the growth equations are correct.
