#!/usr/bin/env python3
"""
Unit test for the species-free (Leg B) config_loader additions.

Validates, against ne.sf_preview.json (the dry-run block from 62b):
  1. The version switch recognizes conus_sf / conus_hybrid and routes to the
     variant JSON.
  2. The species-free runtime block decomposes correctly (intercept, gamma,
     RE tables, species trait effects, hybrid source map).
  3. The runtime trait evaluator reproduces the precomputed trait_effect from
     raw traits + standardization constants + gamma (the core math).
  4. The hybrid per-species source resolver returns leg_a / leg_b correctly.
  5. The full species-free linear predictor sums intercept + trait + REs +
     covariates and returns a finite number.

Dev/test only. Reads the preview file, never the production configs.

Run:
  python3 test_sf_loader.py /path/to/config_loader.py /path/to/ne.sf_preview.json
"""
import sys, os, json, tempfile, shutil, importlib.util

def load_loader(loader_path):
    spec = importlib.util.spec_from_file_location("config_loader", loader_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.FvsConfigLoader

def main():
    loader_path = sys.argv[1] if len(sys.argv) > 1 else "config_loader.py"
    preview     = sys.argv[2] if len(sys.argv) > 2 else "ne.sf_preview.json"
    FvsConfigLoader = load_loader(loader_path)

    # Stage the preview as a normal variant JSON in a temp config dir.
    tmp = tempfile.mkdtemp(prefix="sf_test_")
    cal = os.path.join(tmp, "calibrated"); os.makedirs(cal)
    shutil.copy(preview, os.path.join(cal, "ne.json"))

    fails = []
    def check(name, cond, detail=""):
        print(("PASS" if cond else "FAIL"), "-", name, ("" if cond else f":: {detail}"))
        if not cond: fails.append(name)

    # 1. version switch
    L = FvsConfigLoader("ne", version="conus_sf", config_dir=tmp)
    check("conus_sf is a valid version", "conus_sf" in L.VALID_VERSIONS)
    check("conus_hybrid is a valid version", "conus_hybrid" in L.VALID_VERSIONS)
    check("config_path routes to calibrated/ne.json",
          str(L.config_path).endswith("calibrated/ne.json"), str(L.config_path))
    check("has_conus_sf_block", L.has_conus_sf_block)
    comps = L.conus_sf_components_present()
    check("height_growth present", "height_growth" in comps, str(comps))

    # 2. runtime block decomposition
    rt = L.get_conus_sf_runtime_block("height_growth")
    check("intercept name is a0", rt["intercept_name"] == "a0", rt["intercept_name"])
    check("intercept ~ -1.60", abs(rt["intercept"] - (-1.6019503)) < 1e-3, rt["intercept"])
    check("8 trait gammas", len(rt["gamma"]) == 8, len(rt["gamma"]))
    check("11 L1 REs", len(rt["re_L1"]) == 11, len(rt["re_L1"]))
    check("149 FT REs", len(rt["re_FT"]) == 149, len(rt["re_FT"]))
    check("covariate excludes intercept+sigmas",
          all(not p.startswith(("sigma", "phi")) and p != "a0"
              for p in rt["covariate"]), list(rt["covariate"]))

    # 3. evaluator math: trait_effect = sum_j std_trait[j] * gamma[j] should
    #    reproduce the precomputed trait_effect_mean (std columns are the
    #    imputed + standardized values 61b actually used; raw may carry NA).
    max_err = 0.0
    n_checked = 0
    for spcd, te_pre in list(rt["trait_effect"].items())[:30]:
        std = rt["std_traits"].get(spcd, {})
        eta = sum(std[c] * g for c, g in rt["gamma"].items() if c in std)
        max_err = max(max_err, abs(eta - te_pre))
        n_checked += 1
    check(f"evaluator reproduces precomputed trait_effect ({n_checked} spp, max err {max_err:.2e})",
          max_err < 1e-6, f"max_err={max_err}")
    # and via the public helper
    any_spcd = next(iter(rt["trait_effect"]))
    helper = FvsConfigLoader.sf_trait_effect(rt, any_spcd)
    check("sf_trait_effect helper matches precomputed",
          abs(helper - rt["trait_effect"][any_spcd]) < 1e-9)

    # 4. hybrid source resolver (316=red maple is in the 5 reliable NE species)
    src_316 = L.resolve_species_source("height_growth", 316)
    check("SPCD 316 routes leg_a", src_316 == "leg_a", src_316)
    leg_b_spp = [s for s, v in rt["source_map"].items() if v == "leg_b"]
    check("at least one leg_b species exists", len(leg_b_spp) > 0)
    check("a leg_b species resolves leg_b",
          L.resolve_species_source("height_growth", leg_b_spp[0]) == "leg_b")
    na = sum(1 for v in rt["source_map"].values() if v == "leg_a")
    nb = sum(1 for v in rt["source_map"].values() if v == "leg_b")
    print(f"   hybrid source map: leg_a={na}, leg_b={nb}")

    # 5. full linear predictor for a sample tree
    eco = {"L1": "8", "L2": list(rt["re_L2"])[0], "L3": list(rt["re_L3"])[0],
           "FT": list(rt["re_FT"])[0]}
    covs = {"a1": 2.5, "a2": 2.0, "a3": -0.3, "a5": 5.0, "a6": 0.3,
            "a7": 10.0, "a8": 0.1}
    eta = L.sf_linear_predictor("height_growth", 316, eco, covs)
    check("linear predictor finite", eta == eta and abs(eta) < 1e6, eta)
    print(f"   sample eta(SPCD 316) = {eta:.4f}")

    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if fails:
        print(f"FAILED {len(fails)} check(s): {fails}")
        sys.exit(1)
    print("ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
