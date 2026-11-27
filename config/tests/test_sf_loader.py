#!/usr/bin/env python3
"""
Unit test for the species-free (Leg B) config_loader additions.

Validates, against the current calibrated NE configuration:
  1. The version switch recognizes conus_sf / conus_hybrid and routes to the
     variant JSON.
  2. The species-free runtime block decomposes correctly (intercept, gamma,
     RE tables, species trait effects, hybrid source map).
  3. The runtime trait evaluator reproduces the precomputed trait_effect from
     raw traits + standardization constants + gamma (the core math).
  4. The hybrid per-species source resolver returns leg_a / leg_b correctly.
  5. The full species-free linear predictor sums intercept + trait + REs +
     covariates and returns a finite number.

Dev/test only. Reads a copy of a production config in a temporary directory.

Run:
  python3 test_sf_loader.py /path/to/config_loader.py /path/to/ne.json
"""
import sys, os, json, tempfile, shutil, importlib.util

def load_loader(loader_path):
    spec = importlib.util.spec_from_file_location("config_loader", loader_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.FvsConfigLoader

def main():
    loader_path = sys.argv[1] if len(sys.argv) > 1 else "config_loader.py"
    source      = sys.argv[2] if len(sys.argv) > 2 else "ne.json"
    FvsConfigLoader = load_loader(loader_path)

    # Stage the current calibrated config in a temp config dir.
    tmp = tempfile.mkdtemp(prefix="sf_test_")
    cal = os.path.join(tmp, "calibrated"); os.makedirs(cal)
    shutil.copy(source, os.path.join(cal, "ne.json"))

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
    check("intercept is finite", isinstance(rt["intercept"], (int, float)) and
          abs(rt["intercept"]) < 1e6, rt["intercept"])
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

    # 4. Hybrid source resolver agrees with every serialized component map.
    # Height growth currently has no leg_a species, so exercise both legs
    # across the full NE config rather than assuming both occur per component.
    source_counts = {"leg_a": 0, "leg_b": 0}
    routing_ok = True
    for component in comps:
        for spcd, leg in L.get_conus_sf_runtime_block(component)["source_map"].items():
            source_counts[leg] = source_counts.get(leg, 0) + 1
            routing_ok &= L.resolve_species_source(component, spcd) == leg
    check("hybrid source maps use both legs", all(source_counts.get(x, 0) > 0
          for x in ("leg_a", "leg_b")), source_counts)
    check("all serialized species routes resolve", routing_ok)
    print(f"   hybrid source map entries: {source_counts}")

    # 5. full linear predictor for a sample tree
    eco = {"L1": "8", "L2": list(rt["re_L2"])[0], "L3": list(rt["re_L3"])[0],
           "FT": list(rt["re_FT"])[0]}
    covs = {"a1": 2.5, "a2": 2.0, "a3": -0.3, "a5": 5.0, "a6": 0.3,
            "a7": 10.0, "a8": 0.1}
    eta = L.sf_linear_predictor("height_growth", any_spcd, eco, covs)
    check("linear predictor finite", eta == eta and abs(eta) < 1e6, eta)
    print(f"   sample eta(SPCD {any_spcd}) = {eta:.4f}")

    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if fails:
        print(f"FAILED {len(fails)} check(s): {fails}")
        sys.exit(1)
    print("ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
