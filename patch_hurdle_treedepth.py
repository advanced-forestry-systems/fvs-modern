"""
Patch ingrowth hurdle and compos drivers to accept --max_treedepth and
--adapt_delta on CLI. Default behavior unchanged.
"""

import re

paths = [
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/35b_fit_ingrowth_hurdle.R",
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/36_fit_ingrowth_species_composition.R",
]

for p in paths:
    s = open(p).read()

    # 1) Find the SMOKE handling block and add CLI args after it
    smoke_marker = "SMOKE         <- has_flag(\"smoke\")" if "SMOKE         <-" in s else "SMOKE     <- has_flag(\"smoke\")"
    needle = smoke_marker
    if needle not in s:
        print(f"NEEDLE NOT FOUND in {p}: {needle!r}")
        continue
    inject = (
        needle
        + "\nMAX_TREEDEPTH <- as.integer(get_arg(\"max_treedepth\", \"10\"))"
        + "\nADAPT_DELTA   <- as.numeric(get_arg(\"adapt_delta\", \"0.9\"))"
        + "\nITER_WARMUP   <- as.integer(get_arg(\"iter_warmup\",   if (has_flag(\"smoke\")) \"50\"   else \"1000\"))"
        + "\nITER_SAMPLING <- as.integer(get_arg(\"iter_sampling\", if (has_flag(\"smoke\")) \"50\"   else \"1000\"))"
        + "\nCHAINS        <- as.integer(get_arg(\"chains\",        if (has_flag(\"smoke\")) \"2\"    else \"4\"))"
    )
    if "MAX_TREEDEPTH <- as.integer" in s:
        print(f"already patched: {p}")
    else:
        s = s.replace(needle, inject, 1)

    # 2) Replace the chains/iter assignment lines with our CLI-driven names
    old_chain = "chains <- if (SMOKE) 2 else 4"
    new_chain = "chains <- CHAINS"
    if old_chain in s:
        s = s.replace(old_chain, new_chain, 1)
    old_iw = "iw     <- if (SMOKE) 50 else 1000"
    new_iw = "iw     <- ITER_WARMUP"
    if old_iw in s:
        s = s.replace(old_iw, new_iw, 1)
    old_is = "is_    <- if (SMOKE) 50 else 1000"
    new_is = "is_    <- ITER_SAMPLING"
    if old_is in s:
        s = s.replace(old_is, new_is, 1)

    # 3) Replace the sample call to use CLI-driven max_treedepth/adapt_delta
    s = re.sub(
        r"adapt_delta = 0\.9, max_treedepth = 10",
        "adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH",
        s,
    )

    open(p, "w").write(s)
    print(f"patched: {p}")
