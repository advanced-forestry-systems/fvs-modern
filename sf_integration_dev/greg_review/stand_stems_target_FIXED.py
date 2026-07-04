def stand_stems_target(bundle, n1, years, top_ht=None, rd=None, ln_qmd=None,
                       L1=None):
    """
    Surviving-stems N2 target from the fitted state-space stem-transition bundle
    (76_fit_stand_stems.R). The fitted model is a CLOGLOG EXPOSURE HAZARD, identical
    in scale to 71_fit_stand_survival.R:

        eta_h  = b0 + b_lnN1*ln(N1) + b_topht*top_ht + b_rd*rd + b_lnqmd*ln_qmd + z_L1
        H_stand = exp(eta_h)                 # stand log-hazard, intercept already on raw log(YEARS)
        S(T)    = exp(-H_stand * years)      # survival over the interval (exposure)
        N2      = N1 * S(T)                   # 0 < N2 <= N1 (survivors only; ingrowth separate)

    Time enters ONLY through the *years multiplication inside S(T); the fit folded
    the centered log(YEARS)-3.9 offset back into the intercept, so there is NO
    ln_years fixed effect. Missing covariate terms contribute 0.
    """
    if not bundle or float(years) <= 0 or float(n1) <= 0:
        return float(n1)
    fx = bundle.get("fixed_effects", {}) or {}

    def _m(name):
        v = fx.get(name)
        if isinstance(v, dict):
            return float(v.get("mean", 0.0))
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0

    eta = _m("Intercept") + _m("ln_n1") * math.log(float(n1))
    if top_ht is not None:
        eta += _m("top_ht") * float(top_ht)
    if rd is not None:
        eta += _m("rd") * float(rd)
    if ln_qmd is not None:
        eta += _m("ln_qmd") * float(ln_qmd)
    if L1 is not None:
        re = bundle.get("re_L1") or {}
        lut = {str(lv): float(m) for lv, m in zip(re.get("level", []),
                                                  re.get("mean", []))}
        eta += lut.get(str(L1), 0.0)
    # cloglog exposure hazard -> survival over the interval; guarantees 0 < N2 <= N1
    H = math.exp(min(eta, 30.0))
    return float(n1) * math.exp(-H * float(years))
