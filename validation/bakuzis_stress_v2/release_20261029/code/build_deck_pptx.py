#!/usr/bin/env python3
"""build_deck_pptx.py: the catalog presentation in the CRSF house look (dark green, gold accent, Calibri).
Usage: python3 build_deck_pptx.py --fig fig --sc <scorecard figure dir> --out deck.pptx
Every number on a slide traces to data/catalog_gates_by_arm.csv, data/stress_by_arm.csv, data/stress_by_variant.csv
or the memo Sections 14 to 16 (validation/bakuzis_stress_v2/).
"""
import argparse, os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

GREEN = RGBColor(0x1a, 0x3d, 0x28); GOLD = RGBColor(0xc5, 0xa5, 0x5a); INK = RGBColor(0x22, 0x22, 0x22); GREY = RGBColor(0x6b, 0x70, 0x6a); WHITE = RGBColor(255, 255, 255)
W, H = Inches(13.333), Inches(7.5)

def rect(sl, x, y, w, h, color):
    s = sl.shapes.add_shape(MSO_SHAPE.RECTANGLE, x, y, w, h); s.fill.solid(); s.fill.fore_color.rgb = color; s.line.fill.background(); s.shadow.inherit = False; return s

def text(sl, x, y, w, h, s, size=18, bold=False, color=INK, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, font="Calibri"):
    tb = sl.shapes.add_textbox(x, y, w, h); tf = tb.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.05); tf.margin_top = tf.margin_bottom = Inches(0.03)
    lines = s if isinstance(s, list) else [s]
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); p.alignment = align
        r = p.add_run(); r.text = ln; r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color; r.font.name = font
        p.space_after = Pt(6)
    return tb

def bullets(sl, x, y, w, h, items, size=17):
    tb = sl.shapes.add_textbox(x, y, w, h); tf = tb.text_frame; tf.word_wrap = True
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        r = p.add_run(); r.text = "▪  " + it; r.font.size = Pt(size); r.font.color.rgb = INK; r.font.name = "Calibri"; p.space_after = Pt(8)
    return tb

def footer(sl, n, note):
    rect(sl, 0, H - Inches(0.32), W, Inches(0.32), GREEN)
    text(sl, Inches(0.4), H - Inches(0.31), Inches(10), Inches(0.3), note, 10, color=WHITE, anchor=MSO_ANCHOR.MIDDLE)
    text(sl, W - Inches(1.2), H - Inches(0.31), Inches(0.9), Inches(0.3), str(n), 10, color=GOLD, align=PP_ALIGN.RIGHT, anchor=MSO_ANCHOR.MIDDLE)

def header(sl, title):
    rect(sl, 0, 0, W, Inches(1.0), GREEN); rect(sl, 0, Inches(1.0), W, Inches(0.06), GOLD)
    text(sl, Inches(0.4), Inches(0.12), W - Inches(0.8), Inches(0.8), title, 26, True, WHITE, anchor=MSO_ANCHOR.MIDDLE)

def notes(sl, s):
    sl.notes_slide.notes_text_frame.text = s

def table(sl, x, y, w, rows, colw=None, size=11):
    nr, nc = len(rows), len(rows[0]); shp = sl.shapes.add_table(nr, nc, x, y, w, Inches(0.3) * nr); t = shp.table
    tblPr = shp._element.graphic.graphicData.tbl.tblPr; tblPr.set("bandRow", "0"); tblPr.set("firstRow", "0")
    if colw:
        for j, cw in enumerate(colw): t.columns[j].width = cw
    for i in range(nr):
        for j in range(nc):
            c = t.cell(i, j); c.fill.solid(); c.fill.fore_color.rgb = WHITE if i else GREEN
            c.margin_top = c.margin_bottom = Inches(0.02); c.vertical_anchor = MSO_ANCHOR.MIDDLE
            p = c.text_frame.paragraphs[0]; p.alignment = PP_ALIGN.CENTER; r = p.add_run(); r.text = rows[i][j]
            r.font.size = Pt(size); r.font.name = "Calibri"; r.font.bold = (i == 0); r.font.color.rgb = WHITE if i == 0 else INK
    return t

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fig", default="fig"); ap.add_argument("--sc", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    prs = Presentation(); prs.slide_width = W; prs.slide_height = H; blank = prs.slide_layouts[6]
    src = "Source: fvs-modern main 1a520f5e, validation/bakuzis_stress_v2/ and config/catalog/ (October 29, 2026)"
    n = 0
    # 1 title
    sl = prs.slides.add_slide(blank); n += 1
    rect(sl, 0, 0, W, H, GREEN); rect(sl, Inches(0.6), Inches(3.55), Inches(2.2), Inches(0.08), GOLD)
    text(sl, Inches(0.6), Inches(1.4), Inches(12), Inches(2.1), "A configuration catalog for fvs-modern across the FVS variants of the conterminous United States", 36, True, WHITE)
    text(sl, Inches(0.6), Inches(3.8), Inches(12), Inches(1.2), ["Bakuzis admissibility gates, a stress test of every verdict, and a selector", "Aaron R. Weiskittel, Center for Research on Sustainable Forests, University of Maine"], 20, color=WHITE)
    text(sl, Inches(0.6), Inches(6.3), Inches(12), Inches(0.5), "October 29, 2026", 16, color=GOLD)
    notes(sl, "Deck accompanies the CRSF technical report of the same date and the Zenodo record. Method: modified Bakuzis matrix after Leary (1997), Ecological Modelling 98:35-46.")
    # 2 question
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "The question a user asks")
    bullets(sl, Inches(0.5), Inches(1.4), Inches(7.6), Inches(5), [
        "fvs-modern runs every FVS variant with native mortality or GOMPIT crown closure survival, and with default or Bayesian calibrated MAP parameters: four arms per variant.",
        "The Bakuzis work so far tested NE and ACD on a synthetic 36 cell design and the NE Paper 2 cohort. Nothing said which arm to use in the other twenty variants.",
        "On October 28 a six argument GOMPSURV interface was found updated in one mortality routine and not the other 21: every GOMPIT result outside NE, CS and LS since July 25 read undefined memory. Fixed, 22 libraries rebuilt.",
        "So: can a user pick survival model and parameter set for their variant on evidence? The catalog is that evidence."], 16)
    rect(sl, Inches(8.5), Inches(1.5), Inches(4.3), Inches(4.6), RGBColor(0xf3, 0xf1, 0xea))
    text(sl, Inches(8.7), Inches(1.7), Inches(3.9), Inches(4.3), ["Admissible means:", "stays inside the FIA fitting envelope for 50 and 100 years,", "self-thins the way Reineke says,", "does not pass the ceiling by killing the stand,", "and is a real alternative to default."], 16, color=GREEN)
    footer(sl, n, src)
    # 3 what was run
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "What was run")
    bullets(sl, Inches(0.5), Inches(1.4), Inches(12.3), Inches(5), [
        "Stand sample: up to 300 FIA plots per variant from the fvs-conus pairs frame (17 variants; ACD takes NE plots in ME, NH and VT, an assumption). Trees expanded by CONDPROP_UNADJ, species outside the variant relabelled to its other softwood or hardwood code.",
        "Year zero gate: engine basal area within 20 percent of FIA BALIVE. Passes 95 to 100 percent of stands in 15 variants, 89 percent LS, 78 CS, 54 SN; 4,485 stands carry the gates.",
        "Four engine arms x 100 years at a 5 year step, no regeneration (NOAUTOES in EM and IE), no posterior draws. Cardinal arrays 15976824 and 15984133 after library build 15967530 with the GOMPSURV fix.",
        "Thresholds: plot level FIA percentiles by variant (BA p50, p99, p99.9; SDI p99.9). Age censored, so the 50 year horizon is primary.",
        "Projector arms (fvs-conus b1, b2 by native and GOMPIT) listed for comparison from the October 27 rescore, not on the common sample."], 16)
    footer(sl, n, src)
    # 4 rule
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "The pre-registered rule (fixed before any run was scored)")
    table(sl, Inches(0.5), Inches(1.4), Inches(12.3), [
        ["Gate", "Test on seed passing stands", "Why it is there"],
        ["Magnitude, BA", "< 5 percent of stands above FIA BA p99.9 by year 50, Wilson 95 percent upper bound < 10 percent", "The envelope the engine was fit to"],
        ["Magnitude, SDI", "Same for the SDI p99.9", "A second axis of the same envelope"],
        ["Reineke", "Flag rate among trajectories above relative density 0.55 under 10 percent, on at least 10 evaluated; fewer is NOT TESTED and inadmissible", "Self-thinning is the strong law"],
        ["Stocking", "Median year 50 BA at or above FIA p50", "Catches an arm that passes the ceilings by killing the stand"],
        ["Distinct", "Calibrated trajectories differ from default", "Caught ACD and IE, whose calibrated block never loaded"],
        ["Choice", "Simplest admissible arm (default before calibrated, native before GOMPIT); NE prefers GOMPIT on the Paper 2 near term cohort", "Simplicity is the only tie break outside NE"]],
        colw=[Inches(1.8), Inches(6.5), Inches(4.0)], size=12)
    text(sl, Inches(0.5), Inches(5.6), Inches(12.3), Inches(1.2), "Alternative Reineke gate (Wilson upper bound < 25 percent on at least 30 evaluated) is carried beside the primary rule in the gate table and flags eight rule sensitive variants. The catalog ships under the primary rule (decision of October 29, 2026) because the alternative admits western arms by declining to test their self-thinning.", 13, color=GREY)
    footer(sl, n, src)
    # 5 red team
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "What the blind red team changed before the catalog shipped")
    bullets(sl, Inches(0.5), Inches(1.4), Inches(12.3), Inches(5), [
        "Plan review, seven confirmed defects folded in: one common stand sample per variant, a 50 year primary horizon with Wilson bounds, an under stocking lower bound, one shared variant assignment, projector GOMPIT confirmed R side.",
        "Output review of the first catalog (PR 120), four defects that each made the gates easier than the README claimed; regenerated rather than patched:",
        "    the under stocking bound used the 100 year peak, which includes year 0 and could not fail; now the median year 50 BA against the p50",
        "    a Reineke gate with under 10 evaluated trajectories passed as NA; now NOT TESTED (this removed TT default GOMPIT, which had been RECOMMENDED on 5 trajectories)",
        "    stands failing the year 0 gate were scored; now excluded, with the Reineke counts on the same stands",
        "    ACD and IE calibrated arms equalled default: the keyword loader raised on their per species BAMAX array; loader fixed, both rerun, a distinctness check added everywhere",
        "The same review explained the seeding failures carried since Section 14: partial conditions expanded per plot acre (the exact 0.75), and an all plot fallback in the runner. Both gone."], 15)
    footer(sl, n, src)
    # 6 results
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "Seven variants have an admissible fvs-modern arm, eleven do not")
    sl.shapes.add_picture(os.path.join(a.fig, "fig_variant_stability.png"), Inches(0.4), Inches(1.25), width=Inches(8.4))
    bullets(sl, Inches(9.0), Inches(1.3), Inches(4.0), Inches(5.5), [
        "RECOMMENDED: ACD, BM, CS, NE, SN (all four arms), LS (calibrated arms), WS (calibrated native).",
        "NOT RECOMMENDED: CA, CI, CR, EC, EM, IE, PN, SO, TT, UT, WC. In eight of them every arm passes magnitude and stocking and fails on Reineke alone.",
        "Magnitude passes almost everywhere: BA p99.9 exceedance at or below 1.7 percent in 14 of 18 variants, SDI never above 4.7 percent, lower bound passes in every arm.",
        "AK and NC: no fitting plots, NOT EVALUATED."], 13)
    text(sl, Inches(0.4), Inches(4.95), Inches(8.4), Inches(0.5), "Bars: probability of at least one admissible engine arm over 500 stand bootstraps. Points: share of 144 alternative rule cells. Data: stress_by_variant.csv.", 11, color=GREY)
    footer(sl, n, src)
    # 7 gate space
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "The two gates that decide: 72 engine arms in one plane")
    sl.shapes.add_picture(os.path.join(a.fig, "fig_gate_space.png"), Inches(0.5), Inches(1.2), height=Inches(5.8))
    bullets(sl, Inches(7.4), Inches(1.4), Inches(5.6), Inches(5.5), [
        "Green symbols are admissible. Nearly every arm sits under the 5 percent magnitude line; the Reineke axis sorts them.",
        "Exceptions on magnitude: EM native (22 and 10 percent) and EM default GOMPIT (8), LS default (14 and 10), WC native (11 and 9), PN native (34 and 32 on 38 stands).",
        "Reineke failures come in two forms: default native arms that do not self-thin (CA, TT, UT, WC, EM, CR: median slopes -0.83 to -1.46, flag rates 27 to 84 percent), and calibrated GOMPIT arms whose median slope is right (-1.59 to -1.90) but dispersed (17 to 41 percent flagged).",
        "Western samples are thin on the self-thinning limb: 18 to 202 density gated trajectories per variant; SO rests on 17 to 18, TT on 7 to 31."], 13)
    footer(sl, n, src + "; catalog_gates_by_arm.csv")
    # 8 self-thinning figure
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "Self-thinning limb by variant, default native arm")
    sl.shapes.add_picture(os.path.join(a.sc, "scorecard_20261029_default_native_selfthinning.png"), Inches(0.9), Inches(1.15), height=Inches(5.9))
    text(sl, Inches(8.7), Inches(1.4), Inches(4.3), Inches(5), ["Dotted line: Reineke slope -1.605 through the maximum SDI reached in the arm.", "", "ACD, BM, CS, NE, SN, WS run parallel to it. CA, TT, UT, WC and EM run flat: stems hold while diameter grows, which is a mortality model that does not respond to density.", "", "The per variant scorecards (validation/bakuzis_stress_v2/catalog/scorecards/) carry this panel for all four arms with the component attribution."], 14, color=INK)
    footer(sl, n, "Source: bakuzis_scorecards.R, Cardinal job 15986836")
    # 9 Bakuzis matrix NE
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "The full matrix for one recommended arm: NE, default GOMPIT")
    sl.shapes.add_picture(os.path.join(a.sc, "scorecard_20261029_default_gompit_NE_bakuzis_matrix.png"), Inches(0.9), Inches(1.15), height=Inches(5.9))
    text(sl, Inches(7.7), Inches(1.4), Inches(5.3), Inches(5), ["Nine relations, lower triangle. Site and density classes are within variant terciles of SICOND and BALIVE on real FIA stands, so the site based relations (site ordering, Eichhorn, Sukachev, yield density) are confounded by stand selection and are reported, not gated.", "", "The engine verdicts rest on Reineke (T3, T3b), the magnitude envelope (T4), plateau (T8), inversions (T9) and the initial step shock (T0).", "", "A pass is a failure to falsify, not a confirmation of accuracy."], 14)
    footer(sl, n, "Source: bakuzis_scorecards.R, Cardinal job 15986836")
    # 10 stress design
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "Stress test after shipping: could the analyst choices have gone another way?")
    table(sl, Inches(0.5), Inches(1.4), Inches(12.3), [
        ["Axis", "Shipped", "Alternatives in the 144 cell grid"],
        ["Horizon", "50 years", "100 years"], ["Exceedance tolerance", "5 percent", "3, 10 percent"], ["Wilson upper bound", "10 percent", "15 percent"],
        ["Basal area ceiling", "p99.9", "p99"], ["Reineke flag tolerance", "10 percent", "15, 25 percent"], ["Minimum evaluated trajectories", "10", "30"]],
        colw=[Inches(4.0), Inches(3.0), Inches(5.3)], size=13)
    bullets(sl, Inches(0.5), Inches(4.4), Inches(12.3), Inches(2.5), [
        "Part B: 500 stand bootstraps per variant under the shipped rule, giving each arm's probability of admissibility and each variant's probability of keeping an admissible engine arm.",
        "No new engine runs; recomputed from the trajectory level data (stress_catalog.R, Cardinal job 15990793).",
        "The recomputed primary rule reproduced all 72 shipped verdicts, an independent check of catalog_gates.R."], 15)
    footer(sl, n, src)
    # 11 stress results
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "No verdict overturned; 8 of 23 admissible and 8 of 49 not admissible arms fragile")
    sl.shapes.add_picture(os.path.join(a.fig, "fig_arm_stability.png"), Inches(0.4), Inches(1.2), height=Inches(5.5))
    table(sl, Inches(6.8), Inches(1.3), Inches(6.2), [
        ["Fragile admissible arm", "Bootstrap", "Grid share"],
        ["ACD calibrated GOMPIT", "0.58", "0.33"], ["ACD default GOMPIT", "0.58", "0.33"], ["CS default GOMPIT", "0.77", "0.50"], ["LS calibrated native", "0.84", "0.25"],
        ["NE calibrated GOMPIT", "0.57", "0.58"], ["NE calibrated native", "0.74", "0.54"], ["NE default GOMPIT", "0.84", "0.33"], ["WS calibrated native", "0.90", "0.50"]],
        colw=[Inches(3.2), Inches(1.5), Inches(1.5)], size=12)
    text(sl, Inches(6.8), Inches(4.3), Inches(6.2), Inches(2.5), ["Recommended variants keep an admissible arm in 100 percent of resamples (ACD, BM, CS, LS, SN), 99 (NE), 90 (WS).", "Fragile not admissible: CI default native 0.51, SO four arms 0.19 to 0.46, CR calibrated native 0.24, TT calibrated GOMPIT 0.25, EC default native 0.16; all among the 13 arms a 25 percent Reineke tolerance admits, seven of eight below even odds.", "Reineke tolerance is the only axis that changes recommended variants; the p99 ceiling admits five arms anywhere; 100 years drops LS calibrated native only."], 12)
    footer(sl, n, "Source: stress_by_arm.csv, stress_by_variant.csv (job 15990793)")
    # 12 selector
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "Using it: fvs2py.catalog.select_config")
    rect(sl, Inches(0.5), Inches(1.4), Inches(6.3), Inches(4.0), RGBColor(0xf3, 0xf1, 0xea))
    text(sl, Inches(0.7), Inches(1.5), Inches(6.0), Inches(3.8), [
        "from fvs2py.catalog import select_config",
        "c = select_config('NE')",
        "c.arm            # 'default_gompit'",
        "c.config_version # passed to FVS()",
        "c.env            # {'FVS_GOMPIT': '1', ...}",
        "c.admissible     # four arms",
        "c.notes          # near term tie break",
        "",
        "select_config('SO', survival='native')",
        "# least bad engine arm, NOT RECOMMENDED,",
        "# rule_sensitive note attached"], 14, font="Consolas")
    bullets(sl, Inches(7.2), Inches(1.4), Inches(5.8), Inches(5), [
        "Reads config/catalog/catalog_recommendation_by_variant.csv and catalog_gates_by_arm.csv; regenerate only through the pipeline, never by hand.",
        "Returns a note when the caller forces an arm the catalog did not admit, when Reineke was not testable, when a calibrated arm is not distinct, or when the verdict is rule sensitive.",
        "No admissible arm: the least bad engine arm with the status attached, falling back from a projector arm where that is least bad overall (TT).",
        "AK and NC: NOT EVALUATED.",
        "Stability labels sit beside the catalog in catalog/stress/; read them before choosing a fragile admissible arm."], 14)
    footer(sl, n, src)
    # 13 open items
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "What the catalog does not cover")
    bullets(sl, Inches(0.5), Inches(1.4), Inches(12.3), Inches(5), [
        "Per variant near term accuracy. Only the NE Paper 2 cohort exists, so outside NE the tie break is simplicity alone and the catalog says so.",
        "Posterior uncertainty: calibrated MAP is a point estimate here; the Paper 1 ensemble gives the band for NE and ACD only.",
        "Regeneration: no arm establishes it, so a 100 year plateau or decline in stems is partly a design choice.",
        "Site based relations are confounded on the unstratified sample; the synthetic 36 cell matrix stays the clean test.",
        "263 stands (256 eastern) fail the year zero gate, nearly all entering at 1.3 to 1.8 times their FIA basal area with verified input; excluded, unexplained, open.",
        "Thresholds are the sample maximum in small variants; the projector b1 GOMPIT arm has no NE or SN rows.",
        "In NE the stable arm is default native (0.98); the shipped GOMPIT preference trades long horizon stability for near term evidence."], 15)
    footer(sl, n, src + "; implementation-notes.md")
    # 14 availability
    sl = prs.slides.add_slide(blank); n += 1; header(sl, "Where everything is")
    bullets(sl, Inches(0.5), Inches(1.4), Inches(12.3), Inches(4), [
        "GitHub: advanced-forestry-systems/fvs-modern, main at 1a520f5e (PRs 118 to 126). config/catalog/ holds the two tables the selector reads; validation/bakuzis_stress_v2/ holds the memo (Sections 1 to 16), the catalog folder, scorecards, stress test and run index with md5 lists.",
        "Zenodo: concept DOI 10.5281/zenodo.23777723 (version 1.0.0 at 10.5281/zenodo.23777724) carries the tables, scripts, scorecard cells, the technical report and this deck. No FIA plot identifiers or coordinates in any file.",
        "Scripts: conus_catalog_runs.py (engine runs), bakuzis_stress_v2.R (scorer), catalog_gates.R (gates), bakuzis_scorecards.R (per variant scorecards), stress_catalog.R (stress test), fvs2py/catalog.py (selector).",
        "Cardinal jobs: build 15967530, runs 15976824 and 15984133, finalize 15984134, scorecards 15986836, stress 15990793."], 15)
    rect(sl, Inches(0.5), Inches(5.3), Inches(12.3), Inches(1.0), RGBColor(0xf3, 0xf1, 0xea))
    text(sl, Inches(0.7), Inches(5.4), Inches(12), Inches(0.8), "Aaron R. Weiskittel, Center for Research on Sustainable Forests, University of Maine. aaron.weiskittel@maine.edu", 14, color=GREEN, anchor=MSO_ANCHOR.MIDDLE)
    footer(sl, n, src)
    prs.save(a.out); print("wrote", a.out, n, "slides")

if __name__ == "__main__":
    main()
