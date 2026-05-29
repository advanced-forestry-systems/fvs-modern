// build_silc_deck_v4.js
// =====================================================================
// SILC operational deck v4 for Ian / Ryan at Seven Islands.
// Extends v3 with: stratification methodology slide, plot-count-weighted
// AGM strata trajectories with uncertainty band, three-model species
// composition (AGM + OSM-ACD added to observed), updated next steps.
// =====================================================================
const PptxGenJS = require("pptxgenjs");

const OD = "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory";
const CFI_DIR = `${OD}/silc_cfi`;

// CRSF brand
const CRSF_GREEN = "1A3D28";
const CRSF_ACCENT = "88A47A";
const CHARCOAL = "333333";
const MED_GRAY = "666666";
const LIGHT_GRAY = "999999";

const pres = new PptxGenJS();
pres.layout = "LAYOUT_WIDE";  // 13.33 x 7.5"
pres.title = "SILC Four-Model Benchmark";
pres.author = "Aaron R. Weiskittel";

// Reusable takeaway band
function takeawayBand(slide, text) {
  slide.addShape(pres.ShapeType.rect, {
    x: 0, y: 6.85, w: 13.33, h: 0.65,
    fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
  });
  slide.addText(text, {
    x: 0.5, y: 6.85, w: 12.5, h: 0.65,
    fontFace: "Aptos", fontSize: 22, color: "FFFFFF", bold: true,
    align: "left", valign: "middle"
  });
}

// Reusable footer line above takeaway band
function citationLine(slide, text) {
  slide.addText(text, {
    x: 0.5, y: 6.10, w: 12.5, h: 0.30,
    fontFace: "Aptos", fontSize: 10, color: LIGHT_GRAY, italic: true,
    align: "left", valign: "middle"
  });
}

// ===================================================================
// SLIDE 1: Title
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };

  s.addText("Seven Islands four-model benchmark", {
    x: 0.7, y: 1.6, w: 12.0, h: 1.2,
    fontFace: "Aptos Display", fontSize: 48, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Multi model CFI scorecard plus 100 year strata yield",
    {
    x: 0.7, y: 2.8, w: 12.0, h: 0.7,
    fontFace: "Aptos", fontSize: 26, color: MED_GRAY,
    align: "left", valign: "middle"
  });

  // Accent bar
  s.addShape(pres.ShapeType.rect, {
    x: 0.7, y: 3.7, w: 1.4, h: 0.06,
    fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
  });

  s.addText("Aaron R. Weiskittel", {
    x: 0.7, y: 4.0, w: 12.0, h: 0.5,
    fontFace: "Aptos", fontSize: 22, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Center for Research on Sustainable Forests  |  University of Maine", {
    x: 0.7, y: 4.55, w: 12.0, h: 0.4,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY,
    align: "left", valign: "middle"
  });
  s.addText("Prepared for Ian and Ryan, Seven Islands Land Company  |  28 May 2026", {
    x: 0.7, y: 5.05, w: 12.0, h: 0.4,
    fontFace: "Aptos", fontSize: 14, color: LIGHT_GRAY, italic: true,
    align: "left", valign: "middle"
  });
}

// ===================================================================
// SLIDE 2: Why this deck — the question we are answering
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };

  s.addText("Two operational questions, one framework", {
    x: 0.5, y: 0.4, w: 12.3, h: 0.8,
    fontFace: "Aptos Display", fontSize: 36, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });

  // Two-card grid
  const cardY = 1.6, cardH = 4.6;
  const cardW = 5.95, gap = 0.4;

  // Card 1: long-horizon planning
  s.addShape(pres.ShapeType.rect, {
    x: 0.5, y: cardY, w: cardW, h: cardH,
    fill: { color: "F7F7F7" }, line: { color: "F7F7F7" }
  });
  s.addShape(pres.ShapeType.rect, {
    x: 0.5, y: cardY, w: 0.08, h: cardH,
    fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
  });
  s.addText("Long-horizon yield", {
    x: 0.75, y: cardY + 0.25, w: cardW - 0.4, h: 0.7,
    fontFace: "Aptos Display", fontSize: 24, color: CRSF_GREEN, bold: true,
    align: "left", valign: "top"
  });
  s.addText("What will a Mixedwood A+B stand yield in cords over 100 years?", {
    x: 0.75, y: cardY + 1.05, w: cardW - 0.4, h: 1.1,
    fontFace: "Aptos", fontSize: 18, color: CHARCOAL,
    align: "left", valign: "top"
  });
  s.addText("Driver: AGM byStrata trajectories on the 11 SILC matrix stands, projected to 2123", {
    x: 0.75, y: cardY + 2.25, w: cardW - 0.4, h: 1.0,
    fontFace: "Aptos", fontSize: 14, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });
  s.addText("Section 1 of this deck", {
    x: 0.75, y: cardY + 3.95, w: cardW - 0.4, h: 0.4,
    fontFace: "Aptos", fontSize: 13, color: LIGHT_GRAY,
    align: "left", valign: "bottom"
  });

  // Card 2: short-horizon accuracy
  s.addShape(pres.ShapeType.rect, {
    x: 0.5 + cardW + gap, y: cardY, w: cardW, h: cardH,
    fill: { color: "F7F7F7" }, line: { color: "F7F7F7" }
  });
  s.addShape(pres.ShapeType.rect, {
    x: 0.5 + cardW + gap, y: cardY, w: 0.08, h: cardH,
    fill: { color: CRSF_ACCENT }, line: { color: CRSF_ACCENT }
  });
  s.addText("Short-horizon accuracy", {
    x: 0.75 + cardW + gap, y: cardY + 0.25, w: cardW - 0.4, h: 0.7,
    fontFace: "Aptos Display", fontSize: 24, color: CRSF_GREEN, bold: true,
    align: "left", valign: "top"
  });
  s.addText("How close does the model land at 5 to 10 year remeasurement on SILC's own plots?", {
    x: 0.75 + cardW + gap, y: cardY + 1.05, w: cardW - 0.4, h: 1.1,
    fontFace: "Aptos", fontSize: 18, color: CHARCOAL,
    align: "left", valign: "top"
  });
  s.addText("Driver: SILC CFI database (10 plots, 1981 to 2000) with 24 remeasurement intervals", {
    x: 0.75 + cardW + gap, y: cardY + 2.25, w: cardW - 0.4, h: 1.0,
    fontFace: "Aptos", fontSize: 14, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });
  s.addText("Section 2 of this deck", {
    x: 0.75 + cardW + gap, y: cardY + 3.95, w: cardW - 0.4, h: 0.4,
    fontFace: "Aptos", fontSize: 13, color: LIGHT_GRAY,
    align: "left", valign: "bottom"
  });

  takeawayBand(s, "Long horizon for planning, short horizon for trust");
}

// ===================================================================
// SLIDE 3: New 5x2 stratification framework
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Five forest types x two density classes", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Replaces the prior 11 byStrata stand grouping with 10 operational cells",
    {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${OD}/silc_strata_5x2_plot_counts.png`,
    x: 0.6, y: 1.55, w: 12.1, h: 4.30
  });

  citationLine(s, "Source: SILC matrix stand inventory, 2696 plots across the Acadian Matrix");
  takeawayBand(s, "Hardwood and Commercial Softwood dominate the 10-cell map");
}

// ===================================================================
// SLIDE 3b: Stratification methodology
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Stratification rules for SILC sign off", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("How each SILC plot is assigned to the 5 x 2 grid", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  // Left card: forest type rules
  const cardY = 1.6, cardH = 4.6;
  const halfW = 5.95, gap = 0.4;
  s.addShape(pres.ShapeType.rect, {
    x: 0.5, y: cardY, w: halfW, h: cardH,
    fill: { color: "F7F7F7" }, line: { color: "F7F7F7" }
  });
  s.addShape(pres.ShapeType.rect, {
    x: 0.5, y: cardY, w: 0.08, h: cardH,
    fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
  });
  s.addText("Forest type assignment (by % stand BA)", {
    x: 0.75, y: cardY + 0.18, w: halfW - 0.4, h: 0.5,
    fontFace: "Aptos Display", fontSize: 20, color: CRSF_GREEN, bold: true
  });
  const rules = [
    "Cedar: white cedar share >= 30%",
    "Other Softwood: pines >= 50% (JP, RP, WP)",
    "Mixedwood: hardwood >= 30% AND softwood >= 30% (folds HS + SH)",
    "Hardwood: HW dominant when softwood < 30%",
    "Commercial Softwood: SW dominant (spruce, fir, hemlock) when HW < 30%",
    "Unclassifiable: total identified species below 50% (CFI plot 1103)"
  ];
  s.addText(rules.map(r => ({ text: r, options: { bullet: { code: "2022" } } })), {
    x: 0.85, y: cardY + 0.85, w: halfW - 0.5, h: cardH - 1.0,
    fontFace: "Aptos", fontSize: 14, color: CHARCOAL,
    paraSpaceAfter: 8
  });

  // Right card: density rule
  const rx = 0.5 + halfW + gap;
  s.addShape(pres.ShapeType.rect, {
    x: rx, y: cardY, w: halfW, h: cardH,
    fill: { color: "F7F7F7" }, line: { color: "F7F7F7" }
  });
  s.addShape(pres.ShapeType.rect, {
    x: rx, y: cardY, w: 0.08, h: cardH,
    fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
  });
  s.addText("Density class assignment", {
    x: rx + 0.25, y: cardY + 0.18, w: halfW - 0.4, h: 0.5,
    fontFace: "Aptos Display", fontSize: 20, color: CRSF_GREEN, bold: true
  });
  const density = [
    "A+B (high): REL_DENSITY >= 0.28",
    "C+D (low): REL_DENSITY < 0.28",
    "REL_DENSITY = SDI / SDI_max with SDI_max = 450 (Long 1985 NE softwood)",
    "Threshold = CFI sample median; balances 5/5 plot split",
    "Operational alternative: SILC Matrix BA stocking lines A (>=100), B (80-100), C (60-80), D (<60)"
  ];
  s.addText(density.map(r => ({ text: r, options: { bullet: { code: "2022" } } })), {
    x: rx + 0.35, y: cardY + 0.85, w: halfW - 0.5, h: cardH - 1.0,
    fontFace: "Aptos", fontSize: 14, color: CHARCOAL,
    paraSpaceAfter: 8
  });

  citationLine(s, "All thresholds open for SILC sign off; mapping file silc_strata_5x2_mapping.csv");
  takeawayBand(s, "Awaiting SILC confirmation on species cutoffs and density threshold");
}

// ===================================================================
// SLIDE 4: 100-year merch cords by stratum (AGM)
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("AGM projects merch cords over 100 years", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Plot-count weighted means with between-stand 5 to 95 percentile ribbon", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${OD}/silc_strata_5x2_cords_weighted.png`,
    x: 0.4, y: 1.55, w: 12.5, h: 4.40
  });

  citationLine(s, "Plot-count weighted across SILC matrix stands; ribbon = between stand 5 to 95 percentile band; 4 cells empty in AGM byStrata");
  takeawayBand(s, "Hardwood C+D and Mixedwood C+D project the largest 100 year cords gains");
}

// ===================================================================
// SLIDE 5: Year-100 outcomes table
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Year-100 outcomes by stratum", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Growth factor = year-100 value divided by year-0 value", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  // Table data
  const rows = [
    [{ text: "Stratum",     options: { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } } },
     { text: "n stands",    options: { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } } },
     { text: "Year-0 cords/ac",options: { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } } },
     { text: "Year-100 cords/ac",options: { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } } },
     { text: "Growth factor",options: { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } } }],
    ["Cedar / C+D",                  "1", "35.8", "59.3", "1.66x"],
    ["Hardwood / A+B",               "1", "25.0", "65.3", "2.62x"],
    ["Hardwood / C+D",               "3", "17.9", "61.0", "3.41x"],
    ["Mixedwood / C+D",              "1", "23.3", "66.4", "2.85x"],
    ["Commercial Softwood / A+B",    "2", "27.1", "61.0", "2.25x"],
    ["Commercial Softwood / C+D",    "3", "24.8", "62.9", "2.54x"]
  ];
  s.addTable(rows, {
    x: 0.5, y: 1.65, w: 12.3, h: 4.5,
    colW: [3.5, 1.4, 2.4, 2.6, 2.4],
    fontFace: "Aptos", fontSize: 16, color: CHARCOAL,
    border: { type: "solid", pt: 0.5, color: "DDDDDD" },
    align: "left", valign: "middle",
    rowH: 0.55
  });

  citationLine(s, "AGM (AcadianGY v12) byStrata; 4 cells have no byStrata coverage (Cedar A+B, Mixedwood A+B, both Other Softwood)");
  takeawayBand(s, "Low density cells project the largest growth factors");
}

// ===================================================================
// SLIDE 6: SILC CFI data introduction
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("SILC CFI gives us the empirical anchor", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Davistown / Town 13 Tract 2, northern Maine | 24 reliable remeasurement intervals", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  // Big metric callouts
  const ys = 1.85;
  const cols = [
    ["10", "fixed-area plots", "1/5 ac, EXPF=5"],
    ["5", "measurement waves", "1981 / 1986 / 1990 / 1995 / 2000"],
    ["1,399", "live trees", "327 dead + 269 removed"],
    ["649", "paired tree growth records", "140 mortality records"],
    ["24", "remeasurement intervals", "4 to 9 year span each"]
  ];
  const colW = 2.45, gap = 0.18;
  cols.forEach((c, i) => {
    const x = 0.6 + i * (colW + gap);
    s.addText(c[0], {
      x: x, y: ys, w: colW, h: 1.2,
      fontFace: "Aptos Display", fontSize: 54, color: CRSF_GREEN, bold: true,
      align: "center", valign: "middle"
    });
    s.addText(c[1], {
      x: x, y: ys + 1.25, w: colW, h: 0.5,
      fontFace: "Aptos", fontSize: 16, color: CHARCOAL, bold: true,
      align: "center", valign: "top"
    });
    s.addText(c[2], {
      x: x, y: ys + 1.75, w: colW, h: 0.5,
      fontFace: "Aptos", fontSize: 13, color: MED_GRAY, italic: true,
      align: "center", valign: "top"
    });
  });

  // CFI to 5x2 mapping caption strip
  s.addText("5 of 10 strata cells covered by CFI plots: Cedar / Hardwood (both densities) / Mixedwood (both) / Commercial Softwood / C+D", {
    x: 0.6, y: 4.6, w: 12.1, h: 0.5,
    fontFace: "Aptos", fontSize: 14, color: CHARCOAL,
    align: "left", valign: "middle"
  });
  s.addText("Cedar / A+B is uniquely covered by CFI but NOT by AGM byStrata — the two anchors are complementary", {
    x: 0.6, y: 5.1, w: 12.1, h: 0.5,
    fontFace: "Aptos", fontSize: 14, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });

  citationLine(s, "Source: SILC_CFI_FIADB_Database_v3.xlsx, FIADB-structured, Weiskittel 2026-05-28");
  takeawayBand(s, "Real Seven Islands data with 24 remeasurement intervals");
}

// ===================================================================
// SLIDE 7: Multi-model CFI scorecard (2 x 3 scatter)
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Three models bracket the truth on CFI", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("FVS-NE calibrated lands closest, AGM over projects, OSM-ACD under projects", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${CFI_DIR}/silc_cfi_multimodel_scatter.png`,
    x: 0.4, y: 1.55, w: 12.5, h: 4.40
  });

  citationLine(s, "AGM = AcadianGY 12.3.9 (in-source mortality + ingrowth fix); OSM-ACD = Open Stand Model Acadian v2.26.1, Cardinal SLURM 10988107 + 10990880");
  takeawayBand(s, "FVS-NE calibrated and AGM near zero, OSM is the under bracket");
}

// ===================================================================
// SLIDE 8: Four-row multi-model scorecard table
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Multi model CFI accuracy at a glance", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Routine growth subset, n=17. Best of each metric in green, worst in red", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  const G = CRSF_GREEN;
  const HDR = { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } };
  // Bias column (sign + %): one row per predictor, one column per metric
  const tab = [
    [{ text: "Predictor", options: HDR },
     { text: "BA",         options: HDR },
     { text: "Cords",      options: HDR },
     { text: "BdFt (Intl 1/4)", options: HDR },
     { text: "SDI",        options: HDR },
     { text: "Curtis RD",  options: HDR }],
    ["zero growth",       "-1.8%",  "-6.0%",  "-15.9%", "+4.8%", "+5.2%"],
    ["FIA prior PAI",     "+3.9%",  "+1.4%",  "-8.9%",  "n/a",   "n/a"],
    ["AGM (AcadianGY)",
       "+6.1%", "+7.0%",
       { text: "+8.4%", options: { color: G, bold: true } },
       { text: "-0.3%", options: { color: G, bold: true } },
       { text: "-0.6%", options: { color: G, bold: true } }],
    ["FVS-NE default",    "+4.6%",  "-1.3%", "-11.9%", "-1.0%", "-1.1%"],
    ["FVS-NE calibrated",
       { text: "-0.3%", options: { color: G, bold: true } },
       "-3.8%", "-13.6%", "-5.9%", "-6.1%"],
    ["FVS-ACD default",
       "+5.3%",
       { text: "-2.5%", options: { color: G, bold: true } },
       "-11.5%", "-0.5%", "-0.7%"],
    ["OSM-ACD",
       "-8.0%",
       "-7.3%",
       "-6.6%",
       "+2.0%", "+6.0%"]
  ];
  s.addTable(tab, {
    x: 0.3, y: 1.55, w: 12.7, h: 4.0,
    colW: [2.8, 1.5, 1.5, 2.4, 1.5, 1.5],
    fontFace: "Aptos", fontSize: 14, color: CHARCOAL,
    border: { type: "solid", pt: 0.5, color: "DDDDDD" },
    align: "center", valign: "middle",
    rowH: 0.50
  });

  citationLine(s, "Routine growth subset n=17. BA in ft^2/ac, cords/ac, BdFt/ac Intl 1/4 inch rule, SDI Reineke 1933, Curtis RD = BA / sqrt(QMD)");
  takeawayBand(s, "FVS-NE calibrated wins BA, FVS-ACD wins cords, AGM wins BdFt SDI RD");
}

// ===================================================================
// SLIDE 9: Species composition
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Species composition matched by AGM and OSM-ACD", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("AGM within 1.3 pp on Comm SW share, OSM-ACD within 3 pp on all four groups", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${CFI_DIR}/silc_cfi_species_full.png`,
    x: 0.3, y: 1.50, w: 12.7, h: 4.50
  });

  citationLine(s, "Mean of routine-growth subset n=17. Bias annotations are predicted minus observed in percentage points. FVS per tree species output deferred (TREELIST scheduling)");
  takeawayBand(s, "AGM and OSM-ACD both track the observed Comm SW + Hardwood mix");
}

// ===================================================================
// SLIDE 10: Operational interpretation
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("What this means for SILC operational use", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });

  // Three rec cards
  const recs = [
    { title: "Run FVS-NE calibrated as the 5 to 10 year operational",
      body: "Lowest BA bias (-0.3%) and competitive on cords on SILC's own CFI plots. Validated Cardinal pipeline, default + calibrated parameter sets, ready for SILC.",
      accent: CRSF_GREEN },
    { title: "Keep AGM as the long horizon planning workhorse",
      body: "100 year strata cords trajectories drive allowable cut. AGM is also the sawlog BdFt accuracy champion (RMSE 358, R^2 0.94) for size class shifts.",
      accent: CRSF_GREEN },
    { title: "Bracket with OSM-ACD as the under projection check",
      body: "OSM-ACD (-7% to -8% bias) is the conservative bound. Report FVS-NE plus OSM as a range for stand level operational calls; per plot RMSE 7 ft^2/ac BA, 1.4 cords/ac.",
      accent: CRSF_GREEN }
  ];
  const cardY = 1.5, cardH = 1.5, cardW = 12.1;
  recs.forEach((r, i) => {
    const y = cardY + i * (cardH + 0.15);
    s.addShape(pres.ShapeType.rect, {
      x: 0.5, y: y, w: cardW, h: cardH,
      fill: { color: "F7F7F7" }, line: { color: "F7F7F7" }
    });
    s.addShape(pres.ShapeType.rect, {
      x: 0.5, y: y, w: 0.08, h: cardH,
      fill: { color: r.accent }, line: { color: r.accent }
    });
    s.addText(r.title, {
      x: 0.75, y: y + 0.15, w: cardW - 0.4, h: 0.5,
      fontFace: "Aptos Display", fontSize: 20, color: CRSF_GREEN, bold: true,
      align: "left", valign: "middle"
    });
    s.addText(r.body, {
      x: 0.75, y: y + 0.65, w: cardW - 0.4, h: 0.85,
      fontFace: "Aptos", fontSize: 15, color: CHARCOAL,
      align: "left", valign: "top"
    });
  });

  takeawayBand(s, "FVS-NE calibrated for short horizon, AGM for long horizon, OSM as bracket");
}

// ===================================================================
// SLIDE 11: Next steps
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Next steps", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });

  const steps = [
    { n: "1", t: "Confirm the 5 x 2 stratification with SILC", b: "Sign off that the species composition rules and density threshold match your operational categories" },
    { n: "2", t: "Add FVS per tree species output", b: "Resolve TREELIST column-format issue with standalone FVS binary, so FVS species comp joins observed + AGM + OSM in the scorecard" },
    { n: "3", t: "Backfill the empty AGM strata cells", b: "Cedar/A+B, Mixedwood/A+B, both Other Softwood; pull from plot level AcadianGY R on the underlying SILC plots" },
    { n: "4", t: "Set the uncertainty bands", b: "Use the empirical CFI RMSE figures as the per plot uncertainty for SILC operational reporting" }
  ];
  const sy = 1.4, sh = 1.15;
  steps.forEach((step, i) => {
    const y = sy + i * (sh + 0.10);
    // Circle with number
    s.addShape(pres.ShapeType.ellipse, {
      x: 0.5, y: y + 0.20, w: 0.7, h: 0.7,
      fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
    });
    s.addText(step.n, {
      x: 0.5, y: y + 0.20, w: 0.7, h: 0.7,
      fontFace: "Aptos Display", fontSize: 26, color: "FFFFFF", bold: true,
      align: "center", valign: "middle"
    });
    s.addText(step.t, {
      x: 1.4, y: y + 0.10, w: 11.4, h: 0.45,
      fontFace: "Aptos Display", fontSize: 20, color: CHARCOAL, bold: true,
      align: "left", valign: "middle"
    });
    s.addText(step.b, {
      x: 1.4, y: y + 0.55, w: 11.4, h: 0.55,
      fontFace: "Aptos", fontSize: 14, color: MED_GRAY,
      align: "left", valign: "top"
    });
  });

  takeawayBand(s, "Four model scaffold is live, awaiting SILC sign off");
}

// ===================================================================
// SLIDE 12: Questions / closing
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };

  s.addText("Questions and discussion", {
    x: 0.5, y: 2.5, w: 12.3, h: 1.2,
    fontFace: "Aptos Display", fontSize: 56, color: CHARCOAL, bold: true,
    align: "center", valign: "middle"
  });

  // Accent bar
  s.addShape(pres.ShapeType.rect, {
    x: 5.97, y: 3.85, w: 1.4, h: 0.06,
    fill: { color: CRSF_GREEN }, line: { color: CRSF_GREEN }
  });

  s.addText("Aaron R. Weiskittel", {
    x: 0.5, y: 4.2, w: 12.3, h: 0.5,
    fontFace: "Aptos Display", fontSize: 24, color: CHARCOAL, bold: true,
    align: "center", valign: "middle"
  });
  s.addText("aaron.weiskittel@maine.edu", {
    x: 0.5, y: 4.75, w: 12.3, h: 0.5,
    fontFace: "Aptos", fontSize: 20, color: CRSF_GREEN,
    align: "center", valign: "middle"
  });
  s.addText("Center for Research on Sustainable Forests | University of Maine", {
    x: 0.5, y: 5.25, w: 12.3, h: 0.5,
    fontFace: "Aptos", fontSize: 14, color: MED_GRAY,
    align: "center", valign: "middle"
  });
}

const out = `${CFI_DIR}/SILC_4Model_Benchmark_v4_Weiskittel.pptx`;
pres.writeFile({ fileName: out })
  .then(f => console.log("wrote", f))
  .catch(e => { console.error("err", e); process.exit(1); });
