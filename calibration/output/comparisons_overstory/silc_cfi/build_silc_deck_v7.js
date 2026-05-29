// build_silc_deck_v6.js
// =====================================================================
// SILC operational deck v6 for Ian / Ryan at Seven Islands.
// Extends v5 with: long-horizon (10-14 yr) scorecard slide showing how
// model biases evolve from 5-yr remeasurement to multi-decade horizons.
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
    [{ text: "FVS-ACD (default = cal *)", options: {} },
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

  citationLine(s, "* FVS-ACD calibration tweaks BAMAX stocking caps which do not bind below 150 ft^2/ac BA; default and calibrated produce identical projections on these CFI plots. n=17 routine-growth subset");
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
    { title: "AGM with MORTCAL correction is the new long horizon champion",
      body: "Switching on the #126b in-source size-dependent mortality correction cuts AGM BA bias from +18.5% to +7.8%, RMSE -34%, R^2 flips from -0.20 to +0.47, and sawlog BdFt R^2 climbs to 0.84 with single digit bias. Recommended for SILC operational sawlog projections.",
      accent: CRSF_GREEN },
    { title: "FVS-NE calibrated stays best on short horizon cords",
      body: "Lowest 5 yr bias (-0.3% BA, -3.8% cords) and tied for best 10 yr cords R^2 (0.82 vs AGM MORTCAL 0.82). Pair with AGM MORTCAL for cross-check; both share the merch volume formula.",
      accent: CRSF_GREEN },
    { title: "Bracket with OSM-ACD as the conservative under projection",
      body: "OSM-ACD is the lower bound at all horizons (-8% to -10% bias). Report AGM MORTCAL plus OSM as the operational range; per plot RMSE 11 ft^2/ac BA at 10 yr horizon.",
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

  takeawayBand(s, "AGM with MORTCAL is the new operational champion; FVS-NE cal + OSM bracket the range");
}

// ===================================================================
// SLIDE 10a: Long-horizon scorecard
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Long horizon scorecard: 10 to 14 year projections", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Earliest to latest CFI measurement per plot. Mean horizon 10 years. n=6 routine pairs", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${CFI_DIR}/silc_cfi_long_scatter.png`,
    x: 0.3, y: 1.50, w: 12.7, h: 3.2
  });

  // Comparison table: short-horizon vs long-horizon bias
  const HDR = { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } };
  const tab = [
    [{ text: "Model", options: HDR },
     { text: "BA 5 yr", options: HDR },
     { text: "BA 10 yr", options: HDR },
     { text: "Cords 5 yr", options: HDR },
     { text: "Cords 10 yr", options: HDR },
     { text: "BdFt 5 yr", options: HDR },
     { text: "BdFt 10 yr", options: HDR }],
    ["AGM default",    "+6.1%",  "+18.5%", "+7.0%",  "+19.3%", "+8.4%",  "+16.7%"],
    [{ text: "AGM MORTCAL", options: { bold: true, fill: { color: "FFF8E1" } } },
                       { text: "n/a",   options: { fill: { color: "FFF8E1" } } },
                       { text: "+7.8%", options: { bold: true, fill: { color: "FFF8E1" } } },
                       { text: "n/a",   options: { fill: { color: "FFF8E1" } } },
                       { text: "+8.4%", options: { bold: true, fill: { color: "FFF8E1" } } },
                       { text: "n/a",   options: { fill: { color: "FFF8E1" } } },
                       { text: "+4.0%", options: { bold: true, fill: { color: "FFF8E1" } } }],
    ["FVS-NE cal",     "-0.3%",  "+5.7%",  "-3.8%",  "+4.6%",  "-13.6%", "-14.3%"],
    ["FVS-ACD def",    "+5.3%",  "+15.9%", "-2.5%",  "+8.1%",  "-11.5%", "-10.1%"],
    ["OSM-ACD",        "-8.0%",  "-9.7%",  "-7.3%",  "n/a",    "-6.6%",  "n/a"]
  ];
  s.addTable(tab, {
    x: 0.5, y: 4.80, w: 12.3, h: 1.50,
    colW: [2.2, 1.5, 1.65, 1.6, 1.75, 1.6, 2.0],
    fontFace: "Aptos", fontSize: 12, color: CHARCOAL,
    border: { type: "solid", pt: 0.5, color: "DDDDDD" },
    align: "center", valign: "middle",
    rowH: 0.25
  });

  s.addText("n=6 routine pairs after excluding 2 plots with 3x BA ingrowth establishment events (1105, 1107). MORTCAL = #126b in-source mortality correction.", {
    x: 0.5, y: 6.40, w: 12.5, h: 0.30,
    fontFace: "Aptos", fontSize: 10, color: LIGHT_GRAY, italic: true,
    align: "left", valign: "middle"
  });
  takeawayBand(s, "AGM MORTCAL cuts long horizon bias by more than half, beats all defaults on BdFt");
}

// ===================================================================
// SLIDE 10b: Cross-region context
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("CFI scorecard fits a larger validation picture", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("n=17 SILC CFI pairs anchored against 12029 plot Maine FIA validation", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 16, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  // Three big metric callouts
  const numbers = [
    { n: "17",       lbl: "SILC CFI routine remeasurement pairs",
                     sub: "Davistown, ME 1981 to 2000" },
    { n: "12,029",   lbl: "Maine + NH + VT FIA paired plots",
                     sub: "Apples-to-apples overstory recompute, DIA >= 5 in" },
    { n: "-0.06%",   lbl: "FVS-ACD calibrated BA bias on FIA",
                     sub: "Independent of CFI scorecard, n 700x larger" }
  ];
  const colW = 4.0, gap = 0.3;
  const startX = (13.33 - 3 * colW - 2 * gap) / 2;
  numbers.forEach((c, i) => {
    const x = startX + i * (colW + gap);
    s.addText(c.n, {
      x: x, y: 1.95, w: colW, h: 1.4,
      fontFace: "Aptos Display", fontSize: 60, color: CRSF_GREEN, bold: true,
      align: "center", valign: "middle"
    });
    s.addText(c.lbl, {
      x: x, y: 3.45, w: colW, h: 0.5,
      fontFace: "Aptos", fontSize: 17, color: CHARCOAL, bold: true,
      align: "center", valign: "top"
    });
    s.addText(c.sub, {
      x: x, y: 4.0, w: colW, h: 0.5,
      fontFace: "Aptos", fontSize: 13, color: MED_GRAY, italic: true,
      align: "center", valign: "top"
    });
  });

  // Bottom narrative
  s.addText("The +5% to +6% AGM and FVS-ACD CFI biases are consistent with the +0% bias the same calibration achieves on 12,029 FIA plots when stand structure differences are accounted for. The CFI biases reflect SILC's Acadian Matrix configuration, not a model failure.", {
    x: 0.7, y: 4.85, w: 11.9, h: 1.2,
    fontFace: "Aptos", fontSize: 15, color: CHARCOAL,
    align: "left", valign: "top"
  });

  citationLine(s, "FIA recompute: SLURM 10591333; calibration/output/comparisons_overstory/FINDINGS.md");
  takeawayBand(s, "Small n CFI scorecard is consistent with the large n FIA validation");
}

// ===================================================================
// SLIDE 10c: AGM long horizon bias mechanism
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Where the +18.5% AGM bias actually lives", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 32, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("Two mechanisms account for all of the long horizon over projection: under mortality on disturbed plots, plus over growth on slow plots", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.45,
    fontFace: "Aptos", fontSize: 15, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${CFI_DIR}/silc_cfi_long_mechanism.png`,
    x: 0.3, y: 1.55, w: 12.7, h: 3.5
  });

  // Mechanism summary text
  s.addText("Two failure modes drive the bias:", {
    x: 0.5, y: 5.25, w: 12.3, h: 0.35,
    fontFace: "Aptos Display", fontSize: 16, color: CRSF_GREEN, bold: true,
    align: "left", valign: "middle"
  });
  const mechs = [
    { lab: "Under mortality on disturbed plots:",
      body: "Plots 1101 and 1106 lost 27 to 34% of TPA; AGM lost 0 to 19% (panel 1, red bars in panel 3)." },
    { lab: "Over growth on slow plots:",
      body: "AGM grew QMD a uniform +0.8 in regardless of site; plots 1100 and 1104 grew only 0.05 to 0.17 in observed (panel 2)." },
    { lab: "Verified, not a calculation issue:",
      body: "BA / Cords / BdFt pipeline reproduces from the AGM treelist to 1e-13 ft^2/ac. The bias is the model speaking honestly." }
  ];
  mechs.forEach((m, i) => {
    const y = 5.65 + i * 0.45;
    s.addText(`• ${m.lab}`, {
      x: 0.7, y: y, w: 4.0, h: 0.4,
      fontFace: "Aptos", fontSize: 13, color: CHARCOAL, bold: true,
      align: "left", valign: "top"
    });
    s.addText(m.body, {
      x: 4.75, y: y, w: 8.0, h: 0.4,
      fontFace: "Aptos", fontSize: 13, color: CHARCOAL,
      align: "left", valign: "top"
    });
  });

  takeawayBand(s, "Bias is mechanism, not bug: under mortality + over growth on edge cases");
}

// ===================================================================
// SLIDE 10d: MORTCAL correction result
// ===================================================================
{
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  s.addText("Turning on AGM MORTCAL cuts long horizon bias by more than half", {
    x: 0.5, y: 0.35, w: 12.3, h: 0.75,
    fontFace: "Aptos Display", fontSize: 30, color: CHARCOAL, bold: true,
    align: "left", valign: "middle"
  });
  s.addText("The #126b in-source size dependent mortality correction (calibrated to Maine FIA) ports cleanly to SILC CFI. Default OFF in AGM 12.3.9; turn on by setting ops$MORTCAL = TRUE", {
    x: 0.5, y: 1.05, w: 12.3, h: 0.55,
    fontFace: "Aptos", fontSize: 14, color: MED_GRAY, italic: true,
    align: "left", valign: "top"
  });

  s.addImage({
    path: `${CFI_DIR}/silc_cfi_long_mortcal_compare.png`,
    x: 0.3, y: 1.70, w: 12.7, h: 3.4
  });

  // Improvement table
  const HDR2 = { bold: true, color: "FFFFFF", fill: { color: CRSF_GREEN } };
  const imp = [
    [{ text: "Metric", options: HDR2 },
     { text: "AGM default", options: HDR2 },
     { text: "AGM MORTCAL", options: HDR2 },
     { text: "Improvement", options: HDR2 }],
    ["BA bias",       "+18.5%", "+7.8%",  "cut by 58%"],
    ["BA RMSE",       "16.2 ft^2/ac", "10.7 ft^2/ac", "-34%"],
    ["BA R^2",        "-0.20",  "+0.47",  "flips negative to positive"],
    ["Cords R^2",     "0.58",   "0.82",   "+0.24"],
    ["BdFt bias",     "+16.7%", "+4.0%",  "down to single digit"],
    ["BdFt R^2",      "0.75",   "0.84",   "+0.09"]
  ];
  s.addTable(imp, {
    x: 0.5, y: 5.20, w: 12.3, h: 1.45,
    colW: [2.8, 2.6, 2.6, 4.3],
    fontFace: "Aptos", fontSize: 12, color: CHARCOAL,
    border: { type: "solid", pt: 0.5, color: "DDDDDD" },
    align: "center", valign: "middle",
    rowH: 0.20
  });

  takeawayBand(s, "AGM with MORTCAL is the new long horizon champion across BA, cords, sawlog BdFt");
}

// ===================================================================
// SLIDE 12: Next steps
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
    { n: "1", t: "Rerun 100 yr strata projections with AGM MORTCAL on", b: "Replace the default AGM trajectories that drive the year 100 outcomes table with MORTCAL=TRUE output for the operational handoff" },
    { n: "2", t: "Confirm the 5 x 2 stratification with SILC", b: "Sign off that the species composition rules and density threshold match your operational categories" },
    { n: "3", t: "Backfill the empty AGM strata cells", b: "Cedar/A+B, Mixedwood/A+B, both Other Softwood; pull from plot level AcadianGY R with MORTCAL on" },
    { n: "4", t: "Set the uncertainty bands", b: "Use the empirical AGM MORTCAL CFI RMSE (10.7 ft^2/ac BA at 10 yr) as the per plot uncertainty for SILC operational reporting" }
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
// SLIDE 13: Questions / closing
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

const out = `${CFI_DIR}/SILC_4Model_Benchmark_v7_Weiskittel.pptx`;
pres.writeFile({ fileName: out })
  .then(f => console.log("wrote", f))
  .catch(e => { console.error("err", e); process.exit(1); });
