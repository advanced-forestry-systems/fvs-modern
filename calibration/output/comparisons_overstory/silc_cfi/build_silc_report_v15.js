// build_silc_report_v8.js
// Companion Word report to SILC_4Model_Benchmark_v8 deck.
const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        ImageRun, AlignmentType, LevelFormat, HeadingLevel, BorderStyle,
        WidthType, ShadingType, PageBreak } = require("docx");

const CFI_DIR = "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi";
const OD = "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory";

// Helpers
const border = { style: BorderStyle.SINGLE, size: 4, color: "BBBBBB" };
const borders = { top: border, bottom: border, left: border, right: border };

const p = (text, opts = {}) => new Paragraph({
  spacing: { after: 120 },
  ...opts,
  children: [new TextRun({ text, ...(opts.run || {}) })],
});

const h1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  spacing: { before: 300, after: 200 },
  children: [new TextRun({ text, bold: true })],
});
const h2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  spacing: { before: 240, after: 160 },
  children: [new TextRun({ text, bold: true })],
});

const bullet = (text) => new Paragraph({
  numbering: { reference: "bullets", level: 0 },
  spacing: { after: 80 },
  children: [new TextRun(text)],
});

const tcell = (text, opts = {}) => new TableCell({
  borders,
  width: { size: opts.w || 1500, type: WidthType.DXA },
  margins: { top: 80, bottom: 80, left: 120, right: 120 },
  shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
  children: [new Paragraph({
    alignment: opts.align || AlignmentType.LEFT,
    children: [new TextRun({ text, bold: opts.bold || false, color: opts.color || undefined })],
  })],
});

const img = (path, w, h) => new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 120, after: 120 },
  children: [new ImageRun({
    type: "png",
    data: fs.readFileSync(path),
    transformation: { width: w, height: h },
    altText: { title: path.split("/").pop(), description: "figure", name: path.split("/").pop() },
  })],
});

const children = [];

// === Title block ===
children.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { before: 0, after: 80 },
  children: [new TextRun({ text: "SILC Four-Model Forest Growth Benchmark", bold: true, size: 40 })],
}));
children.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 80 },
  children: [new TextRun({ text: "Companion Report v15 — v16 refined CSI BGI + merch + FVS-ACD configs", italics: true, size: 26 })],
}));
children.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 320 },
  children: [new TextRun({ text: "Aaron R. Weiskittel, CRSF, University of Maine — May 2026", size: 22 })],
}));

// === Model clarification ===
children.push(h1("Model clarification — AGM and AcadianGY are different models"));
children.push(p("Three models are benchmarked in this report: AcadianGY R (CRSF research implementation, the Hennigar et al. growth and yield model, v12.3.9 with the in source MORTCAL patch), FVS-NE (USDA standalone Acadian variant, default and calibrated), and OSM-ACD (Open Stand Model v2.26.1 Acadian variant)."));
children.push(p("AGM is Seven Islands' operational Acadian Growth Model. It is a different model from AcadianGY and has NOT been run on the SILC CFI data or on the same 11 byStrata stands in this report. Earlier drafts conflated the two names. Every reference to AGM in v8 and v9 of this report referred to AcadianGY 12.3.9. The clearest single next step is for SILC to share an AGM trajectory file on the same 11 byStrata stands so the four model picture in this report becomes a five model picture with the operational anchor included."));

// === Executive summary ===
children.push(h1("Executive summary"));
children.push(p("Three findings drive this update."));
children.push(bullet("(1) The +18.5% AcadianGY bias on the 10 year SILC CFI long horizon scorecard is real model behavior, not a pipeline bug. Triple checking the AcadianGY treelist to BA / Cords / BdFt computation reproduces to 1e-13 ft²/ac, and observed STAND_METRICS BA matches direct tree list recompute to 0.005 ft²/ac. The bias decomposes into two mechanisms: AcadianGY under projects mortality on disturbed CFI plots (1101 and 1106 lost 27 to 34% of their stems, AcadianGY lost 0 to 19%), plus AcadianGY grows DBH at a uniform +0.8 in per decade regardless of site (plots 1100 and 1104 grew +0.05 to 0.17 in observed)."));
children.push(bullet("(2) Switching on the #126b in-source MORTCAL correction in AcadianGY 12.3.9 (default OFF) cuts the 10 year BA bias from +18.5% to +7.8%, RMSE -34%, R² from -0.20 to +0.47. Sawlog BdFt bias drops to +4.0% with R² = 0.84, the best of any model on this metric."));
children.push(bullet("(3) Pushed to 100 years, MORTCAL pulls projected BA from a ceiling around 200 ft²/ac down to a credible carrying capacity of 98 to 108 ft²/ac. Year 100 NetCords lands at 35 to 39 cords/ac across all six byStrata cells, plus CFI plot backfill at Cedar A+B (35.8 cords/ac) and Mixedwood A+B (21.9 cords/ac)."));
children.push(bullet("(4) Full 4 model year 100 view added: FVS-NE default 209 to 291 BA, FVS-NE calibrated 125 to 225, OSM-ACD 172 to 221 (covered cells). AcadianGY MORTCAL still the lowest projection. Refined inputs from CFI v3 lat/long swap BGI 3000 to 3902 and CSI 12 to 15.78 m; minimal effect on relative ranking."));
children.push(p("Recommended operational stack: AcadianGY with MORTCAL on as the long horizon operational anchor (35 to 39 cords/ac year 100); FVS-NE calibrated as the upper bracket (44 to 70 cords/ac); OSM-ACD as additional cross check where coverage exists."));

// === Triple check ===
children.push(h1("Triple check of AcadianGY long horizon calculations"));
children.push(p("The +18.5% AcadianGY bias at the 10 year CFI scorecard prompted a forensic audit. Every step of the pipeline checked out."));
children.push(p("Source code: ~/AcadianGY_12.3.9.r on Cardinal (identical to local v12.3.5 plus the opt-in MORTCAL correction). AcadianGYOneStand is annual (cyclen=1 hardcoded since 2022-09-01). The driver loops PERIOD_YR times, giving the correct horizon."));

children.push(h2("Pipeline verification checks"));
const verifyHdr = ["Check", "Result"];
children.push(new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: [5500, 3860],
  rows: [
    new TableRow({ children: [
      tcell("Check", { w: 5500, bold: true, fill: "1A3D28", color: "FFFFFF" }),
      tcell("Result", { w: 3860, bold: true, fill: "1A3D28", color: "FFFFFF", align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("AcadianGY treelist to BA/Cords/BdFt reproduces results CSV", { w: 5500 }),
      tcell("exact to 1e-13", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("STAND_METRICS BA_CURR vs direct TREE.csv recompute", { w: 5500 }),
      tcell("match to 0.005 ft²/ac", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("AcadianGYOneStand cycle length", { w: 5500 }),
      tcell("annual (cyclen=1 hardcoded)", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("Driver loop count", { w: 5500 }),
      tcell("PERIOD_YR iterations", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("BA formula 0.005454 * DBH² * EXPF", { w: 5500 }),
      tcell("correct", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("Cords formula 0.0025 * DBH² * HT * 0.90 * EXPF / 79", { w: 5500 }),
      tcell("correct", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("BdFt formula 0.01 * DBH² * HT * EXPF for DBH >= 9 in (Intl ¼)", { w: 5500 }),
      tcell("correct", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
    new TableRow({ children: [
      tcell("DBH threshold symmetry (pred vs obs)", { w: 5500 }),
      tcell("both at 4.5 in", { w: 3860, align: AlignmentType.CENTER }),
    ]}),
  ],
}));

children.push(h2("Per-plot bias decomposition"));
children.push(p("Six routine CFI pairs after excluding two establishment plots (1105, 1107 had 3x BA ingrowth). The decomposition shows that AcadianGY grows DBH a uniform +0.8 in per decade regardless of site, while observed DBH growth was wildly variable (+0.05 to +1.04 in). AcadianGY also missed two significant mortality events (1101 and 1106 lost 27 to 34% of stems)."));

const decompHdr = ["PLOT", "obs BA chg", "obs TPA chg", "obs QMD chg", "AcadianGY TPA chg", "AcadianGY QMD chg", "AcadianGY BA bias"];
const decompRows = [
  ["1100", "+5%", "0%", "+0.17", "0%", "+0.80", "+18.6%"],
  ["1101", "-7%", "-27%", "+1.04", "0%", "+0.82", "+29.7%"],
  ["1102", "+16%", "-3%", "+0.88", "0%", "+0.83", "+2.6%"],
  ["1104", "-5%", "-6%", "+0.05", "0%", "+0.83", "+28.6%"],
  ["1106", "-22%", "-34%", "+0.67", "-19%", "+0.81", "+25.7%"],
  ["1109", "+14%", "+13%", "+0.05", "0%", "+0.84", "+4.1%"],
];
children.push(new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: [1080, 1280, 1280, 1280, 1280, 1280, 1880],
  rows: [
    new TableRow({ children: decompHdr.map((t, i) => tcell(t, { w: [1080,1280,1280,1280,1280,1280,1880][i], bold: true, fill: "1A3D28", color: "FFFFFF", align: AlignmentType.CENTER })) }),
    ...decompRows.map(r => new TableRow({ children: r.map((t, i) => tcell(t, { w: [1080,1280,1280,1280,1280,1280,1880][i], align: AlignmentType.CENTER, bold: i === 6 })) })),
  ],
}));

children.push(p("Mechanism figure:", { spacing: { before: 200 } }));
children.push(img(`${CFI_DIR}/silc_cfi_long_mechanism.png`, 600, 270));

// === MORTCAL improvement ===
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(h1("MORTCAL correction cuts 10 year bias by more than half"));
children.push(p("The #126b in-source size dependent mortality correction was originally calibrated to Maine FIA. It is OFF by default in AcadianGY 12.3.9. The SILC CFI test enables it via ops$MORTCAL=TRUE, MORTCAL_INTERVAL=5."));
children.push(p("Side by side scorecard on the n=6 routine pairs (mean horizon 10 yr):"));

const mcHdr = ["Metric", "AcadianGY default", "AcadianGY MORTCAL", "Improvement"];
const mcRows = [
  ["BA bias", "+18.5%", "+7.8%", "cut by 58%"],
  ["BA RMSE", "16.2 ft²/ac", "10.7 ft²/ac", "-34%"],
  ["BA R²", "-0.20", "+0.47", "flips negative to positive"],
  ["Cords bias", "+19.3%", "+8.4%", "cut by 56%"],
  ["Cords R²", "0.58", "0.82", "+0.24"],
  ["BdFt bias", "+16.7%", "+4.0%", "single digit"],
  ["BdFt R²", "0.75", "0.84", "+0.09"],
];
children.push(new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: [2200, 2200, 2200, 2760],
  rows: [
    new TableRow({ children: mcHdr.map((t, i) => tcell(t, { w: [2200,2200,2200,2760][i], bold: true, fill: "1A3D28", color: "FFFFFF", align: AlignmentType.CENTER })) }),
    ...mcRows.map(r => new TableRow({ children: r.map((t, i) => tcell(t, { w: [2200,2200,2200,2760][i], align: AlignmentType.CENTER, bold: i === 3 && t !== "+0.09" && t !== "+0.24" })) })),
  ],
}));

children.push(p("Scatter overlay (green = default, gold = MORTCAL):", { spacing: { before: 200 } }));
children.push(img(`${CFI_DIR}/silc_cfi_long_mortcal_compare.png`, 600, 230));

// === 100 yr trajectories ===
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(h1("100 year trajectories: default vs MORTCAL"));
children.push(p("Cardinal job 11076533 ran AcadianGY 12.3.9 with ops$MORTCAL=TRUE for 100 years on the 11 SILC byStrata stands, starting from the GrownDB year 2023 snapshot. The default AcadianGY trajectory has every stratum approaching ~200 ft²/ac at year 100, an artifact of the under mortality issue. MORTCAL pulls year 100 BA into the 98 to 108 ft²/ac carrying capacity range, with sustainable yield 35 to 39 cords/ac across the populated cells."));
children.push(img(`${OD}/silc_strata_5x2_AGM_MORTCAL_BA.png`, 620, 360));

children.push(h2("Year 100 outcomes (year 2123) — 4 model comparison"));
children.push(p("Side by side year 100 BA in ft²/ac on the same 11 byStrata starting stands plus CFI backfill on 2 of 4 empty cells. AcadianGY MORTCAL anchors the operational floor; FVS-NE calibrated and OSM-ACD cluster as the upper bracket."));
const y100Hdr = ["Stratum", "n", "AcadianGY def", "AcadianGY mc", "FVS-NE def", "FVS-NE cal", "OSM-ACD"];
const y100Rows = [
  ["Cedar / A+B (CFI backfill)", "1", "n/a", "81",  "n/a",  "n/a",  "n/a"],
  ["Cedar / C+D",                "1", "206", "99",  "291",  "218",  "249"],
  ["Hardwood / A+B",             "1", "209", "107", "209",  "144",  "n/a"],
  ["Hardwood / C+D",             "3", "198", "108", "214",  "125",  "174"],
  ["Mixedwood / A+B (CFI backfill)", "3", "n/a", "52", "n/a", "n/a", "n/a"],
  ["Mixedwood / C+D",            "1", "212", "102", "239",  "153",  "183"],
  ["Commercial SW / A+B",        "2", "203", "98",  "274",  "225",  "222"],
  ["Commercial SW / C+D",        "3", "211", "99",  "264",  "186",  "199"],
  ["Other Softwood (no CFI)",   "0",  "—",  "—",   "—",   "—",    "—"],
];
const y100W = [2400, 600, 1050, 1050, 1100, 1100, 2060];
children.push(new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: y100W,
  rows: [
    new TableRow({ children: y100Hdr.map((t, i) => tcell(t, { w: y100W[i], bold: true, fill: "1A3D28", color: "FFFFFF", align: AlignmentType.CENTER })) }),
    ...y100Rows.map((r) => {
      const fill = r[0].includes("backfill") ? "FFF8E1" : (r[0].includes("no CFI") ? "F4F4F4" : undefined);
      const color = r[0].includes("no CFI") ? "888888" : undefined;
      return new TableRow({ children: r.map((t, i) => tcell(t, { w: y100W[i], align: AlignmentType.CENTER, fill, color })) });
    }),
  ],
}));

children.push(h2("Year 100 cords/ac (operational sustained yield)"));
const cHdr = ["Stratum", "n", "AcadianGY MORTCAL", "FVS-NE cal", "OSM-ACD"];
const cRows = [
  ["Cedar / A+B (CFI)",  "1", "35.8", "n/a",  "n/a"],
  ["Cedar / C+D",        "1", "36.9", "66.4", "41.5"],
  ["Hardwood / A+B",     "1", "39.3", "51.1", "n/a"],
  ["Hardwood / C+D",     "3", "38.9", "44.0", "32.4"],
  ["Mixedwood / A+B (CFI)", "3", "21.9", "n/a", "n/a"],
  ["Mixedwood / C+D",    "1", "37.6", "51.3", "34.6"],
  ["Commercial SW / A+B","2", "35.9", "69.7", "40.0"],
  ["Commercial SW / C+D","3", "36.2", "58.9", "36.2"],
];
const cW = [2700, 600, 2020, 2020, 2020];
children.push(new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: cW,
  rows: [
    new TableRow({ children: cHdr.map((t, i) => tcell(t, { w: cW[i], bold: true, fill: "1A3D28", color: "FFFFFF", align: AlignmentType.CENTER })) }),
    ...cRows.map((r) => {
      const fill = r[0].includes("CFI") ? "FFF8E1" : undefined;
      return new TableRow({ children: r.map((t, i) => tcell(t, { w: cW[i], align: AlignmentType.CENTER, fill })) });
    }),
  ],
}));

// === Refined inputs ===
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(h1("Refined input parameters from SILC CFI v3 lat/long"));
children.push(p("v3 of the SILC CFI database includes approximate lat/long (46.4628 N, 68.4253 W, Davistown). Used to sample regional rasters and replace conservative defaults."));
const refHdr = ["Parameter", "Default", "Refined", "Source", "Effect"];
const refRows = [
  ["BGI (OSM)", "3000", "3902", "ME_BGI_V1.tif sample", "OSM productivity input"],
  ["CSI (AcadianGY/FVS-ACD)", "12.0 m", "15.78 m", "CSI_2030.tif sample", "Was already 14.33 m for byStrata"],
  ["FVS-NE site index", "42 ft", "35 ft", "Empirical from CFI dom heights", "Year 100 BA effect <0.3%"],
  ["AcadianGY MORTCAL", "off", "TRUE", "#126b in-source patch", "Cuts 10 yr bias from +18.5% to +7.8%"],
];
const refW = [1700, 1300, 1300, 2500, 2560];
children.push(new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: refW,
  rows: [
    new TableRow({ children: refHdr.map((t, i) => tcell(t, { w: refW[i], bold: true, fill: "1A3D28", color: "FFFFFF", align: AlignmentType.CENTER })) }),
    ...refRows.map((r) => new TableRow({ children: r.map((t, i) => tcell(t, { w: refW[i], align: AlignmentType.CENTER })) })),
  ],
}));
children.push(p("FVS-NE Acadian variant year 100 BA changed less than 0.3% per cell when SI dropped from 42 to 35. The variant uses tree level dDBH equations not strongly driven by stand level SI. The OSM BGI swap (3000 to 3902) pushed OSM productivity up slightly; OSM remained convergent."));

// === Recommendations ===
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(h1("Operational recommendations"));
children.push(bullet("AcadianGY with MORTCAL on is the new long horizon champion. Use it for SILC operational sawlog projections (BdFt R² = 0.84) and Cords (R² = 0.82). Set ops$MORTCAL = TRUE, MORTCAL_INTERVAL = 5 in any AcadianGY call."));
children.push(bullet("FVS-NE calibrated stays best on the short horizon BA bias (-0.3% at 5 yr). Pair with AcadianGY MORTCAL as a cross check on the operational 5 yr plans."));
children.push(bullet("OSM-ACD is the conservative under projection (-9.7% BA bias at 10 yr). Report AcadianGY MORTCAL plus OSM as the operational range; per plot RMSE 10.7 ft²/ac BA at 10 yr horizon is the empirical uncertainty band."));
children.push(bullet("Year 100 cords/ac estimates: 35 to 39 across six byStrata cells and CFI backfill Cedar A+B; Mixedwood A+B noticeably lower (22) reflecting heavier mortality in the high density mixed pine / hardwood stands."));

children.push(h1("Open questions for SILC sign off"));
children.push(p("1. Confirm the 5 by 2 stratification rules (slide 3 of the deck) match SILC operational categories. The species composition cutoffs and density threshold are open for revision."));
children.push(p("2. Validate the year 100 BA target of 98 to 108 ft²/ac against SILC's expected carrying capacity for managed Acadian stands."));
children.push(p("3. Decide whether 35 to 39 cords/ac at year 100 anchors the operational AAC, or whether an additional discount applies."));

children.push(h1("File index"));
children.push(p("All artifacts in /home/aweiskittel/Documents/Claude/fvs-modern/silc_cfi_v8_deliverables/. Key files:"));
children.push(bullet("SILC_4Model_Benchmark_v8_Weiskittel.pdf — 17 slide operational deck"));
children.push(bullet("silc_strata_5x2_AGM_MORTCAL_trajectories.csv — strata level 100 yr trajectories"));
children.push(bullet("silc_strata_5x2_year100_mortcal_vs_default.csv — year 2123 outcomes table"));
children.push(bullet("silc_cfi_backfill_100yr_mortcal_trajectories.csv — Cedar A+B and Mixedwood A+B from CFI plots"));
children.push(bullet("silc_cfi_long_mechanism.png and silc_cfi_long_mortcal_compare.png — supporting figures"));
children.push(bullet("draft_email_to_ian_ryan.md — handoff email"));

// === Compose document ===
const doc = new Document({
  styles: {
    default: { document: { run: { font: "Arial", size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 32, bold: true, font: "Arial", color: "1A3D28" },
        paragraph: { spacing: { before: 240, after: 180 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: "Arial", color: "1A3D28" },
        paragraph: { spacing: { before: 200, after: 140 }, outlineLevel: 1 } },
    ],
  },
  numbering: {
    config: [{
      reference: "bullets",
      levels: [{ level: 0, format: LevelFormat.BULLET, text: "•",
        alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 720, hanging: 360 } } } }],
    }],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 },
      },
    },
    children,
  }],
});

const out = `${CFI_DIR}/SILC_4Model_Benchmark_Report_v15_Weiskittel.docx`;
Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(out, buf);
  console.log("wrote", out);
}).catch(err => { console.error(err); process.exit(1); });
