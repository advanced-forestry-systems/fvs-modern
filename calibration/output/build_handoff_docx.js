const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Footer, AlignmentType, LevelFormat, HeadingLevel, BorderStyle,
  WidthType, ShadingType, PageNumber, ExternalHyperlink } = require("docx");
const fs = require("fs");

const BLUE = "1F4E79", GREY = "595959", LINE = "CCCCCC";
const border = { style: BorderStyle.SINGLE, size: 1, color: LINE };
const borders = { top: border, bottom: border, left: border, right: border };
const cmarg = { top: 80, bottom: 80, left: 120, right: 120 };

function p(runs) {
  if (typeof runs === "string") runs = [new TextRun(runs)];
  return new Paragraph({ spacing: { after: 120 }, children: runs });
}
function bullet(text) {
  return new Paragraph({ numbering: { reference: "b", level: 0 },
    spacing: { after: 40 }, children: [new TextRun(text)] });
}
function numd(text) {
  return new Paragraph({ numbering: { reference: "n", level: 0 },
    spacing: { after: 40 }, children: [new TextRun(text)] });
}
function h(text, level) { return new Paragraph({ heading: level, children: [new TextRun(text)] }); }
function cell(text, w, fill, bold) {
  return new TableCell({ borders, width: { size: w, type: WidthType.DXA }, margins: cmarg,
    shading: fill ? { fill, type: ShadingType.CLEAR } : undefined,
    children: [new Paragraph({ children: [new TextRun({ text, bold: !!bold, size: 20 })] })] });
}
function table(headers, rows, widths) {
  const hr = new TableRow({ tableHeader: true, children: headers.map((t, i) => cell(t, widths[i], "D5E8F0", true)) });
  const body = rows.map(r => new TableRow({ children: r.map((t, i) => cell(t, widths[i], null, false)) }));
  return new Table({ width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA }, columnWidths: widths, rows: [hr, ...body] });
}

const doc = new Document({
  styles: {
    default: { document: { run: { font: "Arial", size: 22 } } },
    paragraphStyles: [
      { id: "Title", name: "Title", basedOn: "Normal", next: "Normal", run: { size: 40, bold: true, font: "Arial", color: BLUE }, paragraph: { spacing: { after: 120 } } },
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 28, bold: true, font: "Arial", color: BLUE }, paragraph: { spacing: { before: 260, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 24, bold: true, font: "Arial", color: GREY }, paragraph: { spacing: { before: 160, after: 80 }, outlineLevel: 1 } },
    ],
  },
  numbering: { config: [
    { reference: "b", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] },
    { reference: "n", levels: [{ level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] },
  ] },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [
      new TextRun({ text: "FVS × PERSEUS handoff  —  page ", size: 16, color: GREY }),
      new TextRun({ children: [PageNumber.CURRENT], size: 16, color: GREY })] })] }) },
    children: [
      new Paragraph({ style: "Title", children: [new TextRun("CONUS Forest Vegetation Simulator → PERSEUS")] }),
      p([new TextRun({ text: "Integration handoff — calibrated FVS forest-carbon projections on the PERSEUS Forest Intelligence dashboard", italics: true, color: GREY })]),
      p([new TextRun({ text: "Aaron Weiskittel, University of Maine (Center for Research on Sustainable Forests). June 2026.", size: 20, color: GREY })]),
      p([new TextRun({ text: "Live: ", bold: true }),
        new ExternalHyperlink({ link: "https://holoros.github.io/perseus-forest-intelligence/", children: [new TextRun({ text: "holoros.github.io/perseus-forest-intelligence", style: "Hyperlink" })] }),
        new TextRun("   •   Data DOI: "),
        new ExternalHyperlink({ link: "https://doi.org/10.5281/zenodo.20555666", children: [new TextRun({ text: "10.5281/zenodo.20555666", style: "Hyperlink" })] })]),

      h("Executive summary", HeadingLevel.HEADING_1),
      p("A modernized, Bayesian-calibrated Forest Vegetation Simulator was run over the entire USDA FIA database (49 states, every plot, to 2125) and wired into the PERSEUS dashboard as three growth/mortality engines (default, calibrated, gompit), each under four management scenarios anchored to FIA carbon."),
      p("The most consequential moment was a QA catch before publication: the DataMart tree lists silently under-expanded overstory stem counts by about 6.5x, making the first campaign's biomass roughly 6x too light. That was diagnosed against the raw FIA TREE table, fixed, and the whole campaign re-run. On the corrected projections we added a data-driven harvest-plus-disturbance managed scenario with intensive management confined to FIADB plantations; a CONUS TreeMap spatially-explicit layer and a FIADB-vs-TreeMap multi-scale comparison; bootstrap trend breakdowns by landowner, ecoregion, state, and forest type; and a Bayesian posterior-draw parameter-uncertainty band."),
      p([new TextRun("Headline result: "), new TextRun({ text: "structural (engine-choice) uncertainty of 30 to 60% dominates parameter uncertainty of 0 to 18%", bold: true }), new TextRun(" — the mortality model matters far more than the calibrated coefficients. Everything is live, passes the project API-integrity check and a dedicated invariant stress test (0 violations), is archived on Zenodo with a DOI, and is reproducible from raw FIA via one pipeline script.")]),

      h("What is live on the dashboard", HeadingLevel.HEADING_1),
      table(["Engine", "Character"],
        [["default", "Native (Dixon/VARMRT) mortality; over-accumulates with no harvest"],
         ["calibrated", "Bayesian-calibrated growth; 13 to 32% below default; carries a posterior parameter band on 34+ states"],
         ["gompit", "Johnson national density-dependent mortality; caps and gently declines in late succession"]],
        [2200, 7160]),
      p(""),
      p([new TextRun({ text: "Scenarios per engine: ", bold: true }), new TextRun("reserve (no harvest); managed (conservation), the lightest; managed (harvest), the realistic case with intensive clearcut on plantations and extensive partial harvest on natural stands; managed (intensive), an all-intensive bound; plus the annual harvest carbon flux. Metrics: live aboveground carbon (Tg C) and dry biomass (Tg). All 49 forested states.")]),

      h("Key findings", HeadingLevel.HEADING_1),
      bullet("Treeinit expansion bug (the big one): the DataMart tree lists under-expand overstory ~6.5x, concentrated in eastern variants; caught by a QA check before publication; fixed against the raw FIA TREE table (TPA_UNADJ + heights)."),
      bullet("gompit caps carbon: density-dependent national mortality plateaus and gently declines late-succession biomass — the most realistic ceiling of the three engines."),
      bullet("Plantation-confined intensive management: intensive harvest only on FIADB plantations (STDORGCD=1, 10.2% of CONUS forest). Where plantations are common (GA 28%, OR 20%) the realistic managed path drops below the light scenario; where forests are natural (ME 2.4%) they coincide."),
      bullet("Uncertainty hierarchy: parameter uncertainty (posterior draws) is 0 to 18% at 2075, far inside the 30 to 60% structural spread between engines. The calibrated engine is nonetheless informative everywhere (13 to 32% below default)."),
      bullet("FIADB vs TreeMap across scales: the area-expansion choice is negligible at CONUS (~2%) but grows with resolution (state 0.36 to 1.39, forest type 0.23 to 1.96). FVS-on-TreeMap CONUS carbon cross-validates to within 9% of TreeMap's own imputed carbon."),

      h("The pipeline (raw FIA to dashboard)", HeadingLevel.HEADING_1),
      p([new TextRun({ text: "run_fvs_perseus_pipeline.sh", bold: true }), new TextRun(" chains every stage. The treeinit TPA + height fix is folded in at stage 1 so future campaigns start from correct tree lists.")]),
      numd("treeinit_fix_v2.py — restore TPA_UNADJ + complete heights from the raw FIA TREE table. Always first."),
      numd("build_plantation_flag.py + build_state_harvest_rates.R — plantation flag and conus_hcs harvest/disturbance rates sampled at plot locations."),
      numd("Campaign SLURM arrays — every FIA plot through its regional FVS variant, three engine arms, 100-yr projection, biomass via NSBE."),
      numd("fvs_perseus_aggregate.py — per-state density series with plot-percentile band."),
      numd("fvs_managed_v2.py — plantation-aware managed scenarios."),
      numd("fvs_posterior_uncertainty.py — Bayesian posterior-draw parameter CI (SLURM array over state x variant x draw)."),
      numd("Merge + ribbon + stress + push — on a clean origin/main checkout; fvs_dashboard_stress.py must report 0 and check_api_integrity.py must PASS before any push."),

      h("Deliverables", HeadingLevel.HEADING_1),
      table(["Deliverable", "Status / locator"],
        [["Dashboard engines", "3 engines x 4 scenarios x 49 states, live, integrity-checked"],
         ["Zenodo dataset", "10.5281/zenodo.20555666 (v1.0.0 .20555667), 21 files, CC-BY-4.0"],
         ["Reproducible pipeline", "run_fvs_perseus_pipeline.sh"],
         ["Trend products", "owner / ecoregion / state / forest-type, bootstrap CI"],
         ["FIADB-vs-TreeMap", "multi-scale comparison + CONUS spatial layer"],
         ["Uncertainty", "posterior parameter band (34+ states) + structural engine spread"],
         ["Documentation", "this handoff + 8 findings docs + roadmap (calibration/output/)"]],
        [2700, 6660]),

      h("Coexistence with the master DB-to-API pipeline", HeadingLevel.HEADING_1),
      p("The three FVS engines are direct-injected into the dashboard series JSON (not produced by the master database-to-API pipeline). They have survived the master regenerations to date because that regen is additive by model id, but a full series rebuild from the database would not recreate them — after any such rebuild, re-run pipeline stage 7 (merge + ribbon + stress). Both gates must pass before a push."),

      h("Remaining roadmap (none blocking)", HeadingLevel.HEADING_1),
      bullet("Dashboard UI cut to expose the landowner / ecoregion / forest-type breakdowns (data + CIs already computed)."),
      bullet("Harvested-wood-products carry-over: route the harvest carbon flux into an HWP pool for total-system carbon."),
      bullet("Stand-age / rotation-aware harvest: couple the owner-rotation logic so harvest responds to maturity."),
      bullet("Resolve the disturbance layer's temporal basis (currently annualized over an assumed 20-year window)."),
      bullet("Posterior CI for the multi-variant western states and the remaining sparse-forest states (full 49-state coverage)."),
      bullet("Climate-sensitive variant (current runs are climate-static, by agreed scope)."),
    ],
  }],
});

Packer.toBuffer(doc).then(b => { fs.writeFileSync(process.argv[2], b); console.log("wrote", process.argv[2], b.length, "bytes"); });
