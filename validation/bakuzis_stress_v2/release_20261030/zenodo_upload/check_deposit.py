#!/usr/bin/env python3
"""check_deposit.py: the six pre-publish gates for the fvs-modern catalog Zenodo record. Exit non-zero on any failure.
Run from the zenodo_upload folder: python3 check_deposit.py [--skip-freshness]
Gate 1 lead file; 2 source freshness (mtime); 3 restricted data audit (no override); 4 SHA-256 of every staged file
(written to SHA256SUMS.txt); 5 producers verified by reading the scripts; 6 metadata.
"""
import os, re, sys, json, hashlib, csv, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--skip-freshness", action="store_true"); a = ap.parse_args()
fails = []
def fail(m): fails.append(m); print("FAIL", m)
def ok(m): print("ok  ", m)

# manifest
man = []
for ln in open(os.path.join(HERE, "files_to_upload.txt"), encoding="utf8"):
    ln = ln.split("#", 1)[0].strip()
    if ln: man.append(os.path.normpath(os.path.join(HERE, ln)))
missing = [p for p in man if not os.path.isfile(p)]
if missing: fail("manifest paths missing: " + ", ".join(missing))
names = [os.path.basename(p) for p in man]
if len(set(names)) != len(names): fail("duplicate basenames in manifest")
meta = json.load(open(os.path.join(HERE, "zenodo_metadata.json")))["metadata"]

# gate 1
lead = names[0] if names else ""
if lead.startswith("00_") and lead.endswith(".png") and meta.get("default_preview") == lead and sorted(names)[0] == lead: ok(f"gate 1 lead file {lead}")
else: fail("gate 1 lead file is not first, not 00_ prefixed, or not default_preview")

# gate 2 freshness
def mt(p): return os.path.getmtime(p)
R = os.path.join(HERE, "..")
pairs = [("report/2026-10-29_fvs-modern-catalog-report_DRAFT.docx", "report/2026-10-29_fvs-modern-catalog-report_DRAFT.md"),
         ("report/2026-10-29_fvs-modern-catalog-report_DRAFT.pdf", "report/2026-10-29_fvs-modern-catalog-report_DRAFT.docx"),
         ("deck/2026-10-29_fvs-modern-catalog-deck_DRAFT.pptx", "code/build_deck_pptx.py"),
         ("deck/2026-10-29_fvs-modern-catalog-deck_DRAFT.pdf", "deck/2026-10-29_fvs-modern-catalog-deck_DRAFT.pptx")]
for f in ["00_fvs-modern_catalog_graphical_abstract.png", "fig_variant_stability.png", "fig_gate_space.png", "fig_arm_stability.png"]:
    pairs += [(f"fig/{f}", "code/release_figures.py"), (f"fig/{f}", "data/stress_by_arm.csv"), (f"fig/{f}", "data/catalog_gates_by_arm.csv")]
if a.skip_freshness: print("gate 2 skipped by flag")
else:
    bad = [(x, y) for x, y in pairs if mt(os.path.join(R, x)) < mt(os.path.join(R, y))]
    if bad: fail("gate 2 stale renders: " + "; ".join(f"{x} older than {y}" for x, y in bad))
    else: ok(f"gate 2 freshness, {len(pairs)} pairs")

# gate 3 restricted data audit, every file, every row
HDR = re.compile(r"^(lat|lon|long|latitude|longitude|lat_dd|lon_dd|plot_lat|plot_lon|lat_wgs84|lon_wgs84|x_coord|y_coord|point_x|point_y|geometry|wkt|utm_e|utm_n|lat_actual|lon_actual|actual_lat|actual_lon|true_lat|true_lon|lat_fuzz|lon_fuzz|plt_cn|cn)$", re.I)
TOK = re.compile(r"_SERVERSIDE|FIA\.xy|per_plot_embedding|plot_coords|LAT_ACTUAL|LON_ACTUAL|ACTUAL_LAT|ACTUAL_LON|TRUE_LAT|TRUE_LON|LAT_FUZZ|LON_FUZZ", re.I)
PAIR = re.compile(r"(?<![\w.])(2[4-9]|[3-4]\d)\.\d{3,}(?![\w.])\s*[, ]\s*(?<![\w.])-?(6[5-9]|[7-9]\d|1[0-2]\d)\.\d{3,}(?![\w.])")
CTX = re.compile(r"lat|lon|coord", re.I)
TEXT_EXT = {".csv", ".md", ".R", ".py", ".txt", ".json", ".cff", ".slurm"}
scanned = 0; unscannable = []; chars = 0; waived = []
for p in man:
    ext = os.path.splitext(p)[1]
    if ext not in TEXT_EXT: unscannable.append(os.path.basename(p)); continue
    scanned += 1; txt = open(p, encoding="utf8", errors="replace").read(); chars += len(txt)
    if ext == ".csv":
        hdr = next(csv.reader([txt.splitlines()[0]]))
        hits = [h for h in hdr if HDR.match(h.strip())]
        if hits: fail(f"gate 3 coordinate or CN column in {os.path.basename(p)}: {hits}")
        for i, row in enumerate(txt.splitlines()):
            if TOK.search(row): fail(f"gate 3 restricted token in {os.path.basename(p)} row {i}")
    else:
        for m in TOK.finditer(txt):
            if ext in (".py", ".R", ".md"): waived.append(f"{os.path.basename(p)}: provenance token {m.group(0)}")
            else: fail(f"gate 3 restricted token {m.group(0)} in {os.path.basename(p)}")
        for m in PAIR.finditer(txt):
            win = txt[max(0, m.start() - 40): m.end() + 40]
            if CTX.search(win): fail(f"gate 3 coordinate pair with context in {os.path.basename(p)}: {m.group(0)}")
print(f"gate 3 audit {time.strftime('%Y-%m-%d %H:%M')}: {scanned} text files scanned, {chars} characters, {len(unscannable)} not text-scannable: {unscannable}; waived provenance hits: {waived or 'none'}")

# gate 4 checksums
with open(os.path.join(HERE, "SHA256SUMS.txt"), "w") as fh:
    for p in man:
        h = hashlib.sha256(open(p, "rb").read()).hexdigest(); fh.write(f"{h}  {os.path.basename(p)}\n")
ok(f"gate 4 SHA-256 written for {len(man)} files at {time.strftime('%Y-%m-%d %H:%M:%S')}")

# gate 5 producers, verified by reading the scripts
WRITE = re.compile(r"(?:fwrite|write\.csv|to_csv|write_csv|writeLines|savefig)\s*\(", re.I)
declared = {"catalog_seedgate_v6_noid.csv": "seed gate rows written by conus_catalog_runs.py (runs_v5 default rows and runs_v6 calibrated rows, concatenated by concat_catalog_v6.R in finalize_catalog_v7.slurm) with the plot identifier column removed before deposit",
            "fia_variant_thresholds_plot_acd.csv": "compiled reference table: plot level percentiles of the fvs-conus fitting data, script not deposited (restricted input)",
            "scorecard_20261029_default_native_selfthinning.png": "bakuzis_scorecards.R", "scorecard_20261029_default_gompit_NE_bakuzis_matrix.png": "bakuzis_scorecards.R"}
code = {os.path.basename(p): open(p, encoding="utf8", errors="replace").read() for p in man if os.path.splitext(p)[1] in (".R", ".py")}
for p in man:
    b = os.path.basename(p)
    if os.path.splitext(b)[1] not in (".csv", ".png") or b == "data_dictionary.csv": continue
    prod = [s for s, src in code.items() if WRITE.search(src) and b in src]
    if prod: ok(f"gate 5 {b} <- {prod}")
    elif b in declared: print(f"decl gate 5 {b}: {declared[b]}")
    else: fail(f"gate 5 no producer for {b}")

# data dictionary covers every column
dd = {}
for r in csv.DictReader(open(os.path.join(HERE, "data_dictionary.csv"))): dd.setdefault(r["file"], set()).add(r["column"])
for p in man:
    b = os.path.basename(p)
    if b.endswith(".csv") and b in dd:
        hdr = set(next(csv.reader(open(p, encoding="utf8"))))
        if hdr != dd[b]: fail(f"data dictionary mismatch for {b}: missing {hdr - dd[b]}, extra {dd[b] - hdr}")

# gate 6 metadata
req = ["upload_type", "title", "version", "publication_date", "creators", "description", "license", "keywords"]
for k in req:
    if not meta.get(k): fail(f"gate 6 missing {k}")
for c in meta["creators"]:
    if ", " not in c["name"]: fail("gate 6 creator not Family, Given")
if any("XXXX" in json.dumps(r) or "NNNN" in json.dumps(r) for r in meta.get("related_identifiers", [])): fail("gate 6 placeholder DOI")
if any(k != k.strip() for k in meta["keywords"]): fail("gate 6 whitespace in keyword")
if meta.get("default_preview") not in names: fail("gate 6 default_preview not in manifest")
print("gate 6 related identifiers:", [r["identifier"] for r in meta.get("related_identifiers", [])])
print("\nRESULT:", "PASS" if not fails else f"{len(fails)} FAILURE(S)")
sys.exit(1 if fails else 0)
