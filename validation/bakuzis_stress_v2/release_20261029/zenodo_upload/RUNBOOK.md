# Runbook: fvs-modern configuration catalog record (October 29, 2026)

1. Build: `python3 code/release_figures.py --data data --out fig`; `python3 code/build_report_docx.py --md report/<md> --out report/<docx>`; `python3 code/build_deck_pptx.py --fig fig --sc <scorecard figures> --out deck/<pptx>`; LibreOffice headless to PDF for both. Re-render, then gate, then upload, nothing in between.
2. Gates: `cd zenodo_upload && python3 check_deposit.py` (writes SHA256SUMS.txt; exit non-zero blocks the upload). Gate 3 sweeps every text file for coordinate or plot CN columns, restricted tokens and contextual coordinate pairs; .png, .pdf, .docx and .pptx are reported as not text-scannable. Gate 5 reads the scripts for the write call that names each table; `fia_variant_thresholds_plot_acd.csv` is declared as a compiled reference table whose producing script reads the restricted pairs file and is not deposited.
3. Draft upload from firebreather (token at `~/.config/zenodo/token`, search order per zenodo-deposit; never printed): `setsid nohup python3 upload_to_zenodo.py --token-file ~/.config/zenodo/token --metadata zenodo_metadata.json --files-list files_to_upload.txt > upload.log 2>&1 &`, then poll `upload.log`. No `--publish`: the record stays a draft.
4. Read the reserved DOI from the log, re-run `release_figures.py --doi <doi>` so the graphical abstract carries it, add the DOI to CITATION.cff and README, re-run the gates, re-upload the changed files (PUT by name overwrites in a draft).
5. Publishing is Aaron's action in the Zenodo web interface after review of the draft.
6. Backfill the concept DOI: fvs-modern `validation/bakuzis_stress_v2/release_20261029/README.md`, the memo's data availability line, Paper 1 and Paper 2 data availability statements.

Restricted data audit result is printed by `check_deposit.py` at every run and copied here after the final run.
