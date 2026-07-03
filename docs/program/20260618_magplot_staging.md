# Staging MAGPlot to unblock ACD and AK calibration
2026-06-18. The Canadian NFI / MAGPlot data package must be on Cardinal before ACD and AK calibrate. The
assistant cannot fetch URLs in its environment, so run this once on Cardinal:

    mkdir -p /fs/scratch/PUOM0008/crsfaaron/MAGPlot && cd /fs/scratch/PUOM0008/crsfaaron/MAGPlot
    wget -O MAGPlot_latest_data_package.zip "https://ca.nfis.org/fss/fss?command=retrieveByName&fileName=MAGPlot_latest_data_package.zip&fileNameSpace=magplot&format=xml&promptToSave=true"
    unzip -o MAGPlot_latest_data_package.zip

Then the turnkey ingester (diagnostics_2026-06-16/magplot_ingest.py) inspects the tables, detects the
site and tree tables and key columns, filters to the jurisdictions that map to the variants (NB/NS/PE to
Acadian, coastal BC to the AK variant), and reports remeasurement candidates and a schema manifest:

    module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
    python3 ~/fvs-modern/diagnostics_2026-06-16/magplot_ingest.py /fs/scratch/PUOM0008/crsfaaron/MAGPlot

From the manifest the pair-builder is finalized (t1/t2 per site, live DBH), MAGPlot species are crosswalked
to FVS species via the data dictionary, and the per-variant calibration runs as for the US variants. OGL
Canada licence (attribution). NL and SK are not in the open release.
