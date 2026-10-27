## wrapper: run the v2 stress test on the fvs-modern posterior-median Bakuzis ensemble (real FIA stands, NE and ACD)
system2("Rscript", c("bakuzis_stress_v2.R",
  "--input=repos/fvs-modern/calibration/output/bakuzis/bakuzis_uncertainty_summary_long.csv",
  "--source=uncert", "--out=out_uncert", "--sdimax=sdimax_ne_engine_tph.csv", "--label=fvsmodern_fiaensemble_20261019"),
  stdout = "out_uncert/run_report.txt", stderr = "out_uncert/run_report.txt")
