#!/usr/bin/env Rscript
# Biological-plausibility gate for the CORRECTED GOMPSURV (issue #75 fix).
# Unlike the self-consistency check, this asserts the predictions are sane:
# a vigorous open-grown tree survives at a high annual rate, and survival
# declines under competition for the large majority of species.
fo <- read.csv("surv2_fortran.csv")
n <- nrow(fo)
med_b <- median(fo$s_benign_ann); med_s <- median(fo$s_stress_ann)
q05_b <- quantile(fo$s_benign_ann, 0.05)
n_implaus  <- sum(fo$s_benign_ann < 0.90)
n_wrongdir <- sum(fo$s_stress_ann > fo$s_benign_ann)
cat(sprintf("species: %d\n", n))
cat(sprintf("benign annual survival: median %.4f  q05 %.4f\n", med_b, q05_b))
cat(sprintf("stress annual survival: median %.4f\n", med_s))
cat(sprintf("benign < 0.90 (implausible): %d / %d (%.0f%%)\n", n_implaus, n, 100*n_implaus/n))
cat(sprintf("stress > benign (wrong direction): %d / %d (%.0f%%)\n", n_wrongdir, n, 100*n_wrongdir/n))
pass <- (med_b > 0.95) && (n_implaus < 0.15*n) && (n_wrongdir < 0.15*n)
cat(if (pass) "PASS: corrected GOMPSURV is biologically plausible (issue #75 resolved at unit level)\n"
    else "FAIL: still implausible -- investigate\n")
