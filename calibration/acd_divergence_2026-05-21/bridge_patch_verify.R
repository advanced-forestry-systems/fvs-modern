## bridge_patch_verify.R
##
## Quantify the error of the legacy bridge behaviour (applying the FVS basal
## area increment multiplier baimult directly to the linear Acadian diameter
## increment) and verify the proposed fix (convert through basal area).
##
## Legacy  : dDBH_new = dDBH * m
## Proposed: dDBH_new = sqrt(DBH^2 + m * ((DBH+dDBH)^2 - DBH^2)) - DBH
##           which makes the realized basal area increment exactly m times the
##           unscaled basal area increment, by construction.
##
## Runs in base R. No packages, no FVS.

ba_factor_legacy <- function(DBH, dDBH, m) {
  d2 <- (DBH + m * dDBH)^2 - DBH^2          # realized squared diameter increment
  d1 <- (DBH + dDBH)^2 - DBH^2              # unscaled
  d2 / d1                                   # realized basal area scaling
}

ddbh_proposed <- function(DBH, dDBH, m) {
  sqrt(DBH^2 + m * ((DBH + dDBH)^2 - DBH^2)) - DBH
}
ba_factor_proposed <- function(DBH, dDBH, m) {
  dn <- ddbh_proposed(DBH, dDBH, m)
  ((DBH + dn)^2 - DBH^2) / ((DBH + dDBH)^2 - DBH^2)
}

## grid over realistic Acadian sizes and annual increments, two multipliers
grid <- expand.grid(DBH = c(5, 10, 20, 30, 45, 60),
                    dDBH = c(0.1, 0.3, 0.5, 0.8),
                    m = c(0.80, 1.20))
grid$legacy_BA_factor   <- with(grid, ba_factor_legacy(DBH, dDBH, m))
grid$legacy_error_pct   <- 100 * (grid$legacy_BA_factor - grid$m) / grid$m
grid$proposed_BA_factor <- with(grid, ba_factor_proposed(DBH, dDBH, m))

cat("Intended basal area scaling is m. Legacy applies m to the linear diameter\n")
cat("increment, so its realized basal area scaling drifts from m, most for fast\n")
cat("growing small trees. Proposed conversion hits m exactly.\n\n")

## worst cases for each multiplier
for (mm in c(0.80, 1.20)) {
  sub <- grid[grid$m == mm, ]
  sub <- sub[order(-abs(sub$legacy_error_pct)), ]
  cat(sprintf("--- baimult m = %.2f : largest legacy errors ---\n", mm))
  print(head(data.frame(DBH = sub$DBH, dDBH = sub$dDBH,
                        legacy_BA = round(sub$legacy_BA_factor, 3),
                        legacy_err_pct = round(sub$legacy_error_pct, 1),
                        proposed_BA = round(sub$proposed_BA_factor, 3)), 4),
        row.names = FALSE)
  cat("\n")
}

cat(sprintf("max absolute legacy error across grid = %.1f%% of the intended factor\n",
            max(abs(grid$legacy_error_pct))))

## verify proposed is exact
stopifnot(max(abs(grid$proposed_BA_factor - grid$m)) < 1e-9)
cat("PASS: proposed conversion reproduces the intended basal area factor exactly.\n")
