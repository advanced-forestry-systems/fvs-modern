base <- "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/variants/ne"

cat("===== DG SAMPLES =====\n")
dg <- readRDS(file.path(base, "diameter_growth_samples.rds"))
cat("Class:", class(dg), "\n")
if (is.list(dg)) {
  cat("Names:", paste(names(dg), collapse=", "), "\n")
  for (nm in names(dg)) {
    obj <- dg[[nm]]
    if (is.numeric(obj) && length(obj) <= 20) {
      cat(sprintf("  %s [%s, len=%d]: %s\n", nm, class(obj), length(obj),
                  paste(round(head(obj, 10), 5), collapse=", ")))
    } else if (is.numeric(obj)) {
      cat(sprintf("  %s [%s, len=%d]: range=%.4f to %.4f, first5=%s\n", nm, class(obj), length(obj),
                  min(obj, na.rm=TRUE), max(obj, na.rm=TRUE),
                  paste(round(head(obj, 5), 5), collapse=", ")))
    } else if (is.character(obj)) {
      cat(sprintf("  %s [character, len=%d]: %s\n", nm, length(obj),
                  paste(head(obj, 15), collapse=", ")))
    } else if (is.data.frame(obj)) {
      cat(sprintf("  %s [data.frame, %d x %d]\n", nm, nrow(obj), ncol(obj)))
      cat("    Columns:", paste(names(obj), collapse=", "), "\n")
      print(head(obj, 3))
    } else if (is.matrix(obj)) {
      cat(sprintf("  %s [matrix, %d x %d]\n", nm, nrow(obj), ncol(obj)))
      cat("    Colnames:", paste(head(colnames(obj), 10), collapse=", "), "\n")
      cat("    Rownames:", paste(head(rownames(obj), 10), collapse=", "), "\n")
    } else {
      cat(sprintf("  %s [%s]\n", nm, paste(class(obj), collapse="/")))
    }
  }
} else {
  str(dg, max.level=2)
}

cat("\n===== MORTALITY SAMPLES =====\n")
m <- readRDS(file.path(base, "mortality_samples.rds"))
cat("Class:", class(m), "\n")
if (is.list(m)) {
  cat("Names:", paste(names(m), collapse=", "), "\n")
  for (nm in names(m)) {
    obj <- m[[nm]]
    if (is.numeric(obj) && length(obj) <= 50) {
      cat(sprintf("  %s [%s, len=%d]: %s\n", nm, class(obj), length(obj),
                  paste(round(head(obj, 10), 5), collapse=", ")))
      if (length(obj) > 1 && !is.null(names(obj))) {
        cat("    Names:", paste(head(names(obj), 15), collapse=", "), "\n")
      }
    } else if (is.numeric(obj)) {
      cat(sprintf("  %s [%s, len=%d]: range=%.4f to %.4f\n", nm, class(obj), length(obj),
                  min(obj, na.rm=TRUE), max(obj, na.rm=TRUE)))
    } else if (is.character(obj)) {
      cat(sprintf("  %s [character, len=%d]: %s\n", nm, length(obj),
                  paste(head(obj, 15), collapse=", ")))
    } else {
      cat(sprintf("  %s [%s]\n", nm, paste(class(obj), collapse="/")))
    }
  }
}

cat("\n===== HD SAMPLES =====\n")
hd <- readRDS(file.path(base, "height_diameter_samples.rds"))
cat("Class:", class(hd), "\n")
if (is.list(hd)) {
  cat("Names:", paste(names(hd), collapse=", "), "\n")
  for (nm in names(hd)) {
    obj <- hd[[nm]]
    if (is.numeric(obj) && length(obj) <= 50) {
      cat(sprintf("  %s [%s, len=%d]: %s\n", nm, class(obj), length(obj),
                  paste(round(head(obj, 10), 5), collapse=", ")))
      if (length(obj) > 1 && !is.null(names(obj))) {
        cat("    Names:", paste(head(names(obj), 15), collapse=", "), "\n")
      }
    } else if (is.numeric(obj)) {
      cat(sprintf("  %s [%s, len=%d]: range=%.4f to %.4f\n", nm, class(obj), length(obj),
                  min(obj, na.rm=TRUE), max(obj, na.rm=TRUE)))
    } else {
      cat(sprintf("  %s [%s]\n", nm, paste(class(obj), collapse="/")))
    }
  }
}

cat("\n===== STANDARDIZATION PARAMS =====\n")
sp <- readRDS(file.path(base, "standardization_params.rds"))
cat("Class:", class(sp), "\n")
if (is.data.frame(sp)) {
  print(sp)
} else if (is.list(sp)) {
  str(sp, max.level=2)
  for (nm in names(sp)) {
    if (is.data.frame(sp[[nm]])) {
      cat(nm, ":\n")
      print(sp[[nm]])
    }
  }
}

cat("\n===== CR SAMPLES =====\n")
cr <- readRDS(file.path(base, "crown_ratio_samples.rds"))
cat("Class:", class(cr), "\n")
cat("Names:", paste(names(cr), collapse=", "), "\n")

cat("\n===== TRAINING DATA COLUMNS =====\n")
data_base <- "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data"
for (f in c("diameter_growth.csv", "mortality.csv", "height_diameter.csv")) {
  fp <- file.path(data_base, "ne", f)
  if (file.exists(fp)) {
    d <- read.csv(fp, nrows=2)
    cat(f, "columns:", paste(names(d), collapse=", "), "\n")
  } else {
    # Try uppercase
    fp2 <- file.path(data_base, "NE", f)
    if (file.exists(fp2)) {
      d <- read.csv(fp2, nrows=2)
      cat(f, "(NE) columns:", paste(names(d), collapse=", "), "\n")
    } else {
      cat(f, "- NOT FOUND\n")
    }
  }
}

# Check what data dirs exist
cat("\nData subdirs:", paste(list.files(data_base), collapse=", "), "\n")

cat("\nDone.\n")
