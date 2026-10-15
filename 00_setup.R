# =============================================================================
# 00_setup.R
# -----------------------------------------------------------------------------
# Auto-locates the project folder (path/to/Ferro-aging) regardless of
# where RStudio's working directory currently points, switches to it, creates
# the standard folder layout, and prints a quick status report.
#
# Usage in RStudio (your working dir is path/to):
#   open this file from path/to/Ferro-aging, then run it. Every other
#   step script (step01 ... step09) expects the working directory to be this
#   project folder; source() them after 00_setup.R, or simply run 00_master.R.
#
# Alternative: Session -> Set Working Directory -> To Source File Location.
# =============================================================================

# ---- locate this script's own directory --------------------------------
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) > 0) {
    dirname(normalizePath(sub("^--file=", "", f[1])))
  } else {
    dirname(normalizePath(sys.frames()[[1]]$ofile))
  }
}, error = function(e) getwd())

PROJECT_DIR <- script_dir
setwd(PROJECT_DIR)

cat("== 00_setup ==============================================\n")
cat("Project folder :", PROJECT_DIR, "\n")
cat("Working dir now:", getwd(), "\n")

# ---- standard folders ---------------------------------------------------
dir.create("analysis_out/data", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis_out/data/GEO", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis_out/data/GSE149614", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis_out/hcc_out", showWarnings = FALSE, recursive = TRUE)

cat("Folder layout created under analysis_out/:\n")
cat("  data/            <- put all downloaded files here\n")
cat("  data/GEO/        <- step06 caches GEO downloads here\n")
cat("  data/GSE149614/  <- step08 single-cell files here\n")
cat("  hcc_out/         <- step05-09 results\n")
cat("==========================================================\n")
