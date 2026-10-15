# =============================================================================
# 00_master.R  — one-click runner for the whole HCC pipeline
# -----------------------------------------------------------------------------
# Opens in RStudio from path/to/Ferro-aging. It re-locates the working
# directory automatically, then sources the pipeline steps in order. Steps
# that need downloaded data print instructions and exit gracefully instead of
# crashing, so the master script always runs to the end.
#
#   step01 : primate FARG tissue profile       (your xlsx, no downloads)
#   step02 : monkey -> human orthologs         (optional network)
#   step04 : TCGA-LIHC data prep               (needs Xena downloads)
#   step05 : FARGs signature construction
#   step06 : external validation
#   step07 : immune / drug / TMB
#   step08 : single-cell + virtual knockout   (needs GSE149614 files)
#   step09 : robustness suite                  (bootstrap / comparison /
#                                               subgroups / GSEA / forest)
# =============================================================================

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) > 0) {
    dirname(normalizePath(sub("^--file=", "", f[1])))
  } else {
    dirname(normalizePath(sys.frames()[[1]]$ofile))
  }
}, error = function(e) getwd())
setwd(script_dir)
cat("== 00_master: HCC ferro-aging pipeline ==================\n")
cat("Working dir:", getwd(), "\n")

STEPS <- c("step01_farg_profile.R", "step02_human_map.R",
           "step04_hcc_prep.R", "step05_hcc_signature.R",
           "step06_hcc_validate.R", "step07_hcc_immune_drug.R",
           "step08_hcc_singlecell.R", "step09_hcc_robustness.R",
           "step10_hcc_extras.R", "step12_hcc_aging.R")

t_start <- Sys.time()
for (s in STEPS) {
  if (!file.exists(s)) {
    cat("\n[skip]", s, "- file not found\n")
    next
  }
  cat("\n##########", s, "##########\n")
  ok <- tryCatch({ source(s); TRUE }, error = function(e) {
    cat("!! ", s, "ERROR:", conditionMessage(e), "\n"); FALSE
  })
  if (!ok) {
    cat("!! continuing with next step ...\n")
    if (s == "step04_hcc_prep.R")
      cat("!! NOTE: steps 05-09 depend on hcc_data.rds; without it they run\n",
          "!!       SELF-TEST on synthetic data and write to hcc_out/self_test/.\n",
          "!!       Place the Xena files and re-run to get real results.\n")
  }
}
cat("\n== 00_master finished in",
    round(difftime(Sys.time(), t_start, units = "mins"), 1), "min ==\n")
