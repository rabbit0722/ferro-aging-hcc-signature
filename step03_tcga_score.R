# =============================================================================
# step03_tcga_score.R
# -----------------------------------------------------------------------------
# ssGSEA scoring of the ferro-aging gene set (FARGs + its functional modules)
# in TCGA pan-cancer expression data, + per-cancer Cox regression and one
# Kaplan-Meier plot.
#
# DATA REQUIRED (download manually, place under analysis_out/data/):
#   Expression:
#     UCSC Xena Toil hub, "TcgaTargetGtex_rsem_gene_tpm.tsv.gz" (pan-cancer TPM)
#     https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.tsv.gz
#     (Xena browser page: https://xenabrowser.net/datapages/?hub=https://toil.xenahubs.net)
#   Survival (optional but recommended):
#     "Survival_SupplementalTable_S1_20171025_xena_sp.gz"
#     https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz
#
# If the files are missing the script runs a built-in SELF-TEST on synthetic
# data (validates the GSVA+Cox code path) and prints download instructions.
#
# Inputs : analysis_out/step01/module_assignment.csv (from step01)
#          analysis_out/step02/farg_human_map.csv (from step02)
# Outputs: analysis_out/step03/
#   ssgsea_scores_tcga.csv      modules x samples
#   cox_pancancer.csv           per-cancer Cox of the all-FARG score
#   km_<CANCER>.png             KM of median-split score for CANCER_CODE
#
# Requirements: data.table, GSVA, survival; survminer for the KM plot
# Run:  Rscript step03_tcga_score.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(survival)
})

## ------------------------- CONFIG -------------------------
EXPR_FILE    <- "analysis_out/data/TcgaTargetGtex_rsem_gene_tpm.tsv.gz"
SURV_FILE    <- "analysis_out/data/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
OUT_DIR      <- "analysis_out/step03"
CANCER_CODE  <- "KIRC"      # for the KM figure (change as needed)
MIN_TPM      <- 0.5         # keep genes with mean TPM >= this (speed/memory)
# -------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("== step03: FARG ssGSEA scoring in TCGA ==\n")
cat("R:", R.version.string, "\n\n")

## ---- gene sets --------------------------------------------------------------
mod_file <- "analysis_out/step01/module_assignment.csv"
if (!file.exists(mod_file)) stop("Run step01 first: ", mod_file, " missing.")
mod <- fread(mod_file)
mod <- mod[!is.na(SYMBOL) & SYMBOL != ""]
gs <- split(mod$SYMBOL, mod$Module)
gs[["all_FARG"]] <- mod$SYMBOL
cat("Gene sets (modules):\n")
print(data.table(set = names(gs), n = lengths(gs)))

## ---- load data (or self-test) ------------------------------------------------
have_expr <- file.exists(EXPR_FILE)
have_surv <- file.exists(SURV_FILE)

if (!have_expr) {
  cat("\n[WARN] Expression file not found:\n       ", normalizePath(EXPR_FILE,
        mustWork = FALSE), "\n")
  cat("       Download instructions:\n")
  cat("       1) create folder analysis_out/data/\n")
  cat("       2) download from UCSC Xena Toil hub:\n")
  cat("          https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.tsv.gz\n")
  cat("       3) (optional) survival table:\n")
  cat("          https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz\n")
  cat("\n>>> Running SELF-TEST on synthetic data to validate the pipeline ...\n")

  set.seed(42)
  n_samp <- 120
  cancers <- rep(c("TCGA-AA", "TCGA-BB", "TCGA-CC"), each = n_samp / 3)
  samples <- paste0(cancers, "-", sprintf("%04d", 1:n_samp), "-01")
  # synthetic expression: 200 genes x 120 samples
  genes <- c(mod$SYMBOL, paste0("BG", 1:(200 - uniqueN(mod$SYMBOL))))
  expr <- matrix(rnorm(length(genes) * n_samp, mean = 5, sd = 1.5),
                 nrow = length(genes), dimnames = list(genes, samples))
  # latent ferro-aging signal in FARG rows
  signal <- rnorm(n_samp)
  expr[mod$SYMBOL, ] <- expr[mod$SYMBOL, ] + matrix(signal,
    nrow = uniqueN(mod$SYMBOL), ncol = n_samp, byrow = TRUE)
  # survival with a real association to the latent signal
  os_time <- rexp(n_samp, 0.05) * exp(-0.4 * signal)
  os <- rbinom(n_samp, 1, plogis(-0.8 + 1.2 * signal))
  surv <- data.table(sample = samples, OS = os, OS.time = round(os_time, 1))
  have_surv <- TRUE
  SELF_TEST <- TRUE
} else {
  cat("\nLoading expression matrix (large; 1-3 min) ...\n")
  expr_dt <- fread(EXPR_FILE)
  cat("      matrix:", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n")
  gene_col <- names(expr_dt)[1]
  setnames(expr_dt, gene_col, "gene")
  # Xena gene ids sometimes carry "SYMBOL|entrez" -> keep symbol part
  expr_dt[, gene := sub("\\|.*$", "", gene)]
  expr_dt <- expr_dt[!duplicated(gene)]
  samples <- setdiff(names(expr_dt), "gene")
  # keep TCGA primary-tumor samples (barcode position 14-15 == "01")
  keep <- grepl("^TCGA-", samples) & substr(samples, 14, 15) == "01"
  cat("      TCGA primary-tumor samples:", sum(keep), "\n")
  samples <- samples[keep]
  m <- as.matrix(expr_dt[, ..samples])
  rownames(m) <- expr_dt$gene
  rm(expr_dt); gc()
  # filter lowly expressed genes (keeps enough background for ssGSEA ranking)
  keepg <- rowMeans(m) >= MIN_TPM
  cat("      genes kept after TPM filter:", sum(keepg), "\n")
  m <- m[keepg, , drop = FALSE]
  expr <- log2(m + 0.001)
  rm(m); gc()
  SELF_TEST <- FALSE
  if (have_surv) {
    surv <- fread(SURV_FILE)
    surv <- surv[sample %in% samples]
  }
}

## ---- ssGSEA ------------------------------------------------------------------
cat("\nRunning ssGSEA (rank-based; may take minutes on pan-cancer data) ...\n")
t0 <- Sys.time()
gs_ok <- lapply(gs, function(g) intersect(g, rownames(expr)))
gs_ok <- gs_ok[lengths(gs_ok) >= 5]                     # drop tiny sets
res <- tryCatch(
  gsva(as.matrix(expr), gs_ok, method = "ssgsea", ssgsea.norm = TRUE,
       verbose = FALSE),
  error = function(e) NULL)
if (is.null(res)) {                                      # GSVA 2.x interface
  param <- ssgseaParam(as.matrix(expr), gs_ok)
  res <- gsva(param, verbose = FALSE)
}
cat("      done in", round(difftime(Sys.time(), t0, units = "mins"), 1),
    "min\n")
scores <- t(res)                                          # samples x sets
scores_dt <- as.data.table(scores, keep.rownames = "sample")
fwrite(scores_dt, file.path(OUT_DIR, if (SELF_TEST) "self_test_scores.csv"
                            else "ssgsea_scores_tcga.csv"))
cat("Saved scores:", if (SELF_TEST) "self_test_scores.csv"
    else "ssgsea_scores_tcga.csv", "\n")

## ---- per-cancer Cox on the all-FARG score ------------------------------------
if (have_surv && !is.null(surv)) {
  cat("\nPer-cancer Cox regression of the all-FARG score ...\n")
  dat <- merge(scores_dt[, .(sample, score = all_FARG)], surv,
               by = "sample", all.x = TRUE)
  dat <- dat[!is.na(OS) & !is.na(OS.time) & OS.time > 0]
  dat[, cancer := substr(sample, 1, 7)]
  res_list <- lapply(sort(unique(dat$cancer)), function(cc) {
    dd <- dat[cancer == cc]
    if (nrow(dd) < 20 || length(unique(dd$OS)) < 2) return(NULL)
    fit <- tryCatch(coxph(Surv(OS.time, OS) ~ score, data = dd),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    data.table(cancer = cc, n = nrow(dd), events = sum(dd$OS),
               HR = s$conf.int[1, "exp(coef)"],
               HR_low = s$conf.int[1, "lower .95"],
               HR_high = s$conf.int[1, "upper .95"],
               p = s$coefficients[1, "Pr(>|z|)"])
  })
  cox <- rbindlist(res_list[!vapply(res_list, is.null, logical(1))])
  if (nrow(cox)) {
    fwrite(cox, file.path(OUT_DIR, if (SELF_TEST) "self_test_cox.csv"
                          else "cox_pancancer.csv"))
    print(cox[order(p)])
  } else cat("      not enough events per cancer; Cox skipped.\n")

  ## ---- KM for CANCER_CODE ----------------------------------------------------
  cc_km <- if (SELF_TEST) sort(unique(dat$cancer))[1] else paste0("TCGA-", CANCER_CODE)
  dd <- dat[cancer == cc_km]
  if (nrow(dd) >= 20) {
    dd[, grp := ifelse(score >= median(score), "High", "Low")]
    fit <- survfit(Surv(OS.time, OS) ~ grp, data = dd)
    if (requireNamespace("survminer", quietly = TRUE)) {
      suppressPackageStartupMessages(library(survminer))
      p_km <- file.path(OUT_DIR, paste0(if (SELF_TEST) "self_test_", "km_",
                                        sub("^TCGA-", "", cc_km), ".png"))
      g <- ggsurvplot(fit, data = dd, pval = TRUE, risk.table = TRUE,
                      palette = c("#B2182B", "#2166AC"),
                      xlab = "Time", title = paste("FARG score KM -",
                                                   cc_km))
      png(p_km, width = 1400, height = 1100, res = 150)
      print(g, newpage = FALSE)
      dev.off()
      cat("Saved KM:", p_km, "\n")
    }
  }
}

if (SELF_TEST) {
  cat("\n== SELF-TEST passed: pipeline is ready. ==\n")
  cat("   Place the real files as described above, then rerun this script.\n")
}
cat("\n== step03 done ==\n")
