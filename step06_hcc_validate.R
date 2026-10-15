# =============================================================================
# step06_hcc_validate.R
# -----------------------------------------------------------------------------
# External validation of the FARGs signature built in step05.
# Applies the frozen coefficients to independent cohorts:
#   GSE14520 (Affymetrix, n=242, survival)   <- auto-download via GEOquery
#   GSE16757 (Illumina,   n=100, survival)   <- matrix + clinical supplement
#
# For each cohort: map probes to symbols -> z-score within cohort ->
# risk = sum(coef * z) -> KM (median split) + timeROC + Cox HR per SD.
# Finally a forest plot combining training + validation HRs.
#
# Input : analysis_out/hcc_out/signature_genes.csv   (from step05)
#         analysis_out/hcc_out/risk_training.csv     (from step05)
# If the signature file or GEO download is unavailable the script runs a
# SELF-TEST on synthetic validation cohorts (same code path).
#
# Output: analysis_out/hcc_out/validate/
#   risk_<cohort>.csv, km_<cohort>.png, auc_<cohort>.png, forest_hr.png
#
# Requirements: data.table, GEOquery, survival, survminer, timeROC, ggplot2
# Run:  Rscript step06_hcc_validate.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

SIG_FILE    <- "analysis_out/hcc_out/signature_genes.csv"
RISK_TRAIN  <- "analysis_out/hcc_out/risk_training.csv"
OUT_DIR     <- "analysis_out/hcc_out/validate"
GEO_CACHE   <- "analysis_out/data/GEO"
SEED        <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GEO_CACHE, showWarnings = FALSE, recursive = TRUE)
cat("== step06: external validation of the FARGs signature ==\n")

## ---- signature coefficients ---------------------------------------------------
if (file.exists(SIG_FILE)) {
  sig <- fread(SIG_FILE)
  sig <- sig[!is.na(gene) & !is.na(coef)]
  cat("Signature loaded from step05:", paste(sig$gene, collapse = ", "), "\n")
} else {
  cat("[WARN] signature_genes.csv not found -> demo coefficients (SELF-TEST).\n")
  sig <- data.table(gene = c("ALOX15", "ATG3", "BAP1"),
                    coef = c(0.44, 0.37, 0.51))
}
coefs <- setNames(sig$coef, sig$gene)

## ---- helper: probe -> symbol mapping -------------------------------------------
probe2symbol <- function(expr, fdat) {
  # case-insensitive detection of the symbol column (GEO annot files use
  # "Gene Symbol", "Gene symbol", "Symbol" ... inconsistently)
  cols <- colnames(fdat)
  cand <- cols[grepl("gene.?symbol|^symbol$", cols, ignore.case = TRUE)]
  if (length(cand) == 0) {
    cand <- intersect(c("gene_assignment", "ILMN_Gene"), cols)
    if (length(cand) == 0) stop("no symbol column found in fData: ",
                                paste(cols, collapse = ", "))
  }
  col <- cand[1]
  sym <- fdat[[col]]
  if (col == "gene_assignment") {
    sym <- sub(".*//\\s*([A-Za-z0-9\\.\\-]+)\\s*//.*", "\\1", sym)
  }
  rownames(expr) <- sym
  ok <- !is.na(sym) & sym != "" & sym != "---"
  expr <- expr[ok, , drop = FALSE]
  # collapse multiple probes per symbol: keep the highest-mean row
  if (anyDuplicated(rownames(expr))) {
    keep <- !duplicated(rownames(expr), fromLast =
                          ave(rowMeans(expr), rownames(expr),
                              FUN = function(x) rank(x)) == 1)
    expr <- expr[order(rowMeans(expr), decreasing = TRUE), , drop = FALSE]
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  }
  expr
}

## ---- helper: auto-detect survival columns --------------------------------------
detect_surv <- function(pd, time_col = NA, status_col = NA, time_unit = "auto") {
  cols <- colnames(pd)
  if (is.na(time_col)) {
    cand <- grep("month|day|time|survival|follow|year", cols, ignore.case = TRUE,
                 value = TRUE)
    cand <- cand[!grepl("status|event|alive|dead", cand, ignore.case = TRUE)]
    if (length(cand) == 0) stop("cannot auto-detect time column; set time_col manually. Columns: ",
                                paste(cols, collapse = ", "))
    time_col <- cand[1]
  }
  if (is.na(status_col)) {
    cand <- grep("status|event|alive|dead|censor", cols, ignore.case = TRUE,
                 value = TRUE)
    cand <- cand[!grepl("month|day|time|year", cand, ignore.case = TRUE)]
    status_col <- if (length(cand)) cand[1] else NA_character_
  }
  t <- suppressWarnings(as.numeric(as.character(pd[[time_col]])))
  s <- if (!is.na(status_col)) as.character(pd[[status_col]]) else rep("1", nrow(pd))
  # normalize event indicator
  s_num <- suppressWarnings(as.numeric(s))
  event <- ifelse(!is.na(s_num), s_num, ifelse(grepl("dead|1", s, ignore.case = TRUE), 1, 0))
  event <- pmin(pmax(event, 0), 1)
  if (time_unit == "auto") {
    time_unit <- ifelse(grepl("month", time_col, ignore.case = TRUE),
                        "months", ifelse(grepl("year", time_col, ignore.case = TRUE),
                                         "years", "days"))
  }
  mult <- switch(time_unit, months = 30.44, years = 365.25, days = 1)
  data.table(time = t * mult, status = event)
}

# per-cohort survival column overrides (platform-specific pData column names)
SURV_OVERRIDES <- list()

## ---- helper: run validation for one cohort --------------------------------------
validate_cohort <- function(name, expr, surv, prefix = "") {
  cat("\n---- cohort:", name, "----\n")
  common <- intersect(names(coefs), rownames(expr))
  if (length(common) < 1) {
    cat("      no signature genes present; skipped.\n"); return(NULL)
  }
  cat("      signature genes present:", paste(common, collapse = ", "), "\n")
  z <- t(scale(t(expr[common, , drop = FALSE])))
  risk <- as.numeric(t(z) %*% coefs[common])
  d <- data.table(sample = colnames(expr), risk = risk)
  d <- cbind(d, surv)
  d <- d[!is.na(time) & time > 0 & !is.na(status)]
  d[, group := ifelse(risk >= median(risk), "High", "Low")]
  fwrite(d, file.path(OUT_DIR, paste0(prefix, "risk_", name, ".csv")))
  cat("      n =", nrow(d), "| events =", sum(d$status), "\n")

  # Cox per SD (for forest)
  fitc <- tryCatch(coxph(Surv(time, status) ~ scale(risk), data = d),
                   error = function(e) NULL)
  hr_row <- NULL
  if (!is.null(fitc)) {
    s <- summary(fitc)
    hr_row <- data.table(cohort = name, n = nrow(d),
                         HR = s$conf.int[1, 1], HR_low = s$conf.int[1, 3],
                         HR_high = s$conf.int[1, 4], p = s$coefficients[1, 5])
  }
  # KM
  if (requireNamespace("survminer", quietly = TRUE) &&
      sum(d$status) >= 5 && all(table(d$group) >= 5)) {
    suppressPackageStartupMessages(library(survminer))
    fit_km <- survfit(Surv(time, status) ~ group, data = d)
    p <- file.path(OUT_DIR, paste0(prefix, "km_", name, ".png"))
    g <- ggsurvplot(fit_km, data = d, pval = TRUE, risk.table = TRUE,
                    palette = c("#B2182B", "#2166AC"), xlab = "Time (days)",
                    title = name)
    png(p, width = 1500, height = 1200, res = 150)
    print(g, newpage = FALSE); dev.off()
    cat("      saved", basename(p), "\n")
  }
  # timeROC
  if (requireNamespace("timeROC", quietly = TRUE) && sum(d$status) >= 10) {
    suppressPackageStartupMessages(library(timeROC))
    med <- median(d$time)
    times <- c(1, 3, 5) * 365
    times <- times[times <= max(d$time)]
    if (length(times) >= 2) {
      tr <- timeROC(T = d$time, delta = d$status, marker = d$risk,
                    cause = 1, times = times, iid = FALSE)
      p <- file.path(OUT_DIR, paste0(prefix, "auc_", name, ".png"))
      png(p, width = 1200, height = 1000, res = 150)
      plotAUCcurve(tr, FP = 2, add = FALSE)
      title(main = paste0(name, ": time-dependent AUC"))
      legend("bottomright", legend = paste(names(tr$AUC), "AUC =",
                                           round(tr$AUC, 3)), lty = 1, bty = "n")
      dev.off()
      cat("      saved", basename(p), "time-dependent AUC:",
          paste(round(tr$AUC, 3), collapse = " / "), "\n")
    }
  }
  hr_row
}

## ---- training HR for the forest ---------------------------------------------------
cat("\n[1/4] Training cohort HR (per SD of risk) ...\n")
hr_list <- list()
if (file.exists(RISK_TRAIN)) {
  rt <- fread(RISK_TRAIN)
  rt <- rt[!is.na(OS.time) & OS.time > 0 & !is.na(OS)]
  fitc <- coxph(Surv(OS.time, OS) ~ scale(risk), data = rt)
  s <- summary(fitc)
  hr_list[[1]] <- data.table(cohort = "Training (TCGA-LIHC)", n = nrow(rt),
                             HR = s$conf.int[1, 1], HR_low = s$conf.int[1, 3],
                             HR_high = s$conf.int[1, 4], p = s$coefficients[1, 5])
} else {
  cat("      risk_training.csv missing; forest will show validation only.\n")
}

## ---- validation cohorts -------------------------------------------------------------
cat("\n[2/4] External cohorts ...\n")
# GEOquery downloads with max(getOption("timeout"), 120)s; the ftp mirror is
# slow (~20 MB series matrix), so raise the cap to 10 min.
options(timeout = 600)
geo_ids <- c("GSE14520", "GSE16757")
for (gid in geo_ids) {
  res <- tryCatch({
    if (requireNamespace("GEOquery", quietly = TRUE)) {
      suppressPackageStartupMessages(library(GEOquery))
      # 1) prefer locally cached platform matrix files (fully offline; avoids
      #    GEOquery trying to fetch every platform of multi-platform series)
      local <- list.files(GEO_CACHE,
                          pattern = paste0("^", gid, ".*_series_matrix\\.txt\\.gz$"),
                          full.names = TRUE)
      if (length(local) >= 1) {
        gse_list <- lapply(local, function(lf)
          tryCatch(getGEO(filename = lf, AnnotGPL = TRUE),
                   error = function(e) NULL))
        gse_list <- gse_list[!vapply(gse_list, is.null, logical(1))]
        if (length(gse_list) == 0) stop("local matrix parsing failed")
        # use the platform matrix with the most samples (e.g. GPL3921 for GSE14520)
        pick <- which.max(vapply(gse_list, function(g) ncol(exprs(g)), numeric(1)))
        gse <- gse_list[[pick]]
        cat("      using local matrix:", basename(local[pick]), "\n")
      } else {
        gse <- getGEO(gid, destdir = GEO_CACHE, AnnotGPL = TRUE,
                      getGPL = TRUE)[[1]]
      }
      list(expr = probe2symbol(exprs(gse), fData(gse)), pd = pData(gse))
    } else stop("GEOquery not installed")
  }, error = function(e) {
    cat("      GEO download failed for", gid, ":", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(res)) next
  # GSE14520: survival lives in the GEO supplementary file (pData has none)
  if (gid == "GSE14520") {
    supp <- file.path(GEO_CACHE, "GSE14520_Extra_Supplement.txt.gz")
    if (!file.exists(supp)) {
      cat("      GSE14520 supplement missing; download it first:\n")
      cat("      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz\n")
      cat("      -> save to analysis_out/data/GEO/ then rerun step06.\n")
      next
    }
    sup <- fread(supp)
    nms <- names(sup)
    gsm_col <- nms[grepl("GSM", nms, ignore.case = TRUE)][1]
    st_col  <- nms[grepl("survival.*status|^status", nms, ignore.case = TRUE)][1]
    mo_col  <- nms[grepl("survival.*month", nms, ignore.case = TRUE)][1]
    tis_col <- nms[grepl("tissue", nms, ignore.case = TRUE)][1]
    if (is.na(gsm_col) || is.na(st_col) || is.na(mo_col)) {
      cat("      supplement columns not recognized:", paste(nms, collapse = ", "), "\n")
      next
    }
    sup <- sup[!is.na(get(gsm_col))]
    if (!is.na(tis_col)) sup <- sup[grepl("tumor", get(tis_col), ignore.case = TRUE)]
    mm <- match(colnames(res$expr), as.character(sup[[gsm_col]]))
    surv <- data.table(time = as.numeric(sup[[mo_col]][mm]) * 30.44,
                       status = suppressWarnings(as.numeric(sup[[st_col]][mm])))
    # keep row order aligned with expr columns; validate_cohort filters NAs
    cat("      GSE14520 supplement survival merged (n =", sum(!is.na(surv$time)), ")\n")
  } else if (gid == "GSE16757") {
    # survival lives in the clinical supplement (columns auto-detected)
    supp <- file.path(GEO_CACHE, "GSE16757_Clinical.data.csv.gz")
    if (!file.exists(supp)) {
      cat("      GSE16757 clinical file missing; download it first:\n")
      cat("      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16757/suppl/GSE16757_Clinical.data.csv.gz\n")
      cat("      -> save to analysis_out/data/GEO/ then rerun step06.\n")
      next
    }
    sup <- fread(supp)
    nms <- names(sup)
    gsm_col <- nms[grepl("GSM|GEO_ID", nms, ignore.case = TRUE)][1]
    st_col  <- nms[grepl("status|alive|death|event|censor|vital", nms,
                         ignore.case = TRUE)][1]
    tm_col  <- nms[grepl("month|time|day|year|follow|survival|OS\\.|RFS\\.",
                         nms, ignore.case = TRUE)][1]
    if (is.na(gsm_col) || is.na(st_col) || is.na(tm_col)) {
      cat("      GSE16757 clinical columns not recognized:",
          paste(nms, collapse = ", "), "\n")
      next
    }
    sup <- sup[!is.na(get(gsm_col))]
    mm <- match(colnames(res$expr), as.character(sup[[gsm_col]]))
    unit <- if (grepl("month|\\.m$", tm_col, ignore.case = TRUE)) 30.44
            else if (grepl("year", tm_col, ignore.case = TRUE)) 365.25
            else 1
    tm <- as.numeric(sup[[tm_col]][mm]) * unit
    st <- as.character(sup[[st_col]][mm])
    st_num <- suppressWarnings(as.numeric(st))
    event <- ifelse(!is.na(st_num), pmin(pmax(st_num, 0), 1),
                    ifelse(grepl("dead|deceased|recurr|^1|:1", st,
                                 ignore.case = TRUE), 1, 0))
    surv <- data.table(time = tm, status = event)
    cat("      GSE16757 clinical survival merged (n =", sum(!is.na(surv$time)), ")\n")
  } else {
    ov <- SURV_OVERRIDES[[gid]]
    surv <- tryCatch(
      detect_surv(res$pd,
                  time_col = if (!is.null(ov)) ov$time_col else NA,
                  status_col = if (!is.null(ov)) ov$status_col else NA,
                  time_unit = if (!is.null(ov)) ov$time_unit else "auto"),
      error = function(e) {
        cat("      ", conditionMessage(e), "\n"); NULL })
    if (is.null(surv)) next
  }
  hr_list[[length(hr_list) + 1]] <-
    validate_cohort(gid, res$expr, surv)
}

## ---- self-test fallback --------------------------------------------------------------
if (length(hr_list) < 2) {
  cat("\n[3/4] Fewer than 2 cohorts processed -> SELF-TEST with synthetic cohorts.\n")
  set.seed(SEED)
  for (nm in c("SYN-cohort-A", "SYN-cohort-B")) {
    n <- 120
    genes <- c(names(coefs), paste0("BG", 1:250))
    m <- matrix(rnorm(length(genes) * n, 6, 1.2), nrow = length(genes),
                dimnames = list(genes, paste0(nm, "-", 1:n)))
    signal <- rnorm(n)
    m[names(coefs), ] <- m[names(coefs), ] + matrix(signal, length(coefs), n,
                                                    byrow = TRUE)
    time <- round(rexp(n, 0.0015) * exp(-0.6 * signal)); time <- pmin(time, 5000)
    status <- rbinom(n, 1, plogis(-0.8 + 1.1 * signal))
    hr_list[[length(hr_list) + 1]] <-
      validate_cohort(nm, m, data.table(time = time, status = status),
                      prefix = "self_test_")
  }
}

## ---- forest plot -----------------------------------------------------------------------
cat("\n[4/4] Forest plot ...\n")
hr_all <- rbindlist(hr_list, fill = TRUE)
if (nrow(hr_all) >= 2 && requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  hr_all[, cohort := factor(cohort, levels = rev(cohort))]
  g <- ggplot(hr_all, aes(HR, cohort)) +
    geom_errorbar(aes(xmin = HR_low, xmax = HR_high, y = cohort),
                  orientation = "y", width = 0.2, color = "grey40") +
    geom_point(size = 3, color = "#B2182B") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_x_log10() +
    labs(x = "Hazard ratio (per SD of risk)", y = NULL,
         title = "FARGs signature: training + validation") +
    theme_minimal(base_size = 13)
  ggsave(file.path(OUT_DIR, "forest_hr.png"), g, width = 8, height = 4.5,
         dpi = 200)
  cat("      saved forest_hr.png\n")
}
fwrite(hr_all, file.path(OUT_DIR, "validation_hr_summary.csv"))
print(hr_all)
cat("\n== step06 done ==\n")
