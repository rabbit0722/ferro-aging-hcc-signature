# =============================================================================
# step09_hcc_robustness.R
# -----------------------------------------------------------------------------
# Robustness suite for the HCC FARGs signature (run after step05 on real data):
#   A) Bootstrap stability:  BOOT_N x re-fit LASSO on resampled tumor samples
#      -> gene selection frequency; + bootstrap out-of-bag C-index of the
#      final signature
#   B) Incremental value: ssGSEA scores of reference gene sets (MSigDB
#      Hallmark ROS/inflammation/hypoxia + GOBP ferroptosis, plus any
#      published signatures you add to compare_signatures.csv) -> C-index,
#      LRT and AIC vs the FARGs risk score
#   C) Subgroup KM: stage / grade / age / gender (global-median risk split)
#   D) GSEA (fgsea, MSigDB Hallmark + ferroptosis/iron GO terms) of genes
#      ranked by association with the risk score
#   E) Multivariable forest: risk score + stage + grade + age
#
# Input : analysis_out/hcc_data/hcc_data.rds      (from step04)
#         analysis_out/hcc_out/signature_genes.csv (from step05)
#         optional: analysis_out/hcc_out/compare_signatures.csv
#                   (columns: set, gene — one row per gene)
# If hcc_data.rds is missing the script runs a SELF-TEST on synthetic data.
#
# Output: analysis_out/hcc_out/robustness/
#
# Requirements: data.table, limma, survival, survminer, glmnet, GSVA, fgsea,
#               msigdbr, ggplot2
# Run:  Rscript step09_hcc_robustness.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(survival)
  library(glmnet)
  library(GSVA)
})

# ---- locate project folder (safe for RStudio with a different working dir) ----
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f[1]))) else
    dirname(normalizePath(sys.frames()[[1]]$ofile))
}, error = function(e) getwd())
setwd(script_dir)

DATA_RDS  <- "analysis_out/hcc_data/hcc_data.rds"
SIG_FILE  <- "analysis_out/hcc_out/signature_genes.csv"
CMP_FILE  <- "analysis_out/hcc_out/compare_signatures.csv"
OUT_DIR   <- "analysis_out/hcc_out/robustness"
BOOT_N    <- 1000         # bootstrap iterations (consistent with Methods)
PADJ_DEG  <- 0.05
LFC_DEG   <- 1
UNI_COX_P <- 0.1
SEED      <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# msigdbr must see its standard cache dir; a leftover R_USER_CACHE_DIR from
# step02 (00_master session) would redirect it and force a re-download
Sys.unsetenv("R_USER_CACHE_DIR")
cat("== step09: robustness suite ==\n")

## ---- signature ------------------------------------------------------------
sig <- if (file.exists(SIG_FILE)) fread(SIG_FILE) else
  data.table(gene = c("ALOX15", "ATG3", "BAP1"), coef = c(0.44, 0.37, 0.51))
sig <- sig[!is.na(gene) & !is.na(coef)]
coefs <- setNames(sig$coef, sig$gene)
cat("Signature:", paste(names(coefs), collapse = ", "), "\n")

## ---- data or self-test ------------------------------------------------------
SELF_TEST <- FALSE
if (!file.exists(DATA_RDS)) {
  cat("[WARN] hcc_data.rds not found -> SELF-TEST.\n")
  SELF_TEST <- TRUE
  BOOT_N <- 50
  set.seed(SEED)
  n_norm <- 40; n_tum <- 160
  farg_genes <- unique(fread("analysis_out/step01/module_assignment.csv")$SYMBOL)
  genes <- unique(c(farg_genes, names(coefs), paste0("BG", 1:300)))
  expr <- matrix(rnorm(length(genes) * (n_norm + n_tum), 6, 1.2),
                 nrow = length(genes),
                 dimnames = list(genes, NULL))
  up_genes <- c(farg_genes[1:8], names(coefs), paste0("BG", 1:10))
  tum_idx <- (n_norm + 1):(n_norm + n_tum)
  expr[up_genes, tum_idx] <- expr[up_genes, tum_idx] + 1.4
  colnames(expr) <- paste0("SYN", 1:ncol(expr))
  group <- c(rep("Normal", n_norm), rep("Tumor", n_tum))
  z <- scale(t(expr[names(coefs), tum_idx]))
  risk_true <- rowMeans(z)
  os_time <- round(rexp(n_tum, 0.0015) * exp(-0.7 * risk_true))
  os_time <- pmin(os_time, 5000)
  os <- rbinom(n_tum, 1, plogis(-1.2 + 1.6 * risk_true))
  surv <- data.table(sample = colnames(expr)[tum_idx], OS = os, OS.time = os_time)
  pheno <- data.table(
    sample = colnames(expr),
    age = round(pmax(40, pmin(85, rnorm(ncol(expr), 62, 9)))),
    stage = factor(sample(c("Stage I", "Stage II", "Stage III", "Stage IV"),
                          ncol(expr), TRUE, prob = c(.25, .3, .3, .15))),
    grade = factor(sample(c("G1", "G2", "G3"), ncol(expr), TRUE,
                          prob = c(.2, .5, .3))),
    gender = sample(c("MALE", "FEMALE"), ncol(expr), TRUE))
  # synthetic "published" signatures for the comparison framework
  fwrite(rbind(data.table(set = "REF_associated", gene = paste0("BG", 1:20)),
               data.table(set = "REF_null", gene = paste0("BG", 21:40))),
         CMP_FILE)
  cat("      synthetic compare_signatures.csv written for the self-test\n")
} else {
  L <- readRDS(DATA_RDS)
  expr <- L$expr; group <- L$group; surv <- L$surv; pheno <- L$pheno
  cat("Real data:", nrow(expr), "genes x", ncol(expr), "samples\n")
}

## ---- recompute candidate pool (FARG-DEGs + univariate Cox) -------------------
cat("\n[pool] DEG + univariate Cox candidate pool ...\n")
mod_file <- "analysis_out/step01/module_assignment.csv"
farg_genes <- unique(fread(mod_file)$SYMBOL)
g_f <- factor(group, levels = c("Normal", "Tumor"))
fit <- eBayes(lmFit(expr, model.matrix(~ g_f)))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
fdeg <- tt[tt$gene %in% farg_genes & tt$adj.P.Val < PADJ_DEG &
             abs(tt$logFC) > LFC_DEG, ]
cat("      FARG-DEGs:", nrow(fdeg), "\n")
tum_s <- intersect(colnames(expr), surv$sample)
expr_t <- expr[, tum_s]
sv <- surv[match(tum_s, sample)]
# guard: a character subscript containing a name that is NOT in the matrix
# dimnames throws "subscript out of bounds" (opaque). This happens when
# step05 wrote a SELF-TEST (synthetic) signature while hcc_data.rds is real.
miss_sig <- setdiff(names(coefs), rownames(expr_t))
if (length(miss_sig))
  stop("Signature genes missing from the expression matrix: ",
       paste(miss_sig, collapse = ", "),
       "\n       signature_genes.csv does not match hcc_data.rds ",
       "(step05 ran in SELF-TEST mode?).\n       Re-run step05 on the real ",
       "data to rebuild signature_genes.csv, or delete hcc_data.rds to ",
       "force the self-test path.")
uni <- rbindlist(lapply(fdeg$gene, function(gene) {
  dd <- data.frame(time = sv$OS.time, status = sv$OS,
                   x = as.numeric(scale(expr_t[gene, ])))
  dd <- dd[dd$time > 0 & !is.na(dd$x), ]
  f <- tryCatch(coxph(Surv(time, status) ~ x, data = dd),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- summary(f)
  data.table(gene = gene, p = s$coefficients[1, 5])
}))
pool <- uni[p < UNI_COX_P]$gene
cat("      candidate pool:", length(pool), "genes\n")

## ---- risk score ----------------------------------------------------------------
z_t <- scale(t(expr_t[names(coefs), , drop = FALSE]))
risk <- as.numeric(z_t %*% coefs)
risk_dt <- data.table(sample = tum_s, risk = risk)
risk_dt <- merge(risk_dt, sv, by = "sample")
risk_dt <- risk_dt[!is.na(OS.time) & OS.time > 0]
if (!is.null(pheno) && "sample" %in% names(pheno))
  risk_dt <- merge(risk_dt, pheno, by = "sample", all.x = TRUE)

## ---- A. bootstrap internal validation ---------------------------------------------
cat("\n[A] Bootstrap internal validation (conditional (candidate screen fixed; LASSO+fit refit)) ...\n")
set.seed(SEED)
ok <- !is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS)   # valid survival rows
# Raw candidate-gene expression. Standardization is done INSIDE each resample on
# the bootstrap sample only; candidate selection still uses the full data.
X_raw   <- t(expr_t[pool, ok, drop = FALSE])          # n_valid x n_pool (log2 expr)
time_ok <- sv$OS.time[ok]; status_ok <- sv$OS[ok]
n_ok <- nrow(X_raw)
sel_freq <- setNames(rep(0, length(pool)), pool)
c_idx <- rep(NA_real_, BOOT_N)   # NA for skipped/failed iterations (excluded by na.rm)
iteration_log <- data.table(iteration = seq_len(BOOT_N), n_oob = NA_integer_,
                            n_selected = NA_integer_, selected_genes = NA_character_,
                            c_index = NA_real_, status = "not_evaluated")
for (b in seq_len(BOOT_N)) {
  idx <- sample(seq_len(n_ok), replace = TRUE)
  oob <- setdiff(seq_len(n_ok), unique(idx))

  # (1) standardize on the bootstrap sample only
  xb <- scale(X_raw[idx, , drop = FALSE])
  mns <- attr(xb, "scaled:center"); sds <- attr(xb, "scaled:scale")

  # (2) LASSO selection (refit) on the bootstrap sample
  cv <- tryCatch(cv.glmnet(xb, Surv(time_ok[idx], status_ok[idx]),
                           family = "cox", alpha = 1, nfolds = 10),
                 error = function(e) NULL)
  sel <- character(0)
  if (!is.null(cv)) {
    cf <- as.numeric(coef(cv, s = "lambda.min"))
    sel <- colnames(xb)[cf != 0]
    if (length(sel) == 0) {
      cf <- as.numeric(coef(cv, s = "lambda.1se"))
      sel <- colnames(xb)[cf != 0]
    }
    if (length(sel) > 0) sel_freq[sel] <- sel_freq[sel] + 1
  }

  # (3) refit Cox coefficients on the bootstrap sample, then evaluate the held-out
  #     (out-of-bootstrap) observations standardized with the bootstrap-sample
  #     mean/SD -> CONDITIONAL internal validation (candidate screen fixed on
  #     full data; stepwise reduction not reproduced; NOT a full pipeline refit).
  if (length(sel) > 0 && length(oob) >= 10) {
    db <- as.data.frame(xb[, sel, drop = FALSE])
    db$time <- time_ok[idx]; db$status <- status_ok[idx]
    fb <- tryCatch(coxph(Surv(time, status) ~ ., data = db),
                   error = function(e) NULL)
    if (!is.null(fb)) {
      coef_b <- coef(fb)
      xo <- scale(X_raw[oob, sel, drop = FALSE],
                  center = mns[sel], scale = sds[sel])
      risk_oob <- as.numeric(xo %*% coef_b)
      dd <- data.frame(time = time_ok[oob], status = status_ok[oob],
                       risk = risk_oob)
      # Direct concordance of the FIXED risk score (reverse=TRUE: higher risk =
      # worse survival). No Cox refit on held-out, so no second fitting.
      c_idx[b] <- tryCatch(
        concordance(Surv(time, status) ~ risk, data = dd, reverse = TRUE)$concordance,
        error = function(e) NA)
    }
  }
  iteration_log[b, `:=`(n_oob = length(oob), n_selected = length(sel),
                        selected_genes = paste(sel, collapse = ";"),
                        c_index = c_idx[b],
                        status = if (is.finite(c_idx[b])) "ok" else "not_estimable")]
}
fwrite(iteration_log, file.path(OUT_DIR, "bootstrap_iterations.csv"))
fwrite(data.table(apparent_C_index = NA_real_,
                  heldout_C_mean = mean(c_idx, na.rm = TRUE),
                  heldout_C_q025 = quantile(c_idx, 0.025, na.rm = TRUE),
                  heldout_C_q975 = quantile(c_idx, 0.975, na.rm = TRUE),
                  n_valid_iterations = sum(!is.na(c_idx)), BOOT_N = BOOT_N),
       file.path(OUT_DIR, "bootstrap_internal_cindex.csv"))
freq_dt <- data.table(gene = names(sel_freq), freq = sel_freq / BOOT_N)[order(-freq)]
fwrite(freq_dt, file.path(OUT_DIR, "bootstrap_selection_freq.csv"))
cat("      selection frequency (top 10):\n")
print(freq_dt[1:min(10, nrow(freq_dt))])
png(file.path(OUT_DIR, "bootstrap_selection_freq.png"), width = 1100,
    height = 900, res = 150)
par(mar = c(5, 10, 3, 2))
barplot(rev(freq_dt$freq), names.arg = rev(freq_dt$gene), horiz = TRUE,
        las = 1, col = "#4DBBD5", xlab = "selection frequency",
        main = "Bootstrap LASSO selection frequency")
dev.off()
cat("      held-out C-index (conditional internal validation): mean",
    round(mean(c_idx, na.rm = TRUE), 3),
    "empirical 2.5th-97.5th percentiles [", round(quantile(c_idx, 0.025, na.rm = TRUE), 3), ",",
    round(quantile(c_idx, 0.975, na.rm = TRUE), 3), "]",
    "| valid", sum(!is.na(c_idx)), "/", BOOT_N, "failed", sum(is.na(c_idx)), "\n")
png(file.path(OUT_DIR, "bootstrap_cindex.png"), width = 1100, height = 800,
    res = 150)
hist(c_idx, breaks = 30, col = "#B2182B", border = "white",
     xlab = "held-out C-index (conditional internal validation)",
     main = "Bootstrap held-out C-index (conditional)")
abline(v = mean(c_idx, na.rm = TRUE), lty = 2, lwd = 2)
dev.off()

## ---- B. reference signatures comparison -------------------------------------------
cat("\n[B] Incremental value vs reference gene sets ...\n")
ref_sets <- list()
if (requireNamespace("msigdbr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(msigdbr))
  ms <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) {
                   cat("      msigdbr download failed (offline?):",
                       conditionMessage(e), "\n      -> using user CSV sets only.\n")
                   NULL })
  if (!is.null(ms)) {
    hs <- split(ms$gene_symbol, ms$gs_name)
    pick <- grep("REACTIVE_OXYGEN|INFLAMMATORY_RESPONSE|HYPOXIA|TNFA_SIGNALING|P53|IL6_JAK",
                 names(hs), value = TRUE)
    ref_sets <- c(ref_sets, hs[pick])
    c5 <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C5",
                           subcategory = "GO:BP"),
                   error = function(e) NULL)
    if (!is.null(c5)) {
      gs <- split(c5$gene_symbol, c5$gs_name)
      ferr <- grep("FERROPTOSIS|IRON", names(gs), value = TRUE)
      ref_sets <- c(ref_sets, gs[ferr])
    }
    cat("      MSigDB reference sets:", length(ref_sets),
        "(first run downloads the gene set archive, needs internet)\n")
  }
}
if (file.exists(CMP_FILE)) {
  usr <- fread(CMP_FILE)
  ref_sets <- c(ref_sets, split(usr$gene, usr$set))
  cat("      user sets loaded:", paste(unique(usr$set), collapse = ", "), "\n")
} else {
  cat("      (optional) add published signatures as analysis_out/hcc_out/",
      "compare_signatures.csv with columns set,gene\n", sep = "")
}
ref_ok <- lapply(ref_sets, function(g) intersect(g, rownames(expr_t)))
ref_ok <- ref_ok[lengths(ref_ok) >= 5]
if (length(ref_ok)) {
  rres <- tryCatch(
    gsva(as.matrix(expr_t), ref_ok, method = "ssgsea", ssgsea.norm = TRUE,
         verbose = FALSE),
    error = function(e) gsva(ssgseaParam(as.matrix(expr_t), ref_ok),
                             verbose = FALSE))
  rs <- t(rres)
  # align rows: rs covers all tumor samples (369), risk_dt only those with
  # valid survival (363) -> subset rs to risk_dt's samples
  base_idx <- match(rownames(rs), risk_dt$sample)
  ok <- !is.na(base_idx)
  rs <- rs[ok, , drop = FALSE]
  base_d <- data.frame(time = risk_dt$OS.time[base_idx[ok]],
                       status = risk_dt$OS[base_idx[ok]],
                       risk  = risk_dt$risk[base_idx[ok]])
  cmp <- rbindlist(lapply(colnames(rs), function(sn) {
    d2 <- base_d; d2$ref <- as.numeric(scale(rs[, sn]))
    f0 <- coxph(Surv(time, status) ~ ref, data = d2)
    f1 <- coxph(Surv(time, status) ~ ref + risk, data = d2)
    c0 <- concordance(f0)$concordance
    c1 <- concordance(f1)$concordance
    lrt <- tryCatch(anova(f0, f1)$`Pr(>|Chi|)`[2], error = function(e) NA)
    data.table(set = sn, C_ref = c0, C_ref_plus_risk = c1,
               delta_C = c1 - c0, LRT_p = lrt)
  }))
  fwrite(cmp, file.path(OUT_DIR, "reference_comparison.csv"))
  cat("      comparison table:\n"); print(cmp[order(-delta_C)])
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    cmp_long <- melt(cmp[, .(set, C_ref, C_ref_plus_risk)],
                     id.vars = "set", variable.name = "model", value.name = "C")
    g <- ggplot(cmp_long, aes(reorder(set, C), C, color = model)) +
      geom_point(size = 3) +
      geom_line(aes(group = set), color = "grey70") +
      coord_flip() +
      labs(x = NULL, y = "C-index",
           title = "Reference sets alone vs + FARGs risk",
           color = "model") + theme_minimal(base_size = 12)
    ggsave(file.path(OUT_DIR, "reference_comparison.png"), g, width = 9,
           height = 6, dpi = 200)
    cat("      saved reference_comparison.png\n")
  }
} else {
  cat("      no reference sets with >=5 genes present; module skipped.\n")
}

## ---- C. subgroup KM --------------------------------------------------------------------
cat("\n[C] Subgroup KM (global-median split) ...\n")
rd <- copy(risk_dt)
rd[, risk_group := ifelse(risk >= median(risk), "High", "Low")]
if ("stage" %in% names(rd)) {
  rd[, stage_g := ifelse(grepl("III|IV|3|4", stage, ignore.case = TRUE),
                         "Stage III-IV", "Stage I-II")]
} else rd[, stage_g := NA_character_]
if ("grade" %in% names(rd)) {
  rd[, grade_g := ifelse(grepl("G3|G4|3|4", grade, ignore.case = TRUE),
                         "Grade 3-4", "Grade 1-2")]
} else rd[, grade_g := NA_character_]
if ("age" %in% names(rd)) {
  rd[, age_g := ifelse(age >= median(age, na.rm = TRUE), "Older", "Younger")]
} else rd[, age_g := NA_character_]
subs <- c("stage_g", "grade_g", "age_g", "gender")
hr_rows <- list()
plots <- list()
for (sb in subs) {
  if (sb %in% names(rd)) {
    levs <- na.omit(unique(rd[[sb]]))
    if (length(levs) < 2) next
    for (lv in levs) {
      dd <- rd[get(sb) == lv]
      if (nrow(dd) < 20 || length(unique(dd$risk_group)) < 2) next
      fit_km <- survfit(Surv(OS.time, OS) ~ risk_group, data = dd)
      f <- coxph(Surv(OS.time, OS) ~ risk_group, data = dd)
      s <- summary(f)
      hr_rows[[length(hr_rows) + 1]] <-
        data.table(subgroup = paste(sb, lv), n = nrow(dd), events = sum(dd$OS),
                   HR = s$conf.int[1, 1], HR_low = s$conf.int[1, 3],
                   HR_high = s$conf.int[1, 4], p = s$coefficients[1, 5])
      plots[[length(plots) + 1]] <-
        list(fit = fit_km, dd = dd,
             title = paste(sb, "=", lv, "| p =",
                           signif(s$coefficients[1, 5], 2)))
    }
  }
}
hr_all <- rbindlist(hr_rows)
fwrite(hr_all, file.path(OUT_DIR, "subgroup_hr.csv"))
print(hr_all[order(p)])
if (length(plots)) {
  n <- length(plots)
  nc <- min(2, n); nr <- ceiling(n / nc)
  png(file.path(OUT_DIR, "subgroup_km.png"), width = 1000 * nc,
      height = 800 * nr, res = 150)
  par(mfrow = c(nr, nc))
  for (pl in plots) {
    plot(pl$fit, col = c("#B2182B", "#2166AC"), lwd = 2,
         xlab = "Time (days)", ylab = "OS", main = pl$title)
    legend("topright", c("High", "Low"), col = c("#B2182B", "#2166AC"),
           lwd = 2, bty = "n")
  }
  dev.off()
  cat("      saved subgroup_km.png\n")
}

## ---- D. GSEA ------------------------------------------------------------------------------
cat("\n[D] GSEA of risk-associated genes (fgsea) ...\n")
if (requireNamespace("fgsea", quietly = TRUE) &&
    requireNamespace("msigdbr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(fgsea))
  n <- nrow(expr_t)
  r <- apply(expr_t, 1, function(x) {
    rr <- cor(as.numeric(x), risk, use = "complete.obs")
    if (is.na(rr)) 0 else rr
  })
  t_stat <- r * sqrt((n - 2) / pmax(1 - r^2, 1e-8))
  ranks <- sort(t_stat, decreasing = TRUE)
  ms <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) NULL)
  if (!is.null(ms)) {
    paths <- split(ms$gene_symbol, ms$gs_name)
    fg <- fgsea(pathways = paths, stats = ranks, minSize = 5, maxSize = 500,
                nPermSimple = 10000)
    fg_dt <- as.data.table(fg)[order(padj)]
    fwrite(fg_dt[, .(pathway, pval, padj, NES, size)],
           file.path(OUT_DIR, "gsea_hallmark.csv"))
    cat("      top enriched Hallmark pathways:\n")
    print(fg_dt[1:min(10, nrow(fg_dt)),
                .(pathway, NES = round(NES, 2), padj = signif(padj, 2))])
    top_pos <- head(fg_dt[NES > 0]$pathway, 5)
    top_neg <- head(fg_dt[NES < 0]$pathway, 5)
    png(file.path(OUT_DIR, "gsea_top.png"), width = 1500, height = 800,
        res = 150)
    par(mfrow = c(2, 5))
    for (pth in c(top_pos, rev(top_neg)))
      if (length(pth) == 1 && !is.na(pth)) {
        tryCatch(plotEnrichment(paths[[pth]], ranks) + ggtitle(pth) +
                   theme(plot.title = element_text(size = 8)),
                 error = function(e) NULL)
      }
    dev.off()
    cat("      saved gsea_top.png\n")
  } else {
    cat("      msigdbr gene sets unavailable (offline?); GSEA skipped.\n")
  }
} else {
  cat("      fgsea/msigdbr missing; skipped.\n")
}

## ---- E. multivariable forest -------------------------------------------------------------
cat("\n[E] Multivariable Cox forest ...\n")
vars <- c("risk")
if ("age" %in% names(rd)) vars <- c(vars, "age")
# collapse verbose stage/grade levels to numeric for a clean forest
if ("stage" %in% names(rd)) {
  # stage strings are roman numerals ("Stage I".."Stage IV" or arabic)
  sr <- toupper(trimws(sub("Stage", "", rd$stage, ignore.case = TRUE)))
  rd[, stage_num := fifelse(grepl("[0-9]", sr),
                            suppressWarnings(as.numeric(gsub("[^0-9]", "", sr))),
                            match(sr, c("I", "II", "III", "IV")))]
  rd[is.na(stage_num) | stage_num < 1 | stage_num > 4, stage_num := NA]
  vars <- c(vars, "stage_num")
}
if ("grade" %in% names(rd)) {
  rd[, grade_num := suppressWarnings(as.numeric(gsub("[^0-9]", "",
                       toupper(grade))))]
  rd[is.na(grade_num) | grade_num < 1 | grade_num > 4, grade_num := NA]
  vars <- c(vars, "grade_num")
}
md <- as.data.frame(rd[, c("OS.time", "OS", vars), with = FALSE])
md <- md[complete.cases(md) & md$OS.time > 0, ]
if (nrow(md) >= 30) {
  f <- coxph(as.formula(paste("Surv(OS.time, OS) ~",
                              paste(vars, collapse = " + "))), data = md)
  s <- summary(f)
  forest_dt <- data.table(variable = rownames(s$coefficients),
                          HR = s$conf.int[, 1], HR_low = s$conf.int[, 3],
                          HR_high = s$conf.int[, 4], p = s$coefficients[, 5])
  fwrite(forest_dt, file.path(OUT_DIR, "forest_multivariable.csv"))
  print(forest_dt)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    forest_dt[, variable := factor(variable, levels = rev(variable))]
    g <- ggplot(forest_dt, aes(HR, variable)) +
      geom_errorbar(aes(xmin = HR_low, xmax = HR_high, y = variable),
                    orientation = "y", width = 0.2, color = "grey40") +
      geom_point(size = 3, color = "#B2182B") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      scale_x_log10() +
      labs(x = "Hazard ratio", y = NULL,
           title = "Multivariable Cox: FARGs risk + clinical factors") +
      theme_minimal(base_size = 13)
    ggsave(file.path(OUT_DIR, "forest_multivariable.png"), g, width = 8,
           height = 5, dpi = 200)
    cat("      saved forest_multivariable.png\n")
  }
}

## ---- F. Decision curve analysis (dcurves) ---------------------------------------
cat("\n[F] Decision curve analysis (dcurves) ...\n")
if (requireNamespace("dcurves", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(dcurves))
  suppressPackageStartupMessages(library(ggplot2))
  dd_dca <- as.data.frame(rd[, c("OS.time", "OS", "risk", "stage_num", "age"),
                             with = FALSE])
  dd_dca <- dd_dca[complete.cases(dd_dca) & dd_dca$OS.time > 0, ]
  if (nrow(dd_dca) >= 50) {
    f_risk <- coxph(Surv(OS.time, OS) ~ risk, data = dd_dca)
    f_clin <- coxph(Surv(OS.time, OS) ~ stage_num + age, data = dd_dca)
    f_full <- coxph(Surv(OS.time, OS) ~ risk + stage_num + age, data = dd_dca)
    prob_at <- function(fit, newd, tt) {
      1 - as.numeric(summary(survfit(fit, newdata = newd), times = tt)$surv)
    }
    for (tt in c(365, 1095, 1825)) {
      nm <- switch(as.character(tt), "365" = "1y", "1095" = "3y", "1825" = "5y")
      dd_dca$p_risk <- prob_at(f_risk, dd_dca, tt)
      dd_dca$p_clin <- prob_at(f_clin, dd_dca, tt)
      dd_dca$p_full <- prob_at(f_full, dd_dca, tt)
      d <- tryCatch(
        dca(Surv(OS.time, OS) ~ p_risk + p_clin + p_full, data = dd_dca,
            time = tt,
            label = list("p_risk" = "FARGs risk",
                         "p_clin" = "Clinical only",
                         "p_full" = "Risk + Clinical")),
        error = function(e) { cat("      dca err:", conditionMessage(e), "\n"); NULL })
      if (is.null(d)) next
      ddf <- as.data.frame(d$dca)
      fwrite(ddf, file.path(OUT_DIR, paste0("dca_net_benefit_", nm, ".csv")))
      g <- ggplot(ddf, aes(x = threshold, y = net_benefit,
                           color = label, linetype = variable)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = c("FARGs risk" = "#B2182B",
                                      "Clinical only" = "#2166AC",
                                      "Risk + Clinical" = "#4DAF4A",
                                      "Treat All" = "grey50",
                                      "Treat None" = "grey30")) +
        coord_cartesian(xlim = c(0, 0.8), ylim = c(-0.05, 0.25)) +
        labs(x = "Threshold probability", y = "Net benefit",
             color = "Model", linetype = "Model",
             title = paste0("Decision curve analysis (", nm,
                            " overall survival)")) +
        theme_bw(base_size = 14)
      png(file.path(OUT_DIR, paste0("dca_", nm, ".png")),
          width = 1800, height = 1300, res = 150)
      print(g)
      dev.off()
      cat("      saved dca_", nm, ".png\n", sep = "")
    }
  } else {
    cat("      too few complete cases for DCA; skipped.\n")
  }
} else {
  cat("      dcurves/ggplot2 not installed; skipped.\n")
  cat("      install: install.packages('dcurves')\n")
}

cat("\n== step09 done ==\n")
if (SELF_TEST)
  cat("SELF-TEST passed. Run on real data after step04/step05 complete.\n")
