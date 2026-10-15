# =============================================================================
# step05_hcc_signature.R
# -----------------------------------------------------------------------------
# HCC single-cancer core pipeline (training cohort = TCGA-LIHC):
#   DEG (limma) -> ferro-aging DEGs (FARGs ∩ DEGs) -> GO/KEGG enrichment
#   -> univariate Cox -> LASSO -> multivariate Cox (stepwise)
#   -> FARGs risk signature -> KM + timeROC + nomogram + calibration
#
# Input : analysis_out/hcc_data/hcc_data.rds  (built by step04_hcc_prep.R)
#         analysis_out/step02/farg_human_map.csv (fallback: step01 module file)
# If hcc_data.rds is missing the script runs a SELF-TEST on synthetic data.
#
# Output: analysis_out/hcc_out/
#   deg_volcano.png, farg_deg.csv, enrichment.png
#   signature_genes.csv, risk_training.csv
#   km_training.png, auc_training.png, risk_dist.png
#   nomogram.png, calibration.png
#
# Requirements: data.table, limma, survival, survminer, glmnet, timeROC,
#               ggplot2, rms, clusterProfiler, org.Hs.eg.db
# Run:  Rscript step05_hcc_signature.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(survival)
  library(glmnet)
})

OUT_DIR      <- "analysis_out/hcc_out"
DATA_RDS     <- "analysis_out/hcc_data/hcc_data.rds"
PADJ_DEG     <- 0.05
LFC_DEG      <- 1
UNI_COX_P    <- 0.1
SEED         <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("== step05: HCC FARGs signature construction ==\n")
cat("R:", R.version.string, "\n")

## ---- FARG gene list ----------------------------------------------------------
map_file <- "analysis_out/step02/farg_human_map.csv"
if (file.exists(map_file)) {
  farg_genes <- fread(map_file)$symbol_human
} else if (file.exists("analysis_out/step01/module_assignment.csv")) {
  farg_genes <- fread("analysis_out/step01/module_assignment.csv")$SYMBOL
} else {
  stop("FARG list not found: run step01/step02 first.")
}
farg_genes <- unique(farg_genes[!is.na(farg_genes) & farg_genes != ""])
cat("FARG gene list:", length(farg_genes), "genes\n")

## ---- load data or self-test ----------------------------------------------------
SELF_TEST <- FALSE
if (!file.exists(DATA_RDS)) {
  cat("\n[WARN] hcc_data.rds not found -> SELF-TEST on synthetic data.\n")
  cat("       Run step04_hcc_prep.R after downloading Xena files for real data.\n\n")
  SELF_TEST <- TRUE
  set.seed(SEED)
  n_norm <- 40; n_tum <- 160
  bg_genes <- paste0("BG", 1:200)
  all_genes <- c(farg_genes, bg_genes)
  expr <- matrix(rnorm(length(all_genes) * (n_norm + n_tum), 6, 1.2),
                 nrow = length(all_genes),
                 dimnames = list(all_genes, NULL))
  # DEG signal: 8 FARGs + 10 bg up in tumor; 8 bg down in tumor
  up_genes   <- c(farg_genes[1:8], bg_genes[1:10])
  down_genes <- bg_genes[11:18]
  tum_idx <- (n_norm + 1):(n_norm + n_tum)
  expr[up_genes, tum_idx]   <- expr[up_genes, tum_idx] + 1.4
  expr[down_genes, tum_idx] <- expr[down_genes, tum_idx] - 1.2
  colnames(expr) <- paste0("SYN", 1:ncol(expr))
  group <- c(rep("Normal", n_norm), rep("Tumor", n_tum))
  # survival signal driven by 6 of the 8 up FARGs
  sig_genes <- farg_genes[1:6]
  z <- scale(t(expr[sig_genes, tum_idx]))
  risk_true <- rowMeans(z)
  os_time <- round(rexp(n_tum, 0.0015) * exp(-0.7 * risk_true))
  os_time <- pmin(os_time, 5000)
  os <- rbinom(n_tum, 1, plogis(-1.2 + 1.6 * risk_true))
  surv <- data.table(sample = colnames(expr)[tum_idx], OS = os, OS.time = os_time)
  pheno <- data.table(
    sample = colnames(expr),
    age = round(pmax(40, pmin(85, rnorm(ncol(expr), 62, 9)))),
    stage = factor(sample(c("Stage I", "Stage II", "Stage III", "Stage IV"),
                          ncol(expr), replace = TRUE, prob = c(.25, .3, .3, .15))),
    grade = factor(sample(c("G1", "G2", "G3"), ncol(expr), replace = TRUE,
                          prob = c(.2, .5, .3))),
    gender = sample(c("MALE", "FEMALE"), ncol(expr), replace = TRUE))
} else {
  L <- readRDS(DATA_RDS)
  expr <- L$expr; group <- L$group; surv <- L$surv; pheno <- L$pheno
  cat("\nLoaded real data:", nrow(expr), "genes x", ncol(expr), "samples |",
      "Tumor:", sum(group == "Tumor"), "Normal:", sum(group == "Normal"), "\n")
}

# SELF-TEST outputs go to a separate folder so synthetic results can never
# overwrite real results (a self-test signature written over the real
# signature_genes.csv crashes step09 with "subscript out of bounds").
if (SELF_TEST) {
  OUT_DIR <- file.path(OUT_DIR, "self_test")
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cat("[self-test] outputs redirected to", OUT_DIR,
      "; real hcc_out files are left untouched\n\n")
}

## ---- 1. DEG (limma) ------------------------------------------------------------
cat("\n[1/7] Differential expression (limma) ...\n")
g_f <- factor(group, levels = c("Normal", "Tumor"))
design <- model.matrix(~ g_f)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
deg <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < PADJ_DEG &
            abs(tt$logFC) > LFC_DEG, ]
cat("      total DEGs:", nrow(deg), "(up:", sum(deg$logFC > 0),
    "down:", sum(deg$logFC < 0), ")\n")
farg_deg <- deg[deg$gene %in% farg_genes, ]
cat("      FARGs ∩ DEGs:", nrow(farg_deg), "\n")
print(farg_deg[, c("gene", "logFC", "adj.P.Val")])
fwrite(farg_deg[, c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")],
       file.path(OUT_DIR, "farg_deg.csv"))

# volcano
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  tt$sig <- "NS"
  tt$sig[tt$adj.P.Val < PADJ_DEG & tt$logFC > LFC_DEG] <- "Up"
  tt$sig[tt$adj.P.Val < PADJ_DEG & tt$logFC < -LFC_DEG] <- "Down"
  tt$sig[tt$gene %in% farg_deg$gene] <- "FARG"
  g <- ggplot(tt, aes(logFC, -log10(P.Value), color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = c(NS = "grey80", Up = "#E64B35",
                                  Down = "#4DBBD5", FARG = "#3C5488")) +
    labs(title = "TCGA-LIHC DEGs (FARGs highlighted)", x = "log2FC",
         y = "-log10 P") +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  ggsave(file.path(OUT_DIR, "deg_volcano.png"), g, width = 7, height = 6,
         dpi = 200)
  cat("      saved deg_volcano.png\n")
}

## ---- 2. enrichment ---------------------------------------------------------------
cat("\n[2/7] GO/KEGG enrichment of FARG-DEGs ...\n")
if (nrow(farg_deg) >= 5 &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  tryCatch({
    ego <- enrichGO(gene = farg_deg$gene, OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    p1 <- file.path(OUT_DIR, "enrichment_GO.png")
    png(p1, width = 1600, height = 1100, res = 150)
    print(barplot(ego, showCategory = 12) + ggtitle("GO-BP: FARG-DEGs"))
    dev.off()
    cat("      saved", basename(p1), "(", nrow(as.data.frame(ego)), "terms )\n")
  }, error = function(e) cat("      GO enrichment skipped:", conditionMessage(e), "\n"))
} else {
  cat("      too few FARG-DEGs or packages missing; enrichment skipped.\n")
}

## ---- 3. univariate Cox ---------------------------------------------------------
cat("\n[3/7] Univariate Cox on FARG-DEGs (p <", UNI_COX_P, ") ...\n")
tum_s <- intersect(colnames(expr), surv$sample)
expr_t <- expr[, tum_s]
sv <- surv[match(tum_s, sample)]
uni <- lapply(farg_deg$gene, function(gene) {
  dd <- data.frame(time = sv$OS.time, status = sv$OS,
                   x = as.numeric(scale(expr_t[gene, ])))
  dd <- dd[dd$time > 0 & !is.na(dd$x), ]
  fit <- tryCatch(coxph(Surv(time, status) ~ x, data = dd),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.table(gene = gene, HR = s$conf.int[1, 1], p = s$coefficients[1, 5])
})
uni <- rbindlist(uni[!vapply(uni, is.null, logical(1))])
uni <- uni[p < UNI_COX_P]
cat("      genes passing univariate filter:", nrow(uni), "\n")
print(uni[order(p)])

## ---- 4. LASSO + multivariate Cox -----------------------------------------------
cat("\n[4/7] LASSO + stepwise multivariate Cox ...\n")
if (nrow(uni) >= 3) {
  set.seed(SEED)
  ok <- !is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS)   # valid survival
  x <- scale(t(expr_t[uni$gene, ok]))
  y <- Surv(sv$OS.time[ok], sv$OS[ok])
  cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
  coefs <- as.numeric(coef(cv, s = "lambda.min"))       # glmnet 5.x: vector
  sel <- uni$gene[coefs != 0]
  if (length(sel) == 0) {
    coefs <- as.numeric(coef(cv, s = "lambda.1se"))
    sel <- uni$gene[coefs != 0]
    cat("      lambda.min gave no genes; used lambda.1se.\n")
  }
  if (length(sel) == 0) { sel <- uni$gene; cat("      fallback: all uni genes.\n") }
  cat("      LASSO selected:", paste(sel, collapse = ", "), "\n")

  dd <- data.frame(time = sv$OS.time[ok], status = sv$OS[ok],
                   x = x[, sel, drop = FALSE])
  dd <- dd[dd$time > 0, ]
  f_full <- coxph(Surv(time, status) ~ ., data = dd)
  f_final <- step(f_full, direction = "both", trace = 0)
  model_genes <- sub("^x\\.", "", names(coef(f_final)))   # strip data.frame prefix
  coefs_final <- coef(f_final)
  names(coefs_final) <- model_genes
  cat("      final model genes:", paste(model_genes, collapse = ", "), "\n")
  s <- summary(f_final)
  sig_tab <- data.table(gene = model_genes,
                        coef = coefs_final,
                        HR = s$conf.int[, 1],
                        HR_low = s$conf.int[, 3],
                        HR_high = s$conf.int[, 4],
                        p = s$coefficients[, 5])
  fwrite(sig_tab, file.path(OUT_DIR, "signature_genes.csv"))
  print(sig_tab[order(p)])

  # risk score for ALL tumor samples (z-scored per gene)
  z_full <- scale(t(expr_t[model_genes, , drop = FALSE]))
  risk <- as.numeric(z_full %*% coefs_final)
} else {
  cat("      too few genes after univariate filter; using all FARG-DEGs\n")
  ok <- !is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS)
  dd <- data.frame(time = sv$OS.time[ok], status = sv$OS[ok],
                   x = scale(t(expr_t[farg_deg$gene, ok])))
  dd <- dd[dd$time > 0, ]
  f_final <- coxph(Surv(time, status) ~ ., data = dd)
  model_genes <- farg_deg$gene
  z_full <- scale(t(expr_t[model_genes, , drop = FALSE]))
  risk <- as.numeric(z_full %*% coef(f_final))
  sig_tab <- data.table(gene = model_genes, coef = coef(f_final),
                        HR = exp(coef(f_final)), HR_low = NA_real_,
                        HR_high = NA_real_, p = NA_real_)
  fwrite(sig_tab, file.path(OUT_DIR, "signature_genes.csv"))
}

risk_dt <- data.table(sample = tum_s, risk = risk)
risk_dt[, risk_group := ifelse(risk >= median(risk), "High", "Low")]
risk_dt <- merge(risk_dt, sv, by = "sample")
risk_dt <- risk_dt[!is.na(OS.time) & OS.time > 0]
if (!is.null(pheno) && "sample" %in% names(pheno)) {
  risk_dt <- merge(risk_dt, pheno, by = "sample", all.x = TRUE)
}
fwrite(risk_dt, file.path(OUT_DIR, "risk_training.csv"))

## ---- 5. KM + time-dependent AUC --------------------------------------------------
cat("\n[5/7] KM and time-dependent AUC (training) ...\n")
if (requireNamespace("survminer", quietly = TRUE)) {
  suppressPackageStartupMessages(library(survminer))
  fit_km <- survfit(Surv(OS.time, OS) ~ risk_group, data = risk_dt)
  p_km <- file.path(OUT_DIR, "km_training.png")
  g_km <- ggsurvplot(fit_km, data = risk_dt, pval = TRUE, risk.table = TRUE,
                     palette = c("#B2182B", "#2166AC"), xlab = "Time (days)",
                     legend.title = "Risk", title = "Training: TCGA-LIHC")
  png(p_km, width = 1500, height = 1200, res = 150)
  print(g_km, newpage = FALSE); dev.off()
  cat("      saved km_training.png\n")
}
if (requireNamespace("timeROC", quietly = TRUE)) {
  suppressPackageStartupMessages(library(timeROC))
  t_roc <- timeROC(T = risk_dt$OS.time, delta = risk_dt$OS,
                   marker = risk_dt$risk, cause = 1,
                   times = c(1, 3, 5) * 365, iid = FALSE)
  p_roc <- file.path(OUT_DIR, "auc_training.png")
  png(p_roc, width = 1200, height = 1000, res = 150)
  plotAUCcurve(t_roc, FP = 2, add = FALSE, col = c("#4DBBD5", "#7E6148",
                                                   "#B2182B"))
  title(main = "Training (TCGA-LIHC): time-dependent AUC")
  legend("bottomright",
         legend = paste(c("1y", "3y", "5y"), "AUC =", round(t_roc$AUC, 3)),
         col = c("#4DBBD5", "#7E6148", "#B2182B"), lty = 1, lwd = 2, bty = "n")
  cat("      saved auc_training.png | time-dependent AUC(1/3/5y):",
      paste(round(t_roc$AUC, 3), collapse = " / "), "\n")
  dev.off()
}

# risk distribution
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  rp <- risk_dt[order(risk)]
  rp[, rank := 1:.N]
  g_risk <- ggplot(rp, aes(rank, risk, color = risk_group, shape = factor(OS))) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
    labs(title = "Risk score distribution (training)", x = "Sample rank",
         y = "FARGs risk score", color = "Group", shape = "OS") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUT_DIR, "risk_dist.png"), g_risk, width = 8, height = 5,
         dpi = 200)
  cat("      saved risk_dist.png\n")
}

## ---- 6. nomogram + calibration -----------------------------------------------------
cat("\n[6/7] Nomogram + calibration ...\n")
if (requireNamespace("rms", quietly = TRUE)) {
  suppressPackageStartupMessages(library(rms))
  has_age <- "age" %in% names(risk_dt)
  has_stage <- "stage" %in% names(risk_dt)
  nom_vars <- c("risk", if (has_age) "age", if (has_stage) "stage")
  nd <- as.data.frame(risk_dt[, c("OS.time", "OS", nom_vars), with = FALSE])
  nd <- nd[complete.cases(nd) & nd$OS.time > 0, ]
  if (has_stage) nd$stage <- factor(nd$stage)
  tryCatch({
    ddist <- datadist(nd); options(datadist = "ddist")
    f_nom <- cph(as.formula(paste("Surv(OS.time, OS) ~",
                                  paste(nom_vars, collapse = " + "))),
                 data = nd, x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)
    surv <- Survival(f_nom)          # rms: the predicted-survival function
    nom <- nomogram(f_nom,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel = c("1-year OS", "3-year OS", "5-year OS"))
    p_nom <- file.path(OUT_DIR, "nomogram.png")
    png(p_nom, width = 1600, height = 1100, res = 150)
    plot(nom)
    dev.off()
    cat("      saved nomogram.png\n")
    export_cal <- function(cc, uu) {
      pts <- data.table(horizon_days = uu,
                        mean.predicted = as.numeric(cc[, "mean.predicted"]),
                        KM = as.numeric(cc[, "KM"]))
      fwrite(pts, file.path(OUT_DIR, paste0("calibration_points_", uu, "d.csv")))
      cc
    }
    cal_list <- lapply(c(365, 1095), function(uu)
      tryCatch(export_cal(calibrate(f_nom, cmethod = "KM", u = uu, m = 60, B = 200), uu),
               error = function(e) {
                 # m=60 偶发分组不足，降 m 再试一次
                 tryCatch(export_cal(calibrate(f_nom, cmethod = "KM", u = uu,
                                               m = 40, B = 200), uu),
                          error = function(e2) {
                            cat("      calibrate u =", uu, "skipped:",
                                conditionMessage(e2), "\n"); NULL })
               }))
    cal_list <- cal_list[!vapply(cal_list, is.null, logical(1))]
    if (length(cal_list)) {
      p_cal <- file.path(OUT_DIR, "calibration.png")
      png(p_cal, width = 1400, height = 700, res = 150)
      par(mfrow = c(1, length(cal_list)))
      for (cc in cal_list) plot(cc, subtitles = FALSE)
      par(mfrow = c(1, 1))
      dev.off()
      cat("      saved calibration.png (", length(cal_list), "time points)\n",
          sep = "")
    }
  }, error = function(e)
    cat("      nomogram skipped:", conditionMessage(e), "\n"))
}

## ---- 7. summary ----------------------------------------------------------------------
cat("\n[7/7] == step05 summary ===========================================\n")
cat("Mode:", if (SELF_TEST) "SELF-TEST (synthetic)" else "REAL (TCGA-LIHC)", "\n")
cat("DEGs:", nrow(deg), "| FARG-DEGs:", nrow(farg_deg), "\n")
cat("Signature genes:", paste(model_genes, collapse = ", "), "\n")
cat("Files written to:", OUT_DIR, "\n")
cat("==================================================================\n")
