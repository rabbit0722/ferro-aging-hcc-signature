# =============================================================================
# rerun_bootstrap_only.R —— 仅重跑 step09 的 bootstrap 内部验证（条件内部验证（候选基因用全体样本筛选，LASSO+拟合在重采样内重做））
# 目的：绕开 step09 加载 GSVA/Biostrings 时的页文件内存不足问题。
# 只依赖 data.table / survival / glmnet / limma（均轻量）。
# 输出：bootstrap_cindex.png + bootstrap_internal_cindex.csv（修正后的 C-index）
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet); library(limma)
})
OUT_DIR <- "analysis_out/hcc_out/robustness"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BOOT_N <- 1000; UNI_COX_P <- 0.1; PADJ_DEG <- 0.05; LFC_DEG <- 1; SEED <- 42

sig <- fread("analysis_out/hcc_out/signature_genes.csv")
coefs <- setNames(sig$coef, sig$gene)
cat("Signature:", paste(names(coefs), collapse = ", "), "\n")

L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr <- L$expr; group <- L$group; surv <- L$surv
farg_genes <- unique(fread("analysis_out/step01/module_assignment.csv")$SYMBOL)

# candidate pool（与 step09 一致：FARG-DEG + 单因素 Cox P<0.1）
g_f <- factor(group, levels = c("Normal", "Tumor"))
fit <- eBayes(lmFit(expr, model.matrix(~ g_f)))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
fdeg <- tt[tt$gene %in% farg_genes & tt$adj.P.Val < PADJ_DEG &
             abs(tt$logFC) > LFC_DEG, ]
tum_s <- intersect(colnames(expr), surv$sample)
expr_t <- expr[, tum_s]
sv <- surv[match(tum_s, sample)]
uni <- rbindlist(lapply(fdeg$gene, function(gene) {
  dd <- data.frame(time = sv$OS.time, status = sv$OS,
                   x = as.numeric(scale(expr_t[gene, ])))
  dd <- dd[dd$time > 0 & !is.na(dd$x), ]
  f <- tryCatch(coxph(Surv(time, status) ~ x, data = dd), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  data.table(gene = gene, p = summary(f)$coefficients[1, 5])
}))
pool <- uni[p < UNI_COX_P]$gene
cat("candidate pool:", length(pool), "genes:", paste(pool, collapse = ", "), "\n")

# Conditional internal validation: the full-data candidate screen stays fixed;
# LASSO and Cox fitting are repeated, but stepwise reduction is not reproduced.
set.seed(SEED)
ok <- !is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS)
X_raw <- t(expr_t[pool, ok, drop = FALSE])
time_ok <- sv$OS.time[ok]; status_ok <- sv$OS[ok]
n_ok <- nrow(X_raw)
sel_freq <- setNames(rep(0, length(pool)), pool)
c_idx <- rep(NA_real_, BOOT_N)   # NA for skipped/failed iterations (excluded by na.rm)
iteration_log <- data.table(iteration = seq_len(BOOT_N), n_oob = NA_integer_,
                            n_selected = NA_integer_, selected_genes = NA_character_,
                            c_index = NA_real_, status = "not_evaluated")
cat("running", BOOT_N, "bootstrap iterations ...\n")
for (b in seq_len(BOOT_N)) {
  idx <- sample(seq_len(n_ok), replace = TRUE)
  oob <- setdiff(seq_len(n_ok), unique(idx))
  xb <- scale(X_raw[idx, , drop = FALSE])
  mns <- attr(xb, "scaled:center"); sds <- attr(xb, "scaled:scale")
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
  if (length(sel) > 0 && length(oob) >= 10) {
    db <- as.data.frame(xb[, sel, drop = FALSE])
    db$time <- time_ok[idx]; db$status <- status_ok[idx]
    fb <- tryCatch(coxph(Surv(time, status) ~ ., data = db), error = function(e) NULL)
    if (!is.null(fb)) {
      coef_b <- coef(fb)
      xo <- scale(X_raw[oob, sel, drop = FALSE], center = mns[sel], scale = sds[sel])
      risk_oob <- as.numeric(xo %*% coef_b)
      dd <- data.frame(time = time_ok[oob], status = status_ok[oob], risk = risk_oob)
      # direct concordance of the FIXED risk score (reverse=TRUE: higher risk = worse)
      c_idx[b] <- tryCatch(
        concordance(Surv(time, status) ~ risk, data = dd, reverse = TRUE)$concordance,
        error = function(e) NA)
    }
  }
  iteration_log[b, `:=`(n_oob = length(oob), n_selected = length(sel),
                        selected_genes = paste(sel, collapse = ";"),
                        c_index = c_idx[b],
                        status = if (is.finite(c_idx[b])) "ok" else "not_estimable")]
  if (b %% 100 == 0) cat("completed", b, "iterations\n")
}
fwrite(iteration_log, file.path(OUT_DIR, "bootstrap_iterations.csv"))
writeLines(c(paste("seed:", SEED), paste("BOOT_N:", BOOT_N),
             paste("candidate_pool:", paste(pool, collapse = ", ")),
             "validation: fixed full-data candidate screen; repeated LASSO+fit; no stepwise",
             "interval: empirical 2.5th-97.5th percentiles of out-of-bag C-indices; not a confidence interval",
             capture.output(sessionInfo())), file.path(OUT_DIR, "bootstrap_run_info.txt"))

freq_dt <- data.table(gene = names(sel_freq), freq = sel_freq / BOOT_N)[order(-freq)]
fwrite(freq_dt, file.path(OUT_DIR, "bootstrap_selection_freq.csv"))
cat("\nselection frequency:\n"); print(freq_dt)

res <- data.table(
  apparent_C_index = NA_real_,
  heldout_C_mean = mean(c_idx, na.rm = TRUE),
  heldout_C_q025 = quantile(c_idx, 0.025, na.rm = TRUE),
  heldout_C_q975 = quantile(c_idx, 0.975, na.rm = TRUE),
  n_valid_iterations = sum(!is.na(c_idx)),
  BOOT_N = BOOT_N)
fwrite(res, file.path(OUT_DIR, "bootstrap_internal_cindex.csv"))
cat(sprintf("\nConditional held-out C-index: mean %.3f, empirical 2.5th-97.5th percentiles [%.3f, %.3f] (valid iter %d/%d)\n",
            res$heldout_C_mean, res$heldout_C_q025, res$heldout_C_q975,
            res$n_valid_iterations, BOOT_N))

if (!"--no-plot" %in% commandArgs(trailingOnly = TRUE)) {
png(file.path(OUT_DIR, "bootstrap_cindex.png"), width = 1100, height = 800, res = 150)
hist(c_idx, breaks = 30, col = "#B2182B", border = "white",
     xlab = "held-out C-index (conditional (candidate screen fixed; LASSO+fit refit))",
     main = "Bootstrap held-out C-index (conditional (candidate screen fixed; LASSO+fit refit))")
abline(v = mean(c_idx, na.rm = TRUE), lty = 2, lwd = 2)
dev.off()
}
cat("done.\n")
