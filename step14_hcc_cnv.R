# =============================================================================
# step14_hcc_cnv.R —— 签名基因拷贝数变异（CNV）分析
# -----------------------------------------------------------------------------
# 数据：analysis_out/data/TCGA.LIHC.Gistic2_thresholded_by_genes.gz
#       （UCSC Xena TCGA-LIHC GISTIC2 基因水平阈值化 CNV，370 肿瘤样本，0.4 MB）
#       GISTIC 取值：-2 纯合缺失 / -1 杂合缺失 / 0 二倍体 / 1 获得 / 2 扩增
#
# 分析内容（对 HMOX1 / TFRC / SOCS2 / E2F1）：
#   [C] CNV 分布（获得/缺失/扩增/纯合缺失频率）
#   [D] 高/低风险组 CNV 差异（Fisher + Wilcoxon）
#   [E] CNV ↔ 表达 相关（Spearman）
#   [F] CNV 与生存（Cox on GISTIC score + 按 CNV 类别的 KM）
#   [G] 论文草稿
#
# 输出：analysis_out/hcc_out/cnv/
# 运行：setwd("path/to/Ferro-aging"); source("step14_hcc_cnv.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
OUT <- "analysis_out/hcc_out/cnv"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== step14: signature gene CNV ==\n")

GENES <- c("HMOX1", "TFRC", "SOCS2", "E2F1")
CNV_FILE <- "analysis_out/data/TCGA.LIHC.Gistic2_thresholded_by_genes.gz"

## ---- [A] 数据准备 ---------------------------------------------------------------
if (!file.exists(CNV_FILE)) {
  dir.create(dirname(CNV_FILE), showWarnings = FALSE, recursive = TRUE)
  url <- paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/",
                "TCGA.LIHC.sampleMap%2FGistic2_CopyNumber_Gistic2_all_",
                "thresholded.by_genes.gz")
  cat("下载 CNV 数据（约 0.5 MB）...\n")
  try(download.file(url, CNV_FILE, mode = "wb", quiet = TRUE), silent = TRUE)
}
if (!file.exists(CNV_FILE)) {
  cat("!! CNV 文件不存在且下载失败。请手动下载后放到：", CNV_FILE, "\n")
  cat("   URL:", paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/",
                         "TCGA.LIHC.sampleMap%2FGistic2_CopyNumber_Gistic2_all_",
                         "thresholded.by_genes.gz"), "\n")
  stop("missing CNV file")
}

cnv <- fread(CNV_FILE)
setnames(cnv, 1, "gene")
cat("CNV 矩阵:", nrow(cnv), "基因 x", ncol(cnv) - 1, "样本\n")

missing_g <- setdiff(GENES, cnv$gene)
if (length(missing_g))
  cat("!! 以下基因在 CNV 文件中缺失：", paste(missing_g, collapse = ", "), "\n")
GENES <- intersect(GENES, cnv$gene)

sub <- cnv[gene %in% GENES]
cmat <- as.matrix(sub[, -1])
rownames(cmat) <- sub$gene
storage.mode(cmat) <- "numeric"
cnv_samples <- colnames(cmat)
cat("目标基因:", paste(GENES, collapse = ", "), "\n")
cat("样本类型（01=肿瘤）:",
    paste(names(table(substr(cnv_samples, 14, 15))),
          table(substr(cnv_samples, 14, 15)), sep = "=", collapse = " "), "\n")

L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_t <- L$expr[, L$group == "Tumor"]
grp_all <- L$group
rk <- fread("analysis_out/hcc_out/risk_training.csv")

## ---- [B] CNV 分布 ---------------------------------------------------------------
cat("\n[B] CNV class distribution ...\n")
cls_tab <- rbindlist(lapply(GENES, function(g) {
  v <- cmat[g, ]
  n <- sum(!is.na(v))
  data.table(gene = g, n = n,
             freq_amp = sum(v == 2, na.rm = TRUE) / n * 100,
             freq_gain = sum(v == 1, na.rm = TRUE) / n * 100,
             freq_diploid = sum(v == 0, na.rm = TRUE) / n * 100,
             freq_loss = sum(v == -1, na.rm = TRUE) / n * 100,
             freq_homdel = sum(v == -2, na.rm = TRUE) / n * 100,
             pct_altered = sum(abs(v) >= 1, na.rm = TRUE) / n * 100)
}))
print(cls_tab)

## ---- [C] 高/低风险组 CNV 差异 ----------------------------------------------------
cat("\n[C] CNV by risk group ...\n")
grp_tab <- rbindlist(lapply(GENES, function(g) {
  v <- cmat[g, ]
  gg <- rk$risk_group[match(cnv_samples, rk$sample)]
  okk <- !is.na(gg) & !is.na(v)
  altered <- abs(v[okk]) >= 1
  tb <- table(gg[okk], altered)
  fis <- tryCatch(fisher.test(tb)$p.value, error = function(e) NA)
  wt <- tryCatch(wilcox.test(v[okk][gg[okk] == "High"],
                             v[okk][gg[okk] == "Low"])$p.value,
                 error = function(e) NA)
  data.table(gene = g,
             pct_altered_High = mean(altered[gg[okk] == "High"]) * 100,
             pct_altered_Low = mean(altered[gg[okk] == "Low"]) * 100,
             fisher_p_altered = fis,
             median_cnv_High = median(v[okk][gg[okk] == "High"]),
             median_cnv_Low = median(v[okk][gg[okk] == "Low"]),
             wilcox_p_cnv = wt)
}))
grp_tab[, FDR_wilcox := p.adjust(wilcox_p_cnv, "BH")]
print(grp_tab)

## ---- [D] CNV ↔ 表达 相关 ---------------------------------------------------------
cat("\n[D] CNV vs expression (Spearman) ...\n")
de <- rbindlist(lapply(GENES, function(g) {
  if (!g %in% rownames(expr_t)) return(NULL)
  smp <- intersect(cnv_samples, colnames(expr_t))
  x <- cmat[g, smp]; y <- as.numeric(expr_t[g, smp])
  okk <- !is.na(x) & !is.na(y)
  ct <- suppressWarnings(cor.test(x[okk], y[okk], method = "spearman"))
  data.table(gene = g, n_pairs = sum(okk),
             rho_cnv_expr = as.numeric(ct$estimate), p_cnv_expr = ct$p.value)
}))
print(de)

## ---- [E] CNV 与生存 -------------------------------------------------------------
cat("\n[E] CNV and overall survival ...\n")
ds <- rbindlist(lapply(GENES, function(g) {
  smp <- intersect(cnv_samples, rk$sample)
  idx <- match(smp, rk$sample)
  d <- data.table(cnv = cmat[g, smp], time = rk$OS.time[idx], status = rk$OS[idx])
  d <- d[!is.na(cnv) & !is.na(time) & time > 0]
  if (nrow(d) < 50 || sum(d$status) < 10) return(NULL)
  # (1) 有序 GISTIC score 的 Cox（每单位）
  f <- coxph(Surv(time, status) ~ cnv, data = d)
  s <- summary(f)
  # (2) 改变 vs 二倍体
  d[, altered := abs(cnv) >= 1]
  f2 <- coxph(Surv(time, status) ~ altered, data = d)
  s2 <- summary(f2)
  # (3) 三分类 KM
  d[, cls := fifelse(cnv >= 1, "Gain/Amp",
                     fifelse(cnv <= -1, "Loss", "Diploid"))]
  lr <- survdiff(Surv(time, status) ~ cls, data = d)
  lr_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  data.table(gene = g, n = nrow(d), events = sum(d$status),
             HR_per_unit_CNV = s$conf.int["cnv", "exp(coef)"],
             HR_unit_low = s$conf.int["cnv", "lower .95"],
             HR_unit_high = s$conf.int["cnv", "upper .95"],
             p_cox_unit = s$coefficients["cnv", "Pr(>|z|)"],
             HR_altered = s2$conf.int["alteredTRUE", "exp(coef)"],
             p_cox_altered = s2$coefficients["alteredTRUE", "Pr(>|z|)"],
             logrank_p_3class = lr_p)
}))
if (nrow(ds)) {
  ds[, FDR_cox_unit := p.adjust(p_cox_unit, "BH")]
}
print(ds)

## ---- [F] 汇总表 -----------------------------------------------------------------
summ <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
               Filter(Negate(is.null), list(cls_tab, grp_tab, de, ds)))
fwrite(summ, file.path(OUT, "cnv_gene_summary.csv"))

## ---- [G] 每个基因的三联图 ---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  has_ge <- requireNamespace("gridExtra", quietly = TRUE)
  has_sv <- requireNamespace("survminer", quietly = TRUE)
  for (g in GENES) {
    v <- cmat[g, ]
    cls_levels <- c("Homdel(-2)", "Loss(-1)", "Diploid(0)", "Gain(1)", "Amp(2)")
    d1 <- data.table(
      cls = factor(c("Homdel(-2)", "Loss(-1)", "Diploid(0)", "Gain(1)", "Amp(2)"),
                   levels = cls_levels),
      pct = c(mean(v == -2, na.rm = TRUE), mean(v == -1, na.rm = TRUE),
              mean(v == 0, na.rm = TRUE), mean(v == 1, na.rm = TRUE),
              mean(v == 2, na.rm = TRUE)) * 100)
    g1 <- ggplot(d1, aes(cls, pct, fill = cls)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = c("#2166AC", "#92C5DE", "grey80",
                                   "#F4A582", "#B2182B")) +
      labs(x = NULL, y = "Patients (%)", title = g,
           subtitle = "GISTIC2 thresholded CNV (TCGA-LIHC)") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1))
    # CNV vs 表达
    smp <- intersect(cnv_samples, colnames(expr_t))
    d2 <- data.table(cnv = cmat[g, smp], expr = as.numeric(expr_t[g, smp]))
    d2 <- d2[!is.na(cnv) & !is.na(expr)]
    d2[, cls := factor(fifelse(cnv >= 1, "Gain/Amp",
                               fifelse(cnv <= -1, "Loss", "Diploid")),
                       levels = c("Loss", "Diploid", "Gain/Amp"))]
    g2 <- ggplot(d2, aes(cls, expr, fill = cls)) +
      geom_boxplot(outlier.size = 0.3, width = 0.6) +
      scale_fill_manual(values = c(Loss = "#2166AC", Diploid = "grey75",
                                   `Gain/Amp` = "#B2182B")) +
      labs(x = NULL, y = paste0(g, " (log2)"),
           subtitle = sprintf("CNV vs expression rho = %.2f, P = %s",
                              de[gene == g]$rho_cnv_expr,
                              format.pval(de[gene == g]$p_cnv_expr, digits = 3))) +
      theme_bw(base_size = 10) + theme(legend.position = "none")
    # KM
    g3 <- NULL
    if (has_sv) {
      idx <- match(smp, rk$sample)
      d3 <- data.table(cnv = cmat[g, smp], time = rk$OS.time[idx],
                       status = rk$OS[idx])
      d3 <- d3[!is.na(cnv) & !is.na(time) & time > 0]
      d3[, cls := fifelse(cnv >= 1, "Gain/Amp",
                          fifelse(cnv <= -1, "Loss", "Diploid"))]
      d3[, cls := factor(cls, levels = c("Loss", "Diploid", "Gain/Amp"))]
      fit <- survfit(Surv(time, status) ~ cls, data = d3)
      kp <- survminer::ggsurvplot(fit, data = d3, pval = TRUE,
                                  risk.table = FALSE,
                                  palette = c("#2166AC", "grey60", "#B2182B"),
                                  legend.title = "", ggtheme = theme_bw(base_size = 10))
      g3 <- kp$plot +
        labs(subtitle = sprintf(
          "Cox altered vs diploid: HR = %.2f (P = %s); per unit: HR = %.2f (P = %s)",
          ds[gene == g]$HR_altered,
          format.pval(ds[gene == g]$p_cox_altered, digits = 3),
          ds[gene == g]$HR_per_unit_CNV,
          format.pval(ds[gene == g]$p_cox_unit, digits = 3)))
    }
    if (has_ge && !is.null(g3)) {
      png(file.path(OUT, paste0("cnv_panel_", g, ".png")),
          width = 2200, height = 750, res = 150)
      gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
      dev.off()
      cat("  saved cnv_panel_", g, ".png\n", sep = "")
    } else {
      ggsave(file.path(OUT, paste0("cnv_", g, "_distribution.png")), g1,
             width = 4.5, height = 4, dpi = 200)
      ggsave(file.path(OUT, paste0("cnv_", g, "_vs_expr.png")), g2,
             width = 4.5, height = 4, dpi = 200)
      if (!is.null(g3))
        ggsave(file.path(OUT, paste0("cnv_", g, "_KM.png")), g3,
               width = 5, height = 4, dpi = 200)
    }
  }
}

## ---- [H] 论文草稿 ---------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
md <- c(
  "# 签名基因 CNV 分析（step14 自动生成）",
  "",
  sprintf("数据：UCSC Xena TCGA-LIHC GISTIC2 基因水平阈值化 CNV（%d 样本）",
          ncol(cmat)),
  "",
  "## 一、CNV 分布",
  paste0("| 基因 | n | Amplification | Gain | Diploid | Loss | Homdel | 变异比例 |"),
  paste0("|---|---|---|---|---|---|---|---|"),
  apply(cls_tab, 1, function(r) sprintf(
    "| %s | %s | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.1f%% |",
    r[["gene"]], r[["n"]], as.numeric(r[["freq_amp"]]),
    as.numeric(r[["freq_gain"]]), as.numeric(r[["freq_diploid"]]),
    as.numeric(r[["freq_loss"]]), as.numeric(r[["freq_homdel"]]),
    as.numeric(r[["pct_altered"]]))),
  "",
  "## 二、CNV 与表达 / 生存",
  paste0("| 基因 | 变异比例 | 高/低风险组变异率 | Fisher P | rho(CNV↔表达) | P | HR(每单位CNV) | P(Cox) | HR(变异vs二倍体) | P | P(log-rank, 3类) |"),
  paste0("|---|---|---|---|---|---|---|---|---|---|---|"),
  if (nrow(summ)) apply(summ, 1, function(r) sprintf(
    "| %s | %.1f%% | %.1f%% / %.1f%% | %s | %s | %s | %s | %s | %s | %s | %s |",
    r[["gene"]], as.numeric(r[["pct_altered"]]),
    as.numeric(r[["pct_altered_High"]]), as.numeric(r[["pct_altered_Low"]]),
    format.pval(as.numeric(r[["fisher_p_altered"]]), digits = 3),
    fmt(as.numeric(r[["rho_cnv_expr"]])),
    format.pval(as.numeric(r[["p_cnv_expr"]]), digits = 3),
    fmt(as.numeric(r[["HR_per_unit_CNV"]])),
    format.pval(as.numeric(r[["p_cox_unit"]]), digits = 3),
    fmt(as.numeric(r[["HR_altered"]])),
    format.pval(as.numeric(r[["p_cox_altered"]]), digits = 3),
    format.pval(as.numeric(r[["logrank_p_3class"]]), digits = 3))) else "（无）",
  "",
  "## 三、论文 Methods（英文草稿）",
  paste0("> \"Gene-level copy-number data for TCGA-LIHC (GISTIC2 thresholded calls: ",
         "-2 homozygous deletion to +2 high-level amplification) were obtained from ",
         "the UCSC Xena hub. Alteration frequencies were calculated across tumour ",
         "samples; associations with the FARGs risk group were tested by Fisher exact ",
         "and Wilcoxon tests; associations with gene expression were assessed by ",
         "Spearman correlation; and associations with overall survival by Cox ",
         "regression (per unit of GISTIC score) and Kaplan-Meier analysis ",
         "stratified by CNV class, with Benjamini-Hochberg correction across the ",
         "four signature genes.\""),
  "",
  "## 四、论文 Results（英文草稿，按实际数字替换）",
  paste0("> \"E2F1 copy-number gain and amplification were observed in ",
         fmt(cls_tab[gene == "E2F1"]$freq_gain, 1), "% and ",
         fmt(cls_tab[gene == "E2F1"]$freq_amp, 1), "% of tumours, respectively; ",
         "the total alteration frequency, including losses, was ",
         fmt(cls_tab[gene == "E2F1"]$pct_altered, 1), "%. Copy-number scores ",
         "were correlated with E2F1 expression (Spearman rho = ",
         fmt(de[gene == "E2F1"]$rho_cnv_expr), ", P = ",
         format.pval(de[gene == "E2F1"]$p_cnv_expr, digits = 3), "). Survival ",
         "results must be reported separately for per-unit GISTIC, altered-versus-diploid, ",
         "and three-class comparisons, with their actual P values; they do not ",
         "support a blanket survival-association claim.\""),
  "",
  "## 说明",
  "- GISTIC2 阈值化取值：-2 纯合缺失、-1 杂合缺失、0 二倍体、1 获得、2 扩增",
  "- CNV 数据仅含肿瘤样本（GISTIC 基于肿瘤谱计算），故无肿瘤 vs 正常比较",
  "- 关联方向以实际结果为准，不能由风险系数预设CNV必然影响表达或生存；全部变异比例包含丢失，不等同于获得/扩增。",
  "",
  "## 输出文件",
  "- cnv_gene_summary.csv（4 基因的分布/组间/表达/生存全指标）",
  "- cnv_panel_<GENE>.png（CNV 分布 / CNV-表达 / KM 三联图）"
)
writeLines(md, file.path(OUT, "cnv_methods_and_results.md"))

cat("\n== step14 done ==\n")
cat("outputs ->", OUT, "\n")
