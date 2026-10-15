# =============================================================================
# step15_hcc_acsl4.R —— 源文献靶点 ACSL4 与 FARGs 风险轴的关联分析
# -----------------------------------------------------------------------------
# 科学问题：Cell Metab 2026 报道维生素C 通过抑制 ACSL4 缓解 ferro-aging。
#           那么 ACSL4 在 HCC 中是否与本研究的 ferro-aging 风险轴相关？
#           即：维C 的已知靶点是否落在我们签名的生物学通路上。
#
# 数据（全部本地，无需联网）：
#   analysis_out/hcc_data/hcc_data.rds        表达矩阵（11,206 基因 × 419 样本）
#   analysis_out/hcc_out/risk_training.csv    风险分数 + 生存（363 例）
#   analysis_out/step01/module_assignment.csv 95 FARGs 的模块归属
#
# 分析内容：
#   [B] ACSL4 是否属于 95 FARGs / 是否在表达矩阵中
#   [C] ACSL4 在肿瘤 vs 正常肝 的差异
#   [D] ACSL4 ↔ 风险分数（Spearman）+ 高/低风险组比较
#   [E] ACSL4 ↔ 4 个签名基因 + 铁死亡核心基因（GPX4/SLC7A11/ACSL3/HMOX1/TFRC）
#   [F] Legacy manual-module analysis (disabled for manuscript reporting)
#   [G] 预后价值：ACSL4 单基因 Cox/KM、与风险分数联合模型（增量检验 ΔC + LRT）
#   [H] 交互作用：ACSL4 三分位分层下的风险分组 KM
#   [I] 论文草稿
#
# 输出：analysis_out/hcc_out/acsl4/
# 运行：setwd("path/to/Ferro-aging"); source("step15_hcc_acsl4.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
OUT <- "analysis_out/hcc_out/acsl4"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== step15: ACSL4 (vitamin C target) vs FARGs risk axis ==\n")

TARGET <- "ACSL4"
RUN_LEGACY_MODULES <- FALSE  # Historical files remain available for audit only.
SIG <- c("HMOX1", "TFRC", "SOCS2", "E2F1")
FERRO_CORE <- c("GPX4", "SLC7A11", "ACSL3", "ACSL4", "FTH1", "FTL", "SLC40A1")

## ---- [A] 数据 -------------------------------------------------------------------
L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_all <- L$expr
grp_all <- L$group
names(grp_all) <- colnames(expr_all)
rk <- fread("analysis_out/hcc_out/risk_training.csv")
mods <- fread("analysis_out/step01/module_assignment.csv")
farg_all <- unique(mods$SYMBOL)

tum <- names(grp_all)[grp_all == "Tumor"]
nor <- names(grp_all)[grp_all == "Normal"]
cat("表达矩阵:", nrow(expr_all), "基因 |", length(tum), "肿瘤 +", length(nor), "正常\n")

## ---- [B] ACSL4 是否在 FARGs / 表达矩阵中 -----------------------------------------
cat("\n[B] ACSL4 定位 ...\n")
in_farg <- TARGET %in% farg_all
in_expr <- TARGET %in% rownames(expr_all)
cat("  ACSL4 属于 95 个 FARGs:", in_farg, "\n")
if (in_farg && RUN_LEGACY_MODULES)
  cat("  ACSL4 所属模块:", paste(mods[SYMBOL == TARGET]$Module, collapse = ", "), "\n")
cat("  ACSL4 存在于表达矩阵:", in_expr, "\n")
if (!in_expr) stop("ACSL4 不在表达矩阵中，无法继续")

# 顺便看看 ACSL 家族其他成员
fam <- grep("^ACSL[0-9]", rownames(expr_all), value = TRUE)
cat("  表达矩阵中的 ACSL 家族:", paste(fam, collapse = ", "), "\n")

## ---- [C] 肿瘤 vs 正常 ------------------------------------------------------------
cat("\n[C] ACSL4 tumour vs normal ...\n")
tv <- as.numeric(expr_all[TARGET, tum]); nv <- as.numeric(expr_all[TARGET, nor])
wt <- suppressWarnings(wilcox.test(tv, nv))
logfc <- mean(tv) - mean(nv)
cat(sprintf("  肿瘤 %.2f vs 正常 %.2f | log2FC = %+.2f | Wilcoxon P = %s\n",
            mean(tv), mean(nv), logfc, format.pval(wt$p.value, digits = 3)))
c_tab <- data.table(gene = TARGET, mean_tumor = mean(tv), mean_normal = mean(nv),
                    logFC = logfc, p_tumor_vs_normal = wt$p.value)

## ---- [D] ACSL4 ↔ 风险分数 -------------------------------------------------------
cat("\n[D] ACSL4 vs risk score ...\n")
smp <- intersect(rk$sample, tum)
idx <- match(smp, rk$sample)
a <- as.numeric(expr_all[TARGET, smp])
risk <- rk$risk[idx]; rgrp <- rk$risk_group[idx]
ok <- !is.na(a) & !is.na(risk)
ct <- suppressWarnings(cor.test(a[ok], risk[ok], method = "spearman"))
wt2 <- suppressWarnings(wilcox.test(a[ok][rgrp[ok] == "High"],
                                    a[ok][rgrp[ok] == "Low"]))
cat(sprintf("  Spearman rho = %.3f (P = %s)\n", as.numeric(ct$estimate),
            format.pval(ct$p.value, digits = 3)))
cat(sprintf("  高风险 %.2f vs 低风险 %.2f | Wilcoxon P = %s\n",
            mean(a[ok][rgrp[ok] == "High"]), mean(a[ok][rgrp[ok] == "Low"]),
            format.pval(wt2$p.value, digits = 3)))
d_tab <- data.table(gene = TARGET, rho_vs_risk = as.numeric(ct$estimate),
                    p_vs_risk = ct$p.value,
                    mean_High = mean(a[ok][rgrp[ok] == "High"]),
                    mean_Low = mean(a[ok][rgrp[ok] == "Low"]),
                    p_High_vs_Low = wt2$p.value)

## ---- [E] ACSL4 ↔ 签名基因 + 铁死亡核心基因 ----------------------------------------
cat("\n[E] ACSL4 vs signature / ferroptosis genes ...\n")
genes_e <- intersect(c(SIG, FERRO_CORE), rownames(expr_all))
e_tab <- rbindlist(lapply(setdiff(genes_e, TARGET), function(g) {
  y <- as.numeric(expr_all[g, smp])
  okk <- !is.na(y) & !is.na(a)
  cc <- suppressWarnings(cor.test(a[okk], y[okk], method = "spearman"))
  data.table(gene = g, rho_with_ACSL4 = as.numeric(cc$estimate),
             p = cc$p.value,
             in_signature = g %in% SIG)
}))
e_tab[, abs_rho := abs(rho_with_ACSL4)]
setorder(e_tab, -abs_rho)
print(e_tab)

## ---- [F] ACSL4 ↔ 功能基因集评分 ---------------------------------------------
cat("\n[F] ACSL4 vs FARGs module scores (ssGSEA) ...\n")
f_tab <- NULL
if (RUN_LEGACY_MODULES && requireNamespace("GSVA", quietly = TRUE)) {
  mod_list <- split(mods$SYMBOL, mods$Module)
  sets <- lapply(mod_list, function(g) intersect(g, rownames(expr_all)))
  sets <- sets[lengths(sets) >= 5]
  sc <- tryCatch({
    r <- GSVA::gsva(as.matrix(expr_all[, smp]), sets, method = "ssgsea",
                    ssgsea.norm = TRUE, verbose = FALSE)
    t(r)
  }, error = function(e) {
    tryCatch(t(GSVA::gsva(GSVA::ssgseaParam(as.matrix(expr_all[, smp]), sets),
                          verbose = FALSE)), error = function(e2) NULL)
  })
  if (!is.null(sc)) {
    f_tab <- rbindlist(lapply(colnames(sc), function(sn) {
      cc <- suppressWarnings(cor.test(a, as.numeric(sc[, sn]), method = "spearman"))
      data.table(module = sn, rho_with_ACSL4 = as.numeric(cc$estimate),
                 p = cc$p.value)
    }))
    f_tab[, abs_rho := abs(rho_with_ACSL4)]
    setorder(f_tab, -abs_rho)
    print(f_tab)
  }
} else cat("  Legacy manual-module analysis excluded from manuscript scope.\n")

## ---- [G] 预后价值 + 与风险分数的联合模型 ------------------------------------------
cat("\n[G] ACSL4 prognostic value & joint model with risk ...\n")
surv_d <- data.frame(time = rk$OS.time[idx], status = rk$OS[idx],
                     acsl4 = as.numeric(scale(a)), risk = as.numeric(scale(risk)))
surv_d <- surv_d[!is.na(surv_d$time) & surv_d$time > 0 &
                   !is.na(surv_d$status) & !is.na(surv_d$acsl4), ]
cat("  分析样本:", nrow(surv_d), "| 事件:", sum(surv_d$status), "\n")

f_acsl4 <- coxph(Surv(time, status) ~ acsl4, data = surv_d)
f_risk  <- coxph(Surv(time, status) ~ risk, data = surv_d)
f_both  <- coxph(Surv(time, status) ~ risk + acsl4, data = surv_d)

c_acsl4 <- concordance(f_acsl4)$concordance
c_risk  <- concordance(f_risk)$concordance
c_both  <- concordance(f_both)$concordance
lrt <- tryCatch(anova(f_risk, f_both)$`Pr(>|Chi|)`[2], error = function(e) NA)
lrt2 <- tryCatch(anova(f_acsl4, f_both)$`Pr(>|Chi|)`[2], error = function(e) NA)

g_tab <- data.table(
  model = c("ACSL4 only", "Risk only", "Risk + ACSL4"),
  C_index = c(c_acsl4, c_risk, c_both),
  delta_C_vs_risk = c(NA, NA, c_both - c_risk))
print(g_tab)
cat("  在风险分数之上加入 ACSL4：ΔC =", round(c_both - c_risk, 4),
    "| LRT P =", format.pval(lrt, digits = 3), "\n")
cat("  在 ACSL4 之上加入风险分数：ΔC =", round(c_both - c_acsl4, 4),
    "| LRT P =", format.pval(lrt2, digits = 3), "\n")

mv <- data.table(
  variable = rownames(summary(f_both)$coefficients),
  HR = summary(f_both)$conf.int[, "exp(coef)"],
  HR_low = summary(f_both)$conf.int[, "lower .95"],
  HR_high = summary(f_both)$conf.int[, "upper .95"],
  p = summary(f_both)$coefficients[, "Pr(>|z|)"])
print(mv)
fwrite(mv, file.path(OUT, "acsl4_joint_model_cox.csv"))

# ACSL4 中位分组的 KM
surv_d$acsl4_grp <- ifelse(surv_d$acsl4 > median(surv_d$acsl4),
                           "High ACSL4", "Low ACSL4")
lr <- survdiff(Surv(time, status) ~ acsl4_grp, data = surv_d)
cat("  ACSL4 中位分组 log-rank P =",
    format.pval(1 - pchisq(lr$chisq, 1), digits = 3), "\n")

## ---- [H] 交互作用：ACSL4 三分位分层下的风险效应 -------------------------------------
cat("\n[H] risk effect stratified by ACSL4 tertile ...\n")
surv_d$acsl4_t <- cut(surv_d$acsl4, breaks = quantile(surv_d$acsl4,
                                                      probs = c(0, 1/3, 2/3, 1)),
                      labels = c("Low", "Mid", "High"), include.lowest = TRUE)
surv_d$risk_grp <- ifelse(surv_d$risk > median(surv_d$risk), "High risk", "Low risk")
h_tab <- rbindlist(lapply(levels(surv_d$acsl4_t), function(lv) {
  d <- surv_d[surv_d$acsl4_t == lv, ]
  if (nrow(d) < 30 || sum(d$status) < 5) return(NULL)
  f <- coxph(Surv(time, status) ~ (risk_grp == "High risk"), data = d)
  s <- summary(f)
  data.table(ACSL4_tertile = lv, n = nrow(d), events = sum(d$status),
             HR_highRisk = s$conf.int[1, "exp(coef)"],
             HR_low = s$conf.int[1, "lower .95"],
             HR_high = s$conf.int[1, "upper .95"],
             p = s$coefficients[1, "Pr(>|z|)"])
}))
print(h_tab)
int_p <- tryCatch({
  f_int <- coxph(Surv(time, status) ~ risk * acsl4, data = surv_d)
  summary(f_int)$coefficients["risk:acsl4", "Pr(>|z|)"]
}, error = function(e) NA)
cat("  风险 × ACSL4 交互项 P =", format.pval(int_p, digits = 3), "\n")

## ---- [I] 汇总表 + 图 -------------------------------------------------------------
summ <- list(c_tab, d_tab, g_tab)
fwrite(c_tab, file.path(OUT, "acsl4_tumor_vs_normal.csv"))
fwrite(d_tab, file.path(OUT, "acsl4_vs_risk.csv"))
fwrite(e_tab, file.path(OUT, "acsl4_vs_ferroptosis_genes.csv"))
if (!is.null(f_tab)) fwrite(f_tab, file.path(OUT, "acsl4_vs_module_scores.csv"))
fwrite(g_tab, file.path(OUT, "acsl4_model_comparison.csv"))
if (!is.null(h_tab) && nrow(h_tab))
  fwrite(h_tab, file.path(OUT, "acsl4_tertile_risk_effect.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  # (1) ACSL4 vs 风险分数
  d1 <- data.table(acsl4 = a, risk = risk, grp = rgrp)
  d1 <- d1[!is.na(acsl4) & !is.na(risk)]
  g1 <- ggplot(d1, aes(acsl4, risk)) +
    geom_point(aes(color = grp), size = 0.9, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.8) +
    scale_color_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
    labs(x = "ACSL4 expression (log2)", y = "FARGs risk score", color = NULL,
         title = "ACSL4 vs ferro-aging risk score",
         subtitle = sprintf("Spearman rho = %.3f, P = %s",
                            as.numeric(ct$estimate),
                            format.pval(ct$p.value, digits = 3))) +
    theme_bw(base_size = 12)
  ggsave(file.path(OUT, "acsl4_vs_risk_scatter.png"), g1,
         width = 6, height = 5, dpi = 200)

  # (2) 肿瘤 vs 正常
  d2 <- data.table(expr = c(tv, nv),
                   type = factor(c(rep("Tumour", length(tv)),
                                   rep("Normal", length(nv))),
                                 levels = c("Normal", "Tumour")))
  g2 <- ggplot(d2, aes(type, expr, fill = type)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6) +
    scale_fill_manual(values = c(Tumour = "#B2182B", Normal = "#2166AC")) +
    labs(x = NULL, y = "ACSL4 (log2)", title = "ACSL4 in TCGA-LIHC",
         subtitle = sprintf("log2FC = %+.2f, P = %s", logfc,
                            format.pval(wt$p.value, digits = 3))) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  ggsave(file.path(OUT, "acsl4_tumor_vs_normal.png"), g2,
         width = 4.5, height = 4.5, dpi = 200)

  # (3) KM
  if (requireNamespace("survminer", quietly = TRUE)) {
    fit <- survfit(Surv(time, status) ~ acsl4_grp, data = surv_d)
    kp <- survminer::ggsurvplot(fit, data = surv_d, pval = TRUE,
                                risk.table = TRUE,
                                palette = c("#B2182B", "#2166AC"),
                                legend.title = "", ggtheme = theme_bw(base_size = 12))
    png(file.path(OUT, "km_acsl4.png"), width = 1300, height = 1100, res = 150)
    print(kp)
    dev.off()
  }

  # (4) ACSL4 与铁死亡/签名基因的相关
  if (nrow(e_tab)) {
    g4 <- ggplot(e_tab, aes(rho_with_ACSL4, reorder(gene, rho_with_ACSL4),
                            fill = in_signature)) +
      geom_col() + geom_vline(xintercept = 0, color = "grey40") +
      scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "signature gene", `FALSE` = "other")) +
      labs(x = "Spearman rho with ACSL4", y = NULL, fill = NULL,
           title = "ACSL4 vs ferroptosis / signature genes") +
      theme_bw(base_size = 12)
    ggsave(file.path(OUT, "acsl4_gene_correlations.png"), g4,
           width = 6, height = 5, dpi = 200)
  }
  cat("\n  saved figures\n")
}

## ---- [J] 论文草稿 ---------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
md <- c(
  "# ACSL4（维生素C 靶点）与 FARGs 风险轴（step15 自动生成）",
  "",
  sprintf("- ACSL4 是否属于 95 个 FARGs：**%s**", in_farg),
  "- 旧手工分模块相关性不纳入论文论证；原始输出仅作历史审计。",
  sprintf("- 肿瘤 vs 正常：log2FC = %+.2f，P = %s", logfc,
          format.pval(wt$p.value, digits = 3)),
  sprintf("- ACSL4 ↔ 风险分数：Spearman rho = %s，P = %s",
          fmt(as.numeric(ct$estimate)), format.pval(ct$p.value, digits = 3)),
  sprintf("- 高风险 vs 低风险：%.2f vs %.2f，P = %s",
          d_tab$mean_High, d_tab$mean_Low,
          format.pval(d_tab$p_High_vs_Low, digits = 3)),
  "",
  "## 模型比较（C-index）",
  paste0("| 模型 | C-index | 相对仅风险的 ΔC |"),
  paste0("|---|---|---|"),
  apply(g_tab, 1, function(r) sprintf("| %s | %s | %s |", r[["model"]],
                                      fmt(as.numeric(r[["C_index"]])),
                                      fmt(as.numeric(r[["delta_C_vs_risk"]])))),
  sprintf("- 在风险分数之上加入 ACSL4：ΔC = %s，LRT P = %s",
          fmt(c_both - c_risk), format.pval(lrt, digits = 3)),
  "",
  "## 论文 Methods（英文草稿）",
  paste0("> \"ACSL4, the reported target of ascorbate in ferro-aging, was evaluated in ",
         "TCGA-LIHC. Expression differences between tumours and adjacent normal liver ",
         "were tested by Wilcoxon rank-sum tests; associations with the FARGs risk score, ",
         "signature genes and ferroptosis-related genes were assessed by Spearman ",
         "correlation; and prognostic value was evaluated by Cox regression alone and in ",
         "combination with the risk score, with likelihood-ratio tests and comparison of ",
         "concordance indices. Risk-score effects were further examined within ACSL4 ",
         "tertiles, and an interaction term was tested.\""),
  "",
  "## 论文 Results（英文草稿）",
  sprintf(paste0("> \"ACSL4 expression was %s in tumours compared with adjacent normal ",
                 "liver (log2FC = %+.2f, P = %s). ACSL4 correlated %s with the FARGs ",
                 "risk score (Spearman rho = %s, P = %s) and was %s in high-risk tumours ",
                 "(P = %s). In Cox models, ACSL4 alone achieved a C-index of %s and ",
                 "adding it to the risk score changed the C-index by %s ",
                 "(likelihood-ratio P = %s), whereas adding the risk score to ACSL4 ",
                 "changed the C-index by %s (P = %s).\""),
          ifelse(logfc > 0, "higher", "lower"), logfc,
          format.pval(wt$p.value, digits = 3),
          ifelse(as.numeric(ct$estimate) > 0, "positively", "negatively"),
          fmt(as.numeric(ct$estimate)), format.pval(ct$p.value, digits = 3),
          ifelse(d_tab$mean_High > d_tab$mean_Low, "higher", "lower"),
          format.pval(d_tab$p_High_vs_Low, digits = 3),
          fmt(c_acsl4), fmt(c_both - c_risk), format.pval(lrt, digits = 3),
          fmt(c_both - c_acsl4), format.pval(lrt2, digits = 3)),
  "",
  "## 输出文件",
  "- acsl4_tumor_vs_normal.csv / acsl4_vs_risk.csv",
  "- acsl4_vs_ferroptosis_genes.csv",
  "- acsl4_vs_module_scores.csv：历史手工分组结果，不纳入论文论证",
  "- acsl4_model_comparison.csv / acsl4_joint_model_cox.csv / acsl4_tertile_risk_effect.csv",
  "- acsl4_vs_risk_scatter.png / acsl4_tumor_vs_normal.png / km_acsl4.png / acsl4_gene_correlations.png"
)
writeLines(md, file.path(OUT, "acsl4_methods_and_results.md"))

cat("\n== step15 done ==\n")
cat("outputs ->", OUT, "\n")
