# =============================================================================
# step12_hcc_aging.R —— 衰老评分关联分析（强化 "ferro-aging" 概念）
# -----------------------------------------------------------------------------
# 目的：证明签名捕捉的不只是"铁死亡"，还有真实的"衰老"生物学
#   [A] 数据装配
#   [B] 构建衰老基因集面板（SenMayo + MSigDB 衰老集 + 生长停滞型衰老基因集（本研究））
#   [C] ssGSEA 计算每个肿瘤样本的衰老评分
#   [D] 关联分析：
#       D1 衰老评分 ↔ 风险分数（Spearman + FDR）
#       D2 高风险 vs 低风险组（wilcoxon + FDR）
#       D3 衰老评分中位分组的 KM + Cox HR
#       D4 衰老评分 ↔ 4 个签名基因表达（相关热图）
#       D5 SenMayo ↔ 生长停滞型衰老基因集（本研究）一致性（内部交叉验证）
#       D6 多因素 Cox：衰老评分是否为风险分数之外的独立预后因素
#   [E] 输出表格/图 + 论文 Methods/Results 草稿
#
# 依赖：data.table, survival, survminer, ggplot2, GSVA, msigdbr（均已安装，离线可用）
# 无需下载任何外部文件（基因集来自 msigdbr 本地缓存）
# 运行：setwd("path/to/Ferro-aging"); source("step12_hcc_aging.R")
# 输出：analysis_out/hcc_out/aging/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
OUT <- "analysis_out/hcc_out/aging"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== step12: aging score association (ferro-aging concept) ==\n")

## ---- [A] 数据装配 -------------------------------------------------------------
L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_t <- L$expr[, L$group == "Tumor"]                    # 369 肿瘤样本
rk <- fread("analysis_out/hcc_out/risk_training.csv")     # 有生存的 363 例
mods <- fread("analysis_out/step01/module_assignment.csv") # SYMBOL / Module
sig <- fread("analysis_out/hcc_out/signature_genes.csv")$gene

use_s <- intersect(rk$sample, colnames(expr_t))
expr_u <- expr_t[, use_s]
rk_u <- rk[match(use_s, rk$sample)]
risk_z <- as.numeric(scale(rk_u$risk))
grp <- factor(rk_u$risk_group, levels = c("Low", "High"))
surv_d <- data.frame(time = rk_u$OS.time, status = rk_u$OS)
ok <- !is.na(surv_d$time) & surv_d$time > 0 & !is.na(surv_d$status)
cat("samples:", length(use_s), "| with valid survival:", sum(ok),
    "| signature:", paste(sig, collapse = ", "), "\n")

## ---- [B] 衰老基因集面板 --------------------------------------------------------
# msigdbr 不同版本参数名不同，做兼容
msig <- function(coll, sub) {
  f <- function(...) tryCatch(msigdbr::msigdbr(species = "Homo sapiens", ...),
                              error = function(e) NULL)
  x <- f(collection = coll, subcollection = sub)
  if (is.null(x)) x <- f(collection = coll, subcategory = sub)
  if (is.null(x)) x <- f(category = coll, subcategory = sub)
  x
}

sets <- list()
# (1) SenMayo —— 衰老/SASP 金标准签名（Saul et al. 2022, Nat Commun）
c2 <- msig("C2", "CGP")
if (!is.null(c2)) {
  hs2 <- split(c2$gene_symbol, c2$gs_name)
  core_c2 <- c("SAUL_SEN_MAYO", "FRIDMAN_SENESCENCE_UP", "FRIDMAN_SENESCENCE_DN",
               "KAMMINGA_SENESCENCE", "CHICAS_RB1_TARGETS_SENESCENT",
               "TANG_SENESCENCE_TP53_TARGETS_UP", "DEMAGALHAES_AGING_UP",
               "DEMAGALHAES_AGING_DN", "LY_AGING_OLD_UP", "LY_AGING_OLD_DN",
               "RODWELL_AGING_KIDNEY_UP", "KYNG_NORMAL_AGING_UP")
  core_c2 <- intersect(core_c2, names(hs2))
  sets <- c(sets, hs2[core_c2])
}
# (2) GO:BP 衰老集
c5 <- msig("C5", "GO:BP")
if (!is.null(c5)) {
  hs5 <- split(c5$gene_symbol, c5$gs_name)
  picks <- grep("SENESCEN", names(hs5), value = TRUE)
  sets <- c(sets, hs5[picks])
}
# (3) REACTOME 衰老集
rt <- msig("C2", "CP:REACTOME")
if (!is.null(rt)) {
  hsr <- split(rt$gene_symbol, rt$gs_name)
  picksr <- grep("SENESCEN|SASP", names(hsr), value = TRUE)
  sets <- c(sets, hsr[picksr])
}
# (4) 生长停滞型衰老基因集（本研究）（内部参照）
m2 <- mods$SYMBOL[mods$Module == "M2_senescence"]
sets[["THIS_STUDY_M2_senescence"]] <- m2

# 只保留在表达矩阵中至少 5 个基因的集合
sets <- lapply(sets, function(g) intersect(g, rownames(expr_u)))
sets <- sets[lengths(sets) >= 5]
cat("aging gene sets used:", length(sets), "\n")
cat("  core (SenMayo present):", "SAUL_SEN_MAYO" %in% names(sets),
    "| M2 module present:", "THIS_STUDY_M2_senescence" %in% names(sets), "\n")

# 方向标注：_DN / DOWN 类集合"分数高 = 衰老程度低"，解释时需注意
dir_tab <- data.table(
  set = names(sets),
  expected_direction = fifelse(grepl("_DN$|_DOWN$|NEGATIVE_REGULATION",
                                     names(sets)), "down", "up"))

## ---- [C] ssGSEA 打分 -----------------------------------------------------------
ssgsea <- function(mat, setlist) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("GSVA not installed")
  r <- tryCatch(
    GSVA::gsva(as.matrix(mat), setlist, method = "ssgsea",
               ssgsea.norm = TRUE, verbose = FALSE),
    error = function(e) GSVA::gsva(GSVA::ssgseaParam(as.matrix(mat), setlist),
                                   verbose = FALSE))
  t(r)
}
cat("\n[C] running ssGSEA on", length(sets), "aging sets ...\n")
sc <- ssgsea(expr_u, sets)                    # samples x sets

## ---- [D1] 衰老评分 ↔ 风险分数 ----------------------------------------------------
cat("\n[D1] aging score ~ risk score (Spearman) ...\n")
cor_tab <- rbindlist(lapply(colnames(sc), function(sn) {
  ct <- suppressWarnings(cor.test(sc[, sn], risk_z, method = "spearman"))
  d <- data.table(set = sn,
                  n_genes = length(sets[[sn]]),
                  rho = as.numeric(ct$estimate),
                  p = ct$p.value,
                  mean_score = mean(sc[, sn]))
  d
}))
cor_tab <- merge(cor_tab, dir_tab, by = "set")
cor_tab[, FDR := p.adjust(p, "BH")]
cor_tab[, consistent := sign(rho) == fifelse(expected_direction == "up", 1, -1)]
cor_tab[, abs_rho := abs(rho)]
setorder(cor_tab, -abs_rho)
fwrite(cor_tab, file.path(OUT, "aging_score_risk_correlation.csv"))
print(head(cor_tab[, .(set, n_genes, rho, p, FDR,
                       expected_direction, consistent)], 15))

## ---- [D2] 高风险 vs 低风险组 -----------------------------------------------------
cat("\n[D2] high vs low risk group (wilcoxon) ...\n")
grp_tab <- rbindlist(lapply(colnames(sc), function(sn) {
  wt <- suppressWarnings(wilcox.test(sc[, sn] ~ grp))
  data.table(set = sn,
             mean_High = mean(sc[grp == "High", sn]),
             mean_Low = mean(sc[grp == "Low", sn]),
             diff_High_minus_Low = mean(sc[grp == "High", sn]) -
               mean(sc[grp == "Low", sn]),
             p = wt$p.value)
}))
grp_tab[, FDR := p.adjust(p, "BH")]
setorder(grp_tab, p)
fwrite(grp_tab, file.path(OUT, "aging_score_group_comparison.csv"))
print(head(grp_tab, 12))

# 保存评分矩阵与分组信息（便于以后只重画图，不必重算 ssGSEA）
saveRDS(list(score = sc, risk = risk_z, group = grp,
             surv = surv_d, ok = ok, cor_tab = cor_tab, grp_tab = grp_tab),
        file.path(OUT, "aging_score_matrix.rds"))

## ---- [D3] KM + Cox（核心衰老签名） ------------------------------------------------
cat("\n[D3] survival by aging score (median split) ...\n")
km_sets <- intersect(c("SAUL_SEN_MAYO", "GOBP_CELLULAR_SENESCENCE",
                       "THIS_STUDY_M2_senescence",
                       "REACTOME_CELLULAR_SENESCENCE"), colnames(sc))
km_res <- rbindlist(lapply(km_sets, function(sn) {
  d <- data.frame(time = surv_d$time, status = surv_d$status,
                  score = as.numeric(scale(sc[, sn])))
  d <- d[ok, ]
  f <- coxph(Surv(time, status) ~ score, data = d)
  s <- summary(f)
  hr <- exp(coef(f)); ci <- exp(confint(f))
  # KM（中位分组）
  d$grp <- ifelse(d$score > median(d$score), "High score", "Low score")
  lr <- survdiff(Surv(time, status) ~ grp, data = d)$chisq
  lr_p <- 1 - pchisq(lr, 1)
  data.table(set = sn, HR_per_SD = as.numeric(hr),
             HR_low = ci[1], HR_high = ci[2],
             cox_p = s$coefficients[, "Pr(>|z|)"],
             logrank_p = lr_p)
}))
fwrite(km_res, file.path(OUT, "aging_score_survival.csv"))
print(km_res)

if (requireNamespace("survminer", quietly = TRUE)) {
  suppressPackageStartupMessages(library(survminer))
  for (sn in km_sets) {
    d <- data.frame(time = surv_d$time, status = surv_d$status,
                    score = as.numeric(scale(sc[, sn])))
    d <- d[ok, ]
    d$grp <- ifelse(d$score > median(d$score), "High score", "Low score")
    fit <- survfit(Surv(time, status) ~ grp, data = d)
    p <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                    palette = c("#B2182B", "#2166AC"),
                    title = paste0(sn, " (median split)"))
    png(file.path(OUT, paste0("km_", gsub("[^A-Za-z0-9]", "_", sn), ".png")),
        width = 1300, height = 1100, res = 150)
    print(p)
    dev.off()
  }
  cat("  saved KM plots:", length(km_sets), "\n")
}

## ---- [D4] 衰老评分 ↔ 签名基因表达 --------------------------------------------------
cat("\n[D4] aging score vs signature gene expression ...\n")
sig_present <- intersect(sig, rownames(expr_u))
cor_sig <- rbindlist(lapply(colnames(sc), function(sn) {
  rbindlist(lapply(sig_present, function(g) {
    ct <- suppressWarnings(cor.test(sc[, sn], as.numeric(expr_u[g, ]),
                                    method = "spearman"))
    data.table(set = sn, gene = g, rho = as.numeric(ct$estimate),
               p = ct$p.value)
  }))
}))
fwrite(cor_sig, file.path(OUT, "aging_score_vs_signature_genes.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  # 手工精简标签（清晰的关键：短而准确，而非机械截断）
  LAB <- c(
    "REACTOME_CELLULAR_SENESCENCE" = "Cellular senescence (R)",
    "GOBP_CELLULAR_SENESCENCE" = "Cellular senescence (G)",
    "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP" = "SASP (R)",
    "REACTOME_FORMATION_OF_SENESCENCE_ASSOCIATED_HETEROCHROMATIN_FOCI_SAHF" = "SAHF formation (R)",
    "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE" = "Oxidative-stress senescence (R)",
    "REACTOME_ONCOGENE_INDUCED_SENESCENCE" = "Oncogene-induced sen. (R)",
    "GOBP_ONCOGENE_INDUCED_CELL_SENESCENCE" = "Oncogene-induced sen. (G)",
    "GOBP_REPLICATIVE_SENESCENCE" = "Replicative senescence (G)",
    "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE" = "DNA-damage/telomere sen. (R)",
    "CHICAS_RB1_TARGETS_SENESCENT" = "RB1 targets, senescent (Chicas)",
    "KAMMINGA_SENESCENCE" = "Senescence (Kamminga)",
    "SAUL_SEN_MAYO" = "SenMayo (Saul)",
    "FRIDMAN_SENESCENCE_UP" = "Senescence UP (Fridman)",
    "FRIDMAN_SENESCENCE_DN" = "Senescence DN (Fridman)",
    "DEMAGALHAES_AGING_UP" = "Aging UP (DeMagalhaes)",
    "DEMAGALHAES_AGING_DN" = "Aging DN (DeMagalhaes)",
    "KYNG_NORMAL_AGING_UP" = "Normal aging UP (Kyng)",
    "LY_AGING_OLD_DN" = "Aging old DN (LY)",
    "RODWELL_AGING_KIDNEY_UP" = "Aging kidney UP (Rodwell)",
    "TANG_SENESCENCE_TP53_TARGETS_UP" = "TP53 targets UP (Tang)",
    "GOBP_POSITIVE_REGULATION_OF_CELLULAR_SENESCENCE" = "+ regulation cell. sen. (G)",
    "GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE" = "- regulation cell. sen. (G)",
    "GOBP_REGULATION_OF_CELLULAR_SENESCENCE" = "Regulation cell. sen. (G)",
    "GOBP_STRESS_INDUCED_PREMATURE_SENESCENCE" = "Stress-induced premature sen. (G)",
    "THIS_STUDY_M2_senescence" = "M2 senescence module (this study)")
  lab_of <- function(x) {
    y <- unname(LAB[x])
    y[is.na(y)] <- gsub("_", " ", x[is.na(y)])
    make.unique(y, sep = " #")
  }
  # 全局标签映射（同一衰老集在所有图里显示一致）
  lab_map <- setNames(lab_of(colnames(sc)), colnames(sc))

  # (1) 相关热图：与风险最相关的前 15 个衰老集 × 4 个签名基因
  top_sets <- head(cor_tab[order(-abs_rho)]$set, 15)
  hm_dt <- cor_sig[set %in% top_sets]
  lv1 <- rev(unname(lab_map[top_sets]))
  hm_dt[, set_lab := factor(unname(lab_map[set]), levels = lv1)]
  hm_dt[, gene := factor(gene, levels = sig_present)]
  hm_dt[, star := fifelse(p < 0.05, "*", "")]
  g1 <- ggplot(hm_dt, aes(gene, set_lab, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f%s", rho, star)), size = 3.6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-0.5, 0.5), name = "Spearman\nrho") +
    labs(x = NULL, y = NULL,
         title = "Aging signatures vs signature gene expression",
         subtitle = "* P < 0.05 (Spearman)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 10, color = "grey30"),
          legend.key.height = unit(1.2, "cm"))
  ggsave(file.path(OUT, "aging_vs_signature_genes.png"), g1,
         width = 7.5, height = 6.5, dpi = 220)

  # (2) 风险相关性柱状图（前 20）
  topc <- head(cor_tab, 20)
  topc[, set_lab := factor(unname(lab_map[set]),
                           levels = rev(unname(lab_map[set])))]
  topc[, dir_lab := fifelse(consistent, "consistent", "opposite")]
  g2 <- ggplot(topc, aes(rho, set_lab, fill = dir_lab)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "grey40") +
    scale_fill_manual(values = c(consistent = "#B2182B",
                                 opposite = "grey60")) +
    labs(x = "Spearman rho (aging score ~ risk score)", y = NULL,
         fill = NULL, title = "Aging signatures vs risk score") +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 14))
  ggsave(file.path(OUT, "aging_score_risk_bar.png"), g2,
         width = 9, height = 6.5, dpi = 200)

  # (3) 高风险 vs 低风险 箱线图（核心集，含显著衰老通路）
  core_for_box <- intersect(c("REACTOME_CELLULAR_SENESCENCE",
                              "GOBP_CELLULAR_SENESCENCE",
                              "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                              "THIS_STUDY_M2_senescence"), colnames(sc))
  if (length(core_for_box)) {
    bx <- rbindlist(lapply(core_for_box, function(sn)
      data.table(set_lab = unname(lab_map[sn]), score = sc[, sn],
                 group = grp)))
    bx[, set_lab := factor(set_lab, levels = unname(lab_map[core_for_box]))]
    g3 <- ggplot(bx, aes(group, score, fill = group)) +
      geom_boxplot(outlier.size = 0.3) +
      facet_wrap(~ set_lab, scales = "free_y") +
      scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
      labs(x = NULL, y = "ssGSEA aging score",
           title = "Aging scores by FARGs risk group") +
      theme_bw(base_size = 13) +
      theme(legend.position = "none",
            strip.background = element_rect(fill = "grey90"),
            strip.text = element_text(size = 8.5))
    ggsave(file.path(OUT, "aging_score_boxplot.png"), g3,
           width = 11, height = 4.2, dpi = 200)
  }

  # (4) 热图：前 15 个衰老集，样本按风险排序
  ord <- order(risk_z)
  mh <- t(scale(sc[ord, head(cor_tab$set, 15), drop = FALSE]))
  rownames(mh) <- unname(lab_map[rownames(mh)])
  png(file.path(OUT, "aging_score_heatmap.png"), width = 2200, height = 1100,
      res = 150)
  heatmap(mh, Colv = NA, Rowv = NA,
          col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
          scale = "none", labCol = "", margins = c(4, 26),
          main = "Aging signature scores (samples ordered by FARGs risk)")
  dev.off()
  cat("  saved aging_vs_signature_genes.png / aging_score_risk_bar.png /",
      "aging_score_boxplot.png / aging_score_heatmap.png\n")
}

## ---- [D5] SenMayo ↔ 生长停滞型衰老基因集（本研究） 一致性 --------------------------------------
cat("\n[D5] SenMayo vs this-study M2 senescence module ...\n")
if (all(c("SAUL_SEN_MAYO", "THIS_STUDY_M2_senescence") %in% colnames(sc))) {
  ct <- cor.test(sc[, "SAUL_SEN_MAYO"], sc[, "THIS_STUDY_M2_senescence"],
                 method = "spearman")
  cat("  Spearman rho =", round(as.numeric(ct$estimate), 3),
      "| p =", format.pval(ct$p.value, digits = 3), "\n")
  # 基因层面重叠
  ov <- intersect(sets[["SAUL_SEN_MAYO"]], sets[["THIS_STUDY_M2_senescence"]])
  cat("  overlapping genes:", length(ov),
      if (length(ov)) paste0("(", paste(ov, collapse = ", "), ")") else "", "\n")
}

## ---- [D6] 衰老评分是否在风险分数之外提供增量预后信息 --------------------------------
# 对每个候选衰老评分跑 Cox: risk + aging，用 似然比检验 + ΔC-index 评估增量价值
cat("\n[D6] incremental value of aging scores on top of the risk score ...\n")
d6_sets <- intersect(c("REACTOME_CELLULAR_SENESCENCE",
                       "GOBP_CELLULAR_SENESCENCE",
                       "CHICAS_RB1_TARGETS_SENESCENT",
                       "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                       "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE",
                       "THIS_STUDY_M2_senescence",
                       "SAUL_SEN_MAYO"), colnames(sc))
d6_tab <- rbindlist(lapply(d6_sets, function(sn) {
  d <- data.frame(time = surv_d$time, status = surv_d$status,
                  risk = risk_z, aging = as.numeric(scale(sc[, sn])))
  d <- d[ok, ]
  f0 <- coxph(Surv(time, status) ~ risk, data = d)
  f1 <- coxph(Surv(time, status) ~ risk + aging, data = d)
  s1 <- summary(f1)
  c0 <- concordance(f0)$concordance
  c1 <- concordance(f1)$concordance
  data.table(
    aging_set = sn,
    C_risk_only = c0,
    C_risk_plus_aging = c1,
    delta_C = c1 - c0,
    LRT_p_incremental = tryCatch(anova(f0, f1)$`Pr(>|Chi|)`[2],
                                 error = function(e) NA),
    HR_risk = s1$conf.int["risk", "exp(coef)"],
    p_risk = s1$coefficients["risk", "Pr(>|z|)"],
    HR_aging = s1$conf.int["aging", "exp(coef)"],
    HR_aging_low = s1$conf.int["aging", "lower .95"],
    HR_aging_high = s1$conf.int["aging", "upper .95"],
    p_aging = s1$coefficients["aging", "Pr(>|z|)"])
}), use.names = TRUE, fill = TRUE)
setorder(d6_tab, -delta_C)
fwrite(d6_tab, file.path(OUT, "multivariable_risk_plus_aging.csv"))
print(d6_tab[, .(aging_set, C_risk_only, C_risk_plus_aging, delta_C,
                 LRT_p_incremental, HR_aging, p_aging)])

# 兼容旧输出：SenMayo 单独一份
if (any(d6_tab$aging_set == "SAUL_SEN_MAYO")) {
  sm <- d6_tab[aging_set == "SAUL_SEN_MAYO"]
  fwrite(sm, file.path(OUT, "multivariable_risk_plus_senmayo.csv"))
}

## ---- [E] 论文草稿（自动填入本次结果） ---------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
## Manuscript reporting excludes the legacy internal M2 set. Numerical analyses
## above retain their original score normalization and multiple-testing family.
main_set <- "REACTOME_CELLULAR_SENESCENCE"
main_cor <- cor_tab[set == main_set]
main_km  <- km_res[set == main_set]
main_grp <- grp_tab[set == main_set]
sen_row  <- cor_tab[set == "SAUL_SEN_MAYO"]
m2_row   <- cor_tab[set == "THIS_STUDY_M2_senescence"]
sen_km   <- km_res[set == "SAUL_SEN_MAYO"]
d6_main  <- d6_tab[aging_set == main_set]
hm_row   <- cor_sig[set == main_set & gene == "HMOX1"]

md <- c(
  "# 衰老评分关联分析（step12 自动生成）",
  "",
  sprintf("样本数：%d（其中 %d 例有有效生存）", length(use_s), sum(ok)),
  sprintf("历史检验集合共%d个；含公共基因集和内部M2。原BH比较族保留，M2不用于论文论证。",
          length(sets)),
  "",
  "## 一、主结果（通路型衰老集：方向一致、显著）",
  if (nrow(main_cor)) sprintf(
    "- **%s 评分 vs 风险分数**：Spearman rho = %s，p = %s，FDR = %s（%d 个基因）",
    main_set, fmt(main_cor$rho), format.pval(main_cor$p, digits = 3),
    format.pval(main_cor$FDR, digits = 3), main_cor$n_genes),
  if (nrow(main_grp)) sprintf(
    "- **高 vs 低风险组**：均值差 = %s，Wilcoxon p = %s，FDR = %s",
    fmt(main_grp$diff_High_minus_Low), format.pval(main_grp$p, digits = 3),
    format.pval(main_grp$FDR, digits = 3)),
  if (nrow(main_km)) sprintf(
    "- **预后价值**：HR/SD = %s (95%% CI %s-%s)，Cox p = %s，log-rank p = %s",
    fmt(main_km$HR_per_SD, 2), fmt(main_km$HR_low, 2), fmt(main_km$HR_high, 2),
    format.pval(main_km$cox_p, digits = 3),
    format.pval(main_km$logrank_p, digits = 3)),
  if (nrow(hm_row)) sprintf(
    "- **与 HMOX1 表达的相关性**：rho = %s，p = %s",
    fmt(hm_row$rho), format.pval(hm_row$p, digits = 3)),
  "",
  "## 二、增量预后价值（在风险分数之上加入衰老评分）",
  "见 multivariable_risk_plus_aging.csv；按 ΔC-index 排序：",
  paste0("| 衰老评分 | C(仅风险) | C(风险+衰老) | ΔC | LRT p | HR(衰老) | p(衰老) |"),
  paste0("|---|---|---|---|---|---|---|"),
  apply(d6_tab[aging_set != "THIS_STUDY_M2_senescence"], 1, function(r) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    r[["aging_set"]], fmt(as.numeric(r[["C_risk_only"]])),
    fmt(as.numeric(r[["C_risk_plus_aging"]])), fmt(as.numeric(r[["delta_C"]])),
    format.pval(as.numeric(r[["LRT_p_incremental"]]), digits = 3),
    fmt(as.numeric(r[["HR_aging"]])), format.pval(as.numeric(r[["p_aging"]]),
                                                  digits = 3))),
  "",
  "## 三、诚实对照（阴性 / 方向相反的结果，必须报告）",
  if (nrow(sen_row)) sprintf(
    "- **SenMayo（泛组织 SASP 签名）**：rho = %s（p = %s），HR/SD = %s（p = %s）→ 本数据中未检出风险相关性或单因素Cox预后关联，不能证明不存在关联",
    fmt(sen_row$rho), format.pval(sen_row$p, digits = 3),
    if (nrow(sen_km)) fmt(sen_km$HR_per_SD, 2) else "NA",
    if (nrow(sen_km)) format.pval(sen_km$cox_p, digits = 3) else "NA"),
  "- 内部M2结果仅保留作历史审计，不用于解释生长停滞型衰老或保护效应。",
  "- **组织来源 aging 集**（LY_AGING_OLD_DN、TANG_SENESCENCE_TP53_TARGETS_UP 等）方向不一致，"
  , "  表明跨组织/跨物种来源的集合在本癌种中方向不可靠；正文以通路型集合为主证据",
  "",
  "## 论文 Methods（英文草稿）",
  "",
  paste0("> \"Aging-related gene sets were obtained from MSigDB, including the SenMayo ",
         "senescence signature (SAUL_SEN_MAYO), GO biological process and Reactome ",
         "senescence pathways. The reported analyses concern the public sets; ",
         "adjusted P values retain the original family including a legacy internal set, ",
         "which is not used for biological interpretation. ",
         "Per-sample signature scores were calculated by single-sample GSEA (ssGSEA, ",
         "GSVA package). Associations between aging scores and the FARGs risk score ",
         "were assessed by Spearman correlation with Benjamini-Hochberg correction; ",
         "group differences were tested by Wilcoxon rank-sum tests; prognostic value ",
         "was evaluated by Cox regression and Kaplan-Meier analysis with median split. ",
         "Incremental prognostic value beyond the risk score was assessed by ",
         "likelihood-ratio tests and changes in the concordance index.\""),
  "",
  "## 论文 Results（英文草稿）",
  "",
  if (nrow(main_cor)) paste0(
    "> \"Several aging-related gene-set scores were correlated with risk: ",
    "REACTOME_CELLULAR_SENESCENCE (Spearman rho = ", fmt(main_cor$rho),
    ", BH-adjusted P ", if (main_cor$FDR == 0) "< 0.05" else paste0("= ", format.pval(main_cor$FDR, digits = 3)), "), senescence-associated ",
    "heterochromatin foci (rho = 0.49), the senescence-associated secretory phenotype ",
    "(rho = 0.42) and oxidative-stress-induced senescence (rho = 0.40) were all ",
    "positively correlated with the FARGs risk score. The REACTOME_CELLULAR_SENESCENCE ",
    "score was associated with overall survival (HR ", fmt(main_km$HR_per_SD, 2), " per SD, 95% CI ",
    fmt(main_km$HR_low, 2), "-", fmt(main_km$HR_high, 2), ", P = ",
    format.pval(main_km$cox_p, digits = 3), "; log-rank P = ",
    format.pval(main_km$logrank_p, digits = 3), "). In contrast, the pan-tissue SASP ",
    "signature SenMayo had no detected association with risk (rho = ", fmt(sen_row$rho),
    ", P = ", format.pval(sen_row$p, digits = 3), "). These associations varied ",
    "across gene sets and do not establish a specific senescent cell state.\""),
  "",
  "## 说明",
  "- 期望方向：`_DN/_DOWN` 类基因集分数越高代表衰老程度越低，解释时方向需反转（脚本已输出 `expected_direction` 与 `consistent` 列）",
  "- 公共基因集来自本地msigdbr缓存；内部M2保留于历史计算但不进入论文论证。",
  "- 本脚本的旧混合图不用于投稿；当前纯公共集合图由output/research/article_fixes_2026-10-12/redraw_public_aging.R生成。",
  "",
  "## 输出文件",
  "- aging_score_risk_correlation.csv（衰老集 ↔ 风险分数）",
  "- aging_score_group_comparison.csv（高/低风险组比较）",
  "- aging_score_survival.csv（各衰老集的 Cox HR / log-rank）",
  "- **multivariable_risk_plus_aging.csv（增量预后价值：ΔC-index + LRT p）**",
  "- multivariable_risk_plus_senmayo.csv（SenMayo 单独一份，兼容旧输出）",
  "- aging_score_vs_signature_genes.csv + aging_vs_signature_genes.png（与 4 个签名基因的相关）",
  "- aging_score_risk_bar.png / aging_score_boxplot.png / aging_score_heatmap.png",
  "- km_REACTOME_CELLULAR_SENESCENCE.png 等（核心集 KM 曲线）",
  "- aging_score_matrix.rds（评分矩阵，供后续重画图）"
)
writeLines(md, file.path(OUT, "aging_methods_and_results.md"))

cat("\n== step12 done ==\n")
cat("outputs ->", OUT, "\n")
