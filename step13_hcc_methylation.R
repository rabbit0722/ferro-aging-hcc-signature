# =============================================================================
# step13_hcc_methylation.R —— 签名基因甲基化分析（本地复现 MethSurv + MEXPRESS + UALCAN）
# -----------------------------------------------------------------------------
# 数据（本地，无需联网）：
#   analysis_out/data/TCGA.LIHC.HumanMethylation450.gz   Xena 450K beta 矩阵（429 样本）
#   analysis_out/data/methylation_probeMap.tsv           探针 -> 基因 映射（hg19）
#   analysis_out/hcc_data/hcc_data.rds                   表达矩阵（含 Tumor/Normal 分组）
#   analysis_out/hcc_out/risk_training.csv               生存数据（363 例）
#
# 分析内容（对 HMOX1 / TFRC / SOCS2 / E2F1）：
#   [D] 肿瘤 vs 正常 甲基化差异        （UALCAN 等价）
#   [E] 甲基化 ↔ 表达 相关             （MEXPRESS 等价）
#   [F] 探针级 Cox + KM 生存分析       （MethSurv 等价）
#   [G] 三者方向一致性判定 + 论文草稿
#
# 输出：analysis_out/hcc_out/methylation/
# 运行：setwd("path/to/Ferro-aging"); source("step13_hcc_methylation.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
OUT <- "analysis_out/hcc_out/methylation"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== step13: signature gene methylation ==\n")

GENES <- c("HMOX1", "TFRC", "SOCS2", "E2F1")

## ---- [A] 载入注释与表达/生存数据 ------------------------------------------------
pmap <- fread("analysis_out/data/methylation_probeMap.tsv")
setnames(pmap, 1, "probe")
pmap <- pmap[, .(probe, gene_map = gene, chrom, start = chromStart, end = chromEnd)]

L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_all <- L$expr
grp_all <- L$group
names(grp_all) <- colnames(expr_all)
rk <- fread("analysis_out/hcc_out/risk_training.csv")

# 可选：Illumina 450K 官方区域注释（若已安装）
has_anno <- requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                             quietly = TRUE)
if (has_anno) {
  suppressPackageStartupMessages(
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  cat("450K annotation loaded:", nrow(ann), "probes\n")
}

# 每个基因的探针集合（probeMap 的 gene 列可能是逗号分隔多基因）
probe_list <- lapply(GENES, function(g)
  pmap[!is.na(gene_map) & grepl(paste0("(^|,)", g, "($|,)"), gene_map)]$probe)
names(probe_list) <- GENES
n_probe <- lengths(probe_list)
cat("probes per gene:", paste(GENES, n_probe, sep = "=", collapse = " | "), "\n")
target <- unique(unlist(probe_list))

## ---- [B] 流式读取甲基化矩阵，只保留目标探针 ---------------------------------------
meth_file <- "analysis_out/data/TCGA.LIHC.HumanMethylation450.gz"
cat("\n[B] streaming methylation matrix (", round(file.size(meth_file) / 1e6),
    "MB ) ...\n", sep = "")
con <- gzfile(meth_file, "r")
hdr <- strsplit(readLines(con, n = 1), "\t")[[1]]
samples_meth <- hdr[-1]
keep_lines <- character(0)
repeat {
  chunk <- readLines(con, n = 50000, warn = FALSE)
  if (length(chunk) == 0) break
  hit <- sub("\t.*$", "", chunk) %in% target
  if (any(hit)) keep_lines <- c(keep_lines, chunk[hit])
}
close(con)
cat("  target probes kept:", length(keep_lines), "/", length(target), "\n")

mat <- do.call(rbind, lapply(strsplit(keep_lines, "\t"),
                             function(x) as.numeric(x[-1])))
rownames(mat) <- sub("\t.*$", "", keep_lines)
colnames(mat) <- samples_meth
cat("  beta matrix:", nrow(mat), "probes x", ncol(mat), "samples\n")

## ---- [C] 探针区域/岛注释 + 样本分型 ----------------------------------------------
beta_dt <- as.data.table(mat, keep.rownames = "probe")

probe_meta <- rbindlist(lapply(GENES, function(g) {
  dt <- data.table(gene = g, probe = probe_list[[g]])
  dt
}))
probe_meta <- merge(probe_meta, pmap, by = "probe", all.x = TRUE)

if (has_anno) {
  # 从官方注释提取该基因对应的区域（TSS1500/TSS200/1stExon/5'UTR/Body/3'UTR）
  get_region <- function(probe, gene) {
    nm <- ann$UCSC_RefGene_Name[probe]; gp <- ann$UCSC_RefGene_Group[probe]
    if (is.na(nm) || is.na(gp)) return(NA_character_)
    nms <- strsplit(nm, ";")[[1]]; gps <- strsplit(gp, ";")[[1]]
    hit <- which(nms == gene)
    if (length(hit) == 0) return(NA_character_)
    paste(unique(gps[hit]), collapse = "/")
  }
  probe_meta[, region := mapply(get_region, probe, gene)]
  probe_meta[, island := ann[probe, "Relation_to_Island"]]
} else {
  probe_meta[, region := NA_character_]
  probe_meta[, island := NA_character_]
}

# TSS 相对位置（hg19 refGene；探针坐标同版本，可安全比较）
ref_file <- "analysis_out/data/hg19_refGene.txt.gz"
if (file.exists(ref_file)) {
  ref <- fread(ref_file, header = FALSE)
  setnames(ref, c("bin", "name", "chrom", "strand", "txStart", "txEnd",
                  "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                  "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  tss <- ref[name2 %in% GENES,
             .(TSS = if (strand[1] == "+") min(txStart) else max(txEnd),
               strand = strand[1], chrom_gene = chrom[1]), by = .(gene = name2)]
  probe_meta <- merge(probe_meta, tss[, .(gene, TSS, strand)], by = "gene",
                      all.x = TRUE)
  # 转录取向的带符号距离：正 = 基因内（TSS 下游），负 = 上游
  probe_meta[, dist_to_TSS := fifelse(strand == "+", start - TSS, TSS - start)]
  probe_meta[, region_class := fifelse(
    is.na(dist_to_TSS), NA_character_,
    fifelse(dist_to_TSS >= -1500 & dist_to_TSS <= 1500, "TSS-proximal",
            fifelse(dist_to_TSS > 1500, "gene body", "upstream")))]
  cat("\nTSS 定位（hg19 refGene）：\n"); print(tss)
  cat("探针区域分布：\n"); print(probe_meta[, .N, by = region_class])
} else {
  probe_meta[, `:=`(TSS = NA_real_, strand = NA_character_,
                    dist_to_TSS = NA_real_, region_class = NA_character_)]
}

# 样本类型：TCGA barcode 第 14-15 位（01=肿瘤，11=癌旁正常）
sample_type <- substr(samples_meth, 14, 15)
tum_s <- samples_meth[sample_type == "01"]
nor_s <- samples_meth[sample_type == "11"]
cat("  tumour samples:", length(tum_s), "| normal samples:", length(nor_s), "\n")

## ---- [D] 肿瘤 vs 正常 甲基化差异（UALCAN 等价）-----------------------------------
cat("\n[D] tumour vs normal methylation ...\n")
dn <- rbindlist(lapply(seq_len(nrow(mat)), function(i) {
  tv <- as.numeric(mat[i, tum_s]); nv <- as.numeric(mat[i, nor_s])
  tv <- tv[!is.na(tv)]; nv <- nv[!is.na(nv)]
  wt <- tryCatch(suppressWarnings(wilcox.test(tv, nv)),
                 error = function(e) NULL)
  data.table(probe = rownames(mat)[i],
             n_tumor = length(tv), n_normal = length(nv),
             beta_tumor = mean(tv), beta_normal = mean(nv),
             delta_beta = mean(tv) - mean(nv),
             p_tumor_vs_normal = if (is.null(wt)) NA_real_ else wt$p.value)
}))
dn[, FDR_tumor_vs_normal := p.adjust(p_tumor_vs_normal, "BH")]

## ---- [E] 甲基化 ↔ 表达 相关（MEXPRESS 等价）--------------------------------------
cat("[E] methylation vs expression ...\n")
de <- rbindlist(lapply(seq_len(nrow(mat)), function(i) {
  pb <- rownames(mat)[i]
  g <- probe_meta$gene[probe_meta$probe == pb][1]
  if (is.na(g) || !g %in% rownames(expr_all)) return(NULL)
  common <- intersect(samples_meth, colnames(expr_all))
  common <- common[substr(common, 14, 15) == "01"]
  x <- as.numeric(mat[pb, common]); y <- as.numeric(expr_all[g, common])
  okk <- !is.na(x) & !is.na(y)
  if (sum(okk) < 30) return(NULL)
  ct <- suppressWarnings(cor.test(x[okk], y[okk], method = "spearman"))
  data.table(probe = pb, n_pairs = sum(okk),
             rho_meth_expr = as.numeric(ct$estimate),
             p_meth_expr = ct$p.value)
}))

## ---- [F] 探针级生存分析（MethSurv 等价）------------------------------------------
cat("[F] probe-level survival (Cox + KM) ...\n")
ds <- rbindlist(lapply(seq_len(nrow(mat)), function(i) {
  pb <- rownames(mat)[i]
  smp <- intersect(samples_meth, rk$sample)
  if (length(smp) < 50) return(NULL)
  b <- as.numeric(mat[pb, smp])
  idx <- match(smp, rk$sample)
  d <- data.table(beta = b, time = rk$OS.time[idx], status = rk$OS[idx])
  d <- d[!is.na(beta) & !is.na(time) & time > 0]
  if (nrow(d) < 50 || sum(d$status) < 10) return(NULL)
  # 连续 beta 的 Cox（按 SD 缩放，HR 解释为每 SD 甲基化）
  d[, beta_z := as.numeric(scale(beta))]
  f <- coxph(Surv(time, status) ~ beta_z, data = d)
  s <- summary(f)
  # 中位分组的 KM log-rank
  d[, grp := ifelse(beta > median(beta), "High", "Low")]
  lr <- survdiff(Surv(time, status) ~ grp, data = d)
  lr_p <- 1 - pchisq(lr$chisq, 1)
  data.table(probe = pb, n = nrow(d), events = sum(d$status),
             HR_per_SD = s$conf.int["beta_z", "exp(coef)"],
             HR_low = s$conf.int["beta_z", "lower .95"],
             HR_high = s$conf.int["beta_z", "upper .95"],
             p_cox = s$coefficients["beta_z", "Pr(>|z|)"],
             logrank_p = lr_p,
             beta_median_high = mean(d$beta[d$grp == "High"]),
             beta_median_low = mean(d$beta[d$grp == "Low"]))
}))

## ---- [G] 汇总合并 + 方向一致性 ----------------------------------------------------
parts <- Filter(Negate(is.null), list(probe_meta, dn, de, ds))
res <- Reduce(function(a, b) merge(a, b, by = "probe", all = TRUE), parts)
res[, FDR_cox := p.adjust(p_cox, "BH"), by = gene]
# 预期方向：甲基化沉默表达（rho<0）→ 风险基因沉默=保护(HR<1)，保护基因沉默=风险(HR>1)
res[, expr_role := fifelse(gene == "SOCS2", "protective", "risk")]
res[, mech := fifelse(!is.na(rho_meth_expr) & rho_meth_expr < 0 &
                        p_meth_expr < 0.05, "silencing", "no/positive")]
res[, expected_HR_side := fifelse(expr_role == "risk", "HR<1", "HR>1")]
res[, survival_consistent := fifelse(
  expr_role == "risk", HR_per_SD < 1, HR_per_SD > 1)]
setorder(res, gene, p_cox)
fwrite(res, file.path(OUT, "methylation_probe_results.csv"))

# 每基因最佳探针（按生存 p 值）
best <- res[!is.na(p_cox), .SD[which.min(p_cox)], by = gene]
fwrite(best, file.path(OUT, "methylation_best_probe_per_gene.csv"))

# 机制自洽探针：TSS 近端 + 甲基化与表达显著负相关（启动子沉默证据）
coh <- res[region_class == "TSS-proximal" & !is.na(p_meth_expr) &
             p_meth_expr < 0.05 & rho_meth_expr < 0]
if (nrow(coh)) {
  setorder(coh, gene, p_meth_expr)
  coh_best <- coh[, .SD[1], by = gene]
} else {
  coh_best <- coh
}
fwrite(coh_best, file.path(OUT, "methylation_promoter_silencing_probes.csv"))

cat("\n=== 每基因最优探针（按 Cox p）===\n")
print(best[, .(gene, probe, region_class, dist_to_TSS, n, events, HR_per_SD,
               HR_low, HR_high, p_cox, FDR_cox, rho_meth_expr, p_meth_expr,
               delta_beta, p_tumor_vs_normal, survival_consistent)])
cat("\n=== 启动子沉默证据探针（TSS 近端 + 甲基化与表达负相关）===\n")
if (nrow(coh_best)) {
  print(coh_best[, .(gene, probe, dist_to_TSS, rho_meth_expr, p_meth_expr,
                     delta_beta, p_tumor_vs_normal, HR_per_SD, p_cox,
                     survival_consistent)])
} else {
  cat("  （无：未见 TSS 近端探针同时满足甲基化-表达显著负相关）\n")
}

## ---- [H] 每个基因的 3 联图（肿瘤vs正常 / 甲基化vs表达 / KM）----------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  has_ge <- requireNamespace("gridExtra", quietly = TRUE)
  has_sv <- requireNamespace("survminer", quietly = TRUE)
  for (g in GENES) {
    # 优先用"启动子沉默证据探针"（TSS 近端 + 甲基化-表达负相关），否则用生存最优探针
    bp <- coh_best[gene == g]
    if (nrow(bp) == 0) bp <- best[gene == g]
    if (nrow(bp) == 0) next
    pb <- bp$probe
    # (1) 肿瘤 vs 正常
    d1 <- data.table(beta = c(mat[pb, tum_s], mat[pb, nor_s]),
                     type = factor(c(rep("Tumour", length(tum_s)),
                                     rep("Normal", length(nor_s))),
                                   levels = c("Normal", "Tumour")))
    d1 <- d1[!is.na(beta)]
    g1 <- ggplot(d1, aes(type, beta, fill = type)) +
      geom_boxplot(outlier.size = 0.3, width = 0.6) +
      scale_fill_manual(values = c(Tumour = "#B2182B", Normal = "#2166AC")) +
      labs(x = NULL, y = paste0("Beta (", pb, ")"),
           title = paste0(g, "  [", pb, "]"),
           subtitle = sprintf("TSS %s bp | Tumour vs normal P = %s",
                              ifelse(is.na(bp$dist_to_TSS), "NA",
                                     as.character(round(bp$dist_to_TSS))),
                              format.pval(bp$p_tumor_vs_normal, digits = 3))) +
      theme_bw(base_size = 11) + theme(legend.position = "none")
    # (2) 甲基化 vs 表达
    common <- intersect(samples_meth, colnames(expr_all))
    common <- common[substr(common, 14, 15) == "01"]
    d2 <- data.table(beta = as.numeric(mat[pb, common]),
                     expr = as.numeric(expr_all[g, common]))
    d2 <- d2[!is.na(beta) & !is.na(expr)]
    g2 <- ggplot(d2, aes(beta, expr)) +
      geom_point(size = 0.7, alpha = 0.5, color = "grey30") +
      geom_smooth(method = "lm", formula = y ~ x, color = "#B2182B",
                  fill = "#B2182B", alpha = 0.15) +
      labs(x = paste0("Beta (", pb, ")"), y = paste0(g, " (log2)"),
           subtitle = sprintf("Methylation vs expression rho = %.2f, P = %s",
                              bp$rho_meth_expr,
                              format.pval(bp$p_meth_expr, digits = 3))) +
      theme_bw(base_size = 11)
    # (3) KM（中位分组）
    smp <- intersect(samples_meth, rk$sample)
    idx <- match(smp, rk$sample)
    d3 <- data.table(beta = as.numeric(mat[pb, smp]), time = rk$OS.time[idx],
                     status = rk$OS[idx])
    d3 <- d3[!is.na(beta) & !is.na(time) & time > 0]
    d3[, grp := ifelse(beta > median(beta), "High methylation",
                       "Low methylation")]
    g3 <- NULL
    if (has_sv) {
      fit <- survfit(Surv(time, status) ~ grp, data = d3)
      kp <- survminer::ggsurvplot(fit, data = d3, pval = TRUE,
                                  risk.table = FALSE,
                                  palette = c("#B2182B", "#2166AC"),
                                  legend.title = "", legend.labs = c("High", "Low"),
                                  ggtheme = theme_bw(base_size = 11))
      g3 <- kp$plot +
        labs(subtitle = sprintf("Cox HR/SD = %.2f (P = %s)", bp$HR_per_SD,
                                format.pval(bp$p_cox, digits = 3)))
    }
    if (has_ge && !is.null(g3)) {
      p <- file.path(OUT, paste0("methylation_panel_", g, ".png"))
      png(p, width = 2200, height = 750, res = 150)
      gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
      dev.off()
      cat("  saved methylation_panel_", g, ".png\n", sep = "")
    } else {
      ggsave(file.path(OUT, paste0("methylation_", g, "_tumor_vs_normal.png")),
             g1, width = 4, height = 4, dpi = 200)
      ggsave(file.path(OUT, paste0("methylation_", g, "_vs_expr.png")),
             g2, width = 4.5, height = 4, dpi = 200)
      if (!is.null(g3))
        ggsave(file.path(OUT, paste0("methylation_", g, "_KM.png")), g3,
               width = 5, height = 4, dpi = 200)
    }
  }
}

## ---- [I] 论文草稿 ---------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
md <- c(
  "# 签名基因甲基化分析（step13 自动生成）",
  "",
  sprintf("数据：TCGA-LIHC 450K 甲基化（%d 探针 × %d 样本）", nrow(mat),
          ncol(mat)),
  sprintf("样本：肿瘤 %d 例，癌旁正常 %d 例", length(tum_s), length(nor_s)),
  "",
  "## 每基因最优探针（按 Cox p）",
  paste0("| 基因 | 探针 | 区域 | TSS 距离 | HR/SD | p(Cox) | FDR | rho(甲-表) | p | Δβ(肿瘤-正常) | 方向一致 |"),
  paste0("|---|---|---|---|---|---|---|---|---|---|---|"),
  apply(best, 1, function(r) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r[["gene"]], r[["probe"]],
    ifelse(is.na(r[["region_class"]]), "-", r[["region_class"]]),
    ifelse(is.na(r[["dist_to_TSS"]]), "-", as.character(round(as.numeric(r[["dist_to_TSS"]])))),
    fmt(as.numeric(r[["HR_per_SD"]])), format.pval(as.numeric(r[["p_cox"]]), digits = 3),
    format.pval(as.numeric(r[["FDR_cox"]]), digits = 3),
    fmt(as.numeric(r[["rho_meth_expr"]])),
    format.pval(as.numeric(r[["p_meth_expr"]]), digits = 3),
    fmt(as.numeric(r[["delta_beta"]])),
    ifelse(as.logical(r[["survival_consistent"]]), "yes", "no"))),
  "",
  "## 启动子近端关联探针（甲基化–表达负相关；非沉默机制验证）",
  if (nrow(coh_best)) paste0("| 基因 | 探针 | TSS 距离 | rho(甲-表) | p | Δβ | HR/SD | p(Cox) |"),
  if (nrow(coh_best)) paste0("|---|---|---|---|---|---|---|---|"),
  if (nrow(coh_best)) apply(coh_best, 1, function(r) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    r[["gene"]], r[["probe"]], round(as.numeric(r[["dist_to_TSS"]])),
    fmt(as.numeric(r[["rho_meth_expr"]])),
    format.pval(as.numeric(r[["p_meth_expr"]]), digits = 3),
    fmt(as.numeric(r[["delta_beta"]])),
    fmt(as.numeric(r[["HR_per_SD"]])),
    format.pval(as.numeric(r[["p_cox"]]), digits = 3))) else "（无）",
  "",
  "## 论文 Methods（英文草稿）",
  paste0("> \"DNA methylation (Illumina HumanMethylation450) and gene expression data ",
         "for TCGA-LIHC were obtained from the UCSC Xena hub. Probe-level analyses ",
         "were performed in R: tumour versus adjacent-normal beta-value differences ",
         "were tested by Wilcoxon rank-sum tests; methylation-expression ",
         "associations were assessed by Spearman correlation in tumour samples; and ",
         "probe-level prognostic value was evaluated by Cox proportional-hazards ",
         "models (per standard deviation of beta value) and Kaplan-Meier analysis ",
         "with median methylation split, with Benjamini-Hochberg correction across ",
         "probes. These analyses reproduce the functionality of MethSurv, MEXPRESS ",
         "and UALCAN.\""),
  "",
  "## 说明",
  "- `region` 列来自 Illumina 官方 450K 注释（TSS1500/TSS200/1stExon/5'UTR/Body/3'UTR）",
  "- `survival_consistent`：风险基因（HMOX1/TFRC/E2F1）甲基化高应为保护（HR<1）；保护基因（SOCS2）甲基化高应为风险（HR>1）",
  "- 历史`mech`列的silencing仅编码甲基化–表达负相关且P<0.05；不代表功能性沉默证据。",
  "",
  "## 输出文件",
  "- methylation_probe_results.csv（全部探针的 4 项分析结果）",
  "- methylation_best_probe_per_gene.csv（每基因最优探针）",
  "- methylation_panel_<GENE>.png（肿瘤vs正常 / 甲基化vs表达 / KM 三联图）"
)
writeLines(md, file.path(OUT, "methylation_methods_and_results.md"))

cat("\n== step13 done ==\n")
cat("outputs ->", OUT, "\n")
