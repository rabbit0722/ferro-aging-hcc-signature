# =============================================================================
# step11_hcc_genomics.R —— 签名基因基因组层分析（可选，需大下载）
# -----------------------------------------------------------------------------
# [A] 拷贝数变异 CNV（TCGAbiolinks 基因水平 CNV，下载 ~10-30 MB）
# [B] 启动子甲基化（450K，下载 ~1 GB + 注释包 ~150 MB；不想下载请用
#     MethSurv 网页替代：https://bio-bits.eu/MethSurv 或
#     https://bioinfo.uth.edu/MethSurv —— 直接查 4 个基因在 LIHC 的甲基化
#     与预后的关系并截图）
#
# 依赖: TCGAbiolinks（你已装）+ SummarizedExperiment；甲基化部分还需要
#       BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# 输出: analysis_out/hcc_out/genomics/
# 运行: setwd("path/to/Ferro-aging"); source("step11_hcc_genomics.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
OUT <- "analysis_out/hcc_out/genomics"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== step11: signature gene CNV + methylation ==\n")

rk <- fread("analysis_out/hcc_out/risk_training.csv")
sig <- fread("analysis_out/hcc_out/signature_genes.csv")$gene
L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_t <- L$expr[, L$group == "Tumor"]

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  cat("TCGAbiolinks not installed. Install: BiocManager::install('TCGAbiolinks')\n")
  quit(save = "no")
}
suppressPackageStartupMessages(library(TCGAbiolinks))

## ---- [A] CNV ------------------------------------------------------------------
cat("\n[A] gene-level CNV (GDC download ~10-30 MB) ...\n")
cnv_ok <- tryCatch({
  query <- GDCquery(project = "TCGA-LIHC",
                    data.category = "Copy Number Variation",
                    data.type = "Gene Level Copy Number",
                    access = "open")
  GDCdownload(query, method = "api")
  cnv <- GDCprepare(query)          # SummarizedExperiment, 基因 x 样本
  m <- as.matrix(assay(cnv))
  # rownames 可能是 gene symbol；若是 Ensembl 则映射（GDCprepare 基因水平
  # CNV 默认返回 gene symbol 行名）
  common <- intersect(sig, rownames(m))
  if (length(common) == 0) {
    cat("  signature genes not found in CNV rownames (", head(rownames(m), 5),
        "...); skipped.\n", sep = "")
    FALSE
  } else {
    cnv_tab <- rbindlist(lapply(common, function(g) {
      v <- m[g, ]
      smp <- substr(colnames(m), 1, 15)
      gain <- sum(v > 0.3, na.rm = TRUE) / sum(!is.na(v)) * 100
      loss <- sum(v < -0.3, na.rm = TRUE) / sum(!is.na(v)) * 100
      grp <- rk$risk_group[match(smp, rk$sample)]
      okk <- !is.na(grp) & !is.na(v)
      p_gain <- tryCatch(wilcox.test(v[okk & grp == "High"],
                                     v[okk & grp == "Low"])$p.value,
                         error = function(e) NA)
      data.table(gene = g, pct_gain = round(gain, 1),
                 pct_loss = round(loss, 1), wilcox_gainloss_p = p_gain)
    }))
    fwrite(cnv_tab, file.path(OUT, "signature_gene_cnv.csv"))
    cat("  CNV summary:\n"); print(cnv_tab)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      suppressPackageStartupMessages(library(ggplot2))
      ml <- melt(cnv_tab[, .(gene, Gain = pct_gain, Loss = -pct_loss)],
                 id.vars = "gene", variable.name = "type", value.name = "pct")
      g <- ggplot(ml, aes(gene, pct, fill = type)) +
        geom_col() +
        scale_fill_manual(values = c(Gain = "#B2182B", Loss = "#2166AC")) +
        labs(x = NULL, y = "Frequency (%)",
             title = "Signature gene CNV in TCGA-LIHC") +
        theme_bw(base_size = 13)
      ggsave(file.path(OUT, "signature_gene_cnv.png"), g,
             width = 6, height = 4.5, dpi = 200)
    }
    TRUE
  }
}, error = function(e) {
  cat("  CNV failed:", conditionMessage(e), "\n"); FALSE
})

## ---- [B] Methylation (450K, ~1 GB download) --------------------------------------
cat("\n[B] promoter methylation (WARNING: ~1 GB download) ...\n")
meth_ok <- tryCatch({
  query2 <- GDCquery(project = "TCGA-LIHC",
                     data.category = "DNA Methylation",
                     platform = "Illumina Human Methylation 450",
                     data.type = "Methylation Beta Value",
                     access = "open")
  GDCdownload(query2, method = "api")
  meth <- GDCprepare(query2)        # 探针 x 样本 (beta)
  beta <- as.matrix(assay(meth))
  if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                        quietly = TRUE)) {
    cat("  需要注释包，先安装:\n")
    cat("  BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19')\n")
    cat("  （或跳过本模块，用 MethSurv 网页查 4 基因的甲基化-预后图）\n")
    FALSE
  } else {
    suppressPackageStartupMessages(
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    met_tab <- rbindlist(lapply(sig, function(g) {
      # 启动子区探针 (TSS1500/TSS200/5'UTR/1stExon)
      in_gene <- grepl(paste0("(^|;)", g, "($|;)"), ann$UCSC_RefGene_Name)
      in_pro  <- grepl("TSS1500|TSS200|5'UTR|1stExon", ann$UCSC_RefGene_Group)
      probes <- rownames(ann)[in_gene & in_pro]
      if (length(probes) == 0)
        return(data.table(gene = g, n_probes = 0L))
      pb <- beta[probes, , drop = FALSE]
      # 样本水平均值 beta；样本 id 取前 15 位匹配表达
      mbeta <- colMeans(pb, na.rm = TRUE)
      smp <- substr(names(mbeta), 1, 15)
      # 相同样本取均值
      agg <- tapply(mbeta, smp, mean, na.rm = TRUE)
      smp2 <- names(agg)
      ge <- as.numeric(expr_t[g, smp2[smp2 %in% colnames(expr_t)]])
      # 甲基化-表达 spearman
      rho <- tryCatch(cor.test(agg[smp2 %in% colnames(expr_t)],
                               as.numeric(expr_t[g, smp2[smp2 %in%
                                                          colnames(expr_t)]]),
                               method = "spearman")$estimate,
                      error = function(e) NA)
      grp <- rk$risk_group[match(smp2, rk$sample)]
      okk <- !is.na(grp)
      pw <- tryCatch(wilcox.test(agg[okk & grp == "High"],
                                 agg[okk & grp == "Low"])$p.value,
                     error = function(e) NA)
      # 甲基化中位分组的生存 HR
      srv <- rk[match(smp2, rk$sample), .(OS.time, OS)]
      d4 <- cbind(srv, beta_m = as.numeric(agg))
      d4 <- d4[!is.na(OS.time) & OS.time > 0 & !is.na(beta_m)]
      d4[, grp_m := ifelse(beta_m > median(beta_m), "HighMeth", "LowMeth")]
      hr <- tryCatch({
        f <- coxph(Surv(OS.time, OS) ~ (grp_m == "HighMeth"), data = d4)
        exp(coef(f))
      }, error = function(e) NA)
      data.table(gene = g, n_probes = length(probes),
                 cor_meth_expr = as.numeric(rho),
                 wilcox_HighLow_p = pw, KM_HR_highMeth = hr)
    }))
    fwrite(met_tab, file.path(OUT, "signature_gene_methylation.csv"))
    cat("  methylation summary:\n"); print(met_tab)
    TRUE
  }
}, error = function(e) {
  cat("  methylation failed:", conditionMessage(e), "\n"); FALSE
})

cat("\n== step11 done ==\n")
