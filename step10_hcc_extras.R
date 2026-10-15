# =============================================================================
# step10_hcc_extras.R —— 补强分析（快，~5 分钟，无需大下载）
# -----------------------------------------------------------------------------
# [B] GSVA Hallmark 通路活性与风险分数的相关性热图
# [C] 4 个签名基因的突变频率（高低风险组比较，需要 step07 的 MAF）
# [D] CMap(clue.io) 查询输入文件导出（高低风险组差异基因列表）
#
# 输入 : analysis_out/hcc_data/hcc_data.rds
#        analysis_out/hcc_out/risk_training.csv / signature_genes.csv
#        MAF 文件（与 step07 相同来源，自动查找）
# 输出 : analysis_out/hcc_out/extras/
# 运行 : setwd("path/to/Ferro-aging"); source("step10_hcc_extras.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
OUT <- "analysis_out/hcc_out/extras"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== step10: extras (GSVA cor / mutation / CMap) ==\n")

## ---- 数据加载 -----------------------------------------------------------------
L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_t <- L$expr[, L$group == "Tumor"]                 # 369 肿瘤样本
rk <- fread("analysis_out/hcc_out/risk_training.csv")
sig <- fread("analysis_out/hcc_out/signature_genes.csv")$gene
mods <- fread("analysis_out/step01/module_assignment.csv")
mods <- mods[, .(SYMBOL, Module)]
mod_list <- split(mods$SYMBOL, mods$Module)
farg_all <- unique(mods$SYMBOL)

use_s <- intersect(rk$sample, colnames(expr_t))
expr_u <- expr_t[, use_s]
rk_u <- rk[match(use_s, rk$sample)]
risk_z <- as.numeric(scale(rk_u$risk))
surv_d <- data.frame(time = rk_u$OS.time, status = rk_u$OS)
ok_s <- !is.na(surv_d$time) & surv_d$time > 0 & !is.na(surv_d$status)
surv_d <- surv_d[ok_s, ]
risk_z2 <- risk_z[ok_s]
cat("samples:", nrow(surv_d), "| signature:", paste(sig, collapse = ", "), "\n")

ssgsea <- function(mat, sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("GSVA not installed")
  r <- tryCatch(
    GSVA::gsva(as.matrix(mat), sets, method = "ssgsea",
               ssgsea.norm = TRUE, verbose = FALSE),
    error = function(e) GSVA::gsva(GSVA::ssgseaParam(as.matrix(mat), sets),
                                   verbose = FALSE))
  t(r)
}

## [A] 模块贡献分解 —— 已移除（不再使用四模块划分，见修订记录）
## ---- [B] GSVA Hallmark 通路相关热图 ---------------------------------------------
cat("\n[B] GSVA Hallmark pathway correlation with risk ...\n")
h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) NULL)
if (!is.null(h) && nrow(h) > 0) {
  hs <- lapply(split(h$gene_symbol, h$gs_name),
               function(g) intersect(g, rownames(expr_u)))
  hs <- hs[lengths(hs) >= 15]
  hsc <- tryCatch(ssgsea(expr_u, hs), error = function(e) NULL)
  if (!is.null(hsc)) {
    hsc <- hsc[ok_s, , drop = FALSE]
    cor_tab <- rbindlist(lapply(colnames(hsc), function(sn) {
      ct <- cor.test(as.numeric(hsc[, sn]), risk_z2, method = "spearman")
      data.table(pathway = sub("^HALLMARK_", "", sn),
                 rho = as.numeric(ct$estimate), p = ct$p.value)
    }))
    cor_tab[, FDR := p.adjust(p, "BH")]
    setorder(cor_tab, -rho)
    fwrite(cor_tab, file.path(OUT, "gsva_hallmark_risk_cor.csv"))
    cat("  top pathways (rho>0):\n")
    print(head(cor_tab, 10))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      top <- cor_tab[1:min(20, nrow(cor_tab))]
      top[, pathway := factor(pathway, levels = rev(pathway))]
      g2 <- ggplot(top, aes(rho, pathway, fill = rho)) +
        geom_col() +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
        labs(x = "Spearman rho (pathway score ~ risk)",
             y = NULL, title = "Hallmark pathway activity vs FARGs risk") +
        theme_bw(base_size = 12) + theme(legend.position = "none")
      ggsave(file.path(OUT, "gsva_hallmark_risk_bar.png"), g2,
             width = 8, height = 6, dpi = 200)
    }
    # 热图: 相关通路 x 样本（按风险排序）
    keep <- cor_tab[order(-abs(rho))][1:min(15, nrow(cor_tab))]$pathway
    hm <- hsc[, paste0("HALLMARK_", keep), drop = FALSE]
    ord <- order(risk_z2)
    mh <- t(scale(hm[ord, ]))
    png(file.path(OUT, "gsva_hallmark_heatmap.png"), width = 2000,
        height = 1200, res = 150)
    heatmap(mh, Colv = NA, Rowv = NA, col = colorRampPalette(
      c("#2166AC", "white", "#B2182B"))(50), scale = "none",
      labCol = "", main = "Hallmark pathway scores (samples ordered by risk)")
    dev.off()
    cat("  saved gsva_hallmark_risk_bar.png / gsva_hallmark_heatmap.png\n")
  }
} else {
  cat("  msigdbr unavailable (offline?); module [B] skipped.\n")
}

## ---- [C] 签名基因突变频率（高低风险组） -------------------------------------------
cat("\n[C] signature gene mutation frequency ...\n")
find_maf <- function() {
  cand_dirs <- c("analysis_out/data", "data", ".", "..", "../data")
  prefixes <- c("TCGA.LIHC.mutect", "TCGA.LIHC.mutect2_snv",
                "mc3.v0.2.8.PUBLIC.maf.gz", "mc3")
  for (d in cand_dirs) {
    if (!dir.exists(d)) next
    fs <- list.files(d, recursive = TRUE, full.names = TRUE)
    for (pr in prefixes) {
      hit <- fs[grepl(pr, basename(fs), ignore.case = TRUE) &
                  grepl("\\.(tsv|txt|maf)(\\.gz)?$", fs)]
      if (length(hit)) return(hit[1])
    }
  }
  NA_character_
}
MAF_FILE <- find_maf()
if (!is.na(MAF_FILE) && file.exists(MAF_FILE)) {
  maf <- fread(MAF_FILE)
  setnames(maf, names(maf), sub("^#", "", names(maf)))
  if (!"Tumor_Sample_Barcode" %in% names(maf) && "sampleId" %in% names(maf))
    setnames(maf, "sampleId", "Tumor_Sample_Barcode")
  if (!"Hugo_Symbol" %in% names(maf) && "hugoGeneSymbol" %in% names(maf))
    setnames(maf, "hugoGeneSymbol", "Hugo_Symbol")
  if (!"Variant_Classification" %in% names(maf) && "mutationType" %in% names(maf))
    setnames(maf, "mutationType", "Variant_Classification")
  maf[, s := substr(Tumor_Sample_Barcode, 1, 15)]
  non_syn <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
               "Translation_Start_Site", "Nonstop_Mutation")
  maf <- maf[Variant_Classification %in% non_syn & Hugo_Symbol %in% sig]
  grp <- rk$risk_group[match(maf$s, rk$sample)]
  maf[, group := grp]
  maf <- maf[!is.na(group)]
  mut <- rbindlist(lapply(sig, function(g) {
    m_h <- uniqueN(maf$s[maf$Hugo_Symbol == g & maf$group == "High"])
    m_l <- uniqueN(maf$s[maf$Hugo_Symbol == g & maf$group == "Low"])
    n_h <- sum(rk_u$risk_group == "High"); n_l <- sum(rk_u$risk_group == "Low")
    pv <- tryCatch(fisher.test(matrix(c(m_h, n_h - m_h, m_l, n_l - m_l), 2))$p.value,
                   error = function(e) NA)
    data.table(gene = g, mut_High = m_h, pct_High = m_h / n_h * 100,
               mut_Low = m_l, pct_Low = m_l / n_l * 100, fisher_p = pv)
  }))
  fwrite(mut, file.path(OUT, "signature_gene_mutation.csv"))
  cat("  mutation frequency by group:\n"); print(mut)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ml <- melt(mut[, .(gene, High = pct_High, Low = pct_Low)],
               id.vars = "gene", variable.name = "group", value.name = "pct")
    g3 <- ggplot(ml, aes(gene, pct, fill = group)) +
      geom_col(position = "dodge", width = 0.6) +
      scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
      labs(x = NULL, y = "Mutation frequency (%)",
           title = "Signature gene mutations by risk group") +
      theme_bw(base_size = 13)
    ggsave(file.path(OUT, "signature_gene_mutation.png"), g3,
           width = 6, height = 4.5, dpi = 200)
  }
} else {
  cat("  MAF not found; module [C] skipped (needs step07 MAF).\n")
}

## ---- [D] CMap(clue.io) 输入文件导出 ---------------------------------------------
cat("\n[D] CMap input lists (high vs low risk DEGs) ...\n")
if (requireNamespace("limma", quietly = TRUE)) {
  suppressPackageStartupMessages(library(limma))
  grp <- factor(rk_u$risk_group, levels = c("Low", "High"))
  design <- model.matrix(~ grp)
  fit <- lmFit(expr_u, design)
  fit <- eBayes(fit)
  tt <- as.data.table(topTable(fit, coef = 2, number = Inf, sort.by = "p"))
  tt[, gene := rownames(topTable(fit, coef = 2, number = Inf, sort.by = "p"))]
  up <- tt[logFC > 0 & adj.P.Val < 0.05][order(-logFC)]
  dn <- tt[logFC < 0 & adj.P.Val < 0.05][order(logFC)]
  fwrite(up[, .(gene, logFC)], file.path(OUT, "cmap_up_genes.csv"))
  fwrite(dn[, .(gene, logFC)], file.path(OUT, "cmap_down_genes.csv"))
  writeLines(up$gene[1:min(150, nrow(up))],
             file.path(OUT, "cmap_up_genes.txt"))
  writeLines(dn$gene[1:min(150, nrow(dn))],
             file.path(OUT, "cmap_down_genes.txt"))
  cat("  up:", nrow(up), "genes | down:", nrow(dn), "genes\n")
  cat("  上传 clue.io -> Query -> 粘贴 cmap_up_genes.txt + cmap_down_genes.txt\n")
} else {
  cat("  limma not installed; module [D] skipped.\n")
}

cat("\n== step10 done ==\n")
