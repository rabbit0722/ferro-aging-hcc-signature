# =============================================================================
# export_source_tables.R —— 导出中间源表（全基因 DEG 表）
# 补齐 stage2 台账的"火山图全量 tt 未归档"缺口。
# DCA 逐时间表与校准点需在原 step09/step05 中补导出（此处不做近似重算）。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(limma)})
OUT <- "analysis_out/hcc_out"

L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr <- L$expr; group <- L$group
g_f <- factor(group, levels = c("Normal", "Tumor"))
fit <- eBayes(lmFit(expr, model.matrix(~ g_f)))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt <- as.data.table(tt, keep.rownames = "gene")
fwrite(tt, file.path(OUT, "full_gene_deg.csv"))
cat("full_gene_deg.csv:", nrow(tt), "genes | 列:", paste(names(tt), collapse = ", "), "\n")
cat("（列：gene, logFC, AveExpr, t, P.Value, adj.P.Val, B）\n")
