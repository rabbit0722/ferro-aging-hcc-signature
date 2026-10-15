# =============================================================================
# coexpression_clustering.R —— 检验四模块是否被共表达数据支持
# 用 TCGA-LIHC 肿瘤样本对 95 个 FARGs 做层次聚类，看是否自然形成
# "铁死亡/衰老/炎症"类簇，并与手工 M1–M4 标签比对（调整兰德指数 ARI）。
# 输出：analysis_out/step01/coexpression_*
# =============================================================================
suppressPackageStartupMessages({library(data.table)})

OUT <- "analysis_out/step01"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 数据 ----
L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_t <- L$expr[, L$group == "Tumor"]
mod <- fread(file.path(OUT, "module_assignment.csv"))
genes <- intersect(mod$SYMBOL, rownames(expr_t))
cat("肿瘤样本:", ncol(expr_t), "| 参与聚类的 FARGs:", length(genes), "\n")
miss <- setdiff(mod$SYMBOL, rownames(expr_t))
if (length(miss)) cat("缺失基因:", paste(miss, collapse=", "), "\n")

X <- t(expr_t[genes, , drop = FALSE])   # 样本 x 基因
C <- cor(X, method = "spearman")         # 基因 x 基因相关
D <- as.dist(1 - C)                      # 距离
hc <- hclust(D, method = "average")

# ---- 调整兰德指数（手动实现，避免额外依赖）----
ARI <- function(a, b) {
  a <- as.factor(a); b <- as.factor(b)
  n <- length(a)
  if (n != length(b)) stop("length mismatch")
  # 配对计数
  tbl <- table(a, b)
  ni <- rowSums(tbl); nj <- colSums(tbl)
  comb2 <- function(x) x * (x - 1) / 2
  sum_tbl <- sum(comb2(tbl))
  sum_ni <- sum(comb2(ni)); sum_nj <- sum(comb2(nj))
  expected <- sum_ni * sum_nj / comb2(n)
  max_index <- (sum_ni + sum_nj) / 2
  if (max_index == expected) return(1)
  (sum_tbl - expected) / (max_index - expected)
}

orig <- mod$Module[match(genes, mod$SYMBOL)]

# ---- 不同 k 的 ARI + 平均轮廓宽度 ----
cat("\n=== 不同聚类数 k 的 ARI（数据簇 vs 手工模块）===\n")
ari_tab <- rbindlist(lapply(2:8, function(k) {
  cl <- cutree(hc, k = k)
  data.table(k = k, ARI = round(ARI(cl, orig), 3))
}))
print(ari_tab)

# ---- 主结果：k = 4（与 M1–M4 对齐）----
k <- 4
cl4 <- cutree(hc, k = k)
ari4 <- ARI(cl4, orig)
cat(sprintf("\n[k=4] ARI = %.3f\n", ari4))
conting <- as.data.frame.matrix(table(DataCluster = cl4,
                                      Module = factor(orig, levels = sort(unique(orig)))))
print(conting)

res <- data.table(SYMBOL = genes,
                  Module_orig = orig,
                  DataCluster_k4 = cl4)
# 每个数据簇里各模块的构成
res[, .N, by = .(DataCluster_k4, Module_orig)][order(DataCluster_k4)]
fwrite(res, file.path(OUT, "coexpression_gene_clusters.csv"))
fwrite(ari_tab, file.path(OUT, "coexpression_ARI_by_k.csv"))
fwrite(as.data.table(conting, keep.rownames = "DataCluster"),
       file.path(OUT, "coexpression_contingency_k4.csv"))

# ---- 4 个签名基因落在哪个数据簇 ----
sig <- c("HMOX1", "TFRC", "SOCS2", "E2F1")
cat("\n=== 4 个签名基因的数据簇归属 ===\n")
print(res[SYMBOL %in% sig])

# ---- 热图 + 树状图 ----
# 行注释：手工模块色 + 数据簇色
module_cols <- c(M1_iron_ferroptosis = "#E64B35", M2_senescence = "#4DBBD5",
                 M3_inflammation = "#00A087", M4_signaling_other = "#B09C85")
cluster_cols <- RColorBrewer::brewer.pal(8, "Set1")[1:k]
names(cluster_cols) <- paste0("C", 1:k)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pheatmap))
  ann <- data.frame(Module = orig,
                    DataCluster = paste0("C", cl4),
                    row.names = genes)
  ann_cols <- list(Module = module_cols,
                   DataCluster = setNames(cluster_cols, paste0("C", 1:k)))
  png(file.path(OUT, "coexpression_heatmap.png"), width = 1800, height = 1600, res = 150)
  pheatmap(C, cluster_rows = hc, cluster_cols = hc,
           show_rownames = TRUE, show_colnames = FALSE,
           fontsize_row = 7, annotation_row = ann,
           annotation_colors = ann_cols, annotation_names_row = TRUE,
           main = "TCGA-LIHC FARGs co-expression (Spearman), hierarchical clustering",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
  cat("\n保存 heatmap\n")
} else {
  cat("\npheatmap 不可用，改存纯树状图\n")
  png(file.path(OUT, "coexpression_dendrogram.png"), width = 1800, height = 900, res = 150)
  plot(hc, labels = genes, cex = 0.5, hang = -1,
       main = "FARGs 共表达层次聚类（平均连接，1-Spearman）")
  rect.hclust(hc, k = k, border = cluster_cols)
  dev.off()
}

cat("\n== 完成 ==\n")
cat("关键指标：ARI(数据簇 vs M1-M4) =", round(ari4, 3), "\n")
cat("解读：ARI 接近 0 = 数据不支持手工模块；接近 1 = 支持；>0.3 才有实质一致\n")
