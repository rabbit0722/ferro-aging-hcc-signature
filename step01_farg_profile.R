# =============================================================================
# step01_farg_profile.R
# -----------------------------------------------------------------------------
# Profile the 95-gene ferro-aging gene set (FARGs) across 68 primate tissues
# using the three sheets of "Ferro-aging-基因列表.xlsx" (Cell Metab 2026 suppl).
#
# Input : Ferro-aging-基因列表.xlsx  (workspace root)
# Output: analysis_out/step01/
#   - farg_tissue_log2fc.csv        95 FARGs x 68 tissues (log2FC, NA = not tested)
#   - farg_tissue_sig.csv           signed -log10(padj) for padj<0.05 (0 otherwise)
#   - tissue_index.csv              per-tissue ferro-aging index
#   - module_assignment.csv         FARG -> functional module (editable)
#   - module_tissue_mean.csv        module x tissue mean log2FC
#   - farg_trajectory_counts.csv    per-FARG trajectory cluster counts
#   - fig01_log2fc_heatmap.png, fig02_module_heatmap.png,
#     fig03_trajectory_stacked.png
#
# Requirements: readxl, data.table, ComplexHeatmap (falls back to pheatmap),
#               ggplot2
# Run:  Rscript step01_farg_profile.R   (from the Ferro-aging workspace)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

## ------------------------- CONFIG -------------------------
XLSX_FILE  <- "Ferro-aging-基因列表.xlsx"
OUT_DIR    <- "analysis_out/step01"
CACHE_RDS  <- file.path(OUT_DIR, "xlsx_cache.rds")
FORCE_REREAD <- FALSE
LFC_CLIP   <- 2.5          # clip log2FC color scale to [-2.5, 2.5]
PADJ_CUT   <- 0.05
HEAT_W     <- 2400         # png width (px)
HEAT_H     <- 1600
# -------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== step01: ferro-aging gene set (FARGs) tissue profile ==\n")
cat("R:", R.version.string, "\n\n")

## ------------------------- 1. Load xlsx (cached) -------------------------
if (file.exists(CACHE_RDS) && !FORCE_REREAD) {
  cat("[1/6] Loading cached data from", CACHE_RDS, "\n")
  L <- readRDS(CACHE_RDS)
  farg <- L$farg; deg <- L$deg; tr <- L$tr
} else {
  cat("[1/6] Reading xlsx (large trend sheet may take 1-3 min) ...\n")
  t0 <- Sys.time()
  # skip title rows only; the next row (real header) is auto-detected as names
  farg <- as.data.table(read_excel(XLSX_FILE, sheet = "Ferroaging gene set",
                                   skip = 1))
  deg  <- as.data.table(read_excel(XLSX_FILE, sheet = "Aging DEGs", skip = 3))
  cat("      ... Aging DEGs read, reading trend sheet ...\n")
  tr   <- as.data.table(read_excel(XLSX_FILE, sheet = "Gene expression trends",
                                   skip = 1))
  cat("      ... done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  saveRDS(list(farg = farg, deg = deg, tr = tr), CACHE_RDS)
  cat("      cached ->", CACHE_RDS, "\n")
}

# safety: force expected column names regardless of header auto-detection
setnames(deg, names(deg), c("Tissue", "Geneid", "SYMBOL", "log2FoldChange", "padj"))
setnames(tr,  names(tr),  c("Tissue", "Geneid", "SYMBOL", "Cluster"))
setnames(farg, names(farg), c("Geneid", "SYMBOL"))
farg <- farg[!is.na(Geneid) & !is.na(SYMBOL)]

cat("      FARGs:", nrow(farg), "| DEG rows:", nrow(deg),
    "| trend rows:", nrow(tr), "\n\n")

## ------------------------- 2. Functional modules -------------------------
# Literature-informed split of the 95 FARGs (EDIT THIS LIST as you like;
# unlisted genes go to "M5_unassigned" with a warning).
MODULES <- list(
  "M1_iron_ferroptosis" = c("ALOX15","LPCAT3","PTGS2","KEAP1","SOD1","TFRC",
                            "HMOX1","LCN2","SLC1A5","MCU","DUOX1","NOX4","MPO",
                            "SAT1","ATG3","LACTB"),
  "M2_senescence"       = c("IGFBP7","ATF3","TXNIP","HMGB1","HBP1","SMURF2",
                            "SETD7","KDM6B","E2F1","E2F3","TBX2","KLF6","SP1",
                            "RUNX3","MEN1","SMARCB1","BRD7","FBXO31"),
  "M3_inflammation"     = c("IL6","IL1B","IFNG","NLRP3","TLR4","TNFAIP3",
                            "TNFAIP1","CXCL10","IRF1","IRF7","IRF9","SOCS1",
                            "SOCS2","S100A8","CD74","CD82","SLAMF8","PADI4",
                            "LGMN"),
  "M4_signaling_other"  = c("HIF1A","YAP1","WWTR1","ZEB1","WNT5A","EPHA2",
                            "EPHA4","ERN1","HERPUD1","MAP3K14","MAPK1","MAPK14",
                            "PRKD1","NUAK2","DYRK1A","PDE4B","FOSL1","EGR1",
                            "NR1D1","NR2F2","GMFB","LIFR","ABCC1","ACVR1B",
                            "BAP1","CAVIN1","CDO1","COX7A1","CTSB","DPP4",
                            "DPEP1","EMP1","EBF3","LOX","PPP2R2B","PTBP1",
                            "RBM3","SNCA","SPATA2","ICA1","BCL6","EDN1")
)
mod_df <- rbindlist(lapply(names(MODULES), function(m)
  data.table(SYMBOL = MODULES[[m]], Module = m)), fill = TRUE)
farg <- merge(farg, mod_df, by = "SYMBOL", all.x = TRUE)
farg[is.na(Module), Module := "M5_unassigned"]
if (any(farg$Module == "M5_unassigned")) {
  cat("WARNING: unassigned genes:", paste(farg[Module == "M5_unassigned"]$SYMBOL,
                                          collapse = ", "), "\n")
}
fwrite(farg, file.path(OUT_DIR, "module_assignment.csv"))
cat("[2/6] Module assignment saved. Counts:\n")
print(farg[, .N, by = Module][order(Module)])

## ------------------------- 3. FARG x tissue matrices -------------------------
# 3a. log2FC matrix: all 95 FARGs x all 68 tissues (NA = not tested in tissue)
fg <- deg[Geneid %in% farg$Geneid]                       # keep FARG rows
setkey(farg, Geneid); setkey(fg, Geneid)
fg <- farg[fg, on = "Geneid", nomatch = 0]               # add SYMBOL/Module
all_tissues <- sort(unique(deg$Tissue))
all_sym     <- farg$SYMBOL

make_mat <- function(wide_dt) {
  m <- as.data.frame(wide_dt)
  rownames(m) <- m$SYMBOL; m$SYMBOL <- NULL
  missing_cols <- setdiff(all_tissues, colnames(m))
  if (length(missing_cols)) m[, missing_cols] <- NA_real_   # add NA cols
  m <- m[, all_tissues, drop = FALSE]
  missing_rows <- setdiff(all_sym, rownames(m))
  if (length(missing_rows)) {                               # add NA rows
    add <- matrix(NA_real_, nrow = length(missing_rows), ncol = ncol(m),
                  dimnames = list(missing_rows, colnames(m)))
    m <- rbind(m, add)
  }
  m <- m[match(all_sym, rownames(m)), , drop = FALSE]
  rownames(m) <- all_sym
  as.matrix(m)
}
lfc_mat <- make_mat(dcast(fg, SYMBOL ~ Tissue, value.var = "log2FoldChange"))

# 3b. signed -log10(padj) matrix (0 when not significant / not tested)
fg[, score := ifelse(padj < PADJ_CUT & !is.na(padj),
                     -log10(padj) * sign(log2FoldChange), 0)]
sc_mat <- make_mat(dcast(fg, SYMBOL ~ Tissue, value.var = "score"))

tissues <- colnames(lfc_mat)
cat("\n[3/6] Matrices built:", nrow(lfc_mat), "genes x", ncol(lfc_mat),
    "tissues\n")
cat("      FARGs with at least one DEG entry:", sum(rowSums(!is.na(lfc_mat)) > 0),
    "/", nrow(lfc_mat), "\n")
miss <- farg[!Geneid %in% unique(deg$Geneid)]
cat("      FARGs absent from DEG table entirely:",
    if (nrow(miss)) paste(miss$SYMBOL, collapse = ", ") else "none", "\n")
cat("      Significant (padj<", PADJ_CUT, ") FARG-tissue pairs:",
    sum(fg$padj < PADJ_CUT, na.rm = TRUE), "\n")

## ------------------------- 4. Tissue ferro-aging index -------------------------
sig_fg <- fg[padj < PADJ_CUT & !is.na(padj)]
idx <- fg[, .(
  n_tested = .N,
  n_sig    = sum(padj < PADJ_CUT, na.rm = TRUE),
  n_up     = sum(log2FoldChange > 0 & padj < PADJ_CUT, na.rm = TRUE),
  n_down   = sum(log2FoldChange < 0 & padj < PADJ_CUT, na.rm = TRUE),
  score_mean   = mean(score, na.rm = TRUE),
  mean_lfc     = mean(log2FoldChange, na.rm = TRUE),
  mean_abs_lfc = mean(abs(log2FoldChange), na.rm = TRUE)
), by = Tissue]
idx[, balance := ifelse(n_sig > 0, (n_up - n_down) / n_sig, 0)]  # -1 .. +1
idx[, confidence := ifelse(n_sig >= 5, "high", "low")]
idx <- idx[order(-score_mean)]
fwrite(idx, file.path(OUT_DIR, "tissue_index.csv"))
cat("\n[4/6] Tissue ferro-aging index (top 10, high-confidence: n_sig>=5):\n")
top10 <- idx[confidence == "high"][1:min(10, sum(idx$confidence == "high"))]
print(top10[, .(Tissue, n_sig, n_up, n_down,
                score_mean = round(score_mean, 3),
                balance = round(balance, 2))])
cat("\n      Bottom 5 (high-confidence):\n")
hi <- idx[confidence == "high"]
print(tail(hi, 5)[, .(Tissue, n_sig, n_up, n_down,
                      score_mean = round(score_mean, 3),
                      balance = round(balance, 2))])

## ------------------------- 5. Figures -------------------------
cat("\n[5/6] Drawing figures ...\n")

# Fig 1: log2FC heatmap (95 x 68), rows split by module
use_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
p1 <- file.path(OUT_DIR, "fig01_log2fc_heatmap.png")
m1 <- lfc_mat
m1[is.na(m1)] <- 0
m1[m1 >  LFC_CLIP] <-  LFC_CLIP
m1[m1 < -LFC_CLIP] <- -LFC_CLIP
if (use_ch) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  mod <- farg$Module[match(rownames(m1), farg$SYMBOL)]
  col_fun <- colorRamp2(c(-LFC_CLIP, 0, LFC_CLIP),
                        c("#2166AC", "#F7F7F7", "#B2182B"))
  ht <- Heatmap(m1, name = "log2FC", col = col_fun,
                row_split = factor(mod), cluster_row_slices = FALSE,
                show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
                show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                column_names_rot = 45, border = TRUE,
                heatmap_legend_param = list(title = "log2FC"))
  png(p1, width = HEAT_W, height = HEAT_H, res = 120)
  draw(ht)
  dev.off()
} else if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(m1, filename = p1, width = 20, height = 13,
                     cluster_cols = TRUE, fontsize_row = 4, fontsize_col = 5,
                     color = colorRampPalette(c("#2166AC", "#F7F7F7",
                                                "#B2182B"))(100))
} else {
  png(p1, width = HEAT_W, height = HEAT_H); heatmap(m1); dev.off()
}
cat("      saved:", p1, "\n")

# Fig 2: module x tissue mean log2FC heatmap (4 x 68)
p2 <- file.path(OUT_DIR, "fig02_module_heatmap.png")
mod_lfc <- merge(fg[, .(Geneid, SYMBOL, Module, Tissue, log2FoldChange)],
                 farg[, .(SYMBOL, Module)], by = c("SYMBOL", "Module"))
mm <- dcast(mod_lfc, Module ~ Tissue, value.var = "log2FoldChange",
            fun.aggregate = mean)
mm_mat <- as.matrix(mm[, -1]); rownames(mm_mat) <- mm$Module
mm_mat[is.na(mm_mat)] <- 0
fwrite(as.data.table(mm_mat, keep.rownames = "Module"),
       file.path(OUT_DIR, "module_tissue_mean.csv"))
if (use_ch) {
  col_fun2 <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
  ht2 <- Heatmap(mm_mat, name = "mean log2FC", col = col_fun2,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 7),
                 column_names_rot = 45, border = TRUE)
  png(p2, width = 2400, height = 700, res = 120)
  draw(ht2); dev.off()
} else {
  pheatmap::pheatmap(mm_mat, filename = p2, width = 18, height = 5,
                     fontsize_row = 10, fontsize_col = 6,
                     color = colorRampPalette(c("#2166AC", "#F7F7F7",
                                                "#B2182B"))(100))
}
cat("      saved:", p2, "\n")

# Fig 3: trajectory cluster composition (per module, fraction)
p3 <- file.path(OUT_DIR, "fig03_trajectory_stacked.png")
trg <- tr[Geneid %in% farg$Geneid]
trg <- merge(trg, farg[, .(Geneid, SYMBOL, Module)], by = "Geneid")
tc <- trg[, .N, by = .(Module, Cluster)]
tc[, frac := N / sum(N), by = Module]
fwrite(dcast(tc, Module ~ Cluster, value.var = "N", fill = 0),
       file.path(OUT_DIR, "farg_trajectory_counts.csv"))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  g <- ggplot(tc, aes(x = Module, y = frac, fill = Cluster)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = c("Cluster U" = "#E64B35", "Cluster D" = "#4DBBD5",
                                 "Cluster UD" = "#7E6148", "Cluster DU" = "#B09C85")) +
    labs(x = NULL, y = "Fraction of FARG trajectory assignments",
         title = "Aging trajectory clusters of the 95 FARGs (68 tissues)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(p3, g, width = 7, height = 5, dpi = 200)
} else {
  png(p3, width = 1400, height = 900, res = 150)
  tab <- dcast(tc, Module ~ Cluster, value.var = "N", fill = 0)
  m <- as.matrix(tab[, -1]); rownames(m) <- tab$Module
  barplot(t(m), legend.text = colnames(m), col = 1:4,
          args.legend = list(x = "topright", bty = "n"))
  dev.off()
}
cat("      saved:", p3, "\n")

## ------------------------- 6. Save matrices + summary -------------------------
fwrite(as.data.table(lfc_mat, keep.rownames = "SYMBOL"),
       file.path(OUT_DIR, "farg_tissue_log2fc.csv"))
fwrite(as.data.table(sc_mat, keep.rownames = "SYMBOL"),
       file.path(OUT_DIR, "farg_tissue_sig.csv"))
cat("\n[6/6] Matrices saved to", OUT_DIR, "\n")

cat("\n== step01 summary ============================================\n")
cat("Top 5 FARGs by number of tissues with padj<0.05:\n")
topg <- sig_fg[, .N, by = SYMBOL][order(-N)][1:5]
print(topg)
cat("Top 5 tissues by mean |log2FC| of significant FARGs:\n")
print(idx[order(-mean_abs_lfc)][1:5, .(Tissue, n_sig,
                                       mean_abs_lfc = round(mean_abs_lfc, 3))])
cat("Trajectory clusters of FARGs overall:\n")
print(trg[, .N, by = Cluster][order(-N)])
cat("===============================================================\n")
