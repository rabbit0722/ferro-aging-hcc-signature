# =============================================================================
# step08_hcc_singlecell.R
# -----------------------------------------------------------------------------
# Single-cell localization of the FARGs signature + virtual knockout of the
# core gene, using HCC scRNA-seq (GSE149614, Lu et al. 2021) as the default.
#
# Analyses:
#   1) cell-type proportions
#   2) marker QC: top wilcox markers per cell type (verify annotation)
#   3) dotplot/heatmap: signature genes across cell types (mean + pct)
#   4) ferro-aging module scores per cell type (ssGSEA): iron/ferroptosis,
#      senescence, inflammation, signaling + all_FARG; malignant vs others
#   5) virtual knockout with scTenifoldKnk (same tool as the DPEP1 paper):
#      KO the core signature gene in malignant cells -> differentially
#      regulated genes + GO enrichment + plotKO
#
# Data (local, in analysis_out/data/GSE149614/), three accepted layouts:
#   A) GSE149614 processed files (recommended; keep original names):
#      - GSE149614_HCC.scRNAseq.S71915.count.txt.gz   (158 MB, cells in columns)
#      - GSE149614_HCC.metadata.updated.txt.gz        (cell annotations)
#      from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/
#   B) 10x per-sample folders (barcodes/features/matrix, from GEO RAW.tar)
#   C) none -> SELF-TEST on synthetic data
# If metadata has no celltype column, Seurat is used to cluster and annotate.
#
# Download page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614
# Supplementary files: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/
#
# Output: analysis_out/hcc_out/singlecell/
#
# Requirements: data.table, Matrix, GSVA, ComplexHeatmap, scTenifoldKnk,
#               clusterProfiler, org.Hs.eg.db (Seurat optional, for annotation)
# Run:  Rscript step08_hcc_singlecell.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(GSVA)
})

OUT_DIR      <- "analysis_out/hcc_out/singlecell"
# data directory: try several common locations for the GSE149614 files
DATA_DIRS    <- c("analysis_out/data/GSE149614", "analysis_out/GSE149614",
                  "GSE149614", "data/GSE149614", "../GSE149614")
DATA_DIR     <- DATA_DIRS[file.exists(DATA_DIRS)][1]
if (is.na(DATA_DIR)) DATA_DIR <- DATA_DIRS[1]
MOD_FILE     <- "analysis_out/step01/module_assignment.csv"
SIG_FILE     <- "analysis_out/hcc_out/signature_genes.csv"
VKO_TARGET   <- NA_character_     # NA -> first gene of the signature
MALIG_PAT    <- "malignant|tumor|hepat|HCC"   # celltype regex
MAX_CELLS    <- 12000             # downsample cap (memory safety for 72k cells)
VKO_MAX_CELLS<- 3000              # cap for the knockout network
VKO_N_GENES  <- 800               # panel size for the knockout network
SEED         <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("== step08: single-cell localization + virtual knockout ==\n")

## ---- signature + modules ---------------------------------------------------------
sig <- if (file.exists(SIG_FILE)) fread(SIG_FILE) else
  data.table(gene = c("ALOX15", "ATG3", "BAP1"), coef = c(0.44, 0.37, 0.51))
sig_genes <- sig$gene
if (is.na(VKO_TARGET)) VKO_TARGET <- sig_genes[1]
cat("Signature genes:", paste(sig_genes, collapse = ", "),
    "| VKO target:", VKO_TARGET, "\n")
mods <- if (file.exists(MOD_FILE)) split(fread(MOD_FILE)$SYMBOL,
                                         fread(MOD_FILE)$Module) else
  list(all_FARG = sig_genes)
mods[["all_FARG"]] <- unique(unlist(mods))
farg_all <- unique(unlist(mods))          # used by marker QC in both modes

## ---- data loading ---------------------------------------------------------------
load_sc <- function() {
  A_cnt <- list.files(DATA_DIR, pattern = "count.*\\.txt(\\.gz)?$",
                      full.names = TRUE, ignore.case = TRUE)
  A_meta<- list.files(DATA_DIR, pattern = "metadata.*\\.txt(\\.gz)?$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(A_cnt) > 0 && length(A_meta) > 0) {
    cat("[data] GSE149614 count matrix + metadata found:\n",
        "      ", basename(A_cnt[1]), "\n      ", basename(A_meta[1]), "\n")
    # Memory-safe loading: sample cells DURING fread (select=), because the
    # full matrix is 24k genes x 71,915 cells (~14 GB as a dense matrix).
    hdr <- names(fread(A_cnt[1], nrows = 0))
    if (length(hdr) > 50000) {
      # genes in rows, cells in columns (GSE149614 layout): select cell cols
      cell_cols <- setdiff(hdr, hdr[1])
      if (length(cell_cols) > MAX_CELLS) {
        set.seed(SEED); cell_cols <- sample(cell_cols, MAX_CELLS)
      }
      m <- fread(A_cnt[1], select = c(1L, match(cell_cols, hdr)))
      names(m) <- c(hdr[1], cell_cols)
      gn <- m[[1]]
      counts <- as.matrix(m[, ..cell_cols]); rownames(counts) <- gn
    } else {
      # cells in rows: read fully, downsample rows, transpose
      m <- fread(A_cnt[1])
      if (nrow(m) > MAX_CELLS) {
        set.seed(SEED); m <- m[sample(.N, MAX_CELLS)]
      }
      m <- as.data.table(transpose(m, keep.names = names(m)[1]))
      gn <- m[[1]]; m[[1]] <- NULL
      counts <- as.matrix(m); rownames(counts) <- gn
    }
    rm(m); gc()
    cat("      cells kept after downsampling:", ncol(counts), "\n")
    meta <- fread(A_meta[1])
    cat("      metadata columns:", paste(names(meta), collapse = ", "), "\n")
    cell_col <- names(meta)[1]
    ct_cand <- grep("celltype|cell_type|major|annotation|cluster|type",
                    names(meta), value = TRUE, ignore.case = TRUE)
    ct_col <- if (length(ct_cand)) ct_cand[1] else NA_character_
    extra_cols <- intersect(c("site", "sample", "patient", "stage"),
                            names(meta)[-1])
    if (is.na(ct_col)) {
      meta <- data.table(cell = meta[[cell_col]])
      cat("      no celltype column detected; Seurat annotation will be used.\n")
    } else {
      keep <- c(cell_col, ct_col, extra_cols)
      meta <- meta[, ..keep]
      setnames(meta, c("cell", "celltype", extra_cols))
    }
    cat("      kept metadata columns:",
        paste(c("cell", "celltype", extra_cols), collapse = ", "), "\n")
    # reconcile cell ids (metadata may use full barcodes or stripped ids)
    # IMPORTANT: align meta ROWS to the downsampled counts COLUMNS (reorder,
    # not just filter) so celltype labels stay 1:1 with the kept cells.
    meta_keep <- c("cell", "celltype", extra_cols)
    cnt_cells <- colnames(counts)
    m2 <- match(cnt_cells, meta$cell)          # counts col -> meta row (NA if missing)
    if (sum(!is.na(m2)) < 100) {
      meta[, cell2 := sub("-\\d+$", "", cell)]
      m2 <- match(sub("-\\d+$", "", cnt_cells), meta$cell2)
    }
    keep <- !is.na(m2)
    counts <- counts[, keep, drop = FALSE]     # drop cells without metadata
    meta <- meta[m2[keep], ..meta_keep]        # align meta to counts column order
    return(list(counts = counts, meta = meta))
  }
  mtx_files <- list.files(DATA_DIR, pattern = "matrix\\.mtx(\\.gz)?$",
                          recursive = TRUE, full.names = TRUE)
  if (length(mtx_files) > 0) {
    cat("[data] 10x per-sample files found (", length(mtx_files), "samples).\n")
    smats <- list(); metas <- list()
    for (mf in mtx_files) {
      d <- dirname(mf)
      bf <- list.files(d, pattern = "barcodes", full.names = TRUE)[1]
      ff <- list.files(d, pattern = "features|genes", full.names = TRUE)[1]
      if (is.na(bf) || is.na(ff)) next
      mm <- readMM(if (grepl("\\.gz$", mf)) gzfile(mf) else mf)
      feat <- fread(if (grepl("\\.gz$", ff)) gzfile(ff) else ff, header = FALSE)
      barc <- fread(if (grepl("\\.gz$", bf)) gzfile(bf) else bf, header = FALSE)
      gn <- feat[[2]]; colnames(mm) <- barc[[1]]; rownames(mm) <- gn
      smats[[length(smats) + 1]] <- mm
      metas[[length(metas) + 1]] <-
        data.table(cell = colnames(mm), sample = basename(d))
    }
    common <- Reduce(intersect, lapply(smats, rownames))
    cat("      common genes:", length(common), "\n")
    counts <- do.call(cbind, lapply(smats, function(x) x[common, , drop = FALSE]))
    meta <- rbindlist(metas)
    return(list(counts = counts, meta = meta))
  }
  NULL
}

SELF_TEST <- FALSE
sc <- load_sc()
if (is.null(sc)) {
  cat("[data] nothing local -> SELF-TEST synthetic single-cell data.\n")
  SELF_TEST <- TRUE
  set.seed(SEED)
  cts <- rep(c("Malignant", "T_cell", "Myeloid", "Fibroblast", "Endothelial"),
             each = 120)
  cells <- paste0("SYNC", 1:length(cts))
  farg_all <- unique(unlist(mods))
  markers <- c("ALB", "AFP", "GPC3", "EPCAM", "CD3D", "CD3E", "NKG7", "GNLY",
               "LYZ", "CD68", "S100A8", "COL1A1", "DCN", "ACTA2", "PECAM1",
               "VWF", "ENG")
  partners <- paste0("P", 1:8); downstream <- paste0("D", 1:5)
  genes <- unique(c(farg_all, markers, partners, downstream, paste0("BG", 1:150)))
  lambda <- matrix(0.5, length(genes), length(cts),
                   dimnames = list(genes, cells))
  # cell-type-specific expression
  type_up <- list(Malignant = c("ALB", "AFP", "GPC3", "EPCAM",
                                "TFRC", "HMOX1", "ALOX15", "LPCAT3", "IL6"),
                  T_cell = c("CD3D", "CD3E", "NKG7", "GNLY", "IFNG", "CXCL10"),
                  Myeloid = c("LYZ", "CD68", "S100A8", "S100A9", "NLRP3"),
                  Fibroblast = c("COL1A1", "DCN", "ACTA2", "LOX"),
                  Endothelial = c("PECAM1", "VWF", "ENG", "EDN1"))
  for (ct in names(type_up)) {
    type_up[[ct]] <- intersect(type_up[[ct]], rownames(lambda))
    idx <- which(cts == ct)
    lambda[type_up[[ct]], idx] <- lambda[type_up[[ct]], idx] + 2
  }
  # co-expression module around the VKO target in malignant cells
  mal_idx <- which(cts == "Malignant")
  latent <- rnorm(length(mal_idx))
  for (g in c(VKO_TARGET, partners, downstream))
    lambda[g, mal_idx] <- lambda[g, mal_idx] + 0.9 * latent
  lambda <- pmax(lambda, 0.05)
  counts <- matrix(rpois(length(lambda), lambda), nrow = nrow(lambda),
                   dimnames = dimnames(lambda))
  meta <- data.table(cell = cells, celltype = cts, sample = "SYN")
  sc <- list(counts = counts, meta = meta)
}

counts <- sc$counts; meta <- sc$meta
cat("Loaded:", nrow(counts), "genes x", ncol(counts), "cells\n")

## ---- celltype column ----------------------------------------------------------------
ct_col <- intersect(c("celltype", "CellType", "cell_type", "annotation",
                      "cluster", "major_type"), names(meta))[1]
if (is.na(ct_col)) {
  cat("[annot] no celltype column -> Seurat clustering + annotation ...\n")
  if (requireNamespace("Seurat", quietly = TRUE)) {
    suppressPackageStartupMessages(library(Seurat))
    use_cells <- colnames(counts)
    if (length(use_cells) > MAX_CELLS) {
      set.seed(SEED); use_cells <- sample(use_cells, MAX_CELLS)
    }
    obj <- CreateSeuratObject(counts[, use_cells])
    obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 2000) |>
      ScaleData() |> RunPCA() |> FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = 0.5) |> RunUMAP(dims = 1:20)
    mk <- list(Malignant = c("ALB", "AFP", "GPC3", "EPCAM"),
               T_NK = c("CD3D", "CD3E", "NKG7", "GNLY"),
               Myeloid = c("LYZ", "CD68", "CD14", "S100A8"),
               B_Plasma = c("CD79A", "MS4A1", "SDC1", "MZB1"),
               Fibroblast = c("COL1A1", "DCN", "ACTA2", "FAP"),
               Endothelial = c("PECAM1", "VWF", "ENG", "CDH5"))
    for (nm in names(mk)) {
      g <- intersect(mk[[nm]], rownames(obj))
      if (length(g)) obj <- AddModuleScore(obj, features = list(g), name = nm)
    }
    sc_cols <- grep("^Malignant|^T_NK|^Myeloid|^B_Plasma|^Fibroblast|^Endothelial",
                    colnames(obj@meta.data), value = TRUE)
    if (length(sc_cols) > 1) {
      assign <- apply(obj@meta.data[, sc_cols], 1,
                      function(r) sub("1$", "", sc_cols[which.max(r)])[1])
      obj$celltype_auto <- assign
    } else {
      obj$celltype_auto <- as.character(Idents(obj))
    }
    p_umap <- file.path(OUT_DIR, "umap_annotation.png")
    png(p_umap, width = 1300, height = 1100, res = 150)
    print(DimPlot(obj, group.by = "celltype_auto", label = TRUE) +
            ggtitle("GSE149614: marker-based annotation"))
    dev.off()
    meta <- data.table(cell = colnames(obj),
                       celltype = obj$celltype_auto,
                       sample = "GSE149614")
    rm(obj); gc()
    ct_col <- "celltype"
  } else {
    stop("no celltype column and Seurat not installed")
  }
}
meta$celltype <- meta[[ct_col]]
meta <- meta[cell %in% colnames(counts)]
cat("Cell types:", paste(names(table(meta$celltype)),
                         table(meta$celltype), collapse = " | "), "\n")

## ---- 1. proportions ------------------------------------------------------------------
cat("\n[1/5] Cell-type proportions ...\n")
prop <- meta[, .N, by = celltype][order(-N)]
p1 <- file.path(OUT_DIR, "celltype_proportions.png")
png(p1, width = 1200, height = 800, res = 150)
par(mar = c(12, 5, 3, 2))
barplot(prop$N / sum(prop$N), names.arg = prop$celltype, las = 2,
        col = "#4DBBD5", ylab = "Fraction", main = "Cell-type proportions")
dev.off()
cat("      saved celltype_proportions.png\n")

## ---- 2. marker QC ----------------------------------------------------------------------
cat("\n[2/5] Marker QC (top wilcox markers per cell type) ...\n")
cpm <- t(t(counts) / pmax(colSums(counts), 1)) * 1e4
logn <- log1p(cpm)
mk_tab <- rbindlist(lapply(unique(meta$celltype), function(ct) {
  in_ct <- meta$celltype == ct
  if (sum(in_ct) < 3 || sum(!in_ct) < 3) return(NULL)
  g <- intersect(c(farg_all, rownames(counts)[1:min(600, nrow(counts))]),
                 rownames(logn))
  pv <- apply(logn[g, , drop = FALSE], 1, function(x)
    suppressWarnings(wilcox.test(x[in_ct], x[!in_ct])$p.value))
  fc <- apply(logn[g, , drop = FALSE], 1, function(x)
    mean(x[in_ct]) - mean(x[!in_ct]))
  tab <- data.table(celltype = ct, gene = g, p = pv, mean_diff = fc)
  tab[order(p)][1:min(8, nrow(tab))]
}))
fwrite(mk_tab, file.path(OUT_DIR, "celltype_markers.csv"))
cat("      saved celltype_markers.csv (check that markers match the annotation)\n")

## ---- 3. signature gene dotplot --------------------------------------------------------
cat("\n[3/5] Signature genes across cell types ...\n")
ct_avg <- lapply(unique(meta$celltype), function(ct) {
  rowMeans(logn[sig_genes, meta$celltype == ct, drop = FALSE])
})
ct_avg <- do.call(cbind, ct_avg); colnames(ct_avg) <- unique(meta$celltype)
ct_pct <- lapply(unique(meta$celltype), function(ct) {
  rowMeans(logn[sig_genes, meta$celltype == ct, drop = FALSE] > 0) * 100
})
ct_pct <- do.call(cbind, ct_pct); colnames(ct_pct) <- unique(meta$celltype)
z <- t(scale(t(ct_avg)))
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  ht1 <- Heatmap(z, name = "mean (z)", row_names_gp = gpar(fontsize = 10),
                 column_names_rot = 45, border = TRUE)
  ht2 <- Heatmap(ct_pct, name = "pct %", col = colorRamp2(c(0, 50, 100),
                 c("white", "#FDDBC7", "#B2182B")),
                 row_names_gp = gpar(fontsize = 10), column_names_rot = 45,
                 border = TRUE)
  png(file.path(OUT_DIR, "dotplot_signature_genes.png"),
      width = 1800, height = 900, res = 150)
  draw(ht1 + ht2)
  dev.off()
  cat("      saved dotplot_signature_genes.png\n")
}
fwrite(as.data.table(z, keep.rownames = "gene"),
       file.path(OUT_DIR, "signature_celltype_mean.csv"))

## ---- 4. module scores --------------------------------------------------------------------
cat("\n[4/5] Ferro-aging module scores per cell type ...\n")
gs_ok <- lapply(mods, function(g) intersect(g, rownames(logn)))
gs_ok <- gs_ok[lengths(gs_ok) >= 3]
res <- tryCatch(
  gsva(as.matrix(logn), gs_ok, method = "ssgsea", ssgsea.norm = TRUE,
       verbose = FALSE),
  error = function(e) gsva(ssgseaParam(as.matrix(logn), gs_ok), verbose = FALSE))
mod_score <- t(res)
# downsize for plotting if huge
if (nrow(mod_score) > 50000) {
  set.seed(SEED)
  keep <- sample(seq_len(nrow(mod_score)), 50000)
  mod_score <- mod_score[keep, , drop = FALSE]
  meta <- meta[keep]
}
ct <- meta$celltype
mod_mean <- sapply(unique(ct), function(x) colMeans(mod_score[ct == x, , drop = FALSE]))
mod_mean <- t(mod_mean)
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  png(file.path(OUT_DIR, "module_score_heatmap.png"), width = 1600, height = 800,
      res = 150)
  draw(Heatmap(t(scale(t(mod_mean))), name = "module score (z)",
               row_names_gp = gpar(fontsize = 9), column_names_rot = 45,
               border = TRUE))
  dev.off()
  cat("      saved module_score_heatmap.png\n")
}
# malignant cells: hepatocyte-like cells; when site info exists, restrict to
# TUMOR sites (metadata contains matched adjacent-tissue cells too)
mal <- grepl(MALIG_PAT, ct, ignore.case = TRUE)
if ("site" %in% names(meta)) {
  site_tumor <- grepl("tumor", meta$site, ignore.case = TRUE)
  cat("      site info used: tumor-site cells =", sum(site_tumor), "\n")
  mal <- mal & site_tumor
}
cat("      malignant (tumor hepatocyte) cells:", sum(mal), "\n")
mod_wt <- lapply(colnames(mod_score), function(m)
  data.table(module = m,
             p = suppressWarnings(wilcox.test(mod_score[mal, m],
                                              mod_score[!mal, m])$p.value),
             mean_mal = mean(mod_score[mal, m]),
             mean_other = mean(mod_score[!mal, m])))
mod_wt <- rbindlist(mod_wt)
fwrite(mod_wt, file.path(OUT_DIR, "module_malignant_vs_other.csv"))
cat("      malignant vs other (wilcox):\n"); print(mod_wt[order(p)])
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  pd <- data.table(celltype = ct, score = mod_score[, "all_FARG"])
  g <- ggplot(pd, aes(reorder(celltype, score, FUN = median), score,
                      fill = celltype)) +
    geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.size = 0.2) +
    labs(x = NULL, y = "all-FARG ssGSEA score",
         title = "Ferro-aging score per cell type") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(OUT_DIR, "module_score_violin.png"), g, width = 10, height = 6,
         dpi = 200)
  cat("      saved module_score_violin.png\n")
}

## ---- 5. virtual knockout ---------------------------------------------------------------
cat("\n[5/5] Virtual knockout of", VKO_TARGET, "(scTenifoldKnk) ...\n")
if (!VKO_TARGET %in% rownames(counts)) {
  cat("      target gene not in the matrix; skipped.\n")
} else {
  vko_mal <- grepl(MALIG_PAT, meta$celltype, ignore.case = TRUE)
  if ("site" %in% names(meta))
    vko_mal <- vko_mal & grepl("tumor", meta$site, ignore.case = TRUE)
  use_cells <- meta$cell[vko_mal]
  if (length(use_cells) > VKO_MAX_CELLS) {
    set.seed(SEED); use_cells <- sample(use_cells, VKO_MAX_CELLS)
  }
  cat("      malignant cells used:", length(use_cells), "\n")
  if (length(use_cells) >= 50) {
    # Panel: union of all FARGs and signature genes plus genes ranked by mean/SD.
    # This is the implemented selection rule, not a high-variance gene filter.
    # Module labels do not enter the union used for the VKO panel.
    m <- logn[, use_cells, drop = FALSE]
    v <- apply(m, 1, function(x) mean(x) / (sd(x) + 1e-6))
    panel <- unique(c(unlist(mods), sig_genes,
                      names(sort(v, decreasing = TRUE))[1:VKO_N_GENES]))
    panel <- intersect(panel, rownames(counts))
    cnt <- counts[panel, use_cells, drop = FALSE]
    cat("      network panel:", nrow(cnt), "genes x", ncol(cnt), "cells\n")
    vko_dir <- file.path(OUT_DIR, "vko")
    dir.create(vko_dir, showWarnings = FALSE, recursive = TRUE)
    oldwd <- getwd(); setwd(vko_dir)
    ok <- tryCatch({
      suppressPackageStartupMessages(library(scTenifoldKnk))
      knk <- scTenifoldKnk(countMatrix = cnt, gKO = VKO_TARGET,
                           qc = FALSE, nCores = 1)
      dr <- as.data.table(knk$diffRegulation)   # already computed inside
      list(knk = knk, dr = dr)
    }, error = function(e) {
      cat("      scTenifoldKnk failed:", conditionMessage(e), "\n"); NULL
    })
    setwd(oldwd)
    if (!is.null(ok)) {
      dr <- ok$dr[order(p.value)]
      fwrite(dr[1:min(50, nrow(dr))],
             file.path(OUT_DIR, "vko_dRegulation_top50.csv"))
      fwrite(dr, file.path(OUT_DIR, "vko_dRegulation_full.csv"))  # 完整 887 基因面板
      cat("      top differentially regulated genes after KO:\n")
      print(head(dr[, .(gene, distance, p.value, p.adj)], 10))
      # deterministic bar chart of top perturbed genes (plotKO is interactive
      # plotly/HTML and does not write to a png device, leaving a 0-byte file)
      topg <- dr[1:min(20, nrow(dr))]
      topg[, gene := factor(gene, levels = rev(gene))]
      p_ko <- file.path(OUT_DIR, "vko_plot.png")
      png(p_ko, width = 1400, height = 900, res = 150)
      par(mar = c(5, 10, 3, 2))
      barplot(rev(setNames(topg$distance, topg$gene)), horiz = TRUE, las = 1,
              col = "#B2182B", border = "white",
              xlab = "differential regulation distance",
              main = paste0(VKO_TARGET, " virtual knockout: top perturbed genes"))
      dev.off()
      cat("      saved vko_plot.png (bar chart of top perturbed genes)\n")
      # Exploratory GO: when fewer than five genes pass FDR < 0.05, use the
      # top 20 ranked genes, including non-significant genes and the KO target.
      # The historical run used the default annotated-gene universe (no universe
      # argument below), not the tested VKO panel; do not claim a mechanism.
      topg <- dr[p.adj < 0.05]$gene
      if (length(topg) < 5) topg <- head(dr$gene, 20)
      topg <- head(unique(topg), 50)
      if (length(topg) >= 5 &&
          requireNamespace("clusterProfiler", quietly = TRUE)) {
        suppressPackageStartupMessages(library(clusterProfiler))
        suppressPackageStartupMessages(library(org.Hs.eg.db))
        ego <- tryCatch(enrichGO(gene = topg, OrgDb = org.Hs.eg.db,
                                 keyType = "SYMBOL", ont = "BP",
                                 pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                        error = function(e) NULL)
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          png(file.path(OUT_DIR, "vko_GO.png"), width = 1600, height = 1100,
              res = 150)
          print(barplot(ego, showCategory = 12))
          dev.off()
          fwrite(as.data.table(as.data.frame(ego)),
                 file.path(OUT_DIR, "vko_GO.csv"))
          cat("      saved vko_GO.png / vko_GO.csv\n")
        }
      }
    } else {
      cat("      virtual knockout unavailable; rerun after checking scTenifoldKnk\n")
      cat("      (devtools::install_github('cailab-tamu/scTenifoldKnk'))\n")
    }
  } else {
    cat("      too few malignant cells for network construction; skipped.\n")
  }
}

cat("\n== step08 done ==\n")
if (SELF_TEST)
  cat("SELF-TEST passed. Real run: place GSE149614 files under", DATA_DIR,
      "and rerun.\n")
