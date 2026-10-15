# =============================================================================
# step07_hcc_immune_drug.R
# -----------------------------------------------------------------------------
# Immune landscape + drug sensitivity + mutation for the HCC FARGs risk groups.
#
# Modules (each runs when its inputs are available; otherwise prints how-to):
#   1) TME ssGSEA         28 immune-cell signatures (IOBR signature_tme,
#                         Charoentong et al. 2017) scored with GSVA -> High vs
#                         Low wilcox test + boxplot grid
#   2) Immune checkpoints expression of checkpoints (PD-L1, PD-1, CTLA4, ...)
#                         -> wilcox + boxplot; export files for the TIDE web tool
#   3) Drug sensitivity   oncoPredict::calcPhenotype with GDSC2 training data
#                         -> per-drug IC50 High vs Low
#   4) Mutation / TMB     GDC LIHC mutect MAF -> TMB per sample + top genes
#                         (oncoplot if maftools installed)
#
# Input : analysis_out/hcc_data/hcc_data.rds        (from step04)
#         analysis_out/hcc_out/risk_training.csv    (from step05)
# Optional: analysis_out/data/GDSC2_Expr.RData, GDSC2_Res.RData (oncoPredict)
#           analysis_out/data/TCGA.LIHC.mutect2_snv.tsv.gz       (Xena GDC hub)
# If hcc_data.rds is missing the script runs a SELF-TEST on synthetic data.
#
# Output: analysis_out/hcc_out/immune/
#   tme_scores_per_sample.csv, tme_wilcox.csv, tme_boxplot.png, checkpoint.csv, checkpoint_boxplot.png
#   tide_upload/ (expr + sample groups for https://tide.dfci.harvard.edu)
#   drug_ic50.csv, drug_boxplot.png (when GDSC files present)
#   tmb.csv, tmb_boxplot.png (when MAF present)
#
# Requirements: data.table, GSVA, ggplot2; optional IOBR, oncoPredict, maftools
# Run:  Rscript step07_hcc_immune_drug.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
})

DATA_RDS    <- "analysis_out/hcc_data/hcc_data.rds"
RISK_FILE   <- "analysis_out/hcc_out/risk_training.csv"
OUT_DIR     <- "analysis_out/hcc_out/immune"
FIND_DIRS   <- c("analysis_out/data", "data", ".", "..", "../data")
find_data <- function(prefix) {
  if (length(prefix) > 1) {                       # try several prefixes in order
    for (p in prefix) {
      hit <- find_data(p)
      if (!is.na(hit)) return(hit)
    }
    return(NA_character_)
  }
  for (d in FIND_DIRS) {
    if (!dir.exists(d)) next
    fs <- list.files(d, recursive = TRUE, full.names = TRUE)
    hit <- fs[grepl(paste0("^", prefix), basename(fs), ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
GDSC_EXPR   <- find_data("GDSC2_Expr")
GDSC_RES    <- find_data("GDSC2_Res")
# MAF: prefer LIHC-specific files; the incomplete mc3 (if present) comes last
MAF_FILE    <- find_data(c("TCGA.LIHC.mutect", "TCGA.LIHC.mutect2_snv",
                           "mc3.v0.2.8.PUBLIC.maf.gz", "mc3"))
SEED        <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("== step07: immune / drug / mutation analysis of FARGs risk groups ==\n")

## ---- load data or self-test ------------------------------------------------------
SELF_TEST <- FALSE
if (!file.exists(DATA_RDS)) {
  cat("[WARN] hcc_data.rds not found -> SELF-TEST on synthetic data.\n\n")
  SELF_TEST <- TRUE
  set.seed(SEED)
  n <- 200
  markers <- c("CD163", "MRC1", "CD68", "S100A8", "S100A9", "COL1A1", "ACTA2",
               "FAP", "CD8A", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "KLRD1",
               "CD4", "FOXP3", "IL2RA", "CD19", "MS4A1", "PECAM1", "VWF",
               "CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1",
               "CD86", "IL1B", "TNF", "CXCL10")
  genes <- unique(c(paste0("FARG", 1:95), paste0("BG", 1:250), markers))
  expr <- matrix(rnorm(length(genes) * n, 6, 1.2), nrow = length(genes),
                 dimnames = list(genes, paste0("SYN", 1:n)))
  # immune signal: "High" group enriched for myeloid/CAF markers, depleted CD8
  hi <- 101:n
  for (g in c("CD163", "MRC1", "CD68", "S100A8", "S100A9", "COL1A1", "ACTA2",
              "FAP", "CD274", "IL1B", "TNF", "CXCL10"))
    expr[g, hi] <- expr[g, hi] + 1.5
  for (g in c("CD8A", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "KLRD1"))
    expr[g, hi] <- expr[g, hi] - 1.2
  risk <- rnorm(n); risk[hi] <- risk[hi] + 0.8
  group <- ifelse(rank(risk) > n / 2, "High", "Low")
  expr_all_tumor <- expr                     # self-test: same matrix
} else {
  L <- readRDS(DATA_RDS)
  expr_all_tumor <- L$expr[, L$group == "Tumor"]   # ALL tumors (for drug module)
  rk <- fread(RISK_FILE)
  group <- rk$risk_group[match(colnames(expr_all_tumor), rk$sample)]
  keep <- !is.na(group)
  expr <- expr_all_tumor[, keep]                   # risk-annotated subset
  group <- group[keep]
  cat("Real data loaded:", nrow(expr), "genes x", ncol(expr),
      "risk-annotated tumor samples (all tumors:",
      ncol(expr_all_tumor), ")\n")
}
cat("Groups:", paste(names(table(group)), table(group), collapse = " | "), "\n")

## ---- 1. TME ssGSEA ----------------------------------------------------------------
cat("\n[1/5] TME immune-cell ssGSEA ...\n")
sig_tme <- NULL
if (requireNamespace("IOBR", quietly = TRUE)) {
  sig_tme <- tryCatch(IOBR:::signature_tme, error = function(e) NULL)
}
if (is.null(sig_tme)) {
  # compact fallback marker sets (literature-adapted)
  sig_tme <- list(
    CD8_T = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"),
    CD4_T = c("CD4", "IL7R", "CD40LG", "CD28", "ICOS"),
    Treg  = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "CCR8"),
    NK    = c("NKG7", "KLRD1", "KLRK1", "NCR1", "GNLY"),
    Bcell = c("CD19", "MS4A1", "CD79A", "CD79B", "BLK"),
    Macrophage = c("CD68", "CD163", "CSF1R", "ITGAM"),
    M1    = c("IL1B", "TNF", "CXCL10", "IRF5", "CD86"),
    M2    = c("MRC1", "CD163", "MSR1", "CCL18", "IL10"),
    DC    = c("ITGAX", "CD1C", "CLEC9A", "LAMP3", "BATF3"),
    Mast  = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
    Neutrophil = c("FCGR3B", "S100A8", "S100A9", "CSF3R", "CXCR2"),
    CAF   = c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "FAP", "PDGFRB"),
    Endothelial = c("PECAM1", "VWF", "CDH5", "ENG", "KDR"))
}
gs_ok <- lapply(sig_tme, function(g) intersect(g, rownames(expr)))
gs_ok <- gs_ok[lengths(gs_ok) >= 3]
cat("      scoring", length(gs_ok), "cell types with ssGSEA ...\n")
res <- tryCatch(
  gsva(as.matrix(expr), gs_ok, method = "ssgsea", ssgsea.norm = TRUE,
       verbose = FALSE),
  error = function(e) gsva(ssgseaParam(as.matrix(expr), gs_ok),
                            verbose = FALSE))
scores <- t(res)
tme <- data.table(sample = colnames(expr), group = group)
for (cc in rownames(res)) tme[[cc]] <- scores[, cc]
# 逐样本 ssGSEA 打分落盘（供独立核验 Fig4A 箱线图）
fwrite(tme, file.path(OUT_DIR, "tme_scores_per_sample.csv"))

tme_long <- melt(tme, id.vars = c("sample", "group"),
                 variable.name = "cell", value.name = "score")
wt <- tme_long[, .(p = wilcox.test(score ~ group)$p.value), by = cell]
wt[, FDR := p.adjust(p, "BH")]
wt <- wt[order(p)]
fwrite(wt, file.path(OUT_DIR, "tme_wilcox.csv"))
cat("      top differential cell types:\n")
print(wt[1:min(8, nrow(wt))])

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  g <- ggplot(tme_long, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6) +
    facet_wrap(~ cell, scales = "free") +
    scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
    labs(x = NULL, y = "ssGSEA score",
         title = "Immune cell infiltration by FARGs risk group") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  ggsave(file.path(OUT_DIR, "tme_boxplot.png"), g, width = 12, height = 9,
         dpi = 200)
  cat("      saved tme_boxplot.png\n")
}

## ---- 2. checkpoints + TIDE export -------------------------------------------------
cat("\n[2/5] Immune checkpoints ...\n")
cps <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1",
         "SIGLEC15", "VSIR", "TNFRSF9", "CD80", "CD86")
cps <- intersect(cps, rownames(expr))
cp_dt <- data.table(sample = colnames(expr), group = group)
for (g in cps) cp_dt[[g]] <- as.numeric(expr[g, ])
cp_long <- melt(cp_dt, id.vars = c("sample", "group"),
                variable.name = "gene", value.name = "expr")
cp_wt <- cp_long[, .(p = wilcox.test(expr ~ group)$p.value), by = gene]
cp_wt[, FDR := p.adjust(p, "BH")]
fwrite(cp_wt, file.path(OUT_DIR, "checkpoint.csv"))
cat("      checkpoint differences:\n"); print(cp_wt[order(p)])
if (requireNamespace("ggplot2", quietly = TRUE)) {
  g2 <- ggplot(cp_long, aes(group, expr, fill = group)) +
    geom_boxplot(outlier.size = 0.3) +
    facet_wrap(~ gene, scales = "free") +
    scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
    labs(x = NULL, y = "log2 expression",
         title = "Immune checkpoint expression by risk group") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(OUT_DIR, "checkpoint_boxplot.png"), g2, width = 10,
         height = 7, dpi = 200)
  cat("      saved checkpoint_boxplot.png\n")
}
# TIDE web upload files
tide_dir <- file.path(OUT_DIR, "tide_upload")
dir.create(tide_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(as.data.table(expr, keep.rownames = "gene"),
       file.path(tide_dir, "expr_tumor.txt"), sep = "\t")
fwrite(data.table(sample = colnames(expr), group = group),
       file.path(tide_dir, "sample_groups.csv"))
writeLines(c(
  "TIDE web tool upload guide",
  "1. Open https://tide.dfci.harvard.edu",
  "2. Upload expr_tumor.txt as the expression matrix (genes x samples)",
  "3. Optional: upload sample_groups.csv to color by risk group",
  "4. Run TIDE; download the prediction table (Dysfunction / Exclusion / TIDE / responder)",
  "5. Compare responder fraction between High and Low groups and add to Figure 5"),
  file.path(tide_dir, "README_tide.txt"))
cat("      TIDE upload files written to", tide_dir, "\n")

## ---- 2b. TIDE in R via tidepy (Python CLI, optional) --------------------------------
# 安装:  pip install tidepy   （清华镜像: pip install -i https://pypi.tuna.tsinghua.edu.cn/simple tidepy）
# 本机实测路径（Python 3.14）；换机器时改这两行
TIDEPY_PYTHON <- "path/to/python.exe"   # 请改为本机 Python 路径
TIDEPY_SCRIPT <- "path/to/tidepy"       # 请改为本机 tidepy 脚本路径
cat("\n[2b/5] TIDE (tidepy CLI) ...\n")
if (file.exists(TIDEPY_PYTHON) && file.exists(TIDEPY_SCRIPT)) {
  tide_out <- file.path(OUT_DIR, "tide_result")
  run_ok <- tryCatch({
    system2(TIDEPY_PYTHON,
            c(TIDEPY_SCRIPT, "-o", tide_out, "-c", "Other", "--ignore_norm",
              file.path(tide_dir, "expr_tumor.txt")),
            stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) {
    cat("      tidepy run failed:", conditionMessage(e), "\n"); FALSE
  })
  if (run_ok) {
    tide_files <- list.files(OUT_DIR, pattern = "^tide_result", full.names = TRUE)
    if (length(tide_files)) {
      td <- fread(tide_files[1])
      cat("      tidepy output:", basename(tide_files[1]),
          "(", nrow(td), "rows )\n      columns:",
          paste(names(td), collapse = ", "), "\n")
      # defensive column detection
      sample_col <- names(td)[1]
      td[[sample_col]] <- as.character(td[[sample_col]])
      td$group <- group[match(td[[sample_col]], colnames(expr))]
      td <- td[!is.na(group)]
      resp_col <- grep("responder", names(td), value = TRUE, ignore.case = TRUE)[1]
      tide_col <- grep("^tide", names(td), value = TRUE, ignore.case = TRUE)[1]
      dysf_col <- grep("dysfunction", names(td), value = TRUE, ignore.case = TRUE)[1]
      excl_col <- grep("exclusion", names(td), value = TRUE, ignore.case = TRUE)[1]
      fwrite(td, file.path(OUT_DIR, "tide_predictions_annotated.csv"))
      if (!is.na(resp_col)) {
        tab <- table(td$group, grepl("true", td[[resp_col]], ignore.case = TRUE))
        ft <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA)
        cat("      Responder fraction by group:\n")
        print(prop.table(tab, 1))
        cat("      fisher p =", ft, "\n")
      }
      for (sc in c(tide_col, dysf_col, excl_col)) {
        if (!is.na(sc)) {
          pv <- tryCatch(wilcox.test(td[[sc]] ~ td$group)$p.value,
                         error = function(e) NA)
          cat("      ", sc, "High vs Low wilcox p =", pv, "\n")
        }
      }
      if (requireNamespace("ggplot2", quietly = TRUE) &&
          !is.na(tide_col)) {
        suppressPackageStartupMessages(library(ggplot2))
        g5 <- ggplot(td, aes(group, .data[[tide_col]], fill = group)) +
          geom_boxplot(outlier.size = 0.3) +
          scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
          labs(x = NULL, y = "TIDE score",
               title = "TIDE score by FARGs risk group") +
          theme_minimal(base_size = 12) + theme(legend.position = "none")
        ggsave(file.path(OUT_DIR, "tide_boxplot.png"), g5,
               width = 6, height = 5, dpi = 200)
        cat("      saved tide_boxplot.png\n")
      }
    } else {
      cat("      tidepy finished but no tide_result* file found; check manually.\n")
    }
  }
} else {
  cat("      tidepy not found at:\n      ", TIDEPY_PYTHON, "\n      ",
      TIDEPY_SCRIPT, "\n", sep = "")
  cat("      Install: pip install tidepy, then update TIDEPY_PYTHON/TIDEPY_SCRIPT\n")
  cat("      at the top of module [2b], or use the web tool (files in",
      tide_dir, ").\n")
}

## ---- 3. drug sensitivity -------------------------------------------------------------
cat("\n[3/5] Drug sensitivity (oncoPredict / GDSC2) ...\n")
if (!is.na(GDSC_EXPR) && !is.na(GDSC_RES) &&
    file.exists(GDSC_EXPR) && file.exists(GDSC_RES)) {
  if (requireNamespace("oncoPredict", quietly = TRUE)) {
    suppressPackageStartupMessages(library(oncoPredict))
    TrainingData  <- readRDS(GDSC_EXPR)     # .rds (genes x cell lines)
    TrainingPtype <- readRDS(GDSC_RES)      # .rds (cell lines x drugs)
    cat("      GDSC2 training:", nrow(TrainingData), "genes x",
        ncol(TrainingData), "cell lines | drugs:", ncol(TrainingPtype), "\n")
    # keep only HCC-relevant drugs (sorafenib/lenvatinib/cytotoxics ...):
    # cutting 198 -> ~25 drugs makes calcPhenotype ~8x faster
    drug_pat <- c("Sorafenib", "Lenvatinib", "Regorafenib", "Cabozantinib",
                  "Cisplatin", "Doxorubicin", "Gemcitabine", "5-Fluorouracil",
                  "Oxaliplatin", "Paclitaxel", "Etoposide", "Irinotecan",
                  "Mitomycin", "Vinblastine", "Vincristine", "Bleomycin",
                  "Camptothecin", "Cytarabine", "Methotrexate",
                  "Epirubicin", "Afatinib", "Erlotinib", "Gefitinib")
    keep_drugs <- unique(unlist(lapply(drug_pat, function(p)
      grep(paste0("^", p), colnames(TrainingPtype), value = TRUE,
           ignore.case = TRUE))))
    if (length(keep_drugs) < 5) keep_drugs <- colnames(TrainingPtype)
    TrainingPtype <- TrainingPtype[, keep_drugs, drop = FALSE]
    cat("      drugs kept for prediction:", length(keep_drugs), "\n")
    # calcPhenotype bug workaround: use ALL tumor samples (369) for the
    # prediction (a specific count, e.g. 363, can trigger an internal
    # result-assembly bug in oncoPredict 1.3.1), subset results afterwards
    test_expr <- as.matrix(expr_all_tumor)
    # keep only genes shared with GDSC training
    common <- intersect(rownames(TrainingData), rownames(test_expr))
    cat("      genes shared with GDSC2:", length(common), "\n")
    cat("      predicting IC50 (may take several minutes) ...\n")
    # oncoPredict 1.3.x notes:
    #  - args: powerTransformPhenotype + removeLowVaringGenesFrom (Varing!)
    #  - 1.3.1 has a result-assembly bug when drug count x sample count is
    #    large (16 drugs x 369 samples failed; 8+8 and 200x16 both worked),
    #    so we predict in TWO drug batches and cbind the results
    fit_drugs <- function(drgs) {
      tryCatch(
        calcPhenotype(trainingExprData = TrainingData[common, ],
                      trainingPtype = TrainingPtype[, drgs, drop = FALSE],
                      testExprData = test_expr[common, ],
                      batchCorrect = "eb", powerTransformPhenotype = TRUE,
                      removeLowVaryingGenes = 0.2,
                      removeLowVaringGenesFrom = "homogenizeData",
                      minNumSamples = 10, printOutput = FALSE,
                      parallel = FALSE),
        error = function(e) {
          cat("      batch failed (", conditionMessage(e),
              "); retrying single-drug loop ...\n", sep = "")
          res_list <- lapply(drgs, function(dg) {
            tryCatch(calcPhenotype(
              trainingExprData = TrainingData[common, ],
              trainingPtype = TrainingPtype[, dg, drop = FALSE],
              testExprData = test_expr[common, ],
              batchCorrect = "eb", powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = "homogenizeData",
              minNumSamples = 10, printOutput = FALSE,
              parallel = FALSE), error = function(e2) NULL)
          })
          res_list <- res_list[!vapply(res_list, is.null, logical(1))]
          do.call(cbind, res_list)
        })
    }
    drug_n <- length(keep_drugs)
    half_idx <- seq_len(ceiling(drug_n / 2))
    ic50 <- if (drug_n == 1) fit_drugs(keep_drugs) else
      cbind(fit_drugs(keep_drugs[half_idx]),
            fit_drugs(keep_drugs[-half_idx]))
    # orientation: prefer samples x drugs; oncoPredict may return either way
    if (length(intersect(colnames(expr), rownames(ic50))) <
        length(intersect(colnames(expr), colnames(ic50)))) ic50 <- t(ic50)
    # keep only the risk-annotated samples
    ic50 <- ic50[intersect(colnames(expr), rownames(ic50)), , drop = FALSE]
    cat("      predicted IC50 matrix:", nrow(ic50), "samples x",
        ncol(ic50), "drugs\n")
    # GDSC drug names carry suffixes (e.g. Sorafenib_2203) -> prefix match
    drug_pat <- c("Sorafenib", "Lenvatinib", "Cisplatin", "Doxorubicin",
                  "Gemcitabine", "5-Fluorouracil", "Oxaliplatin", "Paclitaxel")
    drugs_interest <- unique(unlist(lapply(drug_pat, function(p)
      grep(paste0("^", p), colnames(ic50), value = TRUE, ignore.case = TRUE))))
    if (length(drugs_interest) == 0)
      drugs_interest <- colnames(ic50)[1:min(8, ncol(ic50))]
    cat("      HCC-relevant drugs found:", length(drugs_interest), "\n")
    ddt <- data.table(sample = rownames(ic50), group = group)
    for (dr in colnames(ic50)) ddt[[dr]] <- as.numeric(ic50[, dr])
    ddt <- ddt[!is.na(group)]   # drop samples without risk annotation (NA group)
    dl <- melt(ddt, id.vars = c("sample", "group"),
               variable.name = "drug", value.name = "IC50")
    dl <- dl[!is.na(IC50)]
    # FDR across ALL kept drugs (16), not just the 8 reported -> unified family
    dw <- dl[, .(p = wilcox.test(IC50 ~ group)$p.value,
                 median_High = median(IC50[group == "High"]),
                 median_Low = median(IC50[group == "Low"])), by = drug]
    dw[, FDR := p.adjust(p, method = "BH")]
    dw[, direction := fifelse(median_High > median_Low,
                              "High higher IC50 (more resistant)",
                              "High lower IC50 (more sensitive)")]
    setorder(dw, p)
    fwrite(dw, file.path(OUT_DIR, "drug_ic50.csv"))
    # per-sample IC50 matrix + per-group quantiles (independent reproduction)
    fwrite(ddt, file.path(OUT_DIR, "drug_ic50_per_sample.csv"))
    qsum <- dl[, .(n = .N, min = round(min(IC50), 4),
                   Q1 = round(quantile(IC50, 0.25), 4),
                   median = round(median(IC50), 4),
                   Q3 = round(quantile(IC50, 0.75), 4),
                   max = round(max(IC50), 4)), by = .(drug, group)]
    fwrite(qsum[order(drug, group)], file.path(OUT_DIR, "drug_ic50_quantiles.csv"))
    cat("      drug IC50 differences (raw P + BH-FDR across all kept drugs):\n")
    print(dw)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      g3 <- ggplot(dl, aes(group, IC50, fill = group)) +
        geom_boxplot(outlier.size = 0.3) +
        facet_wrap(~ drug, scales = "free") +
        scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
        labs(x = NULL, y = "Predicted IC50", title = "Drug sensitivity by risk group") +
        theme_minimal(base_size = 11) + theme(legend.position = "top")
      ggsave(file.path(OUT_DIR, "drug_boxplot.png"), g3, width = 10, height = 7,
             dpi = 200)
      cat("      saved drug_boxplot.png\n")
    }
  } else {
    cat("      oncoPredict not installed:\n")
    cat("      devtools::install_github('danieldunn/oncoPredict')\n")
  }
} else {
  cat("      GDSC2 files not found. To enable drug sensitivity:\n")
  cat("      1) open https://osf.io/c6tfx/ (oncoPredict official data)\n")
  cat("      2) download 'GDSC2_Expr (RData)' and 'GDSC2_Res (RData)'\n")
  cat("      3) save both to analysis_out/data/ (keep original filenames)\n")
  if (SELF_TEST) cat("      (skipped in SELF-TEST mode)\n")
}

## ---- 4. mutation / TMB ------------------------------------------------------------------
cat("\n[4/5] Mutation / TMB ...\n")
if (!is.na(MAF_FILE) && file.exists(MAF_FILE)) {
  # read only the columns needed (MC3 has 100+ columns; keeps memory low)
  maf_hdr <- names(fread(MAF_FILE, nrows = 0))
  keep_cols <- intersect(c("Tumor_Sample_Barcode", "Hugo_Symbol",
                           "Variant_Classification", "sampleId",
                           "hugoGeneSymbol", "mutationType"), maf_hdr)
  maf <- fread(MAF_FILE, select = keep_cols)
  cat("      MAF rows:", nrow(maf), "\n")
  # accept GDC MAF columns or cBioPortal-style columns
  if (!"Tumor_Sample_Barcode" %in% names(maf) && "sampleId" %in% names(maf))
    setnames(maf, "sampleId", "Tumor_Sample_Barcode")
  if (!"Hugo_Symbol" %in% names(maf) && "hugoGeneSymbol" %in% names(maf))
    setnames(maf, "hugoGeneSymbol", "Hugo_Symbol")
  if (!"Variant_Classification" %in% names(maf) && "mutationType" %in% names(maf))
    setnames(maf, "mutationType", "Variant_Classification")
  need <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  if (all(need %in% names(maf))) {
    maf[, sample := substr(Tumor_Sample_Barcode, 1, 15)]
    non_syn <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
                 "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del",
                 "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation")
    maf_ns <- maf[Variant_Classification %in% non_syn]
    n_match <- sum(unique(maf_ns$sample) %in% colnames(expr))
    cat("      LIHC samples with mutations:", n_match, "/", ncol(expr), "\n")
    if (n_match < 100) {
      cat("      WARNING: MAF covers too few LIHC samples (<100) -> the file is\n")
      cat("      likely incomplete (e.g. a truncated MC3). Recommended: download\n")
      cat("      the LIHC-specific mutations from cBioPortal (~2 MB, no login):\n")
      cat("        curl.exe -L -o \"analysis_out\\data\\TCGA.LIHC.mutect_cbioportal.tsv\" \"https://www.cbioportal.org/api/molecular-profiles/lihc_tcga_mutations/mutations?sampleListId=lihc_tcga_all&projection=DETAILED&pageSize=1000000\"\n")
      cat("      (this script accepts cBioPortal column names automatically)\n")
    }
    tmb <- maf_ns[, .(tmb = .N), by = sample]          # TMB ~ count / 38 Mb exome
    tmb[, tmb_per_mb := tmb / 38]
    tmb <- merge(tmb, data.table(sample = colnames(expr), group = group),
                 by = "sample", all.y = TRUE)
    tmb[is.na(tmb), tmb := 0][is.na(tmb_per_mb), tmb_per_mb := 0]
    fwrite(tmb, file.path(OUT_DIR, "tmb.csv"))
    cat("      TMB (per Mb) by group:\n")
    print(tmb[, .(median_tmb = median(tmb_per_mb)), by = group])
    cat("      wilcox p =", wilcox.test(tmb_per_mb ~ group, data = tmb)$p.value, "\n")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      g4 <- ggplot(tmb, aes(group, tmb_per_mb, fill = group)) +
        geom_boxplot() +
        scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
        labs(x = NULL, y = "TMB (mut/Mb)", title = "Tumor mutation burden") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
      ggsave(file.path(OUT_DIR, "tmb_boxplot.png"), g4, width = 6, height = 5,
             dpi = 200)
      cat("      saved tmb_boxplot.png\n")
    }
    if (requireNamespace("maftools", quietly = TRUE)) {
      suppressPackageStartupMessages(library(maftools))
      # maftools needs ALL MAF fields (Chromosome etc.), so re-read the full
      # file here (the TMB part above used a 3-column subset to save memory)
      maf_full <- fread(MAF_FILE)
      if (!"Tumor_Sample_Barcode" %in% names(maf_full) &&
          "sampleId" %in% names(maf_full))
        setnames(maf_full, "sampleId", "Tumor_Sample_Barcode")
      if (!"Hugo_Symbol" %in% names(maf_full) &&
          "hugoGeneSymbol" %in% names(maf_full))
        setnames(maf_full, "hugoGeneSymbol", "Hugo_Symbol")
      if (!"Variant_Classification" %in% names(maf_full) &&
          "mutationType" %in% names(maf_full))
        setnames(maf_full, "mutationType", "Variant_Classification")
      maf_full[, Tumor_Sample_Barcode := substr(Tumor_Sample_Barcode, 1, 15)]
      grp <- unique(data.table(Tumor_Sample_Barcode = maf_full$Tumor_Sample_Barcode,
                               Group = group[match(maf_full$Tumor_Sample_Barcode,
                                                   colnames(expr))]))
      grp <- grp[!is.na(Group)]
      maf_obj <- tryCatch(
        read.maf(maf = maf_full[Variant_Classification %in% non_syn],
                 clinicalData = grp),
        error = function(e) {
          cat("      read.maf failed:", conditionMessage(e), "\n"); NULL })
      if (!is.null(maf_obj)) {
        p <- file.path(OUT_DIR, "oncoplot.png")
        png(p, width = 2400, height = 1400, res = 150)
        oncoplot(maf = maf_obj, top = 20, clinicalFeatures = "Group")
        dev.off()
        cat("      saved oncoplot.png\n")
      }
    } else {
      cat("      maftools not installed; oncoplot skipped.\n")
      cat("      install: BiocManager::install('maftools')\n")
    }
  }
} else {
  cat("      MAF file not found. To enable TMB/mutation analysis:\n")
  cat("      1) download https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.mutect2_snv.tsv.gz\n")
  cat("      2) save as analysis_out/data/TCGA.LIHC.mutect2_snv.tsv.gz\n")
  if (SELF_TEST) cat("      (skipped in SELF-TEST mode)\n")
}

cat("\n[5/5] == step07 done ==\n")
if (SELF_TEST) {
  cat("SELF-TEST passed: immune ssGSEA + checkpoint modules are ready.\n")
  cat("Drug and TMB modules will activate once GDSC2 / MAF files are placed\n")
  cat("in analysis_out/data/ (see the instructions printed above).\n")
}
