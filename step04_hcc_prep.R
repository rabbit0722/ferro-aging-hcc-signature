# =============================================================================
# step04_hcc_prep.R
# -----------------------------------------------------------------------------
# Prepare TCGA-LIHC data from the UCSC Xena Toil pan-cancer matrix.
#   - finds input files in: analysis_out/data/ | data/ | project root | parent
#     (filename matching is tolerant: *_tpm.gz / *_tpm.tsv.gz etc.)
#   - LIHC samples: tumors from the Xena survival table ("cancer type
#     abbreviation" == LIHC); adjacent normals matched by patient prefix
#     (TCGA-XX-YYYY-11 whose patient TCGA-XX-YYYY has a LIHC tumor)
#   - ENSEMBL ids -> gene symbols via org.Hs.eg.db (no network needed)
#   - Toil values are already log2(TPM+0.001): NO further log transform
#   - optional LIHC clinicalMatrix adds stage / grade / age / gender
#
# Output: analysis_out/hcc_data/hcc_data.rds
#   list(expr = log2-scale matrix genes x samples, group = Tumor/Normal,
#        surv = OS/OS.time data.table, pheno = stage/grade/age/gender)
#
# Requirements: data.table, org.Hs.eg.db
# Run:  Rscript step04_hcc_prep.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
})

FIND_DIRS <- c("analysis_out/data", "data", ".", "..", "../data")
find_data <- function(prefix) {
  for (d in FIND_DIRS) {
    if (!dir.exists(d)) next
    fs <- list.files(d, full.names = TRUE)
    hit <- fs[grepl(paste0("^", prefix), basename(fs), ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
EXPR_FILE <- find_data("TcgaTargetGtex_rsem_gene_tpm")
SURV_FILE <- find_data("Survival_SupplementalTable_S1_20171025_xena_sp")
CLIN_FILE <- find_data("TCGA\\.LIHC\\.sampleMap_LIHC_clinicalMatrix")
OUT_DIR   <- "analysis_out/hcc_data"
MIN_LOG2  <- 0.5          # keep genes with mean log2(TPM+0.001) >= this

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("== step04: prepare TCGA-LIHC data ==\n")
if (is.na(EXPR_FILE) || is.na(SURV_FILE)) {
  cat("\nMissing input files. Download them first (put them in the project\n")
  cat("folder path/to/Ferro-aging or in analysis_out/data/):\n")
  cat("  1) expression (1.2 GB): https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.tsv.gz\n")
  cat("  2) survival:            https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz\n")
  cat("  3) clinical (optional): https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix\n")
  stop_quietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  stop_quietly()
}
cat("      expression:", EXPR_FILE, "\n")
cat("      survival  :", SURV_FILE, "\n")

## ---- 1. LIHC sample list -----------------------------------------------------
cat("[1/4] LIHC sample selection ...\n")
surv_all <- fread(SURV_FILE)
surv <- surv_all[`cancer type abbreviation` == "LIHC"]
lihc_patients <- unique(substr(surv$sample, 1, 12))
cat("      LIHC tumors with survival:", nrow(surv), "\n")

hdr <- names(fread(EXPR_FILE, nrows = 0))          # header only (fast)
tum_lihc <- intersect(surv$sample, hdr)
norm_lihc <- hdr[grepl("^TCGA-", hdr) &
                   substr(hdr, 14, 15) == "11" &
                   substr(hdr, 1, 12) %in% lihc_patients]
lihc_samples <- c(tum_lihc, norm_lihc)
cat("      tumors in matrix:", length(tum_lihc),
    "| matched adjacent normals:", length(norm_lihc), "\n")
if (length(lihc_samples) < 50) {
  cat("WARNING: very few LIHC samples found; check the survival file format.\n")
}

## ---- 2. read expression (selected columns only) --------------------------------
cat("[2/4] Reading expression (selected columns, 1-2 min) ...\n")
lihc_samples <- unique(lihc_samples)
# select by column POSITION: the Toil matrix can contain duplicate sample
# names, which breaks name-based selection
sel_idx <- c(1L, match(lihc_samples, hdr))
sel_idx <- sel_idx[!is.na(sel_idx)]
expr_dt <- fread(EXPR_FILE, select = sel_idx)
names(expr_dt) <- c("sample", hdr[sel_idx[-1]])
lihc_samples <- hdr[sel_idx[-1]]          # actual columns read
cat("      matrix:", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n")

## ---- 3. ENSG -> SYMBOL ----------------------------------------------------------
cat("[3/4] ENSEMBL -> SYMBOL via org.Hs.eg.db ...\n")
ensg <- sub("\\..*$", "", expr_dt$sample)          # strip .version
mp <- AnnotationDbi::select(org.Hs.eg.db, keys = ensg, keytype = "ENSEMBL",
                            columns = "SYMBOL")
mp <- as.data.table(mp)[!is.na(SYMBOL) & SYMBOL != ""]
mp <- mp[, .SD[1], by = ENSEMBL]                   # one symbol per ENSG
expr_dt[, ensg := ensg]
expr_dt <- merge(expr_dt, mp[, .(ensg = ENSEMBL, SYMBOL)], by = "ensg",
                 all.x = TRUE)
expr_dt <- expr_dt[!is.na(SYMBOL)]
# collapse multi-ENSG symbols: keep the highest-mean row
m <- as.matrix(expr_dt[, ..lihc_samples])
rownames(m) <- expr_dt$SYMBOL
rm(expr_dt); gc()
ord <- order(rowMeans(m), decreasing = TRUE)
m <- m[ord, , drop = FALSE]
m <- m[!duplicated(rownames(m)), , drop = FALSE]
cat("      genes with symbols:", nrow(m), "\n")

## ---- 4. group labels, filter, survival, pheno ------------------------------------
code <- substr(colnames(m), 14, 15)
group <- ifelse(code == "01", "Tumor", ifelse(code == "11", "Normal", NA))
keep_s <- !is.na(group)
m <- m[, keep_s, drop = FALSE]; group <- group[keep_s]
cat("      Tumor:", sum(group == "Tumor"), "| Normal:", sum(group == "Normal"), "\n")
keepg <- rowMeans(m) >= MIN_LOG2
m <- m[keepg, , drop = FALSE]
cat("      genes kept (mean log2 >= ", MIN_LOG2, "):", nrow(m), "\n")
expr <- m                                    # Toil values are ALREADY log2
rm(m); gc()

# keep only primary-tumor samples (the Xena survival table also contains
# adjacent-normal (-11) rows, which must NOT enter survival analyses)
surv_sub <- surv[sample %in% colnames(expr) & substr(sample, 14, 15) == "01",
                 .(sample, OS, OS.time)]
cat("      survival rows (tumor only):", nrow(surv_sub), "\n")
pheno <- data.table(sample = colnames(expr))
if (!is.na(CLIN_FILE) && file.exists(CLIN_FILE)) {
  clin <- fread(CLIN_FILE)
  cols <- names(clin)
  pick <- function(pats) { for (p in pats) if (p %in% cols) return(p); NA_character_ }
  stage_col <- pick(c("ajcc_pathologic_tumor_stage", "pathologic_stage", "tumor_stage"))
  grade_col <- pick(c("neoplasm_histologic_grade", "tumor_grade"))
  age_col   <- pick(c("age_at_initial_pathologic_diagnosis", "age_at_diagnosis"))
  gender_col<- pick("gender")
  sub <- data.table(sample = clin[[1]])
  if (!is.na(stage_col))  sub$stage  <- clin[[stage_col]]
  if (!is.na(grade_col))  sub$grade  <- clin[[grade_col]]
  if (!is.na(age_col))    sub$age    <- clin[[age_col]]
  if (!is.na(gender_col)) sub$gender <- clin[[gender_col]]
  pheno <- merge(pheno, sub, by = "sample", all.x = TRUE)
  cat("      clinical columns used:",
      paste(c(stage_col, grade_col, age_col, gender_col)[!is.na(c(stage_col, grade_col, age_col, gender_col))],
            collapse = ", "), "\n")
} else {
  cat("      clinicalMatrix not found (optional); stage/grade will be NA.\n")
}

saveRDS(list(expr = expr, group = group, surv = surv_sub, pheno = pheno),
        file.path(OUT_DIR, "hcc_data.rds"))
cat("\nSaved:", file.path(OUT_DIR, "hcc_data.rds"), "\n")
cat("expr:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("== step04 done ==\n")
