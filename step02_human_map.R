# =============================================================================
# step02_human_map.R
# -----------------------------------------------------------------------------
# Map the 95 monkey ferro-aging genes (ENSMFAG IDs) to human orthologs.
#
# Strategy (with graceful degradation):
#   1) biomaRt : mfascicularis_gene_ensembl -> hsapiens homolog (preferred;
#                gives human ENSG + symbol + orthology type)
#   2) fallback: use the SYMBOL column of the xlsx (already human symbols)
#                + org.Hs.eg.db to add ENTREZID / ENSEMBL
#   3) last resort: SYMBOL only
#
# Output: analysis_out/step02/farg_human_map.csv
#   columns: geneid_monkey, symbol (xlsx), status, ensg_human, entrez_human,
#            symbol_biomart, orthology_type
#
# Requirements: readxl, data.table; optional biomaRt, org.Hs.eg.db
# Run:  Rscript step02_human_map.R
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

XLSX_FILE <- "Ferro-aging-基因列表.xlsx"
OUT_DIR   <- "analysis_out/step02"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# keep biomaRt's on-disk cache inside the workspace (writable & reproducible);
# remember the old value so we can RESTORE it (otherwise this side effect
# would shift msigdbr's cache dir for step09 in a 00_master run)
OLD_R_USER_CACHE <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
Sys.setenv(R_USER_CACHE_DIR = normalizePath(file.path(OUT_DIR, "R_cache"),
                                            mustWork = FALSE))
on.exit({
  if (is.na(OLD_R_USER_CACHE)) Sys.unsetenv("R_USER_CACHE_DIR")
  else Sys.setenv(R_USER_CACHE_DIR = OLD_R_USER_CACHE)
}, add = TRUE)

cat("== step02: monkey -> human ortholog mapping of the 95 FARGs ==\n")

## ---- read FARG list ---------------------------------------------------------
farg <- as.data.table(read_excel(XLSX_FILE, sheet = "Ferroaging gene set",
                                 skip = 1))
setnames(farg, names(farg), c("Geneid", "SYMBOL"))
farg <- farg[!is.na(Geneid) & !is.na(SYMBOL)]
cat("FARGs read:", nrow(farg), "\n")

map <- data.table(geneid_monkey = farg$Geneid,
                  symbol        = farg$SYMBOL,
                  status        = "pending",
                  ensg_human    = NA_character_,
                  entrez_human  = NA_character_,
                  symbol_biomart = NA_character_,
                  orthology_type = NA_character_)

## ---- strategy 1: biomaRt ----------------------------------------------------
biomart_ok <- FALSE
if (requireNamespace("biomaRt", quietly = TRUE)) {
  suppressPackageStartupMessages(library(biomaRt))
  cat("\n[1/3] Trying biomaRt (needs internet; may take ~1-2 min) ...\n")
  ok <- tryCatch({
    mart <- useEnsembl(biomart = "ensembl",
                       dataset = "mfascicularis_gene_ensembl")
    # note: pin reproducibility with version= e.g. useEnsembl(..., version = 113)
    cat("      connected to Ensembl dataset:", mart@dataset,
        "release:", mart@biomart, "\n")
    bm <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name",
                     "hsapiens_homolog_ensembl_gene",
                     "hsapiens_homolog_associated_gene_name",
                     "hsapiens_homolog_orthology_type"),
      filters = "ensembl_gene_id",
      values = map$geneid_monkey,
      mart = mart)
    TRUE
  }, error = function(e) {
    cat("      biomaRt failed:", conditionMessage(e), "\n")
    FALSE
  })
  if (ok && nrow(bm) > 0) {
    bm <- as.data.table(bm)
    map <- merge(map, bm, by.x = "geneid_monkey", by.y = "ensembl_gene_id",
                 all.x = TRUE)
    map[!is.na(hsapiens_homolog_ensembl_gene), `:=`(
      status = "biomart",
      ensg_human = hsapiens_homolog_ensembl_gene,
      symbol_biomart = hsapiens_homolog_associated_gene_name,
      orthology_type = hsapiens_homolog_orthology_type)]
    map[, c("external_gene_name", "hsapiens_homolog_ensembl_gene",
            "hsapiens_homolog_associated_gene_name",
            "hsapiens_homolog_orthology_type") := NULL]
    biomart_ok <- TRUE
    cat("      biomaRt OK:", sum(map$status == "biomart"), "/",
        nrow(map), "mapped\n")
  } else {
    cat("      biomaRt returned no rows; falling back.\n")
  }
} else {
  cat("\n[1/3] biomaRt not installed; skipping.\n")
}

## ---- strategy 2/3: org.Hs.eg.db / symbol only -------------------------------
if (!biomart_ok) {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    cat("\n[2/3] Falling back to org.Hs.eg.db symbol mapping ...\n")
    # SYMBOL -> ENTREZID (1:1 where possible; take first of multi-hits)
    eg <- AnnotationDbi::select(org.Hs.eg.db, keys = map$symbol,
                                keytype = "SYMBOL",
                                columns = c("ENTREZID", "ENSEMBL"))
    eg <- as.data.table(eg)
    eg <- eg[!is.na(ENTREZID)]
    eg <- eg[, .SD[1], by = SYMBOL]                    # keep first hit
    map <- merge(map, eg, by.x = "symbol", by.y = "SYMBOL", all.x = TRUE)
    map[!is.na(ENTREZID), `:=`(status = "orgdb",
                               entrez_human = as.character(ENTREZID),
                               ensg_human = ENSEMBL)]
    map[, c("ENTREZID", "ENSEMBL") := NULL]
    cat("      org.Hs.eg.db OK:", sum(map$status == "orgdb"), "/",
        nrow(map), "mapped\n")
  } else {
    cat("\n[2/3] org.Hs.eg.db not installed; symbols only.\n")
  }
  map[status == "pending", status := "symbol_only"]
}

## ---- QC + output -------------------------------------------------------------
cat("\n[3/3] Mapping summary by status:\n")
print(map[, .N, by = status][order(status)])
cat("\nGenes with no human ID (check manually):\n")
print(map[status == "symbol_only", .(geneid_monkey, symbol)])

# symbol used downstream = biomart symbol when available, else xlsx symbol
map[, symbol_human := fifelse(!is.na(symbol_biomart), symbol_biomart, symbol)]
n_uniq <- uniqueN(map$symbol_human[!is.na(map$symbol_human)])
cat("\nUnique human symbols for downstream use:", n_uniq, "\n")
if (n_uniq < nrow(map)) {
  cat("NOTE: duplicated symbols after mapping:\n")
  print(map[, .N, by = symbol_human][N > 1])
}

fwrite(map, file.path(OUT_DIR, "farg_human_map.csv"))
cat("\nSaved:", file.path(OUT_DIR, "farg_human_map.csv"), "\n")
cat("== step02 done ==\n")
