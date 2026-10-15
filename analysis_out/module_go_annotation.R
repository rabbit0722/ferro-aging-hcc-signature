# =============================================================================
# module_go_annotation.R —— 生成 95 个 FARGs 的四模块 GO 归类依据表
# 用 org.Hs.eg.db 的 GO 生物学过程(BP)注释，客观支撑 M1–M4 手工分类。
# 输出：analysis_out/step01/module_GO_annotation.csv
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(GO.db)
})

mod <- fread("analysis_out/step01/module_assignment.csv")
stopifnot("Module" %in% names(mod), "SYMBOL" %in% names(mod))
genes <- unique(mod$SYMBOL)
cat("模块分配：", nrow(mod), "行 |", length(genes), "基因\n")

# 每个模块的关键词（GO 术语主题）
KW <- list(
  M1_iron_ferroptosis = c("iron", "ferroptosis", "lipid peroxid", "heme",
                          "glutathione", "reactive oxygen", "oxidative stress",
                          "oxidant", "cellular response to oxidative"),
  M2_senescence = c("senescen", "aging", "cell cycle", "mitotic", "replication",
                    "dna damage", "cell division", "proliferation", "cell growth"),
  M3_inflammation = c("inflamm", "immune", "interleukin", "cytokine", "interferon",
                      "defense response", "nf-kappab", "toll-like", "chemotaxis",
                      "jak-stat", "innate immune", "chemokine"),
  M4_signaling_other = c("signaling", "signal transduction", "hippo", "wnt", "mapk",
                         "er stress", "response to hypoxia", "transmembrane receptor",
                         "intracellular signal", "apoptotic signaling")
)

# SYMBOL -> ENTREZID
s2e <- select(org.Hs.eg.db, keys = genes, columns = "ENTREZID", keytype = "SYMBOL")
s2e <- s2e[!is.na(s2e$ENTREZID), ]
cat("SYMBOL->ENTREZID 映射：", nrow(s2e), "/", length(genes), "\n")

# ENTREZID -> GO BP
go <- select(org.Hs.eg.db, keys = unique(s2e$ENTREZID),
             columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")
go <- go[!is.na(go$GO) & go$ONTOLOGY == "BP", ]
go$term <- unname(Term(go$GO))
go$term[is.na(go$term)] <- ""
go <- go[go$term != "", ]
setDT(go); setDT(s2e)
go2 <- merge(go, s2e[, .(SYMBOL, ENTREZID)], by = "ENTREZID", allow.cartesian = TRUE)
cat("GO BP 注释（去重后按基因）：\n")

res <- rbindlist(lapply(genes, function(g) {
  m <- mod[SYMBOL == g]$Module[1]
  gg <- go2[SYMBOL == g, ]
  gg <- unique(gg[, .(GO, term)])
  gg$low <- tolower(gg$term)
  if (nrow(gg) == 0)
    return(data.table(SYMBOL = g, Module = m, GO_ID = NA,
                      GO_term = "(无 GO BP 注释)", n_support = 0L,
                      also_matches = ""))
  # 每个模块关键词的命中
  hit_list <- lapply(names(KW), function(k) {
    kk <- KW[[k]]
    h <- gg[grepl(paste(kk, collapse = "|"), gg$low), ]
    unique(h$term)
  })
  names(hit_list) <- names(KW)
  own <- hit_list[[m]]
  other <- names(KW)[names(KW) != m & lengths(hit_list)[names(KW) != m] > 0]
  # 若本模块无命中，回退为全部 BP 前 3 条并标记
  if (length(own) == 0) {
    own <- head(gg$term, 3)
    flag <- " (fallback: 无主题命中)"
  } else flag <- ""
  data.table(SYMBOL = g, Module = m,
             GO_ID = paste(gg$GO[match(own, gg$term)], collapse = "; "),
             GO_term = paste0(paste(own, collapse = "; "), flag),
             n_support = length(own),
             also_matches = paste(other, collapse = ", "))
}))
setorder(res, Module, SYMBOL)

fwrite(res, "analysis_out/step01/module_GO_annotation.csv")
cat("\n各模块命中情况：\n")
print(res[, .(n_genes = .N, n_with_hit = sum(n_support > 0 & !grepl("fallback", GO_term))), by = Module])
cat("\n多模块交叉命中的基因（pleiotropy）：\n")
print(res[also_matches != "", .(SYMBOL, Module, also_matches)])
cat("\n前 20 行预览：\n")
print(res[1:20, .(SYMBOL, Module, GO_term)])
cat("\n输出：analysis_out/step01/module_GO_annotation.csv\n")
