# =============================================================================
# rerun_drug_only.R —— 仅重跑 step07 药敏段（免 GSVA，生成 drug_ic50.csv 含 FDR）
# 依赖：oncoPredict 1.3.1 + GDSC2 训练数据；内存约需 4-6 GB
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(oncoPredict); library(ggplot2)
})
OUT_DIR <- "analysis_out/hcc_out/immune"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GDSC_EXPR <- "analysis_out/data/GDSC2_Expr (RMA Normalized and Log Transformed).rds"
GDSC_RES  <- "analysis_out/data/GDSC2_Res.rds"

# ---- 数据 ----
L <- readRDS("analysis_out/hcc_data/hcc_data.rds")
expr_all_tumor <- L$expr[, L$group == "Tumor"]
rk <- fread("analysis_out/hcc_out/risk_training.csv")
group <- rk$risk_group[match(colnames(expr_all_tumor), rk$sample)]
cat("tumor samples:", ncol(expr_all_tumor), "| risk-annotated:",
    sum(!is.na(group)), "\n")

TrainingData  <- readRDS(GDSC_EXPR)
TrainingPtype <- readRDS(GDSC_RES)
cat("GDSC2 training:", nrow(TrainingData), "genes x", ncol(TrainingData),
    "cell lines | drugs:", ncol(TrainingPtype), "\n")

# ---- 保留 HCC 相关药物 ----
drug_pat <- c("Sorafenib", "Lenvatinib", "Regorafenib", "Cabozantinib",
              "Cisplatin", "Doxorubicin", "Gemcitabine", "5-Fluorouracil",
              "Oxaliplatin", "Paclitaxel", "Etoposide", "Irinotecan",
              "Mitomycin", "Vinblastine", "Vincristine", "Bleomycin",
              "Camptothecin", "Cytarabine", "Methotrexate",
              "Epirubicin", "Afatinib", "Erlotinib", "Gefitinib")
keep_drugs <- unique(unlist(lapply(drug_pat, function(p)
  grep(paste0("^", p), colnames(TrainingPtype), value = TRUE, ignore.case = TRUE))))
if (length(keep_drugs) < 5) keep_drugs <- colnames(TrainingPtype)
TrainingPtype <- TrainingPtype[, keep_drugs, drop = FALSE]
cat("drugs kept:", length(keep_drugs), "\n")

test_expr <- as.matrix(expr_all_tumor)
common <- intersect(rownames(TrainingData), rownames(test_expr))
cat("genes shared:", length(common), "| predicting IC50 (several minutes) ...\n")

fit_drugs <- function(drgs) {
  tryCatch(
    calcPhenotype(trainingExprData = TrainingData[common, ],
                  trainingPtype = TrainingPtype[, drgs, drop = FALSE],
                  testExprData = test_expr[common, ],
                  batchCorrect = "eb", powerTransformPhenotype = TRUE,
                  removeLowVaryingGenes = 0.2,
                  removeLowVaringGenesFrom = "homogenizeData",
                  minNumSamples = 10, printOutput = FALSE, parallel = FALSE),
    error = function(e) {
      cat("  batch failed:", conditionMessage(e), "-> single-drug loop\n")
      res_list <- lapply(drgs, function(dg)
        tryCatch(calcPhenotype(
          trainingExprData = TrainingData[common, ],
          trainingPtype = TrainingPtype[, dg, drop = FALSE],
          testExprData = test_expr[common, ],
          batchCorrect = "eb", powerTransformPhenotype = TRUE,
          removeLowVaryingGenes = 0.2,
          removeLowVaringGenesFrom = "homogenizeData",
          minNumSamples = 10, printOutput = FALSE, parallel = FALSE),
          error = function(e2) NULL))
      res_list <- res_list[!vapply(res_list, is.null, logical(1))]
      do.call(cbind, res_list)
    })
}
drug_n <- length(keep_drugs)
half_idx <- seq_len(ceiling(drug_n / 2))
ic50 <- if (drug_n == 1) fit_drugs(keep_drugs) else
  cbind(fit_drugs(keep_drugs[half_idx]), fit_drugs(keep_drugs[-half_idx]))
if (length(intersect(colnames(expr_all_tumor), rownames(ic50))) <
    length(intersect(colnames(expr_all_tumor), colnames(ic50)))) ic50 <- t(ic50)
ic50 <- ic50[intersect(colnames(expr_all_tumor), rownames(ic50)), , drop = FALSE]
cat("predicted IC50:", nrow(ic50), "samples x", ncol(ic50), "drugs\n")

# ---- 8 个 HCC 核心药物 + 全部保留药物，各算 P + FDR ----
drug_pat2 <- c("Sorafenib", "Lenvatinib", "Cisplatin", "Doxorubicin",
               "Gemcitabine", "5-Fluorouracil", "Oxaliplatin", "Paclitaxel")
drugs_interest <- unique(unlist(lapply(drug_pat2, function(p)
  grep(paste0("^", p), colnames(ic50), value = TRUE, ignore.case = TRUE))))
if (length(drugs_interest) == 0) drugs_interest <- colnames(ic50)[1:min(8, ncol(ic50))]

ddt <- data.table(sample = rownames(ic50), group = group)
for (dr in colnames(ic50)) ddt[[dr]] <- as.numeric(ic50[, dr])
ddt <- ddt[!is.na(group)]   # drop samples without risk annotation (NA group)
dl <- melt(ddt, id.vars = c("sample", "group"), variable.name = "drug",
           value.name = "IC50")
dl <- dl[!is.na(IC50)]
dw <- dl[, .(p = wilcox.test(IC50 ~ group)$p.value,
             median_High = median(IC50[group == "High"]),
             median_Low = median(IC50[group == "Low"])), by = drug]
dw[, FDR := p.adjust(p, method = "BH")]
dw[, direction := fifelse(median_High > median_Low,
                          "High higher IC50 (more resistant)",
                          "High lower IC50 (more sensitive)")]
setorder(dw, p)
fwrite(dw, file.path(OUT_DIR, "drug_ic50.csv"))

# 逐样本 IC50 矩阵（sample × drug + group）—— 供独立复现箱线图与核验各组 n
fwrite(ddt, file.path(OUT_DIR, "drug_ic50_per_sample.csv"))
# 分组分位数汇总（每 drug × group：n + 五数概括）—— 供核验箱线分位数
qsum <- dl[, .(n = .N,
               min = round(min(IC50), 4),
               Q1 = round(quantile(IC50, 0.25), 4),
               median = round(median(IC50), 4),
               Q3 = round(quantile(IC50, 0.75), 4),
               max = round(max(IC50), 4)), by = .(drug, group)]
fwrite(qsum[order(drug, group)], file.path(OUT_DIR, "drug_ic50_quantiles.csv"))
cat("      saved drug_ic50_per_sample.csv + drug_ic50_quantiles.csv\n")

cat("\n=== 8 个 HCC 核心药物（原始 P + BH-FDR + 方向）===\n")
print(dw[drug %in% drugs_interest])

g3 <- ggplot(dl[drug %in% drugs_interest], aes(group, IC50, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ drug, scales = "free") +
  scale_fill_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
  labs(x = NULL, y = "Predicted IC50", title = "Drug sensitivity by risk group") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "drug_boxplot.png"), g3, width = 10, height = 7, dpi = 300)
cat("saved drug_ic50.csv + drug_boxplot.png (300 dpi)\n")
