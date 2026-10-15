# A ferro-aging-related four-gene signature for hepatocellular carcinoma

Analysis scripts for the manuscript *"A four-gene ferro-aging-related score is associated with survival in hepatocellular carcinoma"* (Scientific Reports).

## Overview

We derived a four-gene score (HMOX1, TFRC, SOCS2, E2F1) from a published ferro-aging gene set (Liu et al., Cell Metabolism 2026) and evaluated its association with overall survival in TCGA-LIHC (training) and two independent cohorts (GSE14520, GSE16757), followed by multi-omics molecular characterization (CNV, methylation, immune/TMB, predicted drug sensitivity, single-cell, virtual knockout, ACSL4).

## Scripts

| Script | Purpose |
|---|---|
| `00_setup.R` | Install/check required R packages |
| `00_master.R` | Run the full pipeline |
| `step01_farg_profile.R` | Ferro-aging gene-set profiling |
| `step02_human_map.R` | Monkey-to-human gene mapping |
| `step03_tcga_score.R` | TCGA scoring |
| `step04_hcc_prep.R` | HCC data preparation |
| `step05_hcc_signature.R` | Signature construction (DEG → univariate → LASSO → stepwise) |
| `step06_hcc_validate.R` | External cohort validation (GSE14520, GSE16757) |
| `step07_hcc_immune_drug.R` | Immune (TIDE/checkpoint), TMB, drug sensitivity (oncoPredict/GDSC2) |
| `step08_hcc_singlecell.R` | Single-cell (GSE149614) + virtual knockout (scTenifoldKnk) |
| `step09_hcc_robustness.R` | Conditional bootstrap + gene-set comparison + DCA |
| `step10_hcc_extras.R` | Supplementary analyses |
| `step11_hcc_genomics.R` | Genomic features |
| `step12_hcc_aging.R` | Public aging gene-set associations |
| `step13_hcc_methylation.R` | DNA methylation (450K) |
| `step14_hcc_cnv.R` | Copy-number (GISTIC2) |
| `step15_hcc_acsl4.R` | ACSL4 target association |
| `manual_download.ps1` | Download public data (TCGA/GEO) |
| `analysis_out/rerun_bootstrap_only.R` | Lightweight conditional bootstrap (no GSVA) |
| `analysis_out/rerun_drug_only.R` | Lightweight drug-sensitivity re-run |
| `analysis_out/export_source_tables.R` | Export full-gene DEG table |
| `analysis_out/module_go_annotation.R` | GO annotation |
| `analysis_out/coexpression_clustering.R` | Co-expression clustering |
| `analysis_out/hpa_extract.py` | Human Protein Atlas extraction |
| `analysis_out/data/docking/compile_docking_panel.R` | Molecular docking summary |
| `analysis_out/data/string/make_cytoscape_attributes.py` | STRING network attributes |
| `analysis_out/hcc_out/extras/cmap_*.py` | CMap (LINCS L1000) queries |

## Environment

- R 4.4.1 (Bioconductor 3.20)
- Key packages: survival 3.6.4, glmnet 5.0, limma 3.62.2, timeROC 0.4.1, survminer 0.5.2, rms 8.1.1, GSVA 2.0.7, clusterProfiler 4.14.6, oncoPredict 1.3.1, maftools 2.22.0, IOBR 2.2.3, dcurves, scTenifoldKnk

## Notes

- Scripts assume the project directory is the working directory (`setwd("path/to/Ferro-aging")`).
- `step07` uses the `tidepy` Python CLI; set `TIDEPY_PYTHON` / `TIDEPY_SCRIPT` to your local paths.
- `manual_download.ps1` sets `$proj` to the project directory.
- Public data: TCGA-LIHC (UCSC Xena), GSE14520, GSE16757, GSE149614 (GEO), GDSC2 (OSF), MSigDB.
