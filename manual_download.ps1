# =============================================================================
# manual_download.ps1 - download dependency data for the Ferro-aging HCC pipeline
# Usage: powershell -ExecutionPolicy Bypass -File "path/to/Ferro-aging\manual_download.ps1"
# =============================================================================
$ErrorActionPreference = 'Stop'
$proj = 'path/to/Ferro-aging'
$data = Join-Path $proj 'analysis_out\data'

Write-Host ''
Write-Host '=== 1/5 MSigDB 2026.1 (step09 GSEA; ~100-200 MB) ===' -ForegroundColor Cyan
$cache = Join-Path $env:LOCALAPPDATA 'R\cache\R\msigdbr'
New-Item -ItemType Directory -Force -Path $cache | Out-Null
$zip = Join-Path $cache 'msigdb.2026.1.zip'
$sum = Join-Path $cache 'msigdb.2026.1.summary.rds'
if (Test-Path $sum) {
    Write-Host '  already present, skip.'
} else {
    Write-Host '  downloading from Zenodo (may take minutes)...'
    curl.exe -L --retry 3 -o $zip 'https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1'
    if ($LASTEXITCODE -ne 0) { throw 'MSigDB download failed, check network.' }
    Expand-Archive -Path $zip -DestinationPath $cache -Force
    if (Test-Path $sum) { Write-Host '  done; msigdbr will use the local cache.' }
    else { Write-Host '  WARNING: summary.rds not found, check zip layout.' }
}

Write-Host ''
Write-Host '=== 2/5 GEO validation cohorts (step06) ===' -ForegroundColor Cyan
$geo = Join-Path $data 'GEO'
New-Item -ItemType Directory -Force -Path $geo | Out-Null
$geoFiles = @(
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz'; fn='GSE14520-GPL3921_series_matrix.txt.gz'},
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL3nnn/GPL3921/annot/GPL3921.annot.gz'; fn='GPL3921.annot.gz'},
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz'; fn='GSE14520_Extra_Supplement.txt.gz'},
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16757/matrix/GSE16757_series_matrix.txt.gz'; fn='GSE16757_series_matrix.txt.gz'},
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16757/suppl/GSE16757_Clinical.data.csv.gz'; fn='GSE16757_Clinical.data.csv.gz'}
)
foreach ($f in $geoFiles) {
    $dst = Join-Path $geo $f.fn
    if (Test-Path $dst) { Write-Host "  $($f.fn) already present, skip."; continue }
    Write-Host "  downloading $($f.fn) ..."
    curl.exe -L --retry 3 -o $dst $f.url
    if ($LASTEXITCODE -ne 0) { Write-Host "  WARNING: $($f.fn) download failed." }
}

Write-Host ''
Write-Host '=== 3/5 LIHC mutations (step07 TMB) ===' -ForegroundColor Cyan
$maf = Join-Path $data 'TCGA.LIHC.mutect_cbioportal.tsv'
if (Test-Path $maf) { Write-Host '  already present, skip.' }
else {
    Write-Host '  downloading LIHC mutations from cBioPortal API (~2-5 MB)...'
    curl.exe -L --retry 3 -o $maf 'https://www.cbioportal.org/api/molecular-profiles/lihc_tcga_mutations/mutations?sampleListId=lihc_tcga_all&projection=DETAILED&pageSize=1000000'
    if ($LASTEXITCODE -eq 0 -and (Test-Path $maf) -and (Get-Item $maf).Length -gt 1000000) {
        Write-Host '  success (cBioPortal format, step07 compatible).'
    } else {
        Write-Host '  download failed. Alternative (GDC portal):'
        Write-Host '  1) open https://portal.gdc.cancer.gov/ (free registration)'
        Write-Host '  2) search TCGA-LIHC -> Files -> Data Category = Simple Nucleotide Variation -> download MAF'
        Write-Host '  3) put the file in analysis_out\data\ named TCGA.LIHC.mutect*.'
    }
}
Write-Host '  NOTE: the mc3.v0.2.8.PUBLIC.maf(.gz) in the project root is an incomplete download (only ~1,986 samples), do not rely on it.'

Write-Host ''
Write-Host '=== 4/5 GSE149614 single-cell (step08; ~160 MB) ===' -ForegroundColor Cyan
$sc = Join-Path $data 'GSE149614'
New-Item -ItemType Directory -Force -Path $sc | Out-Null
$scFiles = @(
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz'; fn='GSE149614_HCC.scRNAseq.S71915.count.txt.gz'},
    @{url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz'; fn='GSE149614_HCC.metadata.updated.txt.gz'}
)
foreach ($f in $scFiles) {
    $dst = Join-Path $sc $f.fn
    if (Test-Path $dst) { Write-Host "  $($f.fn) already present, skip."; continue }
    Write-Host "  downloading $($f.fn) ..."
    curl.exe -L --retry 3 -o $dst $f.url
    if ($LASTEXITCODE -ne 0) { Write-Host "  WARNING: $($f.fn) download failed." }
}

Write-Host ''
Write-Host '=== 5/5 GDSC2 drug sensitivity (step07; manual via browser) ===' -ForegroundColor Cyan
Write-Host '  1) open https://osf.io/c6tfx/ (free OSF account may be required)'
Write-Host '  2) download two files: GDSC2_Expr (RData) and GDSC2_Res (RData)'
Write-Host '  3) place in analysis_out\data\ (keep original filenames)'

Write-Host ''
Write-Host 'Done. Then in RStudio run:' -ForegroundColor Green
Write-Host '  setwd("path/to/Ferro-aging"); source("00_master.R")'
