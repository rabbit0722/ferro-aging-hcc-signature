# -*- coding: utf-8 -*-
"""HPA automation: extract liver-cancer IHC staining + prognostic for the
4 signature genes from the per-gene XMLs, download representative images."""
import xml.etree.ElementTree as ET
import urllib.request, csv, os, sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"path/to/Ferro-aging"
XML_DIR = os.path.join(BASE, "analysis_out", "data", "hpa")
OUT_DIR = os.path.join(BASE, "analysis_out", "hcc_out", "HPA")
os.makedirs(OUT_DIR, exist_ok=True)

GENES = {"HMOX1": "ENSG00000100292", "TFRC": "ENSG00000072274",
         "SOCS2": "ENSG00000120833", "E2F1": "ENSG00000101412"}
# our mRNA log2FC (tumor vs normal, from step05 farg_deg.csv)
LOGFC = {"HMOX1": -1.08, "TFRC": 1.54, "SOCS2": -1.95, "E2F1": 4.23}

def parse_gene(gene):
    p = os.path.join(XML_DIR, gene + ".xml")
    root = ET.parse(p).getroot()
    out = {"gene": gene, "summary": "", "normal_liver": "",
           "stain_high": 0, "stain_medium": 0, "stain_low": 0,
           "stain_neg": 0, "n_patients": 0,
           "prog_tcga": "", "prog_tcga_p": "", "prog_valid": "", "prog_valid_p": "",
           "img_url": ""}
    for el in root.iter():
        if el.get("assayType") != "cancer":
            continue
        s = el.find("summary")
        # 1) staining summary text
        if s is not None and s.text and "staining" in s.text.lower():
            out["summary"] = s.text.strip()
        for data in el.findall("data"):
            tis = data.find("tissue")
            name = (tis.text or "").strip() if tis is not None else ""
            if name == "Liver cancer":
                tc = data.find("tissueCell")
                if tc is not None:
                    for lv in tc.findall("level"):
                        t = lv.get("type"); c = int(lv.get("count") or 0)
                        if t == "staining":
                            out["n_patients"] += c
                            v = (lv.text or "").strip().lower()
                            if v == "high": out["stain_high"] += c
                            elif v == "medium": out["stain_medium"] += c
                            elif v == "low": out["stain_low"] += c
                            elif v == "not detected": out["stain_neg"] += c
                # representative image from patients
                for pat in data.findall("patient"):
                    if out["img_url"]:
                        break
                    for samp in pat.findall("sample"):
                        for ai in samp.findall("assayImage"):
                            for im in ai.findall("image"):
                                u = im.find("imageUrl")
                                if u is not None and u.text:
                                    out["img_url"] = u.text
                                    break
            if name == "Liver Hepatocellular Carcinoma (TCGA)":
                sa = data.find("survivalAnalysis")
                if sa is not None:
                    out["prog_tcga"] = sa.get("prognosticType", "") + " (" + sa.get("prognostic", "") + ")"
                    out["prog_tcga_p"] = sa.get("pValue", "")
            if name == "Liver Hepatocellular Carcinoma (validation)":
                sa = data.find("survivalAnalysis")
                if sa is not None:
                    out["prog_valid"] = sa.get("prognosticType", "") + " (" + sa.get("prognostic", "") + ")"
                    out["prog_valid_p"] = sa.get("pValue", "")
        # normal liver tissue level (assayType = tissue)
    for el in root.iter():
        if el.get("assayType") == "tissue":
            for data in el.findall("data"):
                tis = data.find("tissue")
                if tis is not None and (tis.text or "").strip() == "Liver":
                    lv = data.find("level")
                    if lv is not None and lv.get("type") == "expression":
                        out["normal_liver"] = (lv.text or "").strip()
    return out

rows = [parse_gene(g) for g in GENES]

# download images
for r in rows:
    if r["img_url"]:
        dst = os.path.join(OUT_DIR, "HPA_" + r["gene"] + "_liver_cancer.jpg")
        try:
            req = urllib.request.Request(r["img_url"],
                                         headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as resp, open(dst, "wb") as f:
                f.write(resp.read())
            r["img_file"] = os.path.basename(dst)
        except Exception as e:
            r["img_file"] = "FAILED: " + str(e)[:60]
    else:
        r["img_file"] = ""

with open(os.path.join(OUT_DIR, "hpa_summary.csv"), "w", newline="",
          encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("=== HPA 蛋白层面验证汇总（肝癌） ===\n")
for r in rows:
    print(f"{r['gene']} | mRNA logFC {LOGFC[r['gene']]:+.2f}")
    print(f"  染色: {r['n_patients']} 例 (high {r['stain_high']}, medium {r['stain_medium']}, "
          f"low {r['stain_low']}, neg {r['stain_neg']}) | 正常肝: {r['normal_liver']}")
    print(f"  总结: {r['summary'][:120]}")
    print(f"  HPA-TCGA 预后: {r['prog_tcga']} (p={r['prog_tcga_p']}) | 验证队列: "
          f"{r['prog_valid']} (p={r['prog_valid_p']})")
    print(f"  图: {r['img_file']}")
    print()
print("saved: analysis_out/hcc_out/HPA/hpa_summary.csv + 4 张代表图")
