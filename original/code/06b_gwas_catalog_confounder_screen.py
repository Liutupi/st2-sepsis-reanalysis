#!/usr/bin/env python3
"""Screen MR instruments against the GWAS Catalog for predefined confounders.

The screen intentionally uses direct, index-variant associations only (no LD
proxies) and is not a claim that a null catalogue result proves instrument
validity.  It produces an auditable record of the query date, raw returned
records, and a prespecified confounder classification for the sST2 UVMR and
joint MVMR instrument sets.
"""

from __future__ import annotations

import csv
import json
import os
import time
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "instrument_validity"
OUT.mkdir(parents=True, exist_ok=True)

GW_FILE = ROOT / "results" / "instruments" / "sst2_gw_r2_0.001_instruments.tsv"
MVMR_FILE = ROOT / "results" / "multivariable_mr" / "harmonized_mvmr_instruments.tsv"

# Potential confounders of circulating sST2 and sepsis susceptibility were
# specified before review of the returned association records.
CONFOUNDER_KEYWORDS = {
    "adiposity_or_metabolic": ["body mass index", "obesity", "waist", "type 2 diabetes", "insulin"],
    "blood_pressure_or_cardiometabolic": ["blood pressure", "hypertension", "cholesterol", "lipid", "triglycer", "apolipoprotein"],
    "kidney_or_liver_function": ["estimated glomerular", "kidney", "renal", "creatinine", "albumin", "alkaline phosphatase", "liver"],
    "smoking_or_alcohol": ["smoking", "tobacco", "alcohol"],
    "socioeconomic_or_behavioral": ["educational", "social deprivation"],
    "hematologic_or_immune": ["red cell", "white blood", "platelet", "lymphocyte", "neutrophil", "blood cell"],
}


def read_rsid(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return sorted({row["rsid"].strip('"') for row in csv.DictReader(handle, delimiter="\t") if row.get("rsid")})


def query_variant(rsid: str) -> list[dict[str, str]]:
    url = f"https://www.ebi.ac.uk/gwas/api/v2/variants/{rsid}/associations/download"
    req = Request(url, headers={"User-Agent": "sST2-MR-methods-audit/1.0"})
    try:
        with urlopen(req, timeout=45) as response:
            payload = response.read().decode("utf-8")
    except (HTTPError, URLError, TimeoutError) as exc:
        return [{"variant": rsid, "query_error": str(exc)}]
    return list(csv.DictReader(payload.splitlines(), delimiter="\t"))


def p_as_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def classify_trait(trait: str) -> list[str]:
    lowered = (trait or "").lower()
    return [name for name, words in CONFOUNDER_KEYWORDS.items() if any(word in lowered for word in words)]


def main() -> None:
    uvrs = read_rsid(GW_FILE)
    mvrs = read_rsid(MVMR_FILE)
    all_rsids = sorted(set(uvrs) | set(mvrs))
    raw_path = OUT / "gwas_catalog_direct_variant_associations.tsv"
    if raw_path.exists() and os.environ.get("SST2_REQUERY_CATALOG") != "1":
        with raw_path.open(newline="", encoding="utf-8") as handle:
            raw_rows = list(csv.DictReader(handle, delimiter="\t"))
    else:
        raw_rows: list[dict[str, str]] = []
        for i, rsid in enumerate(all_rsids, start=1):
            rows = query_variant(rsid)
            for row in rows:
                row["queried_rsid"] = rsid
                raw_rows.append(row)
            if i < len(all_rsids):
                time.sleep(0.15)
        raw_fields = sorted({field for row in raw_rows for field in row})
        with raw_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=raw_fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(raw_rows)


    summaries: list[dict[str, str]] = []
    for rsid in all_rsids:
        direct = [row for row in raw_rows if row.get("queried_rsid") == rsid and not row.get("query_error")]
        relevant = []
        flags: set[str] = set()
        for row in direct:
            p = p_as_float(row.get("pValue", ""))
            trait = row.get("traitName", "")
            if p is not None and p < 5e-8:
                matched = classify_trait(trait)
                if matched:
                    flags.update(matched)
                    relevant.append(f"{trait} (P={row.get('pValue')})")
        summaries.append({
            "rsid": rsid,
            "in_uvmr_genomewide_set": "yes" if rsid in uvrs else "no",
            "in_mvmr_set": "yes" if rsid in mvrs else "no",
            "n_catalog_records": str(len(direct)),
            "potential_confounder_flags": "; ".join(sorted(flags)) if flags else "none identified",
            "flagged_direct_associations": " | ".join(relevant) if relevant else "none at P<5e-8 under prespecified categories",
            "screen_scope": "GWAS Catalog direct index-variant associations; no LD proxies",
            "query_date": str(date.today()),
        })

    summary_path = OUT / "instrument_confounder_screen.tsv"
    fields = list(summaries[0])
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    metadata = {
        "data_source": "NHGRI-EBI GWAS Catalog API v2 direct-variant association export",
        "query_date": str(date.today()),
        "p_value_threshold": "5e-8",
        "proxy_policy": "No LD proxies queried",
        "predefined_categories": CONFOUNDER_KEYWORDS,
        "interpretation": "A null result does not prove absence of pleiotropy; this is a transparent screen for reported direct associations in prespecified confounder categories.",
    }
    (OUT / "instrument_confounder_screen_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Wrote {summary_path}")
    print(f"Wrote {raw_path}")


if __name__ == "__main__":
    main()
