#!/usr/bin/env python3
"""Run prespecified UVMR sensitivity analyses after the direct-variant screen."""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "results" / "instrument_validity" / "instrument_confounder_screen.tsv"
HARMONIZED = ROOT / "results" / "forward_mr" / "all_exact_harmonized_variants.tsv"
OUT = ROOT / "results" / "instrument_validity" / "confounder_screened_uvmr_sensitivity.tsv"


def as_float(value: str) -> float:
    return float(value.strip('"'))


def two_sided_normal_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2))


def chi_square_sf_even_df(value: float, df: int) -> float:
    """Upper-tail chi-square probability for the even degrees used here."""
    if df <= 0 or df % 2:
        raise ValueError("This reproducibility helper requires a positive even df")
    half = value / 2
    return math.exp(-half) * sum(half**j / math.factorial(j) for j in range(df // 2))


def main() -> None:
    with SCREEN.open(newline="", encoding="utf-8") as handle:
        screened = list(csv.DictReader(handle, delimiter="\t"))
    keep = {
        row["rsid"]
        for row in screened
        if row["in_uvmr_genomewide_set"] == "yes"
        and row["potential_confounder_flags"] == "none identified"
    }
    with HARMONIZED.open(newline="", encoding="utf-8") as handle:
        data = list(csv.DictReader(handle, delimiter="\t"))

    rows = []
    for outcome in sorted({row["outcome"].strip('"') for row in data}):
        d = [
            row for row in data
            if row["outcome"].strip('"') == outcome
            and row["rsid"].strip('"') in keep
            and row["align_ok"].strip('"').upper() == "TRUE"
        ]
        bx = [as_float(row["beta_exposure"]) for row in d]
        by = [as_float(row["beta_outcome"]) for row in d]
        sy = [as_float(row["se_outcome"]) for row in d]
        weights = [1 / value**2 for value in sy]
        denom = sum(w * x**2 for w, x in zip(weights, bx))
        estimate = sum(w * x * y for w, x, y in zip(weights, bx, by)) / denom
        se_fixed = math.sqrt(1 / denom)
        q = sum(w * (y - estimate * x) ** 2 for w, x, y in zip(weights, bx, by))
        q_df = len(d) - 1
        q_p = chi_square_sf_even_df(q, q_df)
        se_random = se_fixed * math.sqrt(max(1, q / q_df))
        z = estimate / se_random
        rows.append({
            "analysis": "GWAS-Catalog-screened genome-wide IVW (random effects)",
            "outcome": outcome,
            "nsnp": str(len(d)),
            "retained_rsids": "; ".join(sorted(row["rsid"].strip('"') for row in d)),
            "beta": f"{estimate:.10g}",
            "se": f"{se_random:.10g}",
            "or": f"{math.exp(estimate):.10g}",
            "ci_lower": f"{math.exp(estimate - 1.96 * se_random):.10g}",
            "ci_upper": f"{math.exp(estimate + 1.96 * se_random):.10g}",
            "p": f"{two_sided_normal_p(z):.10g}",
            "q": f"{q:.10g}",
            "q_df": str(q_df),
            "q_p": f"{q_p:.10g}",
            "screen_rule": "Excluded direct GWAS Catalog associations at P<5e-8 within prespecified potential-confounder categories; no proxy screen.",
        })

    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
