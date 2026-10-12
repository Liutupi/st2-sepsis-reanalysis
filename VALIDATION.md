# Offline validation record

Validated locally on 12 October 2026 using R 4.5.1; exact package versions remain in the supplied session-information files. No data downloads or software installations were performed.

| Entry | Check | Result |
| --- | --- | --- |
| `Rscript reproduce_cached.R ld` | Ordinary/generalized fixed/random IVW and weighted PCA against supplied baselines | All 24 numerical checks passed; maximum absolute difference 3.56e-14 |
| `Rscript reproduce_cached.R diagnostics` | Two-SNP estimates/power, full/screened MVMR, expanded covariance stress grid, robust implementation and weights | Seven tables matched; maximum absolute difference 0 |
| `Rscript reproduce_cached.R coloc` | Single-signal prior sensitivities, SuSiE convergence/credible-set summaries, sepsis N/L sensitivities and exposure L-capacity check | Four baseline tables matched; maximum absolute difference 0 |
| Script syntax | 12 R scripts and 2 Python scripts | All parsed successfully |
| Original numerical records | 77 files under `original/results` | Byte-identical to frozen analysis records |
| Signed-LD binary contents | 49×49 and 2914×2914 RDS matrices | Numeric matrices with variant rsID dimnames; no individual identifiers |
| Publication privacy/allowlist check | Text, file types, credentials/private-path patterns, reference-person IDs | No unresolved findings; R S4 slot-access syntax is distinguished from email addresses |

The colocalization check reproduced the documented nonconverged full exposure fit and absence of a sepsis credible set. Reproducibility of these limitations does not validate a causal null or convert non-estimable PP.H4 into zero. Diagnostic warnings about weak outcome association and external LD are expected and retained.

Full source-GWAS reconstruction was not rerun for this public packaging check. Those upstream scripts require separately authorized source data, PLINK, and a reference panel. Public-copy changes are limited to path defaults, explicit environment-variable configuration, external-command argument quoting, portable runners, and documentation; no analysis formula was changed.

Generated verification outputs are excluded from the public tree. `SHA256SUMS.tsv` records the distributed files and must be regenerated if any public file changes.
