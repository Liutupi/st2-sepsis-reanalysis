# Regional colocalization and fine-mapping

This module analyzes 2,914 allele-aligned SCALLOP ST2/FinnGen sepsis variants using single-signal `coloc.abf`, SuSiE-RSS fine-mapping, and `coloc.susie` signal-pair checks. The supplied signed LD matrix was estimated from 503 European 1000 Genomes reference participants; no individual-level genotypes are included.

## Reproduce from the supplied summary statistics and LD

From the repository root, run:

```sh
Rscript revision/multisignal_coloc/01_multisignal_coloc.R
Rscript revision/multisignal_coloc/02_validate_and_summarize.R
```

Defaults are resolved from the script location, so the scripts can also be called by an absolute path from another working directory. Required R packages are `data.table`, `coloc`, `susieR`, and `ggplot2`; recorded execution versions are in the session-information files.

The default `COLOC_USE_PREMATCHED=1` reads only `matched_regional_summary_statistics.tsv` and `signed_LD_EUR_503.rds` as primary inputs. Run script 01 before script 02: script 01 creates the fit, diagnostic, and eigendecomposition files required by script 02. Results are regenerated in this directory. To use a separate existing output directory, copy both primary input files there and set `COLOC_REVISION_OUT` to that directory for both commands. `LD_SNP_order.txt` provides an additional human-readable ordering audit.

## Optional reconstruction of LD

Set `COLOC_USE_PREMATCHED=0` to rebuild LD. This mode requires an available EUR PLINK `.bed/.bim/.fam` reference prefix in `EUR_REFERENCE_PREFIX`, and either an executable path in `PLINK_BIN` or `plink` on `PATH`. The regional input defaults to `original/results/colocalization/harmonized_region_Sepsis.tsv` relative to the repository root; override its parent directory with `SST2_ORIGINAL_DIR` if needed. Reference genotypes are not required for the default cached-input mode. Do not add generated `reference_region.*` individual-level files to the public repository.

## Interpretation and limitations

Fits use 95% credible sets, minimum within-set absolute correlation 0.5, fixed residual variance, and bounded convergence checks. The sepsis data identify no qualifying credible set, so multiple-signal PP.H4 is **not estimable**, not zero. The full exposure fit does not converge, and the common/high-sample-size exposure subset reaches the imposed signal cap at both L=10 and L=20. Failed fits, raw PIPs, and capacity-sensitive candidate sets are audit outputs, not validated biological signals.

The small external reference is neither in-sample SCALLOP LD nor Finnish-specific LD. LD–z inconsistency, variable exposure sample sizes, and weak regional outcome evidence limit fine-mapping and prevent exclusion of a weak shared effect. The recorded LD–z diagnostic parameter is not used to shrink the LD matrix supplied to fine-mapping.
