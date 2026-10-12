# Source-data reconstruction

These scripts reconstruct the historical analysis from original GWAS inputs. They are not necessary for the genotype-free cached checks described in the repository root README. Run them from `original/`, or set `SST2_PROJECT_DIR` explicitly to a separate reconstruction project directory. Existing result files at that location will be regenerated, so use a separate work directory to preserve the distributed baselines.

## Required inputs and paths

Obtain these files separately through their original providers and comply with their access/reuse terms:

- SCALLOP ST2 GWAS `GCST90012040`, GRCh37. Historical assay 151_ST2/Q01638; soluble-isoform exclusivity is not established.
- HERMES heart failure `GCST009541`, GRCh37.
- FinnGen R12 `AB1_OTHER_SEPSIS`, `I9_HEARTFAIL`, `K11_IBD_STRICT`, GRCh38.
- 1000 Genomes Phase 3 EUR PLINK panel, and PLINK 1.9.
- UCSC `hg19ToHg38.over.chain`.
- Comparison workbook from doi:10.1186/s12879-026-12748-2.

Environment variables override conventional paths:

| Variable | Conventional path beneath `SST2_PROJECT_DIR` |
| --- | --- |
| `SST2_SCALLOP_FILE` | `data/GCST90012040.tsv.gz` |
| `SST2_HERMES_FILE` | `data/GCST009541.tsv.gz` |
| `SST2_FINNGEN_SEPSIS_FILE` | `data/finngen_R12_AB1_OTHER_SEPSIS.gz` |
| `SST2_FINNGEN_HF_FILE` | `data/finngen_R12_I9_HEARTFAIL.gz` |
| `SST2_FINNGEN_IBD_FILE` | `data/finngen_R12_K11_IBD_STRICT.gz` |
| `SST2_EUR_BFILE` | `data/1kg_eur/EUR` (prefix without extension) |
| `SST2_CHAIN_FILE` | `data_metadata/hg19ToHg38.over.chain` |
| `SST2_BENCHMARK_SUPPLEMENT` | `data_metadata/benchmark_2026_bmc_infectious_diseases/12879_2026_12748_MOESM2_ESM.xls` |
| `PLINK_BIN` | `plink` on PATH, otherwise `tools/plink` |

No raw inputs, chain, software executable, or third-party workbook is bundled. No script automatically downloads these sources. Shell `gzip` and `awk` are used to subset compressed GWAS files.

## Historical run order

```sh
Rscript code/01_prepare_instruments.R
Rscript code/02_forward_mr.R
Rscript code/03_reproduce_2026_il1rl1_sepsis.R
Rscript code/04_colocalization.R
Rscript code/05_multivariable_mr.R
python3 code/06b_gwas_catalog_confounder_screen.py
python3 code/06c_confounder_screen_sensitivity.py
Rscript code/06d_confounder_screened_mvmr.R
```

The Catalog script reuses its cached table by default. Network requests occur only when that cache is absent or `SST2_REQUERY_CATALOG=1` is explicitly set; the repository's root cached runner does not invoke it. The frozen screen is dated 22 July 2026 and is not a prospectively registered validity test.

Required R packages include `data.table`, `rtracklayer`, `GenomicRanges`, `MendelianRandomization`, `MVMR`, `coloc`, and `readxl`. Exact historical sessions are under `original/results`. Python screen scripts use the standard library.

Outputs retain historical field names for traceability. An `IBD_positive_control` filename denotes an illustrative benchmark, not a validated positive control or novel target discovery. `confounder` flags cannot distinguish vertical/horizontal pleiotropy or prove instrument invalidity. Updated LD-aware, uncertainty, overlap-stress, and estimability analyses are in `revision/`; consult the root README for their limitations and cached reproduction.
