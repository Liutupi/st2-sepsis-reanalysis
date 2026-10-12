# ST2 genetic instruments and sepsis: reproducibility analyses

Research code and selected, machine-readable derived summary data for a source-verified reanalysis of circulating ST2 genetic instruments and sepsis susceptibility. This repository is a correction/reproducibility resource, not a claim of a new biological target or proof that ST2 has no causal effect.

## What is included

- `original/code`: instrument construction, build/allele harmonization, ordinary MR, single-signal colocalization, MVMR, and the frozen direct-variant Catalog screen.
- `original/results`: selected-variant and regional summary statistics, instrument-selection records, harmonization decisions, and historical numerical outputs.
- `revision/ld_mr`: same-49 signed-LD generalized IVW and a base-R cached-input reproduction script.
- `revision/reporting`: two-SNP, power, full/screened MVMR, robust-estimator, and hypothetical covariance diagnostics.
- `revision/multisignal_coloc`: matched regional summaries, signed LD, SuSiE/coloc.susie code, and diagnostics including non-estimable analyses.
- `revision/evidence`: structured source/cohort/assay verification records; the third-party papers and supplements themselves are not redistributed.

No manuscript text, reviewer/editor correspondence, school documents, emails, credentials, individual reference genotypes, full-genome source GWAS, or third-party paper full text is included. Public-reference genotypes were used locally to construct the supplied signed LD, but are not distributed here. Recomputable fitted objects and eigendecompositions are omitted.

## Offline numerical reproduction

Run from the repository root, using R on your `PATH`:

```sh
Rscript reproduce_cached.R ld
Rscript reproduce_cached.R diagnostics
```

The LD check requires base R only. The diagnostics require `data.table`, `MVMR`, `MendelianRandomization`, and `robustbase`. No command above downloads data or installs software. Outputs are written under the ignored `recomputed/` directory; distributed baseline files are not overwritten. The runner checks numerical agreement against the supplied baselines with tolerance `1e-7`.

For the more expensive regional fine-mapping attempt:

```sh
Rscript reproduce_cached.R coloc
```

This requires `data.table`, `coloc`, and `susieR`. It stages the supplied regional summary and LD files, then runs `01_multisignal_coloc.R` followed by `02_validate_and_summarize.R`. It can take several minutes. The documented full exposure fit is expected not to converge; absence of a sepsis credible set means multi-signal PP.H4 is **not estimable, not zero**. Expected model limitations are not software-installation failures. Eigen caches and fitted objects are recreated locally, not included in the public tree.

The direct base-R entry also works from any working directory:

```sh
Rscript revision/ld_mr/recheck_from_cached_summary_and_ld.R
```

That direct entry writes `cached_recheck/` beside its input files. Prefer the root runner to keep regenerated files separate.

## Full source-data reconstruction

See [original/code/README.md](original/code/README.md). Source GWAS, the build-conversion chain, the published comparison workbook, and the PLINK reference panel must be obtained separately through their providers. Set the documented environment variables to files that you are authorized to use. No downloads are performed automatically. Exact original software versions are recorded in the session-information files.

The public copies change path defaults and add portable runners/documentation; analytical formulas and distributed numerical results are retained. Historical filenames/fields such as `positive_control`, `confounder`, `decode`, and `prespecified` are preserved for traceability. They do not establish a validated positive control, a pleiotropic mechanism, true deCODE exposure provenance, or prospective registration. The exposure input is SCALLOP `GCST90012040`, not deCODE; the main interpretation concerns the normalized circulating source-assay signal.

## Scientific limitations

The main cis analysis uses only two variants. Strong exposure F statistics do not establish exclusion restrictions or adequate power for small outcome effects. The direct Catalog screen is post hoc, index-only, lacks LD proxies, and cannot distinguish horizontal pleiotropy, vertical mediation, or LD tagging. Screened models are sensitivities, not validity certification.

The matched-49 LD-adjusted estimate remains modestly positive, with reduced precision. Crossing a P-value threshold is not proof that effects differ; strict reclumping also changes selection/information. A source-study FDR value of 1 is not a probability that an individual finding is false. Finite-panel LD uncertainty is not propagated into the IVW interval.

Multiple-signal colocalization is limited by weak regional sepsis evidence, absent sepsis credible sets, exposure LD–z mismatch/nonconvergence, and a small 503-person broader-European reference rather than in-sample SCALLOP or Finnish-specific LD. Failed exposure PIPs/credible sets are diagnostic candidate outputs, not validated biological signals.

Cohort-name overlap is not measured participant overlap. Exact shared participants, residual trait correlations, and exposure-error covariance are unavailable. The expanded rho grid is a hypothetical stress test, not an empirical overlap estimate or overlap-corrected regression. The historical ST2 assay is not established to be soluble-isoform-exclusive. SCALLOP and UKB-PPP assay scales are not assumed interchangeable; HF liability is on a log-odds scale. Q statistics with different definitions/degrees of freedom must not be conflated.

## Source terms and reuse status

Access and reuse of third-party data remain subject to the original providers' terms, approvals, and attribution requirements. [DATA_SOURCES.tsv](DATA_SOURCES.tsv) identifies the primary sources. This repository does not grant extra rights to SCALLOP, HERMES, FinnGen, 1000 Genomes, the GWAS Catalog, UCSC, or published comparison data. Obtain any required access/permissions independently and do not attempt participant re-identification.

Author-written code is provided under the [MIT License](LICENSE), as approved by the authors. This code license does not apply to third-party data or grant rights beyond the original providers' terms. The supplied derived summary records retain their source attribution and applicable data-use conditions.

`SHA256SUMS.tsv` records the distributed files, sizes, and SHA-256 hashes, excluding itself. Regenerated files are intentionally ignored by Git.
