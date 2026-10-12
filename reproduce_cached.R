#!/usr/bin/env Rscript
# Portable offline runner. Never downloads data or installs packages.
options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)
mode <- if (length(argv)) argv[1] else "ld"
stopifnot(mode %in% c("ld", "diagnostics", "coloc", "all"))
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else normalizePath(getwd())
rscript <- file.path(R.home("bin"), "Rscript")
dir.create(file.path(root, "recomputed"), showWarnings = FALSE)

run <- function(script, args = character()) {
  status <- system2(rscript, vapply(c(file.path(root, script), args), shQuote, ""))
  if (status != 0L) stop("Analysis failed: ", script, "; exit status ", status)
}
stage <- function(source_dir, dest_dir, names) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(file.path(source_dir, names), file.path(dest_dir, names), overwrite = TRUE)
  if (!all(ok)) stop("Could not stage all cached inputs.")
}
check_numeric <- function(reference, observed, keys = character(), tolerance = 1e-7) {
  a <- read.delim(reference, check.names = FALSE)
  b <- read.delim(observed, check.names = FALSE)
  if (length(keys)) {
    key <- function(x) do.call(paste, c(x[keys], sep = "|"))
    stopifnot(!anyDuplicated(key(a)), setequal(key(a), key(b)))
    b <- b[match(key(a), key(b)), , drop = FALSE]
  }
  stopifnot(nrow(a) == nrow(b), identical(names(a), names(b)))
  numeric_columns <- names(a)[vapply(a, is.numeric, logical(1))]
  diffs <- vapply(numeric_columns, function(nm) {
    stopifnot(identical(is.na(a[[nm]]), is.na(b[[nm]])))
    d <- abs(a[[nm]] - b[[nm]])
    if (all(is.na(d))) 0 else max(d, na.rm = TRUE)
  }, numeric(1))
  max_diff <- max(c(0, diffs))
  if (max_diff >= tolerance) stop("Numerical mismatch in ", basename(reference), ": ", max_diff)
  cat("PASS ", basename(reference), ": max absolute difference = ", format(max_diff, digits = 6), "\n", sep = "")
  invisible(max_diff)
}

if (mode %in% c("ld", "all")) {
  src <- file.path(root, "revision/ld_mr")
  dst <- file.path(root, "recomputed/ld_mr")
  stage(src, dst, c("ld_order_harmonized_instruments.tsv", "signed_ld_exposure_effect_allele_aligned.rds",
                    "strict_r2_0.001_variant_ids.txt", "ivw_estimates_fixed_and_random.tsv", "weighted_pca_ivw_sensitivity.tsv"))
  run("revision/ld_mr/recheck_from_cached_summary_and_ld.R", dst)
  a <- read.delim(file.path(dst, "cached_recheck/comparison_to_saved_results.tsv"), check.names = FALSE)
  if (!nrow(a) || !"passed" %in% names(a) || !all(a$passed) ||
      max(a$max_absolute_difference) >= 1e-7) stop("LD comparison failed.")
  cat("PASS cached LD analysis: all saved-result comparisons within 1e-7.\n")
}
if (mode %in% c("diagnostics", "all")) {
  run("revision/reporting/revision_diagnostics.R", root)
  files <- c("cis_recomputed.tsv", "cis_per_variant_complete.tsv", "cis_power_scenarios.tsv",
             "mvmr_recomputed_full_and_screened.tsv", "mvmr_expanded_hypothetical_covariance_grid.tsv",
             "robust_mvmr_exact_implementation.tsv", "robust_mvmr_variant_weights.tsv")
  for (name in files) check_numeric(file.path(root, "revision/reporting", name), file.path(root, "recomputed/reporting", name))
}
if (mode %in% c("coloc", "all")) {
  src <- file.path(root, "revision/multisignal_coloc")
  dst <- file.path(root, "recomputed/multisignal_coloc")
  stage(src, dst, c("matched_regional_summary_statistics.tsv", "signed_LD_EUR_503.rds", "LD_SNP_order.txt"))
  Sys.setenv(COLOC_USE_PREMATCHED = "1", COLOC_REVISION_OUT = dst)
  run("revision/multisignal_coloc/01_multisignal_coloc.R")
  run("revision/multisignal_coloc/02_validate_and_summarize.R")
  check_numeric(file.path(src, "coloc_abf_matched_sets_prior_sensitivity.tsv"),
                file.path(dst, "coloc_abf_matched_sets_prior_sensitivity.tsv"), c("analysis_set", "p12"))
  for (name in c("susie_fit_summary.tsv", "sepsis_rss_specification_sensitivity.tsv", "exposure_L20_capacity_sensitivity.tsv")) {
    check_numeric(file.path(src, name), file.path(dst, name))
  }
  cat("Interpretation remains: multi-signal PP.H4 is not estimable, not zero.\n")
}
cat("Completed mode: ", mode, ". Outputs: recomputed/ (not for Git publication).\n", sep = "")
