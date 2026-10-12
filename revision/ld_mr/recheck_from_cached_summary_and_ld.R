#!/usr/bin/env Rscript
# Portable numerical reproduction using derived summary statistics and signed LD.
# Requires base R only; no PLINK, individual genotypes, internet, or external path.
# Run from any working directory:
#   Rscript /path/to/ld_mr/recheck_from_cached_summary_and_ld.R
# Optional first argument supplies the directory containing the three inputs.
# This checks estimation from the cached LD; it does not regenerate the reference
# LD or independently reperform the original strict clumping. The four strict
# rsIDs are those recorded in the original strict-clumping output.
options(stringsAsFactors = FALSE, digits = 16)
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
input_dir <- normalizePath(if (length(args)) args[1] else script_dir, mustWork = TRUE)
out_dir <- file.path(input_dir, "cached_recheck")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
summary_file <- file.path(input_dir, "ld_order_harmonized_instruments.tsv")
ld_file <- file.path(input_dir, "signed_ld_exposure_effect_allele_aligned.rds")
strict_file <- file.path(input_dir, "strict_r2_0.001_variant_ids.txt")
stopifnot(all(file.exists(c(summary_file, ld_file, strict_file))))
d <- read.delim(summary_file, check.names = FALSE)
R <- readRDS(ld_file)
strict_ids <- trimws(readLines(strict_file))
strict_ids <- strict_ids[nzchar(strict_ids)]
required <- c("rsid", "beta_ukb_protein", "se_ukb_protein", "beta_ukb_sepsis", "se_ukb_sepsis")
stopifnot(all(required %in% names(d)), nrow(d) == 49L, !anyDuplicated(d$rsid),
          length(strict_ids) == 4L, !anyDuplicated(strict_ids), all(strict_ids %in% d$rsid),
          is.matrix(R), !is.null(rownames(R)), !is.null(colnames(R)),
          all(d$rsid %in% rownames(R)), all(d$rsid %in% colnames(R)))
R <- R[d$rsid, d$rsid, drop = FALSE]
stopifnot(!anyNA(R), max(abs(R - t(R))) < 1e-10,
          max(abs(diag(R) - 1)) < 1e-10,
          all(is.finite(as.matrix(d[required[-1]]))),
          all(d$se_ukb_sepsis > 0), all(d$se_ukb_protein > 0))
ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(ev) > 1e-8, max(ev) / min(ev) < 1e8)

ivw <- function(data, rho, analysis, model) {
  bx <- data$beta_ukb_protein; by <- data$beta_ukb_sepsis; sy <- data$se_ukb_sepsis
  omega <- rho * outer(sy, sy)
  denominator <- drop(crossprod(bx, solve(omega, bx)))
  beta <- drop(crossprod(bx, solve(omega, by))) / denominator
  residual <- by - beta * bx
  Q <- drop(crossprod(residual, solve(omega, residual)))
  df <- nrow(data) - 1L
  phi <- max(1, Q / df)
  se_fixed <- sqrt(1 / denominator)
  se <- se_fixed * ifelse(model == "random", sqrt(phi), 1)
  data.frame(analysis = analysis, model = model, nsnp = nrow(data), beta = beta, se = se,
             OR = exp(beta), CI_lower = exp(beta - qnorm(.975) * se),
             CI_upper = exp(beta + qnorm(.975) * se), p = 2 * pnorm(-abs(beta / se)),
             Q = Q, Q_df = df, Q_p = pchisq(Q, df, lower.tail = FALSE),
             phi = phi, se_fixed = se_fixed)
}
strict <- d[d$rsid %in% strict_ids, , drop = FALSE]
strict_R <- R[strict$rsid, strict$rsid, drop = FALSE]
estimates <- list()
for (model in c("fixed", "random")) {
  estimates[[length(estimates) + 1L]] <- ivw(d, diag(nrow(d)), "Published_49_naive_IVW", model)
  estimates[[length(estimates) + 1L]] <- ivw(d, diag(nrow(d)), "LD_available_same_variants_naive_IVW", model)
  estimates[[length(estimates) + 1L]] <- ivw(d, R, "LD_available_same_variants_generalized_IVW", model)
  estimates[[length(estimates) + 1L]] <- ivw(strict, diag(nrow(strict)), "Strict_r2_0.001_4_variant_naive_IVW", model)
  estimates[[length(estimates) + 1L]] <- ivw(strict, strict_R, "Strict_r2_0.001_4_variant_generalized_IVW", model)
}
estimates <- do.call(rbind, estimates)
write.table(estimates, file.path(out_dir, "ivw_estimates_fixed_and_random.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Weighted PCs of Psi = diag(bX/SEY) R diag(bX/SEY), fixed 99% and 99.9% thresholds.
# All variants contribute to the projection; no significance-based PC selection.
psi <- R * outer(d$beta_ukb_protein / d$se_ukb_sepsis,
                 d$beta_ukb_protein / d$se_ukb_sepsis)
eig <- eigen(psi, symmetric = TRUE)
eig$values[eig$values < 0 & abs(eig$values) < 1e-8] <- 0
stopifnot(all(eig$values >= 0))
omega <- R * outer(d$se_ukb_sepsis, d$se_ukb_sepsis)
pca <- list()
for (threshold in c(.99, .999)) {
  K <- which(cumsum(eig$values) / sum(eig$values) >= threshold)[1]
  U <- eig$vectors[, seq_len(K), drop = FALSE]
  bx <- drop(crossprod(U, d$beta_ukb_protein))
  by <- drop(crossprod(U, d$beta_ukb_sepsis))
  covariance <- crossprod(U, omega %*% U)
  denominator <- drop(crossprod(bx, solve(covariance, bx)))
  beta <- drop(crossprod(bx, solve(covariance, by))) / denominator
  residual <- by - beta * bx
  Q <- drop(crossprod(residual, solve(covariance, residual)))
  phi <- max(1, Q / (K - 1L))
  for (model in c("fixed", "random")) {
    se <- sqrt(1 / denominator) * ifelse(model == "random", sqrt(phi), 1)
    pca[[length(pca) + 1L]] <- data.frame(
      threshold = threshold, model = model, source_nsnp = nrow(d), n_pc = K,
      beta = beta, se = se, OR = exp(beta), CI_lower = exp(beta - qnorm(.975) * se),
      CI_upper = exp(beta + qnorm(.975) * se), p = 2 * pnorm(-abs(beta / se)),
      Q = Q, Q_df = K - 1L, Q_p = pchisq(Q, K - 1L, lower.tail = FALSE), phi = phi)
  }
}
pca <- do.call(rbind, pca)
write.table(pca, file.path(out_dir, "weighted_pca_ivw_sensitivity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

check_against_saved <- function(actual, expected_name, keys, metrics) {
  expected_file <- file.path(input_dir, expected_name)
  if (!file.exists(expected_file)) {
    return(data.frame(file = expected_name, metric = "expected_file_absent",
                      max_absolute_difference = NA_real_, passed = NA))
  }
  expected <- read.delim(expected_file, check.names = FALSE)
  stopifnot(nrow(actual) == nrow(expected))
  actual_key <- do.call(paste, c(actual[keys], sep = "|"))
  expected_key <- do.call(paste, c(expected[keys], sep = "|"))
  stopifnot(!anyDuplicated(actual_key), !anyDuplicated(expected_key),
            setequal(actual_key, expected_key))
  expected <- expected[match(actual_key, expected_key), , drop = FALSE]
  checks <- do.call(rbind, lapply(metrics, function(metric) {
    delta <- max(abs(actual[[metric]] - expected[[metric]]))
    data.frame(file = expected_name, metric = metric,
               max_absolute_difference = delta, passed = delta < 1e-7)
  }))
  stopifnot(all(checks$passed))
  checks
}
ivw_metrics <- c("nsnp", "beta", "se", "OR", "CI_lower", "CI_upper", "p", "Q", "Q_df", "Q_p", "phi", "se_fixed")
pca_metrics <- c("source_nsnp", "n_pc", "beta", "se", "OR", "CI_lower", "CI_upper", "p", "Q", "Q_df", "Q_p", "phi")
checks <- rbind(
  check_against_saved(estimates, "ivw_estimates_fixed_and_random.tsv", c("analysis", "model"), ivw_metrics),
  check_against_saved(pca, "weighted_pca_ivw_sensitivity.tsv", c("threshold", "model"), pca_metrics))
write.table(checks, file.path(out_dir, "comparison_to_saved_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
input_checksums <- data.frame(
  file = basename(c(summary_file, ld_file, strict_file)),
  md5 = as.character(tools::md5sum(c(summary_file, ld_file, strict_file))))
write.table(input_checksums, file.path(out_dir, "cached_input_checksums.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
writeLines(c(
  "Portable recheck used only cached SNP summary effects, signed LD, and four recorded strict-clumping rsIDs.",
  "No individual-level genotypes, PLINK, internet, or absolute external input paths are required.",
  "It validates numerical estimation from the supplied LD, not source-genotype LD generation or original clumping.",
  "First-order generalized IVW: Omega=diag(SEY) R diag(SEY); normal Wald intervals; multiplicative random dispersion=max(1,Q/(K-1)).",
  "All full and strict ordinary/generalized fixed/random results and both weighted-PC thresholds are retained.",
  sprintf("Maximum absolute difference from saved numeric results: %.12g; required tolerance <1e-7.", max(checks$max_absolute_difference, na.rm = TRUE)),
  "Limitations of source assays, external LD, and variant selection remain those documented in the manuscript."),
  file.path(out_dir, "README.txt"))
print(estimates)
print(pca)
print(checks)
cat("Portable cached-data numerical recheck completed.\n")
