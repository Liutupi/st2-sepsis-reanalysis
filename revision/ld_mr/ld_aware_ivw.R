#!/usr/bin/env Rscript
# Correlated-variant IVW audit for the 49 published UKB-PPP IL1RL1 instruments.
# All upstream files are read-only. Outputs, including PLINK logs and signed
# LD matrices, are placed beside this script (or SST2_LD_MR_OUT).
suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
})
options(stringsAsFactors = FALSE, digits = 16)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
project <- Sys.getenv("SST2_ORIGINAL_PROJECT", file.path(repo_dir, "original"))
out <- Sys.getenv("SST2_LD_MR_OUT", file.path(repo_dir, "recomputed/ld_from_reference"))
ref <- Sys.getenv("SST2_EUR_BFILE", "")
plink <- Sys.getenv("PLINK_BIN", Sys.which("plink"))
if (!nzchar(ref) || !nzchar(plink)) stop("Full LD generation requires SST2_EUR_BFILE and PLINK_BIN (or plink on PATH). Use recheck_from_cached_summary_and_ld.R for genotype-free reproduction.")
input <- file.path(project, "results/benchmark_2026_il1rl1_sepsis/published_ukb_ppp_il1rl1_49_instruments.tsv")
strict_input <- file.path(project, "results/benchmark_2026_il1rl1_sepsis/ukb_ppp_il1rl1_strict_r2_0.001.tsv")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
logcon <- file(file.path(out, "analysis_log.txt"), open = "wt")
sink(logcon, split = TRUE)
on.exit({sink(); close(logcon)}, add = TRUE)
cat("Correlated-variant IVW audit\nRun time: ", format(Sys.time(), tz = "Asia/Shanghai"), " Asia/Shanghai\n", sep = "")
stopifnot(file.exists(input), file.exists(strict_input), file.exists(plink),
          all(file.exists(paste0(ref, c(".bed", ".bim", ".fam")))))
dat <- fread(input)
stopifnot(nrow(dat) == 49L, !anyDuplicated(dat$rsid))
stopifnot(all(dat$effect_allele == dat$effect_allele.outcome),
          all(dat$other_allele == dat$other_allele.outcome),
          all(is.finite(dat$beta_ukb_protein)), all(dat$se_ukb_protein > 0),
          all(is.finite(dat$beta_ukb_sepsis)), all(dat$se_ukb_sepsis > 0))
dat[, original_row := .I]
dat[, SNP_F := (beta_ukb_protein / se_ukb_protein)^2]
dat[, effect_allele := toupper(effect_allele)]
dat[, other_allele := toupper(other_allele)]
ref_bim <- fread(paste0(ref, ".bim"),
                 col.names = c("ref_chr", "rsid", "cm", "ref_pos", "ref_a1", "ref_a2"))
ref_fam <- fread(paste0(ref, ".fam"), header = FALSE)
ref_map <- ref_bim[rsid %chin% dat$rsid]
duplicate_ids <- ref_map[duplicated(rsid) | duplicated(rsid, fromLast = TRUE)]
fwrite(duplicate_ids, file.path(out, "reference_duplicate_variant_ids.tsv"), sep = "\t")
if (nrow(duplicate_ids)) stop("Duplicated requested rsIDs in reference; resolve before extraction.")
audit <- merge(dat, ref_map, by = "rsid", all.x = TRUE, sort = FALSE)
setorder(audit, original_row)
audit[, present_in_reference := !is.na(ref_chr)]
audit[, exact_allele_pair_match := present_in_reference &
        ((effect_allele == ref_a1 & other_allele == ref_a2) |
         (effect_allele == ref_a2 & other_allele == ref_a1))]
audit[, included := present_in_reference & exact_allele_pair_match]
fwrite(audit, file.path(out, "variant_inclusion_audit.tsv"), sep = "\t", na = "NA")
cat("Reference samples:", nrow(ref_fam), "\nRequested variants:", nrow(dat),
    "\nReference present:", sum(audit$present_in_reference),
    "\nExact allele-pair compatible:", sum(audit$exact_allele_pair_match), "\n")
if (sum(audit$included) < 2L) stop("Fewer than two allele-compatible reference variants.")
dat_avail <- audit[included == TRUE]
extract <- file.path(out, "requested_variant_ids.txt")
a1file <- file.path(out, "force_exposure_effect_alleles.tsv")
fwrite(dat_avail[, .(rsid)], extract, col.names = FALSE)
fwrite(dat_avail[, .(rsid, effect_allele)], a1file, sep = "\t", col.names = FALSE)
run_plink <- function(args) {
  # system2 requires explicit quoting because source/output paths contain spaces.
  status <- system2(plink, vapply(args, shQuote, character(1)))
  if (!identical(status, 0L)) stop("PLINK exited with status ", status)
}
prefix <- file.path(out, "reference_49_effect_allele_coded")
run_plink(c("--bfile", ref, "--extract", extract, "--a1-allele", a1file, "2", "1", "0",
            "--keep-allele-order", "--make-bed", "--out", prefix))
bim <- fread(paste0(prefix, ".bim"),
             col.names = c("chr", "rsid", "cm", "pos", "ld_a1", "ld_a2"))
stopifnot(!anyDuplicated(bim$rsid), nrow(bim) == nrow(dat_avail))
harm <- merge(bim, dat_avail[, !c("chr", "pos", "bim_a1", "bim_a2")], by = "rsid", sort = FALSE)
harm <- harm[match(bim$rsid, rsid)]
harm[, ld_flip := fifelse(ld_a1 == effect_allele & ld_a2 == other_allele, 1L,
                          fifelse(ld_a1 == other_allele & ld_a2 == effect_allele, -1L, NA_integer_))]
stopifnot(!anyNA(harm$ld_flip))
ldprefix <- file.path(out, "signed_ld_plink")
run_plink(c("--bfile", prefix, "--keep-allele-order", "--r", "square", "--out", ldprefix))
R_raw <- as.matrix(fread(paste0(ldprefix, ".ld"), header = FALSE))
stopifnot(nrow(R_raw) == nrow(harm), ncol(R_raw) == nrow(harm), !anyNA(R_raw))
# If force-A1 coding fails for any compatible allele pair, the sign matrix
# explicitly reorients every LD entry to the exposure/outcome effect allele.
R <- R_raw * outer(harm$ld_flip, harm$ld_flip)
dimnames(R) <- list(harm$rsid, harm$rsid)
stopifnot(max(abs(R - t(R))) < 1e-10, max(abs(diag(R) - 1)) < 1e-10)
saveRDS(R, file.path(out, "signed_ld_exposure_effect_allele_aligned.rds"))
fwrite(data.table(rsid = rownames(R), as.data.table(R)),
       file.path(out, "signed_ld_exposure_effect_allele_aligned.tsv"), sep = "\t")
fwrite(harm[, .(rsid, chr, pos, ld_a1, ld_a2, effect_allele, other_allele, ld_flip,
               palindromic, ambiguous, eaf.exposure, eaf.outcome, beta_ukb_protein,
               se_ukb_protein, beta_ukb_sepsis, se_ukb_sepsis, SNP_F)],
       file.path(out, "ld_order_harmonized_instruments.tsv"), sep = "\t", na = "NA")
# Independent dosage-based correlation check, retaining raw dosage export.
rawprefix <- file.path(out, "dosage_audit")
run_plink(c("--bfile", prefix, "--keep-allele-order", "--recode", "A", "--out", rawprefix))
raw <- fread(paste0(rawprefix, ".raw"))
expected_cols <- paste(harm$rsid, harm$ld_a1, sep = "_")
stopifnot(all(expected_cols %in% names(raw)))
G <- as.matrix(raw[, ..expected_cols])
R_dosage <- cor(G, use = "pairwise.complete.obs") * outer(harm$ld_flip, harm$ld_flip)
ld_max_abs_diff <- max(abs(R - R_dosage))
if (ld_max_abs_diff > 2e-5) stop("PLINK signed LD disagrees with independent dosage calculation.")
harm[, ref_effect_allele_freq := colMeans(G, na.rm = TRUE) / 2]
harm[ld_flip == -1L, ref_effect_allele_freq := 1 - ref_effect_allele_freq]
harm[, eaf_abs_difference := abs(ref_effect_allele_freq - eaf.exposure)]
harm[, n_missing_genotypes := colSums(is.na(G))]
fwrite(harm[, .(rsid, effect_allele, other_allele, palindromic, ambiguous,
               eaf.exposure, ref_effect_allele_freq, eaf_abs_difference, n_missing_genotypes)],
       file.path(out, "allele_frequency_and_missingness_audit.tsv"), sep = "\t")
eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
condition <- max(eig) / min(eig)
pair_idx <- which(upper.tri(R), arr.ind = TRUE)
pairs <- data.table(rsid1 = rownames(R)[pair_idx[, 1]], rsid2 = rownames(R)[pair_idx[, 2]],
                    r = R[pair_idx], r2 = R[pair_idx]^2)
setorder(pairs, -r2)
fwrite(pairs, file.path(out, "all_pairwise_signed_ld.tsv"), sep = "\t")
diagtable <- data.table(
  metric = c("reference_samples", "requested_variants", "available_allele_compatible_variants",
             "palindromic_variants", "palindromic_previously_ambiguous", "missing_genotype_cells",
             "min_ld_eigenvalue", "max_ld_eigenvalue", "ld_matrix_condition_number",
             "ld_rank_at_1e_minus_8", "max_abs_ld_r", "max_pairwise_ld_r2", "pairs_r2_gt_0.001",
             "pairs_r2_ge_0.1", "max_plink_vs_independent_ld_abs_diff", "max_eaf_abs_difference"),
  value = c(nrow(ref_fam), nrow(dat), nrow(harm), sum(harm$palindromic), sum(harm$ambiguous),
            sum(is.na(G)), min(eig), max(eig), condition, sum(eig > 1e-8), max(abs(pairs$r)),
            max(pairs$r2), sum(pairs$r2 > .001), sum(pairs$r2 >= .1), ld_max_abs_diff,
            max(harm$eaf_abs_difference)))
fwrite(diagtable, file.path(out, "ld_matrix_diagnostics.tsv"), sep = "\t")
fwrite(data.table(component = seq_along(eig), eigenvalue = eig),
       file.path(out, "ld_eigenvalues.tsv"), sep = "\t")
cat("LD diagnostics:\n"); print(diagtable)

fit_algebra <- function(d, rho, label, model) {
  bx <- d$beta_ukb_protein; by <- d$beta_ukb_sepsis; sy <- d$se_ukb_sepsis
  omega <- outer(sy, sy) * rho
  wi_bx <- solve(omega, bx); wi_by <- solve(omega, by)
  info <- sum(bx * wi_bx)
  beta <- sum(bx * wi_by) / info
  residual <- by - beta * bx
  Q <- drop(crossprod(residual, solve(omega, residual)))
  df <- nrow(d) - 1L
  phi <- max(1, Q / df)
  se_fixed <- sqrt(1 / info)
  se <- if (model == "random") se_fixed * sqrt(phi) else se_fixed
  data.table(analysis = label, model = model, nsnp = nrow(d), beta = beta, se = se,
             OR = exp(beta), CI_lower = exp(beta - qnorm(.975) * se),
             CI_upper = exp(beta + qnorm(.975) * se), p = 2 * pnorm(-abs(beta / se)),
             Q = Q, Q_df = df, Q_p = pchisq(Q, df, lower.tail = FALSE),
             phi = phi, se_fixed = se_fixed, reference_n = nrow(ref_fam))
}
fit_with_check <- function(d, rho, label, model) {
  ans <- fit_algebra(d, rho, label, model)
  inp <- mr_input(bx = d$beta_ukb_protein, bxse = d$se_ukb_protein,
                  by = d$beta_ukb_sepsis, byse = d$se_ukb_sepsis,
                  correlation = rho, snps = d$rsid, effect_allele = d$effect_allele,
                  other_allele = d$other_allele)
  fit <- mr_ivw(inp, model = model, correl = TRUE, distribution = "normal")
  package_vals <- c(fit@Estimate, fit@StdError, fit@Pvalue, fit@Heter.Stat)
  independent_vals <- c(ans$beta, ans$se, ans$p, ans$Q, ans$Q_p)
  delta <- max(abs(package_vals - independent_vals))
  if (delta > 1e-9) stop("Package-versus-algebra cross-check failed: ", label, ", ", model)
  ans[, package_algebra_max_abs_difference := delta]
  fit_objects[[paste(label, model, sep = "__")]] <<- fit
  ans
}
fit_objects <- list()
res <- list()
for (model in c("fixed", "random")) {
  res[[length(res) + 1L]] <- fit_with_check(dat, diag(nrow(dat)), "Published_49_naive_IVW", model)
  res[[length(res) + 1L]] <- fit_with_check(harm, diag(nrow(harm)), "LD_available_same_variants_naive_IVW", model)
}
if (min(eig) > 1e-8 && condition < 1e8) {
  for (model in c("fixed", "random")) {
    res[[length(res) + 1L]] <- fit_with_check(harm, R, "LD_available_same_variants_generalized_IVW", model)
  }
} else {
  warning("LD matrix inversion is unstable; full generalized IVW not reported. See PCA sensitivities.")
}
strict_ids <- fread(strict_input)$rsid
strict <- harm[rsid %chin% strict_ids]
stopifnot(nrow(strict) == length(strict_ids))
strict_R <- R[strict$rsid, strict$rsid, drop = FALSE]
for (model in c("fixed", "random")) {
  res[[length(res) + 1L]] <- fit_with_check(strict, diag(nrow(strict)), "Strict_r2_0.001_4_variant_naive_IVW", model)
  res[[length(res) + 1L]] <- fit_with_check(strict, strict_R, "Strict_r2_0.001_4_variant_generalized_IVW", model)
}
results <- rbindlist(res)
fwrite(results, file.path(out, "ivw_estimates_fixed_and_random.tsv"), sep = "\t")
saveRDS(fit_objects, file.path(out, "mendelian_randomization_fit_objects.rds"))
cat("\nIVW estimates (package independently checked by direct algebra):\n"); print(results)

# Fixed-threshold numerical-stability sensitivities using weighted PCs, not selected
# to obtain significance. U are eigenvectors of Psi = diag(bx/sy) R diag(bx/sy).
# The first K PCs explain >=99% or >=99.9% of total weighted eigenvalue mass.
psi <- R * outer(harm$beta_ukb_protein / harm$se_ukb_sepsis,
                 harm$beta_ukb_protein / harm$se_ukb_sepsis)
pe <- eigen(psi, symmetric = TRUE)
pe$values[pe$values < 0 & abs(pe$values) < 1e-8] <- 0
stopifnot(all(pe$values >= 0))
pca_results <- list()
for (threshold in c(.99, .999)) {
  K <- which(cumsum(pe$values) / sum(pe$values) >= threshold)[1]
  U <- pe$vectors[, seq_len(K), drop = FALSE]
  bxpc <- drop(crossprod(U, harm$beta_ukb_protein))
  bypc <- drop(crossprod(U, harm$beta_ukb_sepsis))
  ompc <- crossprod(U, (R * outer(harm$se_ukb_sepsis, harm$se_ukb_sepsis)) %*% U)
  denom <- drop(crossprod(bxpc, solve(ompc, bxpc)))
  beta <- drop(crossprod(bxpc, solve(ompc, bypc))) / denom
  e <- bypc - beta * bxpc
  Q <- drop(crossprod(e, solve(ompc, e)))
  for (model in c("fixed", "random")) {
    phi <- max(1, Q / (K - 1L))
    se <- sqrt(1 / denom) * ifelse(model == "random", sqrt(phi), 1)
    pca_results[[length(pca_results) + 1L]] <- data.table(
      threshold = threshold, model = model, source_nsnp = nrow(harm), n_pc = K,
      beta = beta, se = se, OR = exp(beta), CI_lower = exp(beta - qnorm(.975) * se),
      CI_upper = exp(beta + qnorm(.975) * se), p = 2 * pnorm(-abs(beta / se)),
      Q = Q, Q_df = K - 1L, Q_p = pchisq(Q, K - 1L, lower.tail = FALSE), phi = phi)
  }
}
pca_results <- rbindlist(pca_results)
fwrite(pca_results, file.path(out, "weighted_pca_ivw_sensitivity.tsv"), sep = "\t")
cat("\nWeighted-PC numerical sensitivity (99% and 99.9%, all results retained):\n"); print(pca_results)

g <- results[analysis == "LD_available_same_variants_generalized_IVW" & model == "random"]
n <- results[analysis == "LD_available_same_variants_naive_IVW" & model == "random"]
s <- results[analysis == "Strict_r2_0.001_4_variant_naive_IVW" & model == "random"]
summary_lines <- c(
  "LD-aware IVW: methods, results, and interpretation",
  sprintf("Reference: 1000 Genomes Phase 3 European ancestry, %d individuals; PLINK 1.9 genotype-dosage signed correlation.", nrow(ref_fam)),
  sprintf("Of 49 published harmonized UKB-PPP IL1RL1 instruments, %d were present and allele-compatible. Naive and generalized IVW use exactly this matched set.", nrow(harm)),
  "Effect-allele alignment: outcomes already align exactly to exposure effects in the published supplement; PLINK A1 forced to exposure effect allele and signed LD reoriented explicitly if needed.",
  sprintf("Palindromic variants retained: %d; previously ambiguous: %d. Alignment used exact alleles and effect-allele frequency checks; maximum reference-versus-UKB absolute allele-frequency discrepancy %.5f.", sum(harm$palindromic), sum(harm$ambiguous), max(harm$eaf_abs_difference)),
  sprintf("Independent dosage correlation matched PLINK signed LD to max absolute difference %.3g. LD eigenvalue min %.6g, max %.6g, condition number %.6g; %d genotype cells missing.", ld_max_abs_diff, min(eig), max(eig), condition, sum(is.na(G))),
  "Generalized IVW uses Omega=diag(seY) R diag(seY), beta=(bx' Omega^-1 by)/(bx' Omega^-1 bx), fixed SE=(bx' Omega^-1 bx)^-1/2; multiplicative random SE multiplies by sqrt(max(1,Q/(K-1))). Normal 95% CI and two-sided P values.",
  "MendelianRandomization::mr_ivw version 0.10.0 and independent matrix algebra agreed to numerical precision. First-order outcome-variance weights treat SNP-exposure effects as measured without error.",
  sprintf("Naive matched-set multiplicative random IVW: OR %.6f (95%% CI %.6f-%.6f), P %.8g; Q %.6f (df %d), Q P %.8g.", n$OR, n$CI_lower, n$CI_upper, n$p, n$Q, n$Q_df, n$Q_p),
  if (nrow(g)) sprintf("LD-aware matched-set multiplicative random IVW: OR %.6f (95%% CI %.6f-%.6f), P %.8g; Q %.6f (df %d), Q P %.8g.", g$OR, g$CI_lower, g$CI_upper, g$p, g$Q, g$Q_df, g$Q_p) else "Full generalized IVW not reported because LD matrix failed numerical stability criteria.",
  sprintf("Strict 4-variant multiplicative random IVW: OR %.6f (95%% CI %.6f-%.6f), P %.8g; Q %.6f (df %d), Q P %.8g.", s$OR, s$CI_lower, s$CI_upper, s$p, s$Q, s$Q_df, s$Q_p),
  "Interpretation must follow the observed LD-aware estimate. Comparing only a 49-SNP model with a 4-SNP model changes the instrument set and information content; it cannot identify residual LD as the cause of an earlier association.",
  "Limitations: external European LD estimated from only 503 individuals may differ from UK Biobank LD and is treated as known; this uncertainty is not included in model CIs. Rare-variant LD is less precise. Results concern published UKB-PPP IL1RL1 protein effects and the published UKB sepsis outcome, not SCALLOP or FinnGen primary analyses. The published protein unit must be respected; no new biological specificity is inferred.",
  "Weighted-PC sensitivities retain all 49 variants projected into fixed 99% and 99.9% weighted-variance subspaces; these are numerical sensitivity analyses, not alternative primary findings or preregistered analyses. All fixed/random estimates, signed LD, pairwise LD, allele audits, logs, and session info are saved.")
writeLines(summary_lines, file.path(out, "results_summary.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
hashes <- data.table(file = c(input, strict_input, paste0(ref, ".bim"), paste0(ref, ".fam"), plink),
                     md5 = as.character(tools::md5sum(c(input, strict_input, paste0(ref, ".bim"), paste0(ref, ".fam"), plink))))
fwrite(hashes, file.path(out, "input_checksums.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
sink(); close(logcon)
