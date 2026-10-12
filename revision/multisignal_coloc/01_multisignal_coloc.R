suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})
options(stringsAsFactors = FALSE)
set.seed(20261012)

# Resolve repository-relative defaults independently of the working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript so its directory can be resolved.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
# Inputs are read-only; generated results are written to the selected directory.
out_dir <- normalizePath(Sys.getenv("COLOC_REVISION_OUT", unset = script_dir), mustWork = TRUE)
original_dir <- Sys.getenv("SST2_ORIGINAL_DIR", unset = file.path(repo_dir, "original"))
input_file <- file.path(original_dir, "results/colocalization/harmonized_region_Sepsis.tsv")
ref_prefix <- Sys.getenv("EUR_REFERENCE_PREFIX", unset = "")
plink <- Sys.getenv("PLINK_BIN", unset = Sys.which("plink"))
ref_out <- file.path(out_dir, "reference_region")
input_mode <- Sys.getenv("COLOC_USE_PREMATCHED", unset = "1")
if (!input_mode %in% c("0", "1")) stop("COLOC_USE_PREMATCHED must be 1 (cached summary/LD) or 0 (rebuild LD).")
use_prematched <- identical(input_mode, "1")
if (use_prematched) {
  x <- fread(file.path(out_dir, "matched_regional_summary_statistics.tsv"))
  R <- readRDS(file.path(out_dir, "signed_LD_EUR_503.rds"))
  reference_n <- 503L
  stopifnot(uniqueN(x$rsid) == nrow(x), identical(rownames(R), x$rsid),
            identical(colnames(R), x$rsid), all(x$pos38 == x$pos_outcome),
            max(abs(R - t(R))) < 1e-12, max(abs(diag(R) - 1)) < 1e-12)
} else {
if (!nzchar(ref_prefix) || !all(file.exists(paste0(ref_prefix, c(".bed", ".bim", ".fam")))))
  stop("Set EUR_REFERENCE_PREFIX to an available EUR .bed/.bim/.fam prefix for COLOC_USE_PREMATCHED=0.")
if (!nzchar(plink) || !file.exists(plink) || file.access(plink, mode = 1L) != 0L)
  stop("Set PLINK_BIN to an executable PLINK path, or make plink available on PATH.")
if (!file.exists(input_file)) stop("Regional summary input is missing: ", input_file)
x <- fread(input_file)
stopifnot(nrow(x) == uniqueN(x$rsid), all(x$pos37 >= 102427962), all(x$pos37 <= 103468497))
fwrite(x[, .(rsid)], file.path(out_dir, "regional_rsids.txt"), col.names = FALSE)

run_plink <- function(args) {
  status <- system2(plink, args = vapply(args, shQuote, ""))
  if (status != 0) stop("PLINK returned ", status)
}
if (!file.exists(paste0(ref_out, ".raw"))) {
  run_plink(c("--bfile", ref_prefix, "--extract", file.path(out_dir, "regional_rsids.txt"),
              "--keep-allele-order", "--make-bed", "--out", ref_out))
  run_plink(c("--bfile", ref_out, "--keep-allele-order", "--recode", "A", "--freq",
              "--out", ref_out))
}
bim <- fread(paste0(ref_out, ".bim"), col.names = c("chr", "rsid", "cm", "pos37_ref", "a1_ref", "a2_ref"))
geno <- fread(paste0(ref_out, ".raw"))
freq <- fread(paste0(ref_out, ".frq"))
stopifnot(nrow(bim) == ncol(geno) - 6L, nrow(bim) == uniqueN(bim$rsid))
# PLINK --recode A dosage is the count of the allele shown in the column header.
# keep-allele-order preserves source BIM A1; verify every header before using it.
stopifnot(identical(names(geno)[-(1:6)], paste0(bim$rsid, "_", bim$a1_ref)))
x <- merge(bim, x, by = "rsid", sort = FALSE)
x <- x[match(bim$rsid, rsid)]
stopifnot(identical(x$rsid, bim$rsid), all(x$pos37_ref == x$pos37))
x[, ld_sign := fifelse(effect_allele == a1_ref & other_allele == a2_ref, 1,
                       fifelse(effect_allele == a2_ref & other_allele == a1_ref, -1, NA_real_))]
stopifnot(!anyNA(x$ld_sign))
G <- as.matrix(geno[, -(1:6)])
storage.mode(G) <- "double"
reference_n <- nrow(G)
missing_counts <- colSums(is.na(G))
for (j in which(missing_counts > 0)) G[is.na(G[, j]), j] <- mean(G[, j], na.rm = TRUE)
reference_a1_freq <- colMeans(G) / 2
x[, eaf_reference := fifelse(ld_sign == 1, reference_a1_freq, 1 - reference_a1_freq)]
x[, maf_reference := pmin(eaf_reference, 1 - eaf_reference)]
x[, ref_exposure_freq_difference := abs(eaf_reference - eaf_exposure)]
x[, ref_outcome_freq_difference := abs(eaf_reference - eaf_outcome)]
x[, n_missing_reference := missing_counts]
x[, z_exposure := beta_exposure / se_exposure]
x[, z_outcome := beta_outcome / se_outcome]
x[, p_exposure := 2 * pnorm(-abs(z_exposure))]
monomorphic <- which(!is.finite(apply(G, 2, var)) | apply(G, 2, var) <= 0)
if (length(monomorphic)) {
  fwrite(x[monomorphic], file.path(out_dir, "excluded_monomorphic_reference.tsv"), sep = "\t")
  x <- x[-monomorphic]; G <- G[, -monomorphic, drop = FALSE]
}
Rraw <- cor(G)
R <- Rraw * outer(x$ld_sign, x$ld_sign)
dimnames(R) <- list(x$rsid, x$rsid)
stopifnot(max(abs(R - t(R))) < 1e-12, max(abs(diag(R) - 1)) < 1e-12,
          all(is.finite(R)), identical(rownames(R), x$rsid))
saveRDS(R, file.path(out_dir, "signed_LD_EUR_503.rds"), compress = "xz")
fwrite(x, file.path(out_dir, "matched_regional_summary_statistics.tsv"), sep = "\t")
writeLines(x$rsid, file.path(out_dir, "LD_SNP_order.txt"))
writeLines(c(
  paste("Input:", input_file), paste("Reference:", ref_prefix),
  paste("Reference individuals:", reference_n), paste("Harmonized input SNPs:", nrow(fread(input_file))),
  paste("Matched polymorphic reference SNPs:", nrow(x)),
  paste("Missing reference dosages (mean-imputed):", sum(missing_counts)),
  paste("LD coding: Pearson correlation of PLINK A1 dosages, multiplied by outer(ld_sign,ld_sign)."),
  paste("All effects and signed LD correspond to SCALLOP effect_allele; direct/swapped matching only."),
  paste("Genome build: GRCh37 in LD / SCALLOP; FinnGen coordinates were lifted and harmonized in the original analysis."),
  paste("Exposure sample size range:", paste(range(x$n_exposure), collapse = "–")),
  paste("Exposure median N:", median(x$n_exposure)),
  paste("Sepsis cases/controls:", 17133, "/", 439048),
  paste("Reference panel is broader European, not an in-sample/Finnish-specific LD panel.")),
  file.path(out_dir, "data_provenance.txt"))
}

# Analysis sets are specified before fitting: the original regional overlap,
# and a sensitivity set avoiding rare variants / marked per-variant N variation.
sets <- list(
  all_matched = rep(TRUE, nrow(x)),
  common_highN = pmin(x$eaf_exposure, 1 - x$eaf_exposure) >= 0.01 &
    pmin(x$eaf_outcome, 1 - x$eaf_outcome) >= 0.01 & x$maf_reference >= 0.01 &
    x$n_exposure >= 0.9 * max(x$n_exposure)
)
cases <- 17133; controls <- 439048
Nout <- cases + controls
Neffout <- 4 / (1/cases + 1/controls)

make_datasets <- function(xx, RR) list(
  exposure = list(beta = xx$beta_exposure, varbeta = xx$se_exposure^2,
                  N = round(median(xx$n_exposure)), type = "quant", sdY = 1,
                  snp = xx$rsid, position = xx$pos37, LD = RR),
  sepsis = list(beta = xx$beta_outcome, varbeta = xx$se_outcome^2,
               N = Nout, s = cases/Nout, type = "cc", snp = xx$rsid,
               position = xx$pos37, LD = RR)
)

extract_cs <- function(fit, trait, analysis_set, xx) {
  if (is.null(fit$sets$cs) || !length(fit$sets$cs)) return(data.table())
  rbindlist(lapply(seq_along(fit$sets$cs), function(k) {
    ids <- fit$sets$cs[[k]]
    effect <- fit$sets$cs_index[k]
    top <- ids[which.max(fit$alpha[effect, ids])]
    data.table(analysis_set, trait, cs = names(fit$sets$cs)[k], effect_index = effect,
      n_snps = length(ids), coverage = fit$sets$coverage[k],
      min_abs_corr = fit$sets$purity$min.abs.corr[k],
      mean_abs_corr = fit$sets$purity$mean.abs.corr[k],
      cs_log10bf = fit$lbf[effect] / log(10), top_snp = xx$rsid[top],
      top_pip = fit$pip[top], snps = paste(xx$rsid[ids], collapse = ";"))
  }), fill = TRUE)
}

all_diag <- list(); all_fit <- list(); all_cs <- list(); all_abf <- list(); all_status <- list(); all_coloc <- list()
for (analysis_set in names(sets)) {
  message("=== ANALYSIS SET: ", analysis_set, " ===")
  keep <- which(sets[[analysis_set]])
  xx <- x[keep]; RR <- R[keep, keep, drop = FALSE]
  eig <- eigen(RR, symmetric = TRUE)
  saveRDS(eig, file.path(out_dir, paste0("LD_eigendecomposition_", analysis_set, ".rds")), compress = FALSE)
  # Only numerical negative eigenvalues are truncated by susieR's r_tol;
  # no arbitrary ridge adjustment is applied to the fine-mapping LD matrix.
  attr(RR, "eigen") <- eig
  d <- make_datasets(xx, RR)
  for (p12 in c(1e-6, 1e-5, 1e-4)) {
    abf <- coloc.abf(d$exposure, d$sepsis, p1 = 1e-4, p2 = 1e-4, p12 = p12)
    s <- abf$summary
    all_abf[[paste(analysis_set, p12)]] <- data.table(analysis_set, nsnps = s[["nsnps"]],
      p1 = 1e-4, p2 = 1e-4, p12,
      PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
      PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]])
  }
  fits <- list()
  for (trait in c("sepsis", "exposure")) {
    dd <- d[[trait]]
    z <- dd$beta/sqrt(dd$varbeta)
    names(z) <- xx$rsid
    # Binary-outcome SuSiE uses logistic-Wald z with effective case-control N.
    # Ordinary coloc.abf above continues to use total N and actual case fraction.
    n_rss <- if (trait == "sepsis") Neffout else dd$N
    lambda <- estimate_s_rss(z, RR, n = n_rss, method = "null-mle")
    kr <- kriging_rss(z, RR, n = n_rss, s = lambda)
    krdt <- as.data.table(kr$conditional_dist)
    krdt[, rsid := xx$rsid]
    fwrite(krdt, file.path(out_dir, paste0("z_LD_diagnostic_", analysis_set, "_", trait, ".tsv")), sep = "\t")
    if (!is.null(kr$plot)) ggplot2::ggsave(file.path(out_dir, paste0("z_LD_diagnostic_", analysis_set, "_", trait, ".png")),
                                         kr$plot, width = 7, height = 6, dpi = 180)
    all_diag[[paste(analysis_set, trait)]] <- data.table(analysis_set, trait, nsnps = nrow(xx),
      n_rss, reference_n = reference_n, lambda_null_mle = lambda,
      min_eigenvalue = min(eig$values), rank_at_1e_8 = sum(eig$values > 1e-8),
      max_abs_z = max(abs(z)), min_p = min(2 * pnorm(-abs(z))),
      eaf_reference_max_difference = if (trait == "sepsis") max(xx$ref_outcome_freq_difference) else max(xx$ref_exposure_freq_difference))
    warnings <- character()
    # runsusie stops and discards the object on nonconvergence when bounded.
    # Use the same susie_rss call plus coloc's annotation so a failed fit is
    # preserved and explicitly excluded from signal-pair inference.
    fit <- withCallingHandlers(
      susie_rss(z = z, R = RR, n = n_rss, L = 10, coverage = 0.95,
                min_abs_corr = 0.5, estimate_residual_variance = FALSE,
                max_iter = 1000, tol = 1e-3, refine = FALSE),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    fit <- getFromNamespace("annotate_susie", "coloc")(fit, xx$rsid, RR)
    fits[[trait]] <- fit
    saveRDS(fit, file.path(out_dir, paste0("susie_", analysis_set, "_", trait, ".rds")), compress = "xz")
    cs <- extract_cs(fit, trait, analysis_set, xx)
    all_cs[[paste(analysis_set, trait)]] <- cs
    all_fit[[paste(analysis_set, trait)]] <- data.table(analysis_set, trait, nsnps = nrow(xx),
      n_rss, L = 10, coverage = 0.95, min_abs_corr = 0.5, estimate_residual_variance = FALSE,
      converged = fit$converged, niter = fit$niter, n_credible_sets = nrow(cs),
      final_elbo = tail(fit$elbo, 1), warnings = paste(unique(warnings), collapse = " | "))
    fwrite(data.table(rsid = xx$rsid, pos37 = xx$pos37, pip = fit$pip),
           file.path(out_dir, paste0("susie_PIP_", analysis_set, "_", trait, ".tsv")), sep = "\t")
    fwrite(rbindlist(all_fit, fill = TRUE), file.path(out_dir, "susie_fit_summary.tsv"), sep = "\t")
    fwrite(rbindlist(all_cs, fill = TRUE), file.path(out_dir, "susie_credible_sets.tsv"), sep = "\t")
    fwrite(rbindlist(all_diag, fill = TRUE), file.path(out_dir, "z_LD_diagnostic_summary.tsv"), sep = "\t")
    fwrite(rbindlist(all_abf, fill = TRUE), file.path(out_dir, "coloc_abf_matched_sets_prior_sensitivity.tsv"), sep = "\t")
  }
  ncs <- vapply(fits, function(f) length(f$sets$cs), 0L)
  converged_all <- all(vapply(fits, function(f) isTRUE(f$converged), FALSE))
  if (all(ncs > 0) && converged_all) {
    for (p12 in c(1e-6, 1e-5, 1e-4)) {
      cc <- coloc.susie(fits$exposure, fits$sepsis, p1 = 1e-4, p2 = 1e-4, p12 = p12)
      sm <- as.data.table(cc$summary)
      sm[, `:=`(analysis_set = analysis_set, p1 = 1e-4, p2 = 1e-4, p12 = p12)]
      all_coloc[[paste(analysis_set, p12)]] <- sm
      saveRDS(cc, file.path(out_dir, paste0("coloc_susie_", analysis_set, "_p12_", p12, ".rds")))
    }
    status <- "Signal pairs estimable; inspect compatibility diagnostics before interpretation."
  } else {
    reason <- character()
    if (any(ncs == 0)) reason <- c(reason, paste("no qualifying 95% credible set for", paste(names(ncs)[ncs == 0], collapse = " and ")))
    if (!converged_all) reason <- c(reason, paste("nonconverged fit for", paste(names(fits)[!vapply(fits, function(f) isTRUE(f$converged), FALSE)], collapse = " and ")))
    status <- paste("Not estimable:", paste(reason, collapse = "; "),
                    "; PP.H4 is not assigned zero. No weakening of credible-set coverage was used.")
  }
  all_status[[analysis_set]] <- data.table(analysis_set, nsnps = nrow(xx),
    exposure_credible_sets = ncs[["exposure"]], sepsis_credible_sets = ncs[["sepsis"]], status)
  fwrite(rbindlist(all_fit, fill = TRUE), file.path(out_dir, "susie_fit_summary.tsv"), sep = "\t")
  fwrite(rbindlist(all_cs, fill = TRUE), file.path(out_dir, "susie_credible_sets.tsv"), sep = "\t")
  fwrite(rbindlist(all_diag, fill = TRUE), file.path(out_dir, "z_LD_diagnostic_summary.tsv"), sep = "\t")
  fwrite(rbindlist(all_abf, fill = TRUE), file.path(out_dir, "coloc_abf_matched_sets_prior_sensitivity.tsv"), sep = "\t")
  fwrite(rbindlist(all_status, fill = TRUE), file.path(out_dir, "coloc_susie_status.tsv"), sep = "\t")
  if (length(all_coloc)) fwrite(rbindlist(all_coloc, fill = TRUE), file.path(out_dir, "coloc_susie_signal_pairs.tsv"), sep = "\t")
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
message("COMPLETE")
print(rbindlist(all_diag, fill = TRUE))
print(rbindlist(all_fit, fill = TRUE))
print(rbindlist(all_status, fill = TRUE))
print(rbindlist(all_abf, fill = TRUE))
