suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(MVMR)
  library(MendelianRandomization)
})

options(stringsAsFactors = FALSE)
set.seed(20260717)

project_dir <- normalizePath(Sys.getenv("SST2_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
sst2_file <- Sys.getenv("SST2_SCALLOP_FILE", unset = file.path(project_dir, "data", "GCST90012040.tsv.gz"))
hermes_file <- Sys.getenv("SST2_HERMES_FILE", unset = file.path(project_dir, "data", "GCST009541.tsv.gz"))
sepsis_file <- Sys.getenv("SST2_FINNGEN_SEPSIS_FILE", unset = file.path(project_dir, "data", "finngen_R12_AB1_OTHER_SEPSIS.gz"))
eur_bfile <- Sys.getenv("SST2_EUR_BFILE", unset = file.path(project_dir, "data", "1kg_eur", "EUR"))
plink_bin <- Sys.getenv("PLINK_BIN", unset = Sys.which("plink"))
if (!nzchar(plink_bin)) plink_bin <- file.path(project_dir, "tools", "plink")
chain_file <- Sys.getenv("SST2_CHAIN_FILE", unset = file.path(project_dir, "data_metadata/hg19ToHg38.over.chain"))
out_dir <- file.path(project_dir, "results/multivariable_mr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sst2_file)) stop("Missing SCALLOP file. Set SST2_SCALLOP_FILE.")
if (!file.exists(hermes_file)) stop("Missing HERMES file. Set SST2_HERMES_FILE.")
if (!file.exists(sepsis_file)) stop("Missing FinnGen sepsis file. Set SST2_FINNGEN_SEPSIS_FILE.")
if (!file.exists(paste0(eur_bfile, ".bed"))) stop("Missing 1000 Genomes EUR prefix. Set SST2_EUR_BFILE.")
if (!file.exists(plink_bin)) stop("Missing PLINK 1.9 executable. Set PLINK_BIN.")
if (!file.exists(chain_file)) stop("Missing hg19-to-hg38 chain file.")

message("Reading genome-wide significant sST2 and heart-failure candidates...")
sst2_candidates <- fread(file.path(project_dir, "results/instruments/sst2_gws_candidates.tsv"))
sst2_candidates <- sst2_candidates[, .(rsid, p_sst2 = pval)]

awk_hf <- 'BEGIN{FS=OFS="\\t"} NR==1 || $9<5e-8 {print}'
hf_cmd <- sprintf("gzip -dc %s | awk %s", shQuote(hermes_file), shQuote(awk_hf))
hf_candidates <- fread(cmd = hf_cmd)
hf_candidates[, A1 := toupper(A1)]
hf_candidates[, A2 := toupper(A2)]

bim <- fread(
  paste0(eur_bfile, ".bim"),
  col.names = c("chr", "rsid", "cm", "pos", "bim_a1", "bim_a2")
)
bim <- bim[rsid != "." & grepl("^rs[0-9]+$", rsid)]
bim[, bim_a1 := toupper(bim_a1)]
bim[, bim_a2 := toupper(bim_a2)]
bim <- unique(bim, by = "rsid")
hf_candidates <- hf_candidates[SNP %chin% bim$rsid]
hf_candidates <- hf_candidates[, .(rsid = SNP, p_hf = p)]

candidate_union <- merge(sst2_candidates, hf_candidates, by = "rsid", all = TRUE)
candidate_union[, p_joint := pmin(p_sst2, p_hf, na.rm = TRUE)]
candidate_union <- candidate_union[is.finite(p_joint) & p_joint < 5e-8]
candidate_union[, p_joint := pmax(p_joint, 1e-300)]

clump_input <- file.path(out_dir, "joint_sst2_hf_clump_input.tsv")
fwrite(candidate_union[, .(SNP = rsid, P = p_joint)], clump_input, sep = "\t")
clump_prefix <- file.path(out_dir, "joint_sst2_hf_r2_0.001")
status <- system2(plink_bin, args = vapply(c(
  "--bfile", eur_bfile, "--clump", clump_input,
  "--clump-snp-field", "SNP", "--clump-field", "P",
  "--clump-p1", "5e-8", "--clump-p2", "5e-8",
  "--clump-r2", "0.001", "--clump-kb", "10000",
  "--out", clump_prefix
), shQuote, character(1)))
if (!identical(status, 0L)) stop("Joint PLINK clumping failed")
clumped <- fread(paste0(clump_prefix, ".clumped"), fill = TRUE)
selected_rsids <- clumped[!is.na(SNP) & SNP != "", SNP]
selected_map <- bim[rsid %chin% selected_rsids]
if (nrow(selected_map) < 3L) stop("Too few jointly clumped variants")

message("Extracting cross-associations for both exposures...")
keys <- paste(selected_map$chr, selected_map$pos, sep = ":")
awk_sst2 <- paste0(
  'BEGIN{FS=OFS="\\t"; n=split(keys,k,","); for(i=1;i<=n;i++) wanted[k[i]]=1} ',
  'NR==1{print; next} {split($1,a,":"); key=a[1] ":" a[2]; if(key in wanted) print}'
)
sst2_cmd <- paste(
  "gzip -dc", shQuote(sst2_file), "| awk -v",
  shQuote(paste0("keys=", paste(keys, collapse = ","))), shQuote(awk_sst2)
)
sst2 <- fread(cmd = sst2_cmd)
sst2[, c("chr_text", "pos_text", "marker_ref_alt") := tstrsplit(MarkerName, ":", fixed = TRUE)]
sst2[, chr := as.integer(chr_text)]
sst2[, pos := as.integer(pos_text)]
sst2[, ea := toupper(Allele1)]
sst2[, oa := toupper(Allele2)]
sst2[, allele_key := paste(pmin(ea, oa), pmax(ea, oa), sep = "/")]
selected_map[, allele_key := paste(pmin(bim_a1, bim_a2), pmax(bim_a1, bim_a2), sep = "/")]
sst2 <- merge(
  sst2, selected_map[, .(rsid, chr, pos, allele_key)],
  by = c("chr", "pos", "allele_key"), all = FALSE, allow.cartesian = TRUE
)
sst2[, beta_sst2 := as.numeric(Effect)]
sst2[, se_sst2 := as.numeric(StdErr)]
sst2[, p_sst2_assoc := suppressWarnings(as.numeric(`P-value`))]
sst2[, eaf_sst2 := as.numeric(Freq1)]
sst2[, n_sst2 := as.numeric(TotalSampleSize)]
setorder(sst2, rsid, se_sst2)
sst2 <- unique(sst2, by = "rsid")

awk_by_rsid <- function(file, rsids) {
  awk <- paste0(
    'BEGIN{FS=OFS="\\t"; n=split(ids,k,","); for(i=1;i<=n;i++) wanted[k[i]]=1} ',
    'NR==1 || ($1 in wanted)'
  )
  cmd <- paste(
    "gzip -dc", shQuote(file), "| awk -v",
    shQuote(paste0("ids=", paste(rsids, collapse = ","))), shQuote(awk)
  )
  fread(cmd = cmd)
}
hf <- awk_by_rsid(hermes_file, sst2$rsid)
hf[, A1 := toupper(A1)]
hf[, A2 := toupper(A2)]

message("Lifting selected variants and extracting FinnGen R12 sepsis associations...")
gr37 <- GRanges(
  seqnames = paste0("chr", sst2$chr),
  ranges = IRanges(start = sst2$pos, width = 1)
)
lifted <- liftOver(gr37, import.chain(chain_file))
keep_lift <- lengths(lifted) == 1L
sst2 <- sst2[keep_lift]
gr38 <- unlist(lifted[keep_lift], use.names = FALSE)
sst2[, chr38 := as.integer(sub("^chr", "", as.character(seqnames(gr38))))]
sst2[, pos38 := start(gr38)]

fg_keys <- paste(sst2$chr38, sst2$pos38, sep = ":")
awk_fg <- paste0(
  'BEGIN{FS=OFS="\\t"; n=split(keys,k,","); for(i=1;i<=n;i++) wanted[k[i]]=1} ',
  'NR==1 || (($1 ":" $2) in wanted)'
)
fg_cmd <- paste(
  "gzip -dc", shQuote(sepsis_file), "| awk -v",
  shQuote(paste0("keys=", paste(fg_keys, collapse = ","))), shQuote(awk_fg)
)
fg <- fread(cmd = fg_cmd)
setnames(fg, "#chrom", "chr38")
fg[, chr38 := as.integer(chr38)]
fg[, pos38 := as.integer(pos)]
fg[, ref := toupper(ref)]
fg[, alt := toupper(alt)]

complement <- function(a) chartr("ACGT", "TGCA", a)
align_effect <- function(ea_ref, oa_ref, eaf_ref, ea_other, oa_other, eaf_other) {
  if (any(is.na(c(ea_ref, oa_ref, eaf_ref, ea_other, oa_other, eaf_other)))) {
    return(list(ok = FALSE, sign = NA_real_, rule = "missing_association", diff = NA_real_))
  }
  candidates <- list()
  add <- function(label, ea, oa) {
    if (ea == ea_other && oa == oa_other) {
      candidates[[length(candidates) + 1L]] <<- list(label = label, sign = 1, eaf = eaf_other)
    }
    if (ea == oa_other && oa == ea_other) {
      candidates[[length(candidates) + 1L]] <<- list(label = label, sign = -1, eaf = 1 - eaf_other)
    }
  }
  add("direct", ea_ref, oa_ref)
  add("strand", complement(ea_ref), complement(oa_ref))
  if (length(candidates) == 0L) return(list(ok = FALSE, sign = NA_real_, rule = "mismatch", diff = NA_real_))
  diffs <- vapply(candidates, function(z) abs(eaf_ref - z$eaf), numeric(1))
  best <- candidates[[which.min(diffs)]]
  if (min(diffs) > 0.20) return(list(ok = FALSE, sign = NA_real_, rule = "frequency_mismatch", diff = min(diffs)))
  list(ok = TRUE, sign = best$sign, rule = best$label, diff = min(diffs))
}

d <- merge(
  sst2[, .(rsid, chr, pos, chr38, pos38, ea, oa, eaf_sst2,
           beta_sst2, se_sst2, p_sst2_assoc, n_sst2)],
  hf[, .(rsid = SNP, hf_a1 = A1, hf_a2 = A2, eaf_hf = as.numeric(freq),
         beta_hf_raw = as.numeric(b), se_hf = as.numeric(se), p_hf_assoc = as.numeric(p), n_hf = as.numeric(N))],
  by = "rsid", all.x = TRUE
)
d <- merge(
  d,
  fg[, .(chr38, pos38, sepsis_ref = ref, sepsis_alt = alt,
         eaf_sepsis_alt = as.numeric(af_alt), beta_sepsis_raw = as.numeric(beta),
         se_sepsis = as.numeric(sebeta), p_sepsis = as.numeric(pval))],
  by = c("chr38", "pos38"), all.x = TRUE
)
d[, c("hf_ok", "hf_sign", "hf_rule", "hf_freq_diff") := {
  z <- align_effect(ea, oa, eaf_sst2, hf_a1, hf_a2, eaf_hf)
  list(z$ok, z$sign, z$rule, z$diff)
}, by = seq_len(nrow(d))]
d[, c("sepsis_ok", "sepsis_sign", "sepsis_rule", "sepsis_freq_diff") := {
  z <- align_effect(ea, oa, eaf_sst2, sepsis_alt, sepsis_ref, eaf_sepsis_alt)
  list(z$ok, z$sign, z$rule, z$diff)
}, by = seq_len(nrow(d))]
d[, beta_hf := beta_hf_raw * hf_sign]
d[, beta_sepsis := beta_sepsis_raw * sepsis_sign]
d <- d[
  hf_ok == TRUE & sepsis_ok == TRUE &
    is.finite(beta_sst2) & is.finite(se_sst2) & se_sst2 > 0 &
    is.finite(beta_hf) & is.finite(se_hf) & se_hf > 0 &
    is.finite(beta_sepsis) & is.finite(se_sepsis) & se_sepsis > 0
]
setorder(d, rsid)
fwrite(d, file.path(out_dir, "harmonized_mvmr_instruments.tsv"), sep = "\t", na = "NA")

if (nrow(d) < 5L) stop("Too few fully harmonized MVMR instruments: ", nrow(d))

message("Estimating direct effects of sST2 and HF liability on sepsis...")
bx <- as.matrix(d[, .(beta_sst2, beta_hf)])
bxse <- as.matrix(d[, .(se_sst2, se_hf)])
mvmr_input <- format_mvmr(
  BXGs = bx, BYG = d$beta_sepsis,
  seBXGs = bxse, seBYG = d$se_sepsis, RSID = d$rsid
)
conditional_f <- strength_mvmr(mvmr_input, gencov = 0)
validity <- pleiotropy_mvmr(mvmr_input, gencov = 0)

# The exact sampling covariance between the two exposure GWAS estimates is not
# available. Test a deliberately broad range of assumed correlations.
rho_grid <- seq(-0.30, 0.30, by = 0.10)
covariance_sensitivity <- rbindlist(lapply(rho_grid, function(rho) {
  pcor <- matrix(c(1, rho, rho, 1), nrow = 2)
  gencov <- phenocov_mvmr(pcor, bxse)
  invisible(capture.output(f_rho <- strength_mvmr(mvmr_input, gencov = gencov)))
  invisible(capture.output(q_rho <- pleiotropy_mvmr(mvmr_input, gencov = gencov)))
  data.table(
    assumed_exposure_error_correlation = rho,
    conditional_F_sST2 = as.numeric(f_rho[1, 1]),
    conditional_F_HF = as.numeric(f_rho[1, 2]),
    validity_Q = q_rho$Qstat, validity_Q_p = q_rho$Qpval
  )
}))
fwrite(
  covariance_sensitivity,
  file.path(out_dir, "mvmr_exposure_covariance_sensitivity.tsv"),
  sep = "\t", na = "NA"
)

w <- 1 / d$se_sepsis^2
X <- bx
y <- d$beta_sepsis
XtWX_inv <- solve(crossprod(X, w * X))
beta_mv <- as.numeric(XtWX_inv %*% crossprod(X, w * y))
resid <- y - as.numeric(X %*% beta_mv)
Q_outcome <- sum(w * resid^2)
Q_outcome_df <- nrow(d) - ncol(X)
phi <- max(1, Q_outcome / Q_outcome_df)
se_mv <- sqrt(diag(XtWX_inv) * phi)
z <- beta_mv / se_mv

results <- data.table(
  exposure = c("Circulating sST2", "Heart-failure liability"),
  outcome = "FinnGen R12 sepsis",
  nsnp = nrow(d), beta = beta_mv, se = se_mv,
  OR = exp(beta_mv), CI_lower = exp(beta_mv - 1.96 * se_mv),
  CI_upper = exp(beta_mv + 1.96 * se_mv), p = 2 * pnorm(-abs(z)),
  conditional_F = as.numeric(conditional_f[1, ]),
  Q = validity$Qstat, Q_df = nrow(d) - ncol(X) - 1L, Q_p = validity$Qpval,
  outcome_residual_Q = Q_outcome, outcome_residual_Q_df = Q_outcome_df,
  mde_or_80pct_alpha_0.05 = exp((qnorm(0.975) + qnorm(0.80)) * se_mv)
)
fwrite(results, file.path(out_dir, "mvmr_direct_effects.tsv"), sep = "\t", na = "NA")

mv_object <- mr_mvinput(
  bx = bx, bxse = bxse, by = d$beta_sepsis, byse = d$se_sepsis,
  exposure = c("Circulating sST2", "Heart-failure liability"),
  outcome = "FinnGen R12 sepsis", snps = d$rsid
)
robust_fit <- mr_mvivw(mv_object, model = "random", robust = TRUE)
robust_results <- data.table(
  exposure = robust_fit@Exposure, outcome = robust_fit@Outcome,
  nsnp = robust_fit@SNPs, method = "Robust random-effects MVMR-IVW",
  beta = robust_fit@Estimate, se = robust_fit@StdError,
  OR = exp(robust_fit@Estimate), CI_lower = exp(robust_fit@CILower),
  CI_upper = exp(robust_fit@CIUpper), p = robust_fit@Pvalue,
  heterogeneity_Q = robust_fit@Heter.Stat[1],
  heterogeneity_p = robust_fit@Heter.Stat[2]
)
fwrite(robust_results, file.path(out_dir, "mvmr_robust_ivw_sensitivity.tsv"), sep = "\t", na = "NA")

# Weighted MVMR-Egger sensitivity analysis; orient all variants to increasing sST2.
eg <- copy(d)
flip <- eg$beta_sst2 < 0
cols_to_flip <- c("beta_sst2", "beta_hf", "beta_sepsis")
eg[flip, (cols_to_flip) := lapply(.SD, function(v) -v), .SDcols = cols_to_flip]
fit_egger <- lm(beta_sepsis ~ beta_sst2 + beta_hf, data = eg, weights = 1 / se_sepsis^2)
egger_coef <- as.data.table(summary(fit_egger)$coefficients, keep.rownames = "term")
setnames(egger_coef, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
         c("beta", "se", "statistic", "p"))
egger_coef[, `:=`(
  OR = fifelse(term == "(Intercept)", NA_real_, exp(beta)),
  CI_lower = fifelse(term == "(Intercept)", NA_real_, exp(beta - 1.96 * se)),
  CI_upper = fifelse(term == "(Intercept)", NA_real_, exp(beta + 1.96 * se))
)]
fwrite(egger_coef, file.path(out_dir, "mvmr_egger_sensitivity.tsv"), sep = "\t", na = "NA")

audit <- data.table(
  metric = c(
    "sST2 GWS candidates represented in 1KG EUR", "HF GWS candidates represented in 1KG EUR",
    "Jointly clumped instruments", "Fully harmonized MVMR instruments",
    "Instrument validity Q", "Instrument validity Q P"
  ),
  value = c(
    nrow(sst2_candidates), nrow(hf_candidates), length(selected_rsids), nrow(d),
    validity$Qstat, validity$Qpval
  )
)
fwrite(audit, file.path(out_dir, "mvmr_audit.tsv"), sep = "\t", na = "NA")

print(results)
print(robust_results)
print(egger_coef)
print(audit)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info_mvmr.txt"))
