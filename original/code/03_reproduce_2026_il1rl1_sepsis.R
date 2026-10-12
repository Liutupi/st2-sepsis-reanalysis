suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

options(stringsAsFactors = FALSE)
set.seed(20260717)

project_dir <- normalizePath(Sys.getenv("SST2_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
supp_file <- Sys.getenv("SST2_BENCHMARK_SUPPLEMENT", unset = file.path(
  project_dir, "data_metadata/benchmark_2026_bmc_infectious_diseases",
  "12879_2026_12748_MOESM2_ESM.xls"
))
sst2_file <- Sys.getenv("SST2_SCALLOP_FILE", unset = file.path(project_dir, "data", "GCST90012040.tsv.gz"))
eur_bfile <- Sys.getenv("SST2_EUR_BFILE", unset = file.path(project_dir, "data", "1kg_eur", "EUR"))
plink_bin <- Sys.getenv("PLINK_BIN", unset = Sys.which("plink"))
if (!nzchar(plink_bin)) plink_bin <- file.path(project_dir, "tools", "plink")
out_dir <- file.path(project_dir, "results/benchmark_2026_il1rl1_sepsis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(supp_file)) stop("Missing published comparison supplement.")
if (!file.exists(sst2_file)) stop("Missing SCALLOP file. Set SST2_SCALLOP_FILE.")
if (!file.exists(paste0(eur_bfile, ".bed"))) stop("Missing 1000 Genomes EUR prefix. Set SST2_EUR_BFILE.")
if (!file.exists(plink_bin)) stop("Missing PLINK 1.9 executable. Set PLINK_BIN.")

message("Reading the published UKB-PPP IL1RL1 instrument-level data...")
ukb <- as.data.table(read_excel(supp_file, sheet = "Sheet3", skip = 1))
ukb <- ukb[`Gene symbol` == "IL1RL1" & mr_keep == TRUE]
stopifnot(nrow(ukb) == 49L)
setnames(
  ukb,
  c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure",
    "se.exposure", "pval.exposure", "beta.outcome", "se.outcome", "pval.outcome"),
  c("rsid", "effect_allele", "other_allele", "beta_ukb_protein", "se_ukb_protein",
    "p_ukb_protein", "beta_ukb_sepsis", "se_ukb_sepsis", "p_ukb_sepsis")
)
ukb[, effect_allele := toupper(effect_allele)]
ukb[, other_allele := toupper(other_allele)]

bim <- fread(
  paste0(eur_bfile, ".bim"),
  col.names = c("chr", "rsid", "cm", "pos", "bim_a1", "bim_a2")
)
bim[, bim_a1 := toupper(bim_a1)]
bim[, bim_a2 := toupper(bim_a2)]
bim <- bim[rsid %chin% ukb$rsid]
setorder(bim, rsid, chr, pos)
bim <- unique(bim, by = "rsid")
ukb <- merge(ukb, bim[, .(rsid, chr, pos, bim_a1, bim_a2)], by = "rsid", all.x = TRUE)
fwrite(ukb, file.path(out_dir, "published_ukb_ppp_il1rl1_49_instruments.tsv"), sep = "\t", na = "NA")

run_strict_clump <- function(dat, label, p_col) {
  inp <- file.path(out_dir, paste0(label, "_clump_input.tsv"))
  fwrite(dat[, .(SNP = rsid, P = pmax(as.numeric(get(p_col)), 1e-300))], inp, sep = "\t")
  prefix <- file.path(out_dir, label)
  args <- c(
    "--bfile", eur_bfile,
    "--clump", inp,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", "1", "--clump-p2", "1",
    "--clump-r2", "0.001", "--clump-kb", "10000",
    "--out", prefix
  )
  status <- system2(plink_bin, args = vapply(args, shQuote, character(1)))
  if (!identical(status, 0L)) stop("PLINK clumping failed for ", label)
  x <- fread(paste0(prefix, ".clumped"), fill = TRUE)
  x <- x[!is.na(SNP) & SNP != ""]
  dat[rsid %chin% x$SNP]
}

ivw <- function(dat, bx, sx, by, sy, label) {
  x <- dat[
    is.finite(get(bx)) & is.finite(get(sx)) & get(sx) > 0 &
      is.finite(get(by)) & is.finite(get(sy)) & get(sy) > 0
  ]
  w <- 1 / x[[sy]]^2
  b <- sum(w * x[[bx]] * x[[by]]) / sum(w * x[[bx]]^2)
  se_fixed <- sqrt(1 / sum(w * x[[bx]]^2))
  q <- sum(w * (x[[by]] - b * x[[bx]])^2)
  df <- nrow(x) - 1L
  phi <- if (df > 0) max(1, q / df) else 1
  se_random <- se_fixed * sqrt(phi)
  z <- b / se_random
  data.table(
    analysis = label,
    nsnp = nrow(x), beta = b, se = se_random,
    OR = exp(b), CI_lower = exp(b - 1.96 * se_random),
    CI_upper = exp(b + 1.96 * se_random), p = 2 * pnorm(-abs(z)),
    Q = q, Q_df = df, Q_p = if (df > 0) pchisq(q, df, lower.tail = FALSE) else NA_real_
  )
}

strict_ukb <- run_strict_clump(ukb[!is.na(chr)], "ukb_ppp_il1rl1_r2_0.001", "p_ukb_protein")

message("Extracting the full SCALLOP Olink sST2 cis region for cross-study comparison...")
cis_start <- 102927962L - 1000000L
cis_end <- 102968497L + 1000000L
awk_program <- sprintf(
  'BEGIN{FS=OFS="\\t"} NR==1{print; next} {split($1,a,":"); if(a[1]==2 && a[2]>=%d && a[2]<=%d) print}',
  cis_start, cis_end
)
cmd <- sprintf("gzip -dc %s | awk %s", shQuote(sst2_file), shQuote(awk_program))
scallop_region <- fread(cmd = cmd)
scallop_region[, c("chr_text", "pos_text", "ref_alt") := tstrsplit(MarkerName, ":", fixed = TRUE)]
scallop_region[, chr := as.integer(chr_text)]
scallop_region[, pos := as.integer(pos_text)]
scallop_region[, a1 := toupper(Allele1)]
scallop_region[, a2 := toupper(Allele2)]
scallop_region[, allele_key := paste(pmin(a1, a2), pmax(a1, a2), sep = "/")]

bim_map <- copy(bim)
bim_map[, allele_key := paste(pmin(bim_a1, bim_a2), pmax(bim_a1, bim_a2), sep = "/")]
scallop_region <- merge(
  scallop_region,
  bim_map[, .(rsid, chr, pos, allele_key)],
  by = c("chr", "pos", "allele_key"), all = FALSE, allow.cartesian = TRUE
)
scallop_region <- unique(scallop_region, by = "rsid")
scallop <- scallop_region[, .(
  rsid,
  scallop_a1 = a1,
  scallop_a2 = a2,
  eaf_scallop = as.numeric(Freq1),
  beta_scallop_raw = as.numeric(Effect),
  se_scallop = as.numeric(StdErr),
  p_scallop = suppressWarnings(as.numeric(`P-value`)),
  n_scallop = as.numeric(TotalSampleSize)
)]

cross <- merge(ukb, scallop, by = "rsid", all.x = TRUE)
cross[, scallop_alignment := fifelse(
  scallop_a1 == effect_allele & scallop_a2 == other_allele, "direct",
  fifelse(scallop_a1 == other_allele & scallop_a2 == effect_allele, "swapped", "mismatch")
)]
cross[, beta_scallop := fifelse(
  scallop_alignment == "direct", beta_scallop_raw,
  fifelse(scallop_alignment == "swapped", -beta_scallop_raw, NA_real_)
)]
cross[, eaf_scallop_aligned := fifelse(
  scallop_alignment == "direct", eaf_scallop,
  fifelse(scallop_alignment == "swapped", 1 - eaf_scallop, NA_real_)
)]
fwrite(cross, file.path(out_dir, "ukb_ppp_vs_scallop_sst2_cross_study.tsv"), sep = "\t", na = "NA")

cross_gws <- cross[!is.na(beta_scallop) & !is.na(p_scallop) & p_scallop < 5e-8]
strict_scallop <- run_strict_clump(cross_gws[!is.na(chr)], "scallop_sst2_overlap_r2_0.001", "p_scallop")

results <- rbindlist(list(
  ivw(ukb, "beta_ukb_protein", "se_ukb_protein", "beta_ukb_sepsis", "se_ukb_sepsis",
      "Published UKB-PPP instruments, r2<0.1 (standard IVW reproduction)"),
  ivw(strict_ukb, "beta_ukb_protein", "se_ukb_protein", "beta_ukb_sepsis", "se_ukb_sepsis",
      "UKB-PPP instruments, strict r2<0.001"),
  ivw(cross_gws, "beta_scallop", "se_scallop", "beta_ukb_sepsis", "se_ukb_sepsis",
      "SCALLOP sST2 effects on overlapping GWS variants, residual LD retained"),
  ivw(strict_scallop, "beta_scallop", "se_scallop", "beta_ukb_sepsis", "se_ukb_sepsis",
      "SCALLOP sST2 effects on overlapping GWS variants, strict r2<0.001")
), fill = TRUE)
fwrite(results, file.path(out_dir, "cross_platform_ivw_results.tsv"), sep = "\t", na = "NA")

cor_stats <- cross[!is.na(beta_scallop), .(
  n_overlap = .N,
  pearson_r = cor(beta_ukb_protein, beta_scallop, method = "pearson"),
  spearman_rho = cor(beta_ukb_protein, beta_scallop, method = "spearman"),
  sign_concordance = mean(sign(beta_ukb_protein) == sign(beta_scallop))
)]
fwrite(cor_stats, file.path(out_dir, "cross_platform_concordance.tsv"), sep = "\t", na = "NA")
fwrite(strict_ukb, file.path(out_dir, "ukb_ppp_il1rl1_strict_r2_0.001.tsv"), sep = "\t", na = "NA")
fwrite(strict_scallop, file.path(out_dir, "scallop_sst2_overlap_strict_r2_0.001.tsv"), sep = "\t", na = "NA")

benchmark <- data.table(
  item = c(
    "Published deCODE IL1RL1-sepsis estimate", "Published UKB-PPP IL1RL1-sepsis estimate",
    "Published deCODE FDR", "Published UKB-PPP FDR",
    "deCODE instrument-level rows present in Supplementary Table S1"
  ),
  value = c("OR 1.0353; P=0.0102; 65 SNPs", "OR 1.0489; P=0.00586; 49 SNPs",
            "1", "1", "No (0 IL1RL1/ST2 rows found)")
)
fwrite(benchmark, file.path(out_dir, "published_claim_audit.tsv"), sep = "\t", na = "NA")

print(results)
print(cor_stats)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
