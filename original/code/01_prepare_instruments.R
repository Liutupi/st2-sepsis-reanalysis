suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(Sys.getenv("SST2_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
sst2_file <- Sys.getenv("SST2_SCALLOP_FILE", unset = file.path(project_dir, "data", "GCST90012040.tsv.gz"))
eur_bfile <- Sys.getenv("SST2_EUR_BFILE", unset = file.path(project_dir, "data", "1kg_eur", "EUR"))
plink_bin <- Sys.getenv("PLINK_BIN", unset = Sys.which("plink"))
if (!nzchar(plink_bin)) plink_bin <- file.path(project_dir, "tools", "plink")
out_dir <- file.path(project_dir, "results/instruments")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sst2_file)) stop("Missing SCALLOP file. Set SST2_SCALLOP_FILE.")
if (!file.exists(paste0(eur_bfile, ".bed"))) stop("Missing 1000 Genomes EUR prefix. Set SST2_EUR_BFILE.")
if (!file.exists(plink_bin)) stop("Missing PLINK 1.9 executable. Set PLINK_BIN.")

message("Reading SCALLOP Olink sST2 summary statistics (GCST90012040)...")
st2 <- fread(
  sst2_file,
  select = c(
    "MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr",
    "P-value", "TotalSampleSize"
  )
)

st2[, c("chr_text", "pos_text", "marker_ref_alt") := tstrsplit(MarkerName, ":", fixed = TRUE)]
st2[, chr := suppressWarnings(as.integer(chr_text))]
st2[, pos := suppressWarnings(as.integer(pos_text))]
st2[, pval := suppressWarnings(as.numeric(`P-value`))]
st2[, effect_allele := toupper(Allele1)]
st2[, other_allele := toupper(Allele2)]
st2[, eaf := as.numeric(Freq1)]
st2[, beta := as.numeric(Effect)]
st2[, se := as.numeric(StdErr)]
st2[, n := as.numeric(TotalSampleSize)]

valid_base <- c("A", "C", "G", "T")
candidates <- st2[
  !is.na(chr) & chr >= 1L & chr <= 22L &
    !is.na(pos) & !is.na(pval) & pval < 5e-8 &
    !is.na(beta) & !is.na(se) & se > 0 &
    !is.na(eaf) & eaf >= 0.01 & eaf <= 0.99 &
    effect_allele %chin% valid_base & other_allele %chin% valid_base
]
candidates[, allele_key := paste(pmin(effect_allele, other_allele),
                                 pmax(effect_allele, other_allele), sep = "/")]

message("Reading the 1000 Genomes EUR reference variant map...")
bim <- fread(
  paste0(eur_bfile, ".bim"),
  col.names = c("chr", "rsid", "cm", "pos", "bim_a1", "bim_a2")
)
bim[, bim_a1 := toupper(bim_a1)]
bim[, bim_a2 := toupper(bim_a2)]
bim[, allele_key := paste(pmin(bim_a1, bim_a2), pmax(bim_a1, bim_a2), sep = "/")]
bim <- bim[rsid != "." & grepl("^rs[0-9]+$", rsid)]

candidates <- merge(
  candidates,
  bim[, .(chr, pos, allele_key, rsid, bim_a1, bim_a2)],
  by = c("chr", "pos", "allele_key"),
  all = FALSE,
  allow.cartesian = TRUE
)
setorder(candidates, pval, -n)
candidates <- unique(candidates, by = "rsid")
candidates[, p_for_clump := fifelse(pval == 0, 1e-300, pval)]

keep_cols <- c(
  "rsid", "MarkerName", "chr", "pos", "effect_allele", "other_allele",
  "eaf", "beta", "se", "pval", "n", "bim_a1", "bim_a2"
)
fwrite(candidates[, ..keep_cols], file.path(out_dir, "sst2_gws_candidates.tsv"), sep = "\t", na = "NA")

run_clump <- function(dat, label, r2, kb = 10000) {
  if (nrow(dat) == 0L) stop("No variants available for ", label)
  clump_input <- file.path(out_dir, paste0(label, "_clump_input.tsv"))
  out_prefix <- file.path(out_dir, label)
  fwrite(dat[, .(SNP = rsid, P = p_for_clump)], clump_input, sep = "\t")

  args <- c(
    "--bfile", eur_bfile,
    "--clump", clump_input,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", "5e-8",
    "--clump-p2", "5e-8",
    "--clump-r2", format(r2, scientific = FALSE),
    "--clump-kb", as.character(kb),
    "--out", out_prefix
  )
  status <- system2(plink_bin, args = vapply(args, shQuote, character(1)))
  if (!identical(status, 0L)) stop("PLINK clumping failed for ", label)

  clumped_file <- paste0(out_prefix, ".clumped")
  if (!file.exists(clumped_file)) stop("Missing PLINK clump output for ", label)
  selected <- fread(clumped_file, fill = TRUE)
  if (!"SNP" %in% names(selected)) stop("No index variants found for ", label)
  selected <- selected[!is.na(SNP) & SNP != ""]
  result <- dat[match(selected$SNP, rsid)]
  result <- result[!is.na(rsid)]
  result[, instrument_set := label]
  result[, clump_r2 := r2]
  result[, clump_kb := kb]
  fwrite(result[, c(keep_cols, "instrument_set", "clump_r2", "clump_kb"), with = FALSE],
         file.path(out_dir, paste0(label, "_instruments.tsv")), sep = "\t", na = "NA")
  result
}

# IL1RL1 GRCh37 coordinates from Ensembl GRCh37 REST lookup on 2026-07-17:
# chr2:102,927,962-102,968,497. The cis window is gene +/- 500 kb.
cis_start <- 102927962L - 500000L
cis_end <- 102968497L + 500000L
cis_candidates <- candidates[chr == 2L & pos >= cis_start & pos <= cis_end]
fwrite(cis_candidates[, ..keep_cols], file.path(out_dir, "sst2_cis_candidates.tsv"), sep = "\t", na = "NA")

gw <- run_clump(candidates, "sst2_gw_r2_0.001", r2 = 0.001)
cis_strict <- run_clump(cis_candidates, "sst2_cis_r2_0.001", r2 = 0.001)
cis_relaxed <- run_clump(cis_candidates, "sst2_cis_r2_0.1", r2 = 0.1)

summary <- rbindlist(list(
  data.table(set = "GWS candidates", n = nrow(candidates)),
  data.table(set = "GW instruments r2<0.001", n = nrow(gw)),
  data.table(set = "Cis candidates", n = nrow(cis_candidates)),
  data.table(set = "Cis instruments r2<0.001", n = nrow(cis_strict)),
  data.table(set = "Cis instruments r2<0.1", n = nrow(cis_relaxed))
))
fwrite(summary, file.path(out_dir, "instrument_selection_summary.tsv"), sep = "\t")
print(summary)

session <- capture.output(sessionInfo())
writeLines(session, file.path(out_dir, "session_info_prepare_instruments.txt"))
