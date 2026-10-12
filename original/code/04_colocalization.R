suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(coloc)
})

options(stringsAsFactors = FALSE)
set.seed(20260717)

project_dir <- normalizePath(Sys.getenv("SST2_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
sst2_file <- Sys.getenv("SST2_SCALLOP_FILE", unset = file.path(project_dir, "data", "GCST90012040.tsv.gz"))
eur_bfile <- Sys.getenv("SST2_EUR_BFILE", unset = file.path(project_dir, "data", "1kg_eur", "EUR"))
chain_file <- Sys.getenv("SST2_CHAIN_FILE", unset = file.path(project_dir, "data_metadata/hg19ToHg38.over.chain"))
out_dir <- file.path(project_dir, "results/colocalization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

finngen_files <- c(
  Sepsis = Sys.getenv("SST2_FINNGEN_SEPSIS_FILE", unset = file.path(project_dir, "data", "finngen_R12_AB1_OTHER_SEPSIS.gz")),
  Heart_failure = Sys.getenv("SST2_FINNGEN_HF_FILE", unset = file.path(project_dir, "data", "finngen_R12_I9_HEARTFAIL.gz")),
  IBD_positive_control = Sys.getenv("SST2_FINNGEN_IBD_FILE", unset = file.path(project_dir, "data", "finngen_R12_K11_IBD_STRICT.gz"))
)
meta <- data.table(
  outcome = names(finngen_files),
  cases = c(17133, 37653, 10960),
  controls = c(439048, 462695, 489388)
)
if (!file.exists(sst2_file)) stop("Missing SCALLOP file. Set SST2_SCALLOP_FILE.")
if (!file.exists(paste0(eur_bfile, ".bed"))) stop("Missing 1000 Genomes EUR prefix. Set SST2_EUR_BFILE.")
if (!file.exists(chain_file)) stop("Missing hg19-to-hg38 chain file.")
if (!all(file.exists(finngen_files))) stop("Missing FinnGen file(s). Set SST2_FINNGEN_SEPSIS_FILE, SST2_FINNGEN_HF_FILE, and SST2_FINNGEN_IBD_FILE.")

# IL1RL1 GRCh37 coordinates: chr2:102,927,962-102,968,497; use gene +/-500 kb.
region_start <- 102927962L - 500000L
region_end <- 102968497L + 500000L

message("Extracting the SCALLOP Olink sST2 IL1RL1 region...")
awk_program <- sprintf(
  'BEGIN{FS=OFS="\\t"} NR==1{print; next} {split($1,a,":"); if(a[1]==2 && a[2]>=%d && a[2]<=%d) print}',
  region_start, region_end
)
cmd <- sprintf("gzip -dc %s | awk %s", shQuote(sst2_file), shQuote(awk_program))
st2 <- fread(cmd = cmd)
st2[, c("chr_text", "pos_text", "marker_ref_alt") := tstrsplit(MarkerName, ":", fixed = TRUE)]
st2[, chr := as.integer(chr_text)]
st2[, pos := as.integer(pos_text)]
st2[, effect_allele := toupper(Allele1)]
st2[, other_allele := toupper(Allele2)]
st2[, allele_key := paste(pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = "/")]
st2[, beta_exposure := as.numeric(Effect)]
st2[, se_exposure := as.numeric(StdErr)]
st2[, eaf_exposure := as.numeric(Freq1)]
st2[, n_exposure := as.numeric(TotalSampleSize)]
st2 <- st2[
  effect_allele %chin% c("A", "C", "G", "T") &
    other_allele %chin% c("A", "C", "G", "T") &
    is.finite(beta_exposure) & is.finite(se_exposure) & se_exposure > 0 &
    is.finite(eaf_exposure) & eaf_exposure > 0 & eaf_exposure < 1
]

message("Annotating exact rsIDs from the 1000 Genomes EUR reference...")
bim <- fread(
  paste0(eur_bfile, ".bim"),
  col.names = c("chr", "rsid", "cm", "pos", "bim_a1", "bim_a2")
)
bim <- bim[chr == 2L & pos >= region_start & pos <= region_end & grepl("^rs[0-9]+$", rsid)]
bim[, bim_a1 := toupper(bim_a1)]
bim[, bim_a2 := toupper(bim_a2)]
bim[, allele_key := paste(pmin(bim_a1, bim_a2), pmax(bim_a1, bim_a2), sep = "/")]
st2 <- merge(
  st2,
  bim[, .(chr, pos, allele_key, rsid)],
  by = c("chr", "pos", "allele_key"), all = FALSE, allow.cartesian = TRUE
)
setorder(st2, rsid, se_exposure)
st2 <- unique(st2, by = "rsid")

message("Lifting the regional coordinates to GRCh38...")
st2[, idx := .I]
gr <- GRanges(seqnames = "chr2", ranges = IRanges(start = st2$pos, width = 1))
mcols(gr)$idx <- seq_len(nrow(st2))
lifted <- liftOver(gr, import.chain(chain_file))
map <- rbindlist(lapply(seq_along(lifted), function(i) {
  if (length(lifted[[i]]) != 1L) return(NULL)
  data.table(idx = i, pos38 = start(lifted[[i]]), chr38 = as.character(seqnames(lifted[[i]])))
}))
st2 <- merge(st2, map, by = "idx", all = FALSE)
st2 <- st2[chr38 == "chr2"]
fwrite(
  st2[, .(rsid, chr37 = chr, pos37 = pos, chr38, pos38, effect_allele,
          other_allele, beta_exposure, se_exposure, eaf_exposure, n_exposure)],
  file.path(out_dir, "sst2_regional_variants.tsv"), sep = "\t", na = "NA"
)

comp <- function(a) chartr("ACGT", "TGCA", a)

harmonize_finngen <- function(path, outcome) {
  lo <- min(st2$pos38) - 1000L
  hi <- max(st2$pos38) + 1000L
  awk <- sprintf(
    'BEGIN{FS=OFS="\\t"} /^#/{print; next} $1==2 && $2>=%d && $2<=%d {print}', lo, hi
  )
  cmd_fg <- sprintf("gzip -dc %s | awk %s", shQuote(path), shQuote(awk))
  x <- fread(cmd = cmd_fg)
  if ("#chrom" %in% names(x)) setnames(x, "#chrom", "chrom")
  x[, rsid := sub("[,;].*$", "", rsids)]
  x <- x[rsid %chin% st2$rsid]
  x[, ref := toupper(ref)]
  x[, alt := toupper(alt)]
  x[, beta_outcome_raw := as.numeric(beta)]
  x[, se_outcome := as.numeric(sebeta)]
  x[, eaf_outcome_raw := as.numeric(af_alt)]
  x <- merge(
    st2[, .(rsid, pos37 = pos, pos38, effect_allele, other_allele,
            beta_exposure, se_exposure, eaf_exposure, n_exposure)],
    x[, .(rsid, pos_outcome = pos, ref, alt, beta_outcome_raw,
          se_outcome, eaf_outcome_raw, p_outcome = as.numeric(pval))],
    by = "rsid", all = FALSE
  )
  x[, alignment := fifelse(
    alt == effect_allele & ref == other_allele, "direct",
    fifelse(ref == effect_allele & alt == other_allele, "swapped",
      fifelse(comp(alt) == effect_allele & comp(ref) == other_allele, "strand_direct",
        fifelse(comp(ref) == effect_allele & comp(alt) == other_allele, "strand_swapped", "mismatch")
      )
    )
  )]
  x[, beta_outcome := fifelse(
    alignment %chin% c("direct", "strand_direct"), beta_outcome_raw,
    fifelse(alignment %chin% c("swapped", "strand_swapped"), -beta_outcome_raw, NA_real_)
  )]
  x[, eaf_outcome := fifelse(
    alignment %chin% c("direct", "strand_direct"), eaf_outcome_raw,
    fifelse(alignment %chin% c("swapped", "strand_swapped"), 1 - eaf_outcome_raw, NA_real_)
  )]
  x[, freq_diff := abs(eaf_exposure - eaf_outcome)]
  x <- x[
    alignment != "mismatch" & is.finite(beta_outcome) & is.finite(se_outcome) &
      se_outcome > 0 & is.finite(freq_diff) & freq_diff <= 0.15
  ]
  setorder(x, rsid, freq_diff)
  x <- unique(x, by = "rsid")
  x[, outcome := outcome]
  x
}

run_coloc <- function(x, outcome, cases, controls) {
  if (nrow(x) < 100L) stop("Insufficient regional overlap for ", outcome, ": ", nrow(x))
  d1 <- list(
    beta = x$beta_exposure,
    varbeta = x$se_exposure^2,
    N = round(median(x$n_exposure, na.rm = TRUE)),
    type = "quant",
    sdY = 1,
    snp = x$rsid,
    position = x$pos37
  )
  d2 <- list(
    beta = x$beta_outcome,
    varbeta = x$se_outcome^2,
    N = cases + controls,
    s = cases / (cases + controls),
    type = "cc",
    snp = x$rsid,
    position = x$pos37
  )
  p12_values <- c(1e-6, 1e-5, 1e-4)
  fits <- lapply(p12_values, function(p12) coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12))
  summary <- rbindlist(lapply(seq_along(fits), function(i) {
    s <- fits[[i]]$summary
    data.table(
      outcome = outcome, nsnp = as.integer(s[["nsnps"]]), p1 = 1e-4, p2 = 1e-4,
      p12 = p12_values[i], PP.H0 = as.numeric(s[["PP.H0.abf"]]),
      PP.H1 = as.numeric(s[["PP.H1.abf"]]), PP.H2 = as.numeric(s[["PP.H2.abf"]]),
      PP.H3 = as.numeric(s[["PP.H3.abf"]]), PP.H4 = as.numeric(s[["PP.H4.abf"]])
    )
  }))
  default_fit <- fits[[which(p12_values == 1e-5)]]
  top <- as.data.table(default_fit$results)
  setorder(top, -SNP.PP.H4)
  top <- top[1:min(.N, 20)]
  top[, outcome := outcome]
  list(summary = summary, top = top)
}

all_summary <- list()
all_top <- list()
for (i in seq_len(nrow(meta))) {
  outcome <- meta$outcome[i]
  message("Harmonizing and colocating: ", outcome)
  x <- harmonize_finngen(finngen_files[[outcome]], outcome)
  fwrite(x, file.path(out_dir, paste0("harmonized_region_", outcome, ".tsv")), sep = "\t", na = "NA")
  fit <- run_coloc(x, outcome, meta$cases[i], meta$controls[i])
  all_summary[[outcome]] <- fit$summary
  all_top[[outcome]] <- fit$top
}

summary <- rbindlist(all_summary, fill = TRUE)
top <- rbindlist(all_top, fill = TRUE)
fwrite(summary, file.path(out_dir, "coloc_prior_sensitivity.tsv"), sep = "\t", na = "NA")
fwrite(top, file.path(out_dir, "coloc_top_variants_default_prior.tsv"), sep = "\t", na = "NA")
print(summary)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info_coloc.txt"))
