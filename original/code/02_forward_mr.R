suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(MendelianRandomization)
})

options(stringsAsFactors = FALSE)
set.seed(20260717)

project_dir <- normalizePath(Sys.getenv("SST2_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
instrument_dir <- file.path(project_dir, "results/instruments")
out_dir <- file.path(project_dir, "results/forward_mr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chain_file <- Sys.getenv("SST2_CHAIN_FILE", unset = file.path(project_dir, "data_metadata/hg19ToHg38.over.chain"))
hermes_file <- Sys.getenv("SST2_HERMES_FILE", unset = file.path(project_dir, "data", "GCST009541.tsv.gz"))
finngen_files <- c(
  Sepsis = Sys.getenv("SST2_FINNGEN_SEPSIS_FILE", unset = file.path(project_dir, "data", "finngen_R12_AB1_OTHER_SEPSIS.gz")),
  Heart_failure = Sys.getenv("SST2_FINNGEN_HF_FILE", unset = file.path(project_dir, "data", "finngen_R12_I9_HEARTFAIL.gz")),
  IBD_positive_control = Sys.getenv("SST2_FINNGEN_IBD_FILE", unset = file.path(project_dir, "data", "finngen_R12_K11_IBD_STRICT.gz"))
)
finngen_metadata <- data.table(
  outcome = names(finngen_files),
  endpoint = c("AB1_OTHER_SEPSIS", "I9_HEARTFAIL", "K11_IBD_STRICT"),
  cases = c(17133L, 37653L, 10960L),
  controls = c(439048L, 462695L, 489388L)
)

if (!file.exists(chain_file)) stop("Missing hg19-to-hg38 chain file.")
if (!file.exists(hermes_file)) stop("Missing HERMES file. Set SST2_HERMES_FILE.")
if (!all(file.exists(finngen_files))) stop("Missing FinnGen file(s). Set SST2_FINNGEN_SEPSIS_FILE, SST2_FINNGEN_HF_FILE, and SST2_FINNGEN_IBD_FILE.")

read_instrument_set <- function(name) {
  path <- file.path(instrument_dir, paste0(name, "_instruments.tsv"))
  x <- fread(path)
  x <- x[!is.na(rsid) & rsid != ""]
  x[, instrument_set := name]
  x
}

gw <- read_instrument_set("sst2_gw_r2_0.001")
cis <- read_instrument_set("sst2_cis_r2_0.001")
instruments <- unique(rbindlist(list(gw, cis), fill = TRUE), by = "rsid")

message("Lifting instrument positions from GRCh37 to GRCh38...")
chain <- import.chain(chain_file)
gr37 <- GRanges(
  seqnames = paste0("chr", instruments$chr),
  ranges = IRanges(start = instruments$pos, width = 1),
  strand = "*"
)
lifted <- liftOver(gr37, chain)
lift_count <- lengths(lifted)
if (any(lift_count != 1L)) {
  stop("Every retained instrument must have exactly one GRCh38 mapping. Problematic rsIDs: ",
       paste(instruments$rsid[lift_count != 1L], collapse = ", "))
}
gr38 <- unlist(lifted, use.names = FALSE)
instruments[, chr38 := as.integer(sub("^chr", "", as.character(seqnames(gr38))))]
instruments[, pos38 := start(gr38)]
fwrite(
  instruments[, .(rsid, chr37 = chr, pos37 = pos, chr38, pos38,
                  effect_allele, other_allele, eaf, beta, se, pval, n)],
  file.path(out_dir, "instrument_liftover_grch37_to_grch38.tsv"),
  sep = "\t", na = "NA"
)

extract_finngen_positions <- function(file, positions, outcome_name) {
  keys <- unique(paste(positions$chr38, positions$pos38, sep = ":"))
  key_text <- paste(keys, collapse = ",")
  awk_program <- paste0(
    "BEGIN{FS=OFS=\"\\t\"; n=split(keys,a,\",\"); ",
    "for(i=1;i<=n;i++) wanted[a[i]]=1} ",
    "NR==1 || (($1 \":\" $2) in wanted)"
  )
  cmd <- paste(
    "gzip -dc", shQuote(file), "|",
    "awk", "-v", shQuote(paste0("keys=", key_text)), shQuote(awk_program)
  )
  x <- fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
  setnames(x, "#chrom", "chr38")
  x[, chr38 := as.integer(chr38)]
  x[, pos := as.integer(pos)]
  x[, outcome := outcome_name]
  x
}

complement <- function(a) {
  chartr("ACGT", "TGCA", a)
}

align_effect <- function(ea_exp, oa_exp, eaf_exp, ea_out, oa_out, eaf_out) {
  ea_exp <- toupper(ea_exp)
  oa_exp <- toupper(oa_exp)
  ea_out <- toupper(ea_out)
  oa_out <- toupper(oa_out)
  candidates <- list()
  add_candidate <- function(label, ea_map, oa_map, strand) {
    if (ea_map == ea_out && oa_map == oa_out) {
      candidates[[length(candidates) + 1L]] <<- list(
        label = label, sign = 1, strand = strand, comparable_eaf = eaf_out
      )
    }
    if (ea_map == oa_out && oa_map == ea_out) {
      candidates[[length(candidates) + 1L]] <<- list(
        label = label, sign = -1, strand = strand, comparable_eaf = 1 - eaf_out
      )
    }
  }
  add_candidate("direct", ea_exp, oa_exp, "+")
  add_candidate("complement", complement(ea_exp), complement(oa_exp), "-")
  if (length(candidates) == 0L) {
    return(list(ok = FALSE, sign = NA_real_, rule = "allele_mismatch", eaf_diff = NA_real_))
  }
  diffs <- vapply(candidates, function(z) abs(eaf_exp - z$comparable_eaf), numeric(1))
  best <- candidates[[which.min(diffs)]]
  if (!is.na(eaf_exp) && !is.na(eaf_out) && min(diffs) > 0.20) {
    return(list(ok = FALSE, sign = NA_real_, rule = "frequency_mismatch", eaf_diff = min(diffs)))
  }
  list(ok = TRUE, sign = best$sign,
       rule = paste(best$label, ifelse(best$sign == 1, "same", "swapped"), sep = "_"),
       eaf_diff = min(diffs))
}

harmonize_finngen <- function(inst, fg, outcome_name) {
  d <- merge(inst, fg, by.x = c("chr38", "pos38"), by.y = c("chr38", "pos"), all.x = TRUE)
  d[, matched_position := !is.na(beta.y)]
  d[, c("align_ok", "align_sign", "align_rule", "eaf_diff") := {
    if (is.na(beta.y)) {
      list(FALSE, NA_real_, "position_not_found", NA_real_)
    } else {
      z <- align_effect(effect_allele, other_allele, eaf, alt, ref, af_alt)
      list(z$ok, z$sign, z$rule, z$eaf_diff)
    }
  }, by = seq_len(nrow(d))]
  d[, beta_outcome := beta.y * align_sign]
  d[, se_outcome := sebeta]
  d[, pval_outcome := pval.y]
  d[, beta_exposure := beta.x]
  d[, se_exposure := se]
  d[, pval_exposure := pval.x]
  d[, n_exposure := n]
  d[, eaf_outcome_effect_allele := fifelse(align_sign == 1, af_alt, 1 - af_alt)]
  d[, outcome := outcome_name]
  d[, data_source := "FinnGen R12"]
  d
}

message("Extracting exact GRCh38 variants from FinnGen endpoints...")
finngen_harmonized <- list()
for (nm in names(finngen_files)) {
  fg <- extract_finngen_positions(finngen_files[[nm]], instruments, nm)
  finngen_harmonized[[nm]] <- harmonize_finngen(instruments, fg, nm)
}

message("Reading HERMES heart-failure summary statistics and matching rsIDs...")
hermes <- fread(
  hermes_file,
  select = c("SNP", "CHR", "BP", "A1", "A2", "freq", "b", "se", "p", "N")
)
hermes <- hermes[SNP %chin% instruments$rsid]
hermes_h <- merge(instruments, hermes, by.x = "rsid", by.y = "SNP", all.x = TRUE)
hermes_h[, c("align_ok", "align_sign", "align_rule", "eaf_diff") := {
  if (is.na(b)) {
    list(FALSE, NA_real_, "rsid_not_found", NA_real_)
  } else {
    z <- align_effect(effect_allele, other_allele, eaf, A1, A2, freq)
    list(z$ok, z$sign, z$rule, z$eaf_diff)
  }
}, by = seq_len(nrow(hermes_h))]
hermes_h[, beta_outcome := b * align_sign]
hermes_h[, se_outcome := se.y]
hermes_h[, pval_outcome := p]
hermes_h[, beta_exposure := beta]
hermes_h[, se_exposure := se.x]
hermes_h[, pval_exposure := pval]
hermes_h[, n_exposure := n]
hermes_h[, eaf_outcome_effect_allele := fifelse(align_sign == 1, freq, 1 - freq)]
hermes_h[, outcome := "Heart_failure_HERMES"]
hermes_h[, data_source := "HERMES 2020"]

all_harmonized <- rbindlist(c(finngen_harmonized, list(hermes_h)), fill = TRUE)
all_harmonized[, f_statistic := (beta_exposure / se_exposure)^2]
all_harmonized[, r2_exposure := 2 * eaf * (1 - eaf) * beta_exposure^2]

harmonized_cols <- c(
  "outcome", "data_source", "rsid", "chr", "pos", "chr38", "pos38",
  "effect_allele", "other_allele", "eaf", "beta_exposure", "se_exposure",
  "pval_exposure", "n_exposure", "f_statistic", "r2_exposure",
  "beta_outcome", "se_outcome", "pval_outcome", "eaf_outcome_effect_allele",
  "align_ok", "align_rule", "eaf_diff"
)
fwrite(all_harmonized[, ..harmonized_cols],
       file.path(out_dir, "all_exact_harmonized_variants.tsv"), sep = "\t", na = "NA")

ivw_core <- function(dat) {
  k <- nrow(dat)
  bx <- dat$beta_exposure
  by <- dat$beta_outcome
  sy <- dat$se_outcome
  denom <- sum(bx^2 / sy^2)
  estimate <- sum(bx * by / sy^2) / denom
  se_fixed <- sqrt(1 / denom)
  q <- sum((by - estimate * bx)^2 / sy^2)
  q_df <- k - 1L
  q_p <- if (q_df > 0) pchisq(q, df = q_df, lower.tail = FALSE) else NA_real_
  rse <- if (q_df > 0) sqrt(max(1, q / q_df)) else 1
  se_random <- se_fixed * rse
  list(k = k, estimate = estimate, se_fixed = se_fixed, se_random = se_random,
       q = q, q_df = q_df, q_p = q_p, rse = rse)
}

result_row <- function(set_name, outcome_name, method, estimate, se, k,
                       q = NA_real_, q_df = NA_integer_, q_p = NA_real_) {
  z <- estimate / se
  data.table(
    instrument_set = set_name,
    outcome = outcome_name,
    method = method,
    nsnp = k,
    beta = estimate,
    se = se,
    p = 2 * pnorm(abs(z), lower.tail = FALSE),
    or = exp(estimate),
    ci_lower = exp(estimate - 1.96 * se),
    ci_upper = exp(estimate + 1.96 * se),
    q = q,
    q_df = q_df,
    q_p = q_p
  )
}

run_mr <- function(dat, set_name, outcome_name) {
  d <- dat[align_ok == TRUE & !is.na(beta_outcome) & !is.na(se_outcome) & se_outcome > 0]
  if (nrow(d) == 0L) return(NULL)
  core <- ivw_core(d)
  out <- rbindlist(list(
    result_row(set_name, outcome_name, "IVW fixed effect", core$estimate, core$se_fixed,
               core$k, core$q, core$q_df, core$q_p),
    result_row(set_name, outcome_name, "IVW multiplicative random effects",
               core$estimate, core$se_random, core$k, core$q, core$q_df, core$q_p)
  ))

  if (nrow(d) >= 3L) {
    input <- mr_input(
      bx = d$beta_exposure,
      bxse = d$se_exposure,
      by = d$beta_outcome,
      byse = d$se_outcome,
      exposure = "Circulating sST2",
      outcome = outcome_name,
      snps = d$rsid
    )
    median_fit <- mr_median(input, weighting = "weighted", iterations = 10000, seed = 20260717)
    out <- rbind(out, result_row(set_name, outcome_name, "Weighted median",
                                 median_fit@Estimate, median_fit@StdError, nrow(d)))
  }
  if (nrow(d) >= 4L) {
    egger_fit <- mr_egger(input, distribution = "normal")
    out <- rbind(out, result_row(set_name, outcome_name, "MR-Egger",
                                 egger_fit@Estimate, egger_fit@StdError.Est, nrow(d)))
    out <- rbind(out, data.table(
      instrument_set = set_name,
      outcome = outcome_name,
      method = "MR-Egger intercept",
      nsnp = nrow(d),
      beta = egger_fit@Intercept,
      se = egger_fit@StdError.Int,
      p = egger_fit@Pvalue.Int,
      or = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      q = NA_real_, q_df = NA_integer_, q_p = NA_real_
    ), fill = TRUE)
  }

  d[, ratio_beta := beta_outcome / beta_exposure]
  d[, ratio_se := se_outcome / abs(beta_exposure)]
  d[, ratio_or := exp(ratio_beta)]
  d[, ratio_ci_lower := exp(ratio_beta - 1.96 * ratio_se)]
  d[, ratio_ci_upper := exp(ratio_beta + 1.96 * ratio_se)]
  fwrite(d[, ..harmonized_cols],
         file.path(out_dir, paste0("harmonized_", set_name, "_", outcome_name, ".tsv")),
         sep = "\t", na = "NA")
  fwrite(d[, .(rsid, ratio_beta, ratio_se, ratio_or, ratio_ci_lower, ratio_ci_upper)],
         file.path(out_dir, paste0("single_variant_", set_name, "_", outcome_name, ".tsv")),
         sep = "\t", na = "NA")

  if (nrow(d) >= 3L) {
    loo <- rbindlist(lapply(seq_len(nrow(d)), function(i) {
      fit <- ivw_core(d[-i])
      result_row(set_name, outcome_name, paste0("Leave-one-out: ", d$rsid[i]),
                 fit$estimate, fit$se_random, fit$k, fit$q, fit$q_df, fit$q_p)
    }))
    fwrite(loo, file.path(out_dir, paste0("leave_one_out_", set_name, "_", outcome_name, ".tsv")),
           sep = "\t", na = "NA")
  }
  out
}

results <- list()
for (set_name in c("sst2_gw_r2_0.001", "sst2_cis_r2_0.001")) {
  set_ids <- if (set_name == "sst2_gw_r2_0.001") gw$rsid else cis$rsid
  for (outcome_name in unique(all_harmonized$outcome)) {
    key <- paste(set_name, outcome_name, sep = "__")
    results[[key]] <- run_mr(
      all_harmonized[rsid %chin% set_ids & outcome == outcome_name],
      set_name,
      outcome_name
    )
  }
}
mr_results <- rbindlist(results, fill = TRUE)
fwrite(mr_results, file.path(out_dir, "forward_mr_results.tsv"), sep = "\t", na = "NA")

instrument_strength <- all_harmonized[
  align_ok == TRUE,
  .(
    nsnp = uniqueN(rsid),
    min_f = min(f_statistic, na.rm = TRUE),
    mean_f = mean(f_statistic, na.rm = TRUE),
    total_r2 = sum(r2_exposure[!duplicated(rsid)], na.rm = TRUE)
  ),
  by = outcome
]
fwrite(instrument_strength, file.path(out_dir, "instrument_strength.tsv"), sep = "\t", na = "NA")

primary <- mr_results[method == "IVW multiplicative random effects"]
primary[, mde_or_80pct_alpha_0.05 := exp((qnorm(0.975) + qnorm(0.80)) * se)]
primary[, mde_or_80pct_alpha_0.025 := exp((qnorm(0.9875) + qnorm(0.80)) * se)]
fwrite(primary, file.path(out_dir, "primary_estimates_and_mde.tsv"), sep = "\t", na = "NA")
fwrite(finngen_metadata, file.path(out_dir, "finngen_endpoint_metadata.tsv"), sep = "\t")

print(primary[, .(instrument_set, outcome, nsnp, or, ci_lower, ci_upper, p, q_p,
                  mde_or_80pct_alpha_0.05)])
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info_forward_mr.txt"))
