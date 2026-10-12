suppressPackageStartupMessages({library(data.table); library(susieR)})
set.seed(20261012)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript so its directory can be resolved.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
p <- normalizePath(Sys.getenv("COLOC_REVISION_OUT", unset = script_dir), mustWork = TRUE)
x <- fread(file.path(p, "matched_regional_summary_statistics.tsv"))
R <- readRDS(file.path(p, "signed_LD_EUR_503.rds"))
stopifnot(identical(rownames(R), x$rsid), identical(colnames(R), x$rsid),
          max(abs(R-t(R))) < 1e-12, max(abs(diag(R)-1)) < 1e-12,
          all(x$pos38 == x$pos_outcome))
diag <- fread(file.path(p, "z_LD_diagnostic_summary.tsv"))
for (i in seq_len(nrow(diag))) {
  z <- fread(file.path(p, paste0("z_LD_diagnostic_", diag$analysis_set[i], "_", diag$trait[i], ".tsv")))
  diag[i, `:=`(max_abs_standardized_conditional_residual = max(abs(z$z_std_diff)),
               max_flip_logLR = max(z$logLR), n_flip_logLR_above_2 = sum(z$logLR > 2))]
}
fwrite(diag, file.path(p, "z_LD_diagnostic_summary.tsv"), sep = "\t")
fit <- fread(file.path(p, "susie_fit_summary.tsv"))
fit[, `:=`(estimate_prior_variance = TRUE, tolerance = 0.001, max_iter = 1000)]
fit[, reaches_L_cap := n_credible_sets >= L]
fit[, inferential_status := fifelse(!converged, "Nonconverged: do not interpret PIPs or credible sets", 
                            fifelse(n_credible_sets == 0, "Converged, no qualifying credible set",
                            fifelse(reaches_L_cap, "Converged but credible sets reach L cap; external-LD diagnostics apply",
                                    "Converged; external-LD diagnostics still apply")))]
fwrite(fit, file.path(p, "susie_fit_summary.tsv"), sep = "\t")
cs <- fread(file.path(p,"susie_credible_sets.tsv"))
cs <- cs[,setdiff(names(cs),c("converged","reaches_L_cap","inferential_status")),with=FALSE]
cs <- merge(cs,fit[,.(analysis_set,trait,converged,reaches_L_cap,inferential_status)],
            by=c("analysis_set","trait"),all.x=TRUE,sort=FALSE)
fwrite(cs,file.path(p,"susie_credible_sets.tsv"),sep="\t")

# Check that the weak binary-outcome signal is not an artifact of choosing
# N_eff rather than total N, or of the cap L=10. These are sensitivity checks,
# not a search for a relaxed credible-set definition.
attr(R, "eigen") <- readRDS(file.path(p, "LD_eigendecomposition_all_matched.rds"))
settings <- data.table(setting = c("total_N_L10", "effective_N_L20"),
                       n = c(17133+439048, 4/(1/17133+1/439048)), L = c(10L,20L))
sens <- rbindlist(lapply(seq_len(nrow(settings)), function(i) {
  warns <- character()
  s <- withCallingHandlers(susie_rss(z = x$z_outcome, R = R, n = settings$n[i], L = settings$L[i],
        coverage = 0.95, min_abs_corr = 0.5, max_iter = 1000, tol = 0.001,
        estimate_residual_variance = FALSE), warning = function(w) {warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
  saveRDS(s, file.path(p, paste0("sepsis_specification_", settings$setting[i], ".rds")), compress = "xz")
  data.table(setting = settings$setting[i], nsnps = nrow(x), N = settings$n[i], L = settings$L[i],
             coverage = 0.95, min_abs_corr = 0.5, converged = s$converged,
             niter = s$niter, credible_sets = length(s$sets$cs), max_pip = max(s$pip),
             warnings = paste(unique(warns), collapse = " | "))
}))
fwrite(sens, file.path(p, "sepsis_rss_specification_sensitivity.tsv"), sep = "\t")

# Explicitly call coloc.susie for the subset where both SuSiE fits converged.
# With no sepsis credible set the package returns nsnps=NA rather than posteriors.
sc <- readRDS(file.path(p,"susie_common_highN_exposure.rds"))
so <- readRDS(file.path(p,"susie_common_highN_sepsis.rds"))
stopifnot(sc$converged,so$converged)
no_pairs <- coloc::coloc.susie(sc,so,p1=1e-4,p2=1e-4,p12=1e-5)
saveRDS(no_pairs,file.path(p,"coloc_susie_common_highN_no_signal_pairs.rds"))
fwrite(as.data.table(no_pairs),file.path(p,"coloc_susie_common_highN_no_signal_pairs.tsv"),sep="\t",na="NA")

# Exposure capacity check: reaching L=10 is not evidence for precisely ten
# biological signals, especially under the observed LD mismatch.
keep <- pmin(x$eaf_exposure,1-x$eaf_exposure)>=0.01 &
 pmin(x$eaf_outcome,1-x$eaf_outcome)>=0.01 & x$maf_reference>=0.01 &
 x$n_exposure>=0.9*max(x$n_exposure)
xx <- x[keep]; RR <- R[keep,keep,drop=FALSE]
attr(RR,"eigen") <- readRDS(file.path(p,"LD_eigendecomposition_common_highN.rds"))
warns <- character()
more <- withCallingHandlers(susie_rss(z=xx$z_exposure,R=RR,n=round(median(xx$n_exposure)),L=20,
 coverage=0.95,min_abs_corr=0.5,max_iter=1000,tol=0.001,estimate_residual_variance=FALSE),
 warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
saveRDS(more,file.path(p,"exposure_common_highN_L20.rds"),compress="xz")
capacity <- data.table(analysis_set="common_highN",trait="exposure",nsnps=nrow(xx),N=round(median(xx$n_exposure)),
 L=20,converged=more$converged,niter=more$niter,credible_sets=length(more$sets$cs),
 reaches_L_cap=length(more$sets$cs)>=20,warnings=paste(unique(warns),collapse=" | "))
fwrite(capacity,file.path(p,"exposure_L20_capacity_sensitivity.tsv"),sep="\t")
abf <- fread(file.path(p, "coloc_abf_matched_sets_prior_sensitivity.tsv"))
status <- fread(file.path(p, "coloc_susie_status.tsv"))
for_supplement <- merge(fit[, .(analysis_set,trait,nsnps,n_rss,converged,niter,n_credible_sets,inferential_status)],
                       diag[,.(analysis_set,trait,lambda_null_mle,max_abs_standardized_conditional_residual,
                               max_abs_z,min_p)], by=c("analysis_set","trait"),sort=FALSE)
for_supplement[, min_p_display:=fifelse(min_p==0,"<1e-300",format(min_p,scientific=TRUE,digits=6))]
fwrite(for_supplement, file.path(p,"SupplementaryTable_multisignal_fine_mapping.tsv"),sep="\t")

summary_text <- c(
 "Multi-signal colocalization: auditable results",
 "",
 "Signed LD was reconstructed from 503 European 1000 Genomes participants for 2,914 regional variants shared by both summary datasets; zero reference genotypes were missing.",
 "All effect alleles, dosage alleles, SNP positions, LD signs, and SNP ordering were checked explicitly.",
 "SCALLOP per-variant N ranged from 6,907 to 20,507.9; the fine-mapping approximation uses median regional N.",
 "Both traits were analysed with SuSiE-RSS, L=10, 95% credible sets, min |r|=0.5, fixed residual variance, tol=0.001, max_iter=1,000, estimated prior variance.",
 "Sepsis uses effective case-control N=65,958.1121002409 (17,133 cases; 439,048 controls).",
 "",
 "Convergence / credible sets:", capture.output(print(fit, digits=10)),
 "", "LD–z diagnostics:", capture.output(print(diag,digits=10)),
 "", "Signal-pair estimability:", capture.output(print(status)),
 "", "Ordinary single-signal coloc (same variant sets), priors p1=p2=1e-4:", capture.output(print(abf,digits=10)),
 "", "Binary-outcome N / L sensitivity:", capture.output(print(sens,digits=10)),
 "", "Exposure L=20 capacity sensitivity:", capture.output(print(capacity,digits=10)),
 "",
 "Interpretation: absence of a sepsis credible set provides no signal pair for multi-signal colocalization. Report PP.H4 as not estimable, not zero. Weak local sepsis association and limitations of external LD prevent exclusion of a shared causal effect.",
 "The nonconverged exposure fit and substantial exposure z–LD mismatch must be disclosed; raw PIPs / candidate sets from a failed fit are audit outputs, not validated causal signals.",
 "The common_highN sensitivity requires MAF>=0.01 in both traits and reference and exposure N>=90% of maximum regional N; no outcome P-value filtering or relaxed coverage was used.",
 "The 503-person broader-European reference is not in-sample SCALLOP LD or Finnish-specific LD; population mismatch and finite-panel uncertainty limit fine-mapping reliability.",
 "",
 "Reproducibility: 01_multisignal_coloc.R, 02_validate_and_summarize.R, matched_regional_summary_statistics.tsv, signed_LD_EUR_503.rds, LD_SNP_order.txt, recorded fit/diagnostic objects, and session_info.txt.",
 "Intermediate individual-level reference_region.* files are not included in this repository; LD eigendecompositions are regenerable, and matched summary statistics plus signed LD are sufficient to rerun modelling."
)
writeLines(summary_text,file.path(p,"results_summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(p,"session_info_validation.txt"))
print(for_supplement,digits=10)
print(sens,digits=10)
