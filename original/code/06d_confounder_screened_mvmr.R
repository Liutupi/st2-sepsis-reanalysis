#!/usr/bin/env Rscript
# Prespecified direct-variant GWAS Catalog screen sensitivity for MVMR.
# This is a weighted two-exposure regression at zero exposure-association
# covariance, matching the primary analysis' main covariance assumption.

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(Sys.getenv("SST2_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
screen <- read.delim(file.path(project_dir, "results/instrument_validity/instrument_confounder_screen.tsv"),
                     check.names = FALSE)
d <- read.delim(file.path(project_dir, "results/multivariable_mr/harmonized_mvmr_instruments.tsv"),
                check.names = FALSE)
out <- file.path(project_dir, "results/instrument_validity/confounder_screened_mvmr_sensitivity.tsv")

keep <- screen$rsid[
  screen$in_mvmr_set == "yes" &
    screen$potential_confounder_flags == "none identified"
]
d <- d[d$rsid %in% keep & d$hf_ok & d$sepsis_ok, ]
if (nrow(d) < 5L) stop("Too few screened MVMR instruments")

X <- as.matrix(d[, c("beta_sst2", "beta_hf")])
y <- d$beta_sepsis
w <- 1 / d$se_sepsis^2
inv <- solve(crossprod(X, w * X))
beta <- as.numeric(inv %*% crossprod(X, w * y))
resid <- y - as.numeric(X %*% beta)
Q <- sum(w * resid^2)
df <- nrow(d) - ncol(X)
phi <- max(1, Q / df)
se <- sqrt(diag(inv) * phi)
z <- beta / se

result <- data.frame(
  exposure = c("Circulating sST2", "Heart-failure liability"),
  outcome = "FinnGen R12 sepsis",
  nsnp = nrow(d),
  retained_rsids = paste(sort(d$rsid), collapse = "; "),
  beta = beta,
  se = se,
  OR = exp(beta),
  CI_lower = exp(beta - 1.96 * se),
  CI_upper = exp(beta + 1.96 * se),
  p = 2 * pnorm(-abs(z)),
  residual_Q = Q,
  residual_Q_df = df,
  residual_Q_p = pchisq(Q, df = df, lower.tail = FALSE),
  screen_rule = "Excluded direct GWAS Catalog associations at P<5e-8 within prespecified potential-confounder categories; no proxy screen.",
  check.names = FALSE
)
write.table(result, out, sep = "\t", row.names = FALSE, quote = FALSE)
print(result)
