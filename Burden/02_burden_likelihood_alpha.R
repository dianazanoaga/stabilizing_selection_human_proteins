# Load libraries

library(data.table)
library(ggplot2)

# -----------------------------
# Paths and files
# -----------------------------

per_prot_dir <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/likelihood_zoltan_for_gene_burden/data/per_protein_0.001"

files <- list.files(
  per_prot_dir,
  pattern = "^gene_burden_.*\\.tsv$",
  full.names = TRUE
)

# -----------------------------
# Likelihoods
# -----------------------------

# Full model: alpha free, sigmaG free
logLr <- function(par, b, se, q) {
  a    <- par[1]
  sigG <- exp(par[2])

  sum(dnorm(
    b,
    mean = 0,
    sd   = sqrt(sigG^2 * (q * (1 - q))^a + se^2),
    log  = TRUE
  ))
}

# Alpha = 0, sigmaG free
logLnr <- function(par, b, se, q) {
  sigG <- exp(par[1])

  sum(dnorm(
    b,
    mean = 0,
    sd   = sqrt(sigG^2 + se^2),
    log  = TRUE
  ))
}

# Alpha = 0, sigmaG = 0
logL_sigG0 <- function(b, se) {
  sum(dnorm(
    b,
    mean = 0,
    sd   = se,
    log  = TRUE
  ))
}

# -----------------------------
# Fit function for one protein
# -----------------------------

#dt <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/likelihood_zoltan_for_gene_burden/data/per_protein_0.001/gene_burden_APOB.tsv")
#dt <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/likelihood_zoltan_for_gene_burden/data/per_protein_0.001/gene_burden_KLK14.tsv")


fit_one_protein <- function(f) {

  dt <- fread(f)

  # filter
  dt <- dt[complete.cases(A1FREQ, BETA, SE)]
  if (nrow(dt) < 5) {
    return(data.frame(
      File  = basename(f),
      Pheno = ifelse(length(unique(dt$PROT)) == 1, unique(dt$PROT), NA),
      Alpha = NA_real_,
      pval_alpha = NA_real_,
      pval_sigG  = NA_real_
    ))
  }

  q  <- as.numeric(dt$A1FREQ)
  b  <- as.numeric(dt$BETA)
  se <- as.numeric(dt$SE)

  # -----------------------------
  # Test 1: Alpha significance
  # -----------------------------

  # Null: alpha = 0, sigmaG free
  par_init0 <- log(sd(b, na.rm = TRUE))

  res0 <- optim(
    par     = par_init0,
    fn      = logLnr,
    method  = "Brent",
    control = list(fnscale = -1),
    b  = b,
    se = se,
    q  = q,
    lower = -100,
    upper = 10
  )

  # Full: alpha free, sigmaG free
  par_init <- c(0, res0$par[1])

  res <- optim(
    par     = par_init,
    fn      = logLr,
    method  = "Nelder-Mead",
    control = list(fnscale = -1),
    b  = b,
    se = se,
    q  = q
  )

  alpha_hat <- res$par[1]
  sigG_hat  <- exp(res$par[2])

  L0 <- res0$value
  L1 <- res$value

  LR_alpha   <- 2 * (L1 - L0)
  pval_alpha <- pchisq(LR_alpha, df = 1, lower.tail = FALSE)

  # ------------------------------------------
  # Test 2: SigmaG significance when alpha = 0
  # ------------------------------------------

  # H0: sigmaG = 0
  L0_sigG <- logL_sigG0(b, se)

  # H1: sigmaG free (alpha = 0)
  par_init_sigG <- log(sd(b, na.rm = TRUE))

  res_sigG <- optim(
    par     = par_init_sigG,
    fn      = logLnr,
    method  = "Brent",
    control = list(fnscale = -1),
    b  = b,
    se = se,
    q  = q,
    lower = -100,
    upper = 10
  )

  L1_sigG <- res_sigG$value
  LR_sigG <- 2 * (L1_sigG - L0_sigG)

  # pval 
  pval_sigG <- pchisq(LR_sigG, df = 1, lower.tail = FALSE)

  sigG_hat_alpha0 <- exp(res_sigG$par[1])

  # -----------------------------
  # Output
  # -----------------------------

  prot_name <- if (length(unique(dt$PROT)) == 1) unique(dt$PROT) else NA

  data.frame(
    File  = basename(f),
    Pheno = prot_name,
    Alpha = alpha_hat,
    sigG_hat_full   = sigG_hat,
    sigG_hat_alpha0 = sigG_hat_alpha0,
    pval_alpha = pval_alpha,
    pval_sigG  = pval_sigG,
    LR_alpha   = LR_alpha,
    LR_sigG    = LR_sigG,
    logSigG    = res$par[2],
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Run over all files
# -----------------------------
res_list <- lapply(files, fit_one_protein)
res_all  <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

# -----------------------------
# Multiple testing correction
# -----------------------------
res_all[, pval_alpha_adj := p.adjust(pval_alpha, method = "BH")]
res_all[, pval_sigG_adj  := p.adjust(pval_sigG,  method = "BH")]

# -----------------------------
# Save results
# -----------------------------
out_path <- file.path(per_prot_dir, "alpha_sigmaG_results_all_proteins_0.001_correction_optim.tsv")
fwrite(res_all, out_path, sep = "\t")

# -----------------------------
# Quick summaries
# -----------------------------
summary(res_all$Alpha)

table(res_all$pval_alpha_adj < 0.05, useNA = "ifany")
table(res_all$pval_sigG_adj  < 0.05, useNA = "ifany")

# -----------------------------
# Plot Alpha distribution
# -----------------------------
ggplot(res_all[!is.na(Alpha)], aes(x = Alpha)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(
    title = "Distribution of estimated Alpha across proteins",
    x = "Alpha",
    y = "Count"
  )

