# 4_Modelling_mortality.R
# Bayesian monthly mortality model for cetacean bycatch using strandings and
# estimated stranding probabilities.
#
# Model structure:
#   y[b, c] ~ Poisson(N_avg[b, c] * str_prb[b, c])
#   log(N_avg[b, c]) = intercept + alpha[c] + beta[b]
#
# where:
#   - alpha[c] is the year-specific effect,
#   - beta[b] is the month-specific effect,
#   - str_prb[b, c] is the estimated probability that a mortality event in that
#     month-year cell is observed as a stranding.

suppressPackageStartupMessages({
  library(tidyverse)
  library(nimble)
  library(coda)
})

# ----------------------------------------------------------------------------
# 1. Read data and filter to the study window
# ----------------------------------------------------------------------------
strand_data <- read.csv("Data/Strandings_flagged_geo_corrected.csv") |>
  filter(flag == FALSE | is.na(flag), Year >= 2017) |>
  filter(Region == "Goa") |>
  filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise"))

strand_prob <- read.csv("Results/str_prob_month_year_g5km.csv") |>
  filter(year >= 2017) |>
  arrange(year, month)

# ----------------------------------------------------------------------------
# 2. Define helper functions
# ----------------------------------------------------------------------------
prepare_stranding_counts <- function(dat, species_name) {
  dat |>
    filter(Species_co == species_name) |>
    group_by(Month, Year) |>
    summarise(Count = n(), .groups = "drop") |>
    filter(!(Year == 2017 & Month %in% 1:6)) |>
    filter(!(Year == 2025 & Month %in% 10:12)) |>
    complete(Month = 1:12, Year = 2017:2025, fill = list(Count = 0)) |>
    arrange(Year, Month) |>
    mutate(Count = replace_na(Count, 0)) |>
    pull(Count) |>
    matrix(nrow = 12, ncol = 9)
}

prepare_stranding_probability <- function(prob_df) {
  prob_df |>
    complete(month = 1:12, year = 2017:2025, fill = list(mean_prob = 0)) |>
    arrange(year, month) |>
    mutate(mean_prob = if_else(mean_prob == 1, 0.99, mean_prob)) |>
    mutate(mean_prob = replace_na(mean_prob, 0)) |>
    pull(mean_prob) |>
    matrix(nrow = 12, ncol = 9)
}

# The stranding probability matrix is constructed so every month-year cell is
# represented, with a simple fill for the final partial year when the last data
# window is incomplete.
fill_final_year_probability <- function(str_prb) {
  str_prb[10:12, ncol(str_prb)] <- apply(str_prb, 1, mean, na.rm = TRUE)[10:12]
  str_prb
}

fit_monthly_mortality_model <- function(
  y,
  str_prb,
  niter = 20000,
  nburnin = 5000,
  nchains = 3,
  thin = 10
) {
  code <- nimbleCode({
    intercept ~ dnorm(0, 1)

    for (c in 1:t) {
      alpha[c] ~ dnorm(0, 1)
    }

    for (b in 1:j) {
      beta[b] ~ dnorm(0, 1)
    }

    for (c in 1:t) {
      for (b in 1:j) {
        log_lambda[b, c] <- intercept + alpha[c] + beta[b]
        N_avg[b, c] <- exp(log_lambda[b, c])
        y[b, c] ~ dpois(N_avg[b, c] * str_prb[b, c])
      }
    }
  })

  constants <- list(j = nrow(y), t = ncol(y), str_prb = str_prb)
  data_list <- list(y = y)
  inits <- list(
    intercept = 0,
    alpha = rep(0, ncol(y)),
    beta = rep(0, nrow(y))
  )

  model <- nimbleMCMC(
    code = code,
    data = data_list,
    constants = constants,
    inits = inits,
    monitors = c("alpha", "beta", "intercept", "N_avg"),
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    thin = thin,
    samplesAsCodaMCMC = TRUE
  )

  model
}

summarise_posterior <- function(mcmc_out) {
  posterior_summary <- summary(mcmc_out)

  out <- data.frame(
    parameter = rownames(posterior_summary$statistics),
    mean = posterior_summary$statistics[, "Mean"],
    sd = posterior_summary$statistics[, "SD"],
    q2_5 = posterior_summary$quantiles[, "2.5%"],
    q97_5 = posterior_summary$quantiles[, "97.5%"],
    row.names = NULL
  )

  out
}

posterior_predictive_check <- function(
  mcmc_out,
  y,
  str_prb,
  n_draws = 500,
  seed = 123
) {
  set.seed(seed)

  mcmc_mat <- as.matrix(mcmc_out)
  j <- nrow(y)
  t <- ncol(y)

  yrep <- replicate(n_draws, {
    theta <- mcmc_mat[sample(1:nrow(mcmc_mat), 1), ]

    log_lambda <- outer(
      theta[paste0("beta[", 1:j, "]")],
      theta[paste0("alpha[", 1:t, "]")],
      function(b, a) theta["intercept"] + a + b
    )

    N_avg <- exp(log_lambda)
    matrix(rpois(j * t, N_avg * str_prb), nrow = j, ncol = t)
  })

  total_diff <- apply(yrep, 3, function(m) sum(m - y, na.rm = TRUE))
  prop_predicted_gt_obs <- length(which(total_diff >= 0)) / length(total_diff)

  list(
    replicated_datasets = yrep,
    total_diff = total_diff,
    prop_predicted_gt_obs = prop_predicted_gt_obs
  )
}

# ----------------------------------------------------------------------------
# 3. Fit the model for each focal species
# ----------------------------------------------------------------------------
model_results <- map(
  c("Indo-Pacific finless porpoise", "Humpback dolphin"),
  function(species_name) {
    y_counts <- strand_data |>
      prepare_stranding_counts(species_name = species_name)

    str_prb <- strand_prob |>
      prepare_stranding_probability() |>
      fill_final_year_probability()

    mcmc_out <- fit_monthly_mortality_model(y_counts, str_prb)

    annual_mortality <- as.matrix(mcmc_out)[, grep(
      "^N_avg\\[",
      colnames(as.matrix(mcmc_out))
    )]

    posterior_summary <- summarise_posterior(mcmc_out)
    ppc <- posterior_predictive_check(mcmc_out, y_counts, str_prb)

    list(
      species = species_name,
      y = y_counts,
      str_prb = str_prb,
      mcmc_out = mcmc_out,
      posterior_summary = posterior_summary,
      annual_mortality = annual_mortality,
      annual_mortality_mean = mean(rowSums(annual_mortality)),
      annual_mortality_sd = sd(rowSums(annual_mortality)),
      annual_mortality_credible_interval = quantile(
        rowSums(annual_mortality),
        c(0.05, 0.975)
      ),
      ppc = ppc
    )
  }
)

names(model_results) <- c("IPFP", "IOHD")

# ----------------------------------------------------------------------------
# 4. Diagnostics and model checking
# ----------------------------------------------------------------------------
for (species_name in names(model_results)) {
  mcmc_out <- model_results[[species_name]]$mcmc_out

  cat("\n##", species_name, "##\n")
  cat("Effective sample size (first 12 monitored parameters):\n")
  print(effectiveSize(mcmc_out)[1:12])

  cat("\nGelman-Rubin PSRF (first 12 monitored parameters):\n")
  print(gelman.diag(mcmc_out, multivariate = FALSE)$psrf[1:12, 1])
}

# ----------------------------------------------------------------------------
# 5. Create a summary plot of monthly mortality estimates
# ----------------------------------------------------------------------------
monthly_mortality <- map_dfr(
  model_results,
  function(result) {
    mcmc_mat <- as.matrix(result$mcmc_out)
    n_avg <- mcmc_mat[, grep("^N_avg\\[", colnames(mcmc_mat))]

    month_means <- apply(n_avg, 2, mean)
    month_q05 <- apply(n_avg, 2, quantile, probs = 0.05, na.rm = TRUE)
    month_q95 <- apply(n_avg, 2, quantile, probs = 0.975, na.rm = TRUE)

    tibble(
      species = result$species,
      month = rep(1:12, times = ncol(result$y)),
      year = rep(2017:2025, each = 12),
      mean = as.numeric(month_means),
      q05 = as.numeric(month_q05),
      q95 = as.numeric(month_q95)
    )
  }
)

monthly_mortality |>
  ggplot(aes(x = as.factor(month), y = mean, fill = species)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_errorbar(
    aes(ymin = q05, ymax = q95),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  theme_minimal() +
  facet_wrap(~species, ncol = 1) +
  labs(
    x = "Month",
    y = "Estimated mortality at-sea\n(mean ± 95% credible interval)",
    fill = "Species"
  )

# ----------------------------------------------------------------------------
# 6. Annual mortality summary
# ----------------------------------------------------------------------------
annual_summary <- map_dfr(
  model_results,
  function(result) {
    tibble(
      species = result$species,
      mean_total = result$annual_mortality_mean,
      sd_total = result$annual_mortality_sd,
      q05 = result$annual_mortality_credible_interval[1],
      q97_5 = result$annual_mortality_credible_interval[2]
    )
  }
)

annual_summary

# ----------------------------------------------------------------------------
# 7. Posterior predictive check summary
# ----------------------------------------------------------------------------
ppc_summary <- map_dfr(
  model_results,
  function(result) {
    tibble(
      species = result$species,
      prop_predicted_gt_obs = result$ppc$prop_predicted_gt_obs
    )
  }
)

ppc_summary

# ----------------------------------------------------------------------------
# 8. Optional output files
# ----------------------------------------------------------------------------
# write.csv(annual_summary, "Results/Model summaries/annual_mortality_summary.csv")
# write.csv(ppc_summary, "Results/Model summaries/posterior_predictive_checks.csv")
