Modelling mortality
================

# Overview

This document estimates the monthly at-sea mortality for the two
cetacean species in Goa using a Poisson model that adjusts observed
strandings for the probability that an animal mortality event is
detected as a stranding. The analysis combines monthly stranding counts
with estimated stranding probabilities for each month-year cell and
interprets the resulting posterior estimates as latent mortality risk.

The workflow is structured in three parts:

1.  prepare monthly stranding counts and stranding-probability matrices,
2.  fit the Bayesian monthly mortality model, and
3.  evaluate posterior summaries and predictive fit.

## 1. Data preparation

The first section imports the cleaned stranding dataset and the monthly
stranding-probability output from the earlier stages of the analysis. It
filters the records to Goa, keeps only valid records from 2017 onward,
and retains the two focal species: Indo-Pacific finless porpoise and
Indian Ocean humpback dolphin.

### Data import and filtering

``` r
library(tidyverse)
library(nimble)
library(coda)

strand_data <- read.csv("Data/Strandings_flagged_geo_corrected.csv") |>
  filter(flag == FALSE | is.na(flag), Year >= 2017) |>
  filter(Region == "Goa") |>
  filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise"))

strand_prob <- read.csv("Results/str_prob_month_year_g5km.csv") |>
  filter(year >= 2017) |>
  arrange(year, month)
```

This step ensures that the model uses only the time window and the
taxonomic scope relevant to the bycatch analysis. It also ensures that
the mortality model is aligned with the estimated stranding
probabilities computed in the previous analytical step.

### Monthly stranding arrays

Monthly strandings are arranged into a 12 x 9 matrix, where rows are
months and columns are years from 2017 to 2025. The code also fills the
partial first and last years so that each month-year cell has an
explicit value for the model.

``` r
prepare_stranding_counts <- function(dat, species_name) {
  matrix(
    dat |>
      filter(Species_co == species_name) |>
      group_by(Month, Year) |>
      summarise(Count = n(), .groups = "drop") |>
      filter(!(Year == 2017 & Month %in% 1:6)) |>
      filter(!(Year == 2025 & Month %in% 10:12)) |>
      complete(Month = 1:12, Year = 2017:2025, fill = list(Count = NA)) |>
      arrange(Year, Month) |>
      mutate(
        Count = if_else(
          (Year > 2017 | (Year == 2017 & Month >= 7)) &
            (Year < 2025 | (Year == 2025 & Month <= 9)) &
            is.na(Count),
          0,
          Count
        )
      ) |>
      arrange(Year, Month) |>
      pull(Count),
    nrow = 12,
    ncol = 9
  )
}

prepare_stranding_probability <- function(prob_df) {
  prob_df |>
    complete(month = 1:12, year = 2017:2025, fill = list(mean_prob = NA)) |>
    arrange(year, month) |>
    mutate(mean_prob = if_else(mean_prob == 1, 0.99, mean_prob)) |>
    pull(mean_prob) |>
    matrix(nrow = 12, ncol = 9)
}

fill_final_year_probability <- function(str_prb) {
  str_prb[10:12, ncol(str_prb)] <- apply(str_prb, 1, mean, na.rm = TRUE)[10:12]
  str_prb
}
```

### Result

The transformed data are ready for the model as month-by-year matrices,
ensuring that the monthly effect and the year effect are identifiable
from a complete seasonal and annual structure.

## 2. Model specification

The mortality model treats observed strandings as a Poisson process with
an offset determined by the monthly stranding probability. The latent
expected mortality is modelled as a log-linear function of a global
intercept, a year effect, and a month effect.

### Bayesian model

``` r
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

  nimbleMCMC(
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
}
```

The Poisson likelihood is appropriate for count-based monthly
strandings, while the stranding probability acts as an observation
process that scales the latent mortality by the probability that a death
becomes a reported stranding. A prior normal(0, 1) is assigned to the
intercept, month effects, and year effects to maintain weakly
informative regularization on the latent scale.

## 3. Model fitting for the two focal species

``` r
model_results <- map(
  c("Indo-Pacific finless porpoise", "Humpback dolphin"),
  function(species_name) {
    y_counts <- prepare_stranding_counts(strand_data, species_name)
    str_prb <- prepare_stranding_probability(strand_prob) |>
      fill_final_year_probability()

    mcmc_out <- fit_monthly_mortality_model(y_counts, str_prb)
    annual_mortality <- as.matrix(mcmc_out)[, grep("^N_avg\\[", colnames(as.matrix(mcmc_out)))]

    list(
      species = species_name,
      y = y_counts,
      str_prb = str_prb,
      mcmc_out = mcmc_out,
      annual_mortality = annual_mortality,
      annual_mortality_mean = mean(rowSums(annual_mortality)),
      annual_mortality_sd = sd(rowSums(annual_mortality)),
      annual_mortality_credible_interval = quantile(
        rowSums(annual_mortality),
        c(0.05, 0.975)
      )
    )
  }
)

names(model_results) <- c("IPFP", "IOHD")
```

### Result

The model is fit separately for each focal species, and the posterior
draws are retained for subsequent summarization and checking. This
allows species-specific estimates of mortality to be compared while
holding the same observation model and month-year structure.

## 4. Diagnostics and convergence checks

The next section inspects mixing and convergence for the MCMC chains.

``` r
for (species_name in names(model_results)) {
  mcmc_out <- model_results[[species_name]]$mcmc_out

  cat("\n##", species_name, "##\n")
  cat("Effective sample size (first 12 monitored parameters):\n")
  print(effectiveSize(mcmc_out)[1:12])

  cat("\nGelman-Rubin PSRF (first 12 monitored parameters):\n")
  print(gelman.diag(mcmc_out, multivariate = FALSE)$psrf[1:12, 1])
}
```

### Result

The effective sample sizes and PSRF values indicate whether the
posterior samples are behaving sensibly. When these statistics are
acceptable, the uncertainty summaries can be interpreted as reliable
estimates of latent monthly mortality.

## 5. Monthly mortality summaries

The posterior draws for the latent mortality surface are summarised as
month-wise estimates and their 95% credible intervals.

``` r
monthly_mortality <- map_dfr(
  model_results,
  function(result) {
    mcmc_mat <- as.matrix(result$mcmc_out)
    n_avg <- mcmc_mat[, grep("^N_avg\\[", colnames(mcmc_mat))]

    month_mean <- apply(n_avg, 2, mean)
    month_q05 <- apply(n_avg, 2, quantile, probs = 0.05, na.rm = TRUE)
    month_q95 <- apply(n_avg, 2, quantile, probs = 0.975, na.rm = TRUE)

    tibble(
      species = result$species,
      month = rep(1:12, times = ncol(result$y)),
      year = rep(2017:2025, each = 12),
      mean = as.numeric(month_mean),
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
```

### Result

This figure provides the estimated monthly mortality pattern for each
species, with uncertainty represented by the 5th and 97.5th posterior
percentiles. The plot is most useful for identifying whether mortality
risk is concentrated in a particular season or comes from a more diffuse
pattern across the year.

## 6. Annual mortality estimates

The annual mortality estimate is derived by summing the latent monthly
mortality values across all months for each posterior draw. The
resulting summary provides a total annual mortality estimate and an
uncertainty interval.

``` r
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
```

### Result

The total annual mortality estimate summarises the expected number of
deaths at sea per year after accounting for imperfect stranding
detection. These numbers should be interpreted as latent mortality
estimates rather than raw stranding counts, because they incorporate the
estimated detection process in the model.

## 7. Posterior predictive checks

Posterior predictive checks evaluate whether the model generates data
that resemble the observed monthly strandings. This is done by
simulating replicated datasets under the posterior and comparing their
total counts to the observed totals.

``` r
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

  tibble(
    prop_predicted_gt_obs = prop_predicted_gt_obs
  )
}

ppc_summary <- map_dfr(
  model_results,
  function(result) {
    posterior_predictive_check(result$mcmc_out, result$y, result$str_prb) |>
      mutate(species = result$species)
  }
)

ppc_summary
```

### Result

The posterior predictive proportion indicates how often the replicated
datasets produce total counts at or above the observed values. Values
close to 0.5 suggest that the model is not systematically over- or
under-predicting total mortality, while more extreme values indicate a
potential mismatch between the model and the observed data.

## 8. Interpretation

The model estimates latent monthly mortality after accounting for
imperfect detection of dead animals as strandings. This is a necessary
adjustment because the observed stranding counts underestimate total
mortality whenever there is a probability that deaths are not recovered
or recorded.

The resulting monthly and annual estimates are therefore useful for
evaluating seasonal mortality risk and for comparing the relative burden
across the two focal species. The interpretation should focus on the
posterior uncertainty and the biological plausibility of the monthly
pattern rather than on the raw strandings alone.
