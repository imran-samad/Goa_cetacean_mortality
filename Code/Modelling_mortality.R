library(tidyverse)
library(nimble)
library(coda)

strand_data = read.csv("Data/2025/Strandings_flagged_geo_corrected.csv")
strand_data = strand_data |> 
  filter(flag == FALSE | is.na(flag), Year >= 2017) |> #good data available only from 2017
  filter(Region == "Goa") |> 
  filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise"))

strand_prb = read.csv("Results/2025/str_prob_month_year_g5km.csv")
# strand_prb = read.csv("Results/2025/str_prob_month_year_g10km.csv")
strand_prb = 
  strand_prb |> filter(year>=2017) |> arrange(year, month)

#Month-wise
y_strand_data = strand_data |> 
  #filter(Species_co == "Humpback dolphin") |> 
  filter(Species_co == "Indo-Pacific finless porpoise") |> 
  group_by(Month, Year) |> 
  summarise(Count = n(), .groups = "drop") |> 
  filter(!(Year == 2017 & Month %in% 1:6)) |> # Since data is only from 7/17 to 9/25)
  filter(!(Year == 2025 & Month %in% 10:12)) |> 
  complete(Month = 1:12, Year = 2017:2025, fill = list(Count = NA)) |> 
  arrange(Year, Month) |> 
  mutate(Count = if_else(
      (Year > 2017 | (Year == 2017 & Month >= 7)) &
      (Year < 2025 | (Year == 2025 & Month <= 9)) &
      is.na(Count),
    0,
    Count
  ))

y_strand_data = array(y_strand_data$Count, c(length(unique(y_strand_data$Month)), length(unique(y_strand_data$Year))))

str_prb = strand_prb |> 
  #filter(!(year == 2017 & month %in% 1:6)) |> # Since data is only from 7/17 to 9/25)
  #filter(!(year == 2025 & month %in% 10:12)) |> 
  complete(month = 1:12, year = 2017:2025, fill = list(mean_prob = NA)) |> 
  arrange(year, month) |> 
  mutate(mean_prob = ifelse(mean_prob == 1, 0.99, mean_prob))
str_prb = array(str_prb$mean_prob, c(length(unique(str_prb$month)), length(unique(str_prb$year))))
str_prb[10:12,9] = apply(str_prb, 1, mean, na.rm=T)[10:12]
  

# Main code
op=data.frame()
## define the model
code = nimbleCode({
  # Priors
  intercept ~ dnorm(0, 1)
  
  for(c in 1:t) { # year
    alpha[c] ~ dnorm(0, 1) # year effect
    #p_disc[c] ~ dunif(0, 1) # prob of reporting varies annually
  }
  
  for(b in 1:j) { # month
    beta[b] ~ dnorm(0, 1) # month effect
  }
  
  for(c in 1:t) {
    for(b in 1:j) {
      log_lambda[b] <- intercept + alpha[c] + beta[b]
      N_avg[b] <- exp(log_lambda[b])
      #N_avg[b] ~ dpois(N[b])
      y[b,c] ~ dpois(N_avg[b] * str_prb[b,c] * 1)
      
      #including p_disc
      # ya[b,c] ~ dpois(N_avg[b] * str_prb[b,c])
      # y[b,c] ~ dpois(N_avg[b] * p_disc[c])
      
      #y[b,c] ~ dbinom(N_avg[b], prob = str_prb[b,c])
      #y[b,c] ~ dnegbin(prob = str_prb[b,c], N_avg[b])
    }
  }
})

y=y_strand_data
constants = list(j=dim(y)[1],t=dim(y)[2],str_prb=str_prb)
data = list(y = y)
inits = list(
  intercept = 0,
  alpha = rep(0, dim(y)[2]),
  beta = rep(0, dim(y)[1])
  #p_disc = rep(0.5, dim(y)[2])
)
monitors=c("alpha","beta","intercept","N_avg")
#monitors=c("alpha","beta","p_disc","intercept","N_avg")

model=nimbleMCMC(code, data=data, constants=constants, inits=inits, monitors=monitors,
                 niter=20000, nburnin=5000, nchains=3, thin=10, samplesAsCodaMCMC=TRUE)

# mcmc.out.iohd = coda::mcmc.list(model)
# mcmc.out.ipfp = coda::mcmc.list(model)

mcmc.out1 = coda::mcmc.list(model)
op_sum1 = summary(mcmc.out.iohd)
op_sum1$statistics[1:12,]
sum(op_sum1$statistics[1:12,1])
sd(op_sum1$statistics[1:12,1])
apply(op_sum1$quantiles[1:12,c(1,5)], 2, sum)

#Diagnostic metrics
effectiveSize(mcmc.out1)[1:12]
gelman.diag(mcmc.out1, multivariate = FALSE)$psrf[1:12,1]
mod_sum_iohd = MCMCvis::MCMCsummary(mcmc.out.iohd)
mod_sum_ipfp = MCMCvis::MCMCsummary(mcmc.out.ipfp)
# op_sum1$statistics |> write.csv("Results/2025/Model summaries/stranding_model_monthly_IOHD.csv")
# op_sum1$statistics |> write.csv("Results/2025/Model summaries/stranding_model_monthly_IPFP.csv")

dev.new()
#mcmc.out.iohd = mcmc.out1
#mcmc.out.ipfp = mcmc.out1
pp = plot(mcmc.out1)

varnames(mcmc.out1)
gelman.diag(mcmc.out1)#[, c("alpha[1]")])
gelman.diag(mcmc.out)

#Mortality trends
op_sum1 = summary(mcmc.out.iohd)
op_sum1 = data.frame(cbind(op_sum1$statistics[1:12,], op_sum1$quantiles[1:12,]), 'Species' = "Indian Ocean humpback dolphin")
op_sum2 = summary(mcmc.out.ipfp)
op_sum2 = data.frame(cbind(op_sum2$statistics[1:12,], op_sum2$quantiles[1:12,]), 'Species' = "Indo-Pacific finless porpoise")
op_sum1 |> rbind(op_sum2) |> mutate(month = rep(1:12,2)) |> 
  ggplot(aes(x = as.factor(month), y = Mean)) +
  geom_bar(stat = "identity", fill = "gray", alpha = 1) +
  geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width = 0.2) +
  #geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),
        strip.text = element_text(size = 14)) +
  facet_wrap(~Species, ncol = 1) +
  #labs(x = "Month", y = "Estimated mortality at-sea\n(mean ± SD)")
  labs(x = "Month", y = "Estimated mortality at-sea\n(mean ± 95% CI)")

ggsave("Results/2025/Graphs/Estimated_mortality_monthly_IOHD_IPFP_95CI.jpg", width = 6, height = 8)

#annual mortality estimates
mcmc_mat <- as.matrix(mcmc.out.ipfp)
N_avg_vals <- (mcmc_mat[, grep("^N_avg\\[", colnames(mcmc_mat))])
mean(rowSums(N_avg_vals))
sd(rowSums(N_avg_vals))
quantile(N_avg_vals, c(0.05,0.975))

####################
#Identifying mass stranding events
library(tidyverse)
library(sf)
library(lubridate)
library(geosphere)
library(igraph)

# Make sure your date columns are properly formatted
strand_data = strand_data |>
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) |> filter(!is.na(lat))

# Convert to sf object (set CRS = WGS84)
strand_sf = st_as_sf(
  strand_data,
  coords = c("lon", "lat"),
  crs = 4326
)

# Transform to a projected CRS in metres for distance calculation
strand_sf = st_transform(strand_sf, 32643)  # UTM zone 43N (covers Goa region)

# Calculate pairwise spatial distances (in metres)
dist_matrix = st_distance(strand_sf)

# Convert to numeric matrix (it’s an "units" object initially)
dist_matrix = as.numeric(dist_matrix)
dist_matrix = matrix(dist_matrix, nrow = nrow(strand_sf))

# Calculate pairwise temporal differences (in days)
time_diff = abs(outer(strand_sf$Date, strand_sf$Date, "-"))

# Identify pairs that meet both criteria
close_pairs = which(dist_matrix < 5000 & time_diff < 2, arr.ind = TRUE)
close_pairs = close_pairs[close_pairs[,1] < close_pairs[,2], ]  # avoid duplicates/self-pairs

# Inspect
close_pairs_df = data.frame(
  id1 = strand_sf$id[close_pairs[,1]],
  id2 = strand_sf$id[close_pairs[,2]],
  Death_code1 = strand_sf$Death_code[close_pairs[,1]],
  Death_code2 = strand_sf$Death_code[close_pairs[,1]],
  Day = strand_sf$Day[close_pairs[,1]],
  Month = strand_sf$Month[close_pairs[,1]],
  Year = strand_sf$Year[close_pairs[,1]],
  dist_m = dist_matrix[close_pairs],
  days_diff = time_diff[close_pairs]
)

close_pairs_df
unique(c(close_pairs_df$id1,close_pairs_df$id2))
#manually clustering strandings based on the above output

# Create graph from pairs
g <- graph_from_data_frame(close_pairs_df[, c("id1", "id2")], directed = FALSE)

# Find connected components (clusters)
components <- components(g)

# Assign cluster ID to each row
close_pairs_df$cluster <- components$membership[as.character(close_pairs_df$id1)]

# Summarize clusters
cluster_summary <- close_pairs_df %>%
  group_by(cluster) %>%
  summarise(
    pairs_in_cluster = n(),
    avg_day = mean(Day),
    avg_month = mean(Month),
    avg_year = mean(Year),
    members = paste(unique(c(id1, id2)), collapse = "-")
  )

print(cluster_summary)

###########Posterior predictive checks##############
library(coda)
library(bayesplot)

# Convert MCMC output to matrix
mcmc_mat <- as.matrix(mcmc.out.ipfp)

# Number of posterior draws to use
n_draws <- 500  
set.seed(123)

# Dimensions
j <- nrow(y)   # months
t <- ncol(y)   # years

# Simulate replicated datasets
yrep <- replicate(n_draws, {
  # sample one posterior draw
  theta <- mcmc_mat[sample(1:nrow(mcmc_mat), 1), ]
  
  # build lambda for each month-year cell
  # alpha is year effect, beta is month effect
  log_lambda <- outer(theta[paste0("beta[",1:j,"]")],
                      theta[paste0("alpha[",1:t,"]")],
                      function(b,a) theta["intercept"] + a + b)
  
  N_avg <- exp(log_lambda)
  
  # simulate strandings using Poisson likelihood
  matrix(rpois(j*t, N_avg * str_prb), nrow = j, ncol = t)
})

# Posterior predictive check: compare observed vs replicated totals
# obs_total <- sum(y, na.rm=T)
# rep_totals <- apply(yrep, 2, sum)
# ppc_pval <- mean(rep_totals >= obs_total)
# cat("Posterior predictive p-value for total strandings:", ppc_pval, "\n")
# 
# # Visualization: overlay observed vs replicated distributions
# ppc_dens_overlay(y = as.vector(y), yrep = t(apply(yrep[1:50, ], 2, as.vector)))

xx=apply(yrep, 3, function(m) sum(m - y, na.rm = TRUE))
length(which(xx>=0))/length(xx) #proportion of replicated datasets where the total predicted stranding count are > observed stranding count
#IOHD:0.37; IPFP: 0.56