library(tidyverse)
library(nimble)

strand_data = read.csv("Data/2025/Strandings_flagged_geo_corrected.csv")
strand_data = strand_data |> filter(!(Region %in% c("Sri Lanka", "Pakistan", "Unknown"))) |> 
  filter(!(Species_co %in% c("Unidentified", "Dugong", ""))) |> 
  mutate(coast = case_when(
    Region %in% c("Gujarat", "Maharashtra", "Goa", "Karnataka", "Kerala") ~ "West coast",
    Region %in% c("Tamil Nadu", "Andhra Pradesh", "Odisha", "West Bengal") ~ "East coast",
    TRUE ~ "Islands"
  )) |> 
  mutate(Species_co = ifelse(Species_co == "Indo-Pacific bottlenose dolphin", "Bottlenose dolphin", Species_co)) |>
  mutate(species_type = case_when(
    Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise", "Indo-Pacific bottlenose dolphin", "Spinner dolphin", "Striped dolphin", "Small-sized delphinid", "Bottlenose dolphin", "Risso's dolphin", "Long-beaked common dolphin", "Rough-toothed dolphin", "Irrawaddy dolphin") ~ "Dolphins",
    Species_co %in% c("Sperm whale", "Large-sized delphinid", "Cuvier's beaked whale", "Longman's beaked whale", "Killer whale", "Dwarf sperm whale", "False killer whale", "Melon-headed whale", "Pygmy sperm whale","Pygmy killer whale") ~ "Large toothed-whales",
    TRUE ~ "Whales"
  ))

strand_data = strand_data |> 
  filter(flag == FALSE | is.na(flag), Year >= 2017)# |> #good data available only from 2017
  #filter(Region == "Goa") |> 
  #filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise"))

strand_prb = read.csv("Results/2025/str_prob_month_year.csv")
strand_prb = 
  strand_prb |> filter(year>=2017) |> arrange(year, month)

#Month-wise
y_strand_data = strand_data |> 
  #filter(Species_co == "Humpback dolphin") |> 
  #filter(Species_co == "Indo-Pacific finless porpoise") |> 
  group_by(Month, Year) |> 
  summarise(Count = n(), .groups = "drop") |> 
  #filter(!(Year == 2017 & Month %in% 1:6)) |> # Since data is only from 7/17 to 9/25)
  filter(!(Year == 2025 & Month %in% 10:12)) |> 
  complete(Month = 1:12, Year = 2017:2025, fill = list(Count = NA)) |> 
  arrange(Year, Month) |> 
  mutate(Count = if_else(
    (Year > 2017 | (Year == 2017 & Month >= 1)) &
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
    p_disc[c] ~ dunif(0, 1) # prob of reporting varies annually
  }
  
  for(b in 1:j) { # month
    beta[b] ~ dnorm(0, 1) # month effect
  }
  
  for(c in 1:t) {
    for(b in 1:j) {
      log_lambda[b] <- intercept + alpha[c] + beta[b]
      N_avg[b] <- exp(log_lambda[b])
      #N_avg[b] ~ dpois(N[b])
      # y[b,c] ~ dpois(N_avg[b] * str_prb[b,c])
      
      #including p_disc
      ya[b,c] ~ dpois(N_avg[b] * str_prb[b,c])
      y[b,c] ~ dpois(N_avg[b] * p_disc[c])
      
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
  beta = rep(0, dim(y)[1]),
  p_disc = rep(0.5, dim(y)[2])
)
#monitors=c("alpha","beta","intercept","N_avg")
monitors=c("alpha","beta","p_disc","intercept","N_avg")

model=nimbleMCMC(code, data=data, constants=constants, inits=inits, monitors=monitors,
                 niter=20000, nburnin=5000, nchains=3, thin=10, samplesAsCodaMCMC=TRUE)

mcmc.out = coda::mcmc.list(model)
op_sum = summary(mcmc.out1)
op_sum$statistics[1:12,]

dev.new()
pp = plot(mcmc.out)

varnames(mcmc.out)
gelman.diag(mcmc.out[, c("alpha[1]")])
gelman.diag(mcmc.out)
