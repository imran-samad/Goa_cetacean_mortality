library(tidyverse)
library(sf)
library(unmarked)

goa_grid = read_sf("Data/2025/Shp files/Goa_coast_survey_area_grid.shp")

dol_dat = read.csv("Data/2025/Population/Dolphins_locations_2312.csv")
dol_dat$PDistance_m = dol_dat$Distance_m*sin((abs(dol_dat$Angle)*(pi/180)))
dol_dat$Start = as.Date(dol_dat$Start, "%d-%m-%Y")

dol_cov = read.csv("Data/2025/Population/Dolphins_covariates_2312.csv")
temp_df = data.frame()
for(i in 1:nrow(dol_dat))
{
  row_loc = which.min(abs(dol_dat$Waypoint_number[i] - dol_cov$Waypoint.number))[1]
  temp_df = rbind(temp_df, dol_cov[row_loc,c("Boat.speed..kmph.", "Number.of.fishing.boats..small.medium", "Number.of.fishing.boats..large", "Number.of.tourist.boats")])
}
dol_dat = cbind(dol_dat, temp_df)
dol_dat[,(ncol(dol_dat)-2):ncol(dol_dat)] <- lapply(
  dol_dat[,(ncol(dol_dat)-2):ncol(dol_dat)],
  function(x) as.numeric(replace(x, is.na(x), 0))
)


dol_dat$grid_id = as.data.frame(st_intersects(goa_grid, st_as_sf(dol_dat, coords = c("X", "Y"), crs = 4326)))[,1]

df = dol_dat |> group_by(grid_id, rep) |> summarise(n = sum(N_best), mean_bs = mean(Boat.speed..kmph.), tot_fb = sum(Number.of.fishing.boats..small.medium, Number.of.fishing.boats..large), tot_tb = sum(Number.of.tourist.boats), mean_size = mean(N_best))

#Fill other combinations
# Get unique combinations of the two columns
combo_df <- expand.grid(
  grid_id = unique(df$grid_id),
  rep = unique(df$rep)
)

# Identify other columns
other_cols <- setdiff(names(df), c("grid_id", "rep"))

# Add them as NA
for (col in other_cols) {
  combo_df[[col]] <- NA
}

combo_df$n = 0

new_combos = anti_join(combo_df, df, by = c("grid_id", "rep"))
df_complete = bind_rows(df, new_combos)

#Run the n-mix model

library(tidyr)
library(dplyr)

# reshape counts
y_wide <- df_complete %>%
  select(grid_id, rep, n) %>%
  pivot_wider(names_from = rep, values_from = n, names_prefix = "rep_")

# site covariates (one row per grid)
site_covs <- df_complete %>%
  group_by(grid_id) %>%
  summarize(mean_bs = mean(mean_bs, na.rm = TRUE),
            tot_fb  = mean(tot_fb, na.rm = TRUE),
            tot_tb  = mean(tot_tb, na.rm = TRUE),
            mean_size = mean(mean_size, na.rm = TRUE))

# convert to data frame for unmarked
y_mat <- as.data.frame(y_wide[,-1])
rownames(y_mat) <- y_wide$grid_id

umf <- unmarkedFramePCount(y = y_mat, siteCovs = site_covs)

mod0 <- pcount(~1 ~1, mixture = "NB", data = umf)
summary(mod0)
mod_zinb <- detect::pcount(~1 ~1, mixture = "ZINB", data = umf)

predict(mod0, type = "state")  # mean λ per grid
predict(mod0, type = "det")    # mean p per visit
