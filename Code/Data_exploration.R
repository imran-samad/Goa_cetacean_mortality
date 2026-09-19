# -----------------------------------------------------------------------------
# Goa cetacean mortality: analysis workflow
#
# This script is organized into three main sections:
#   1. Data cleaning and exploratory summaries for strandings
#   2. Congruence between reviewer estimates and self-reported decomposition/PMI
#   3. Mapping fishing and tourist activity from boat surveys
# -----------------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(sf)
library(purrr)

# -----------------------------------------------------------------------------
# 1. STRANDING DATA: cleaning and exploratory summaries
# -----------------------------------------------------------------------------

dat_str <- read.csv("Data/Strandings_flagged_geo_corrected.csv")

# Keep only unflagged records and data from 2017 onward; restrict to the two
# focal species for this study.
dat_str <- dat_str |>
  filter(flag == FALSE | is.na(flag), Year >= 2017) |>
  filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise"))

# Species composition summary.
dat_str |>
  count(Species_co, name = "n") |>
  arrange(desc(n))

# Total number of strandings by month, faceted by species.
dat_str |>
  ggplot(aes(x = as.factor(Month))) +
  geom_bar() +
  facet_wrap(~Species_co, ncol = 1) +
  labs(x = "Month", y = "Total strandings") +
  theme_minimal()

# Mean monthly counts with uncertainty shown as mean ± SE.
mean_monthly_counts <- dat_str |>
  group_by(Year, Month, Species_co) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(Month, Species_co) |>
  summarise(
    n = mean(count),
    n_se = sd(count) / sqrt(length(count)),
    .groups = "drop"
  )

mean_monthly_counts |>
  ggplot(aes(x = as.factor(Month), y = n, colour = Species_co)) +
  geom_pointrange(
    aes(ymin = n - n_se, ymax = n + n_se),
    position = position_dodge(width = 0.5)
  ) +
  theme_minimal() +
  scale_colour_manual(
    values = c(
      "Humpback dolphin" = "lightsteelblue4",
      "Indo-Pacific finless porpoise" = "gray"
    )
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "bottom"
  ) +
  labs(x = "Month", y = "No. of strandings (mean ± SE)", colour = "Species")
# Figure S3 in the manuscript
# ggsave("Results/Graphs/stranding_counts_A1.jpg", width = 8, height = 5)

# Age class summaries. Compare age structure across species and months.
dat_str |>
  filter(Age.Class != "") |>
  group_by(Year, Month, Species_co, Age.Class) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(Month, Species_co, Age.Class) |>
  summarise(n = mean(count), n_sd = sd(count), .groups = "drop") |>
  ggplot() +
  geom_pointrange(aes(
    x = as.factor(Month),
    y = n,
    ymin = n - n_sd,
    ymax = n + n_sd
  )) +
  facet_grid(cols = vars(Species_co), rows = vars(Age.Class)) +
  theme_minimal()

# Proportion of each age class by month and species.
age_class_counts <- dat_str |>
  filter(Age.Class != "") |>
  group_by(Month, Species_co, Age.Class) |>
  summarise(count = n(), .groups = "drop")

age_class_total <- dat_str |>
  filter(Age.Class != "") |>
  group_by(Month, Species_co) |>
  summarise(total = n(), .groups = "drop")

age_class_prop <- age_class_total |>
  left_join(age_class_counts, by = c("Month", "Species_co")) |>
  mutate(prop = count / total) |>
  group_by(Month, Species_co, Age.Class) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  complete(Month, Species_co, Age.Class, fill = list(prop = NA))

age_class_prop |>
  ggplot() +
  geom_bar(
    aes(x = as.factor(Month), y = prop, fill = Age.Class),
    stat = "identity",
    position = position_dodge2(2)
  ) +
  facet_grid(cols = vars(Species_co), rows = vars(Age.Class)) +
  theme_minimal() +
  labs(x = "Month", y = "Proportion", fill = "Age class")

# Death-code summary across the two focal species.
stranding_death_codes <- dat_str |>
  filter(Death_code != "") |>
  group_by(Species_co, Death_code) |>
  count() |>
  group_by(Species_co) |>
  mutate(percentage = n / sum(n) * 100)

stranding_death_codes

# Example decomposition-code histogram. Replace the x-axis label for PMI ranges
# when plotting alternative panels.
p1 <- dat_str |>
  filter(Death_code %in% c(1, 2, 3, 4)) |>
  mutate(across(c(PMI_min, PMI_max), as.numeric)) |>
  mutate(Death_code = as.factor(Death_code)) |>
  ggplot() +
  geom_histogram(
    aes(x = Death_code, fill = Species_co),
    position = "dodge",
    stat = "count"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Humpback dolphin" = "lightsteelblue4",
      "Indo-Pacific finless porpoise" = "gray"
    )
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) +
  labs(x = "Decomposition code", y = "Count", fill = "Species")

# Additional p2/p3 plots can be created analogously for PMI ranges.
# ggpubr::ggarrange(p1, p2, p3, ncol = 1, common.legend = TRUE, legend = "bottom")
# Figure S4 in manuscript
# ggsave("Results/Graphs/stranding_decompcode_pmi.jpg", width = 8, height = 10)

# -----------------------------------------------------------------------------
# 2. REVIEWER ESTIMATES: expert vs experimental decomposition and PMI
# -----------------------------------------------------------------------------

decomp <- read.csv("Data/Decomposition_estimates.csv")

decomp <- decomp |>
  rowwise() |>
  mutate(
    m_dc = mean(c_across(starts_with("DC")), na.rm = TRUE),
    m_pmil = mean(c_across(starts_with("PMIl")), na.rm = TRUE),
    m_pmiu = mean(c_across(starts_with("PMIu")), na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(
    dc_diff = m_dc - Self_dc,
    pmiu_diff = m_pmiu - Self_PMIu,
    pmil_diff = m_pmil - Self_PMIl,
    dc_diff_perc = dc_diff / Self_dc * 100,
    pmiu_diff_perc = pmiu_diff / Self_PMIu * 100,
    pmil_diff_perc = pmil_diff / Self_PMIl * 100
  )

# Summary of mean and spread in percentage differences between expert estimates
# and the average of the experimental ratings.
mean_decomp_summary <- tibble(
  metric = c("DC", "PMIl", "PMIu"),
  mean_percent_diff = c(
    mean(decomp$dc_diff_perc, na.rm = TRUE),
    mean(decomp$pmil_diff_perc, na.rm = TRUE),
    mean(decomp$pmiu_diff_perc, na.rm = TRUE)
  ),
  sd_percent_diff = c(
    sd(decomp$dc_diff_perc, na.rm = TRUE),
    sd(decomp$pmil_diff_perc, na.rm = TRUE),
    sd(decomp$pmiu_diff_perc, na.rm = TRUE)
  )
)

mean_decomp_summary

# -----------------------------------------------------------------------------
# 3. MAPPING FISHING AND TOURISM FROM BOAT SURVEYS
# -----------------------------------------------------------------------------

goa_grid <- st_read("GIS/Goa_buffer_5km_grid.shp")
goa_coast <- st_read("GIS/Goa_coast.shp")

process_gpx_and_s1 <- function(gpx_folders, s1_files, goa_grid) {
  results <- map2_dfr(gpx_folders, s1_files, function(folder, s1_file) {
    s1 <- read.csv(s1_file)

    waypoints_all <- list.files(
      folder,
      pattern = "Waypoints_.*\\.gpx$",
      full.names = TRUE
    ) |>
      map_dfr(
        ~ st_read(.x, layer = "waypoints", quiet = TRUE) |>
          mutate(source_file = basename(.x))
      ) |>
      mutate(name = as.numeric(name)) |>
      select(name, geometry) |>
      right_join(s1, by = c("name" = "Waypoint.number"))

    waypoints_all$grid_id <- st_intersects(waypoints_all, goa_grid) %>%
      as.integer()

    waypoints_all |>
      group_by(grid_id) |>
      summarise(
        FBs = sum(Number.of.fishing.boats..small.medium, na.rm = TRUE),
        FBl = sum(Number.of.fishing.boats..large, na.rm = TRUE),
        TB = sum(Number.of.tourist.boats, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(FB = FBs + FBl) |>
      select(-FBs, -FBl) |>
      st_drop_geometry() |>
      mutate(grid_id = as.numeric(grid_id)) |>
      right_join(goa_grid, by = c("grid_id" = "id")) |>
      filter(!is.na(TB)) |>
      st_as_sf()
  })

  return(results)
}

# Survey folders and corresponding covariate files.
folders <- c(
  "Data/Surveys/",
  "Data/Surveys/Apr 25/",
  "Data/Surveys/Nov 25/",
  "Data/Surveys/Feb 26/"
)

s1_files <- c(
  "Data/Population/Dolphins_covariates_2312.csv",
  "Data/Population/Dolphins_covariates_2504.csv",
  "Data/Population/Dolphins_covariates_2511.csv",
  "Data/Population/Dolphins_covariates_2602.csv"
)

final_df <- process_gpx_and_s1(folders, s1_files, goa_grid)

# Mean activity values per grid cell.
dat_map <- final_df |>
  group_by(grid_id) |>
  summarise(
    FB_avg = mean(FB, na.rm = TRUE),
    TB_avg = mean(TB, na.rm = TRUE),
    .groups = "drop"
  )

# Use ggarrange rather than patchwork because patchwork can fail with sf objects
# in this environment when stacking/composing multiple maps.
map_TB <- ggplot() +
  geom_sf(data = dat_map, aes(fill = TB_avg), alpha = 0.7) +
  geom_sf(data = goa_coast, fill = "NA", colour = "black") +
  scale_fill_gradientn(
    colours = c("lightsteelblue4", "dodgerblue", "darkblue")
  ) +
  coord_sf() +
  labs(fill = "Average no. of\nTourist Boats") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

map_FN <- ggplot() +
  geom_sf(data = dat_map, aes(fill = FB_avg), alpha = 0.7) +
  geom_sf(data = goa_coast, fill = "NA", colour = "black") +
  scale_fill_gradientn(
    colours = c("lightsteelblue4", "dodgerblue", "darkblue")
  ) +
  coord_sf() +
  labs(fill = "Average no. of\nFishing Gears") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggpubr::ggarrange(map_TB, map_FN, ncol = 2, common.legend = FALSE)
# Figure S6 in the manuscript
# ggsave("Results/Graphs/Human_activity_maps.jpg", width = 12, height = 6)
