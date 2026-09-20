# Load required libraries for data manipulation and spatial analysis

library(sf)
library(parallel)

# Read trajectory data from CSV files

traj2 <- read.csv("Results/Goa_fwd_op_p2.csv")
# Combine trajectory data from both files

rm(traj1, traj2)
# Define coastline cells
coast_cells <- c(
  290,
  338,
  385,
  386,
  433,
  434,
  482,
  435, #Maharashtra
  483,
  484,
  531,
  532,
  533,
  581,
  582,
  630,
  631,
  678,
  632,
  679,
  727,
  728,
  776,
  777,
  778,
  779,
  780,
  827,
  828, #Goa
  875,
  876,
  877,
  924,
  534,
  726,
  729, #Goa
  925,
  972,
  973,
  974,
  927,
  975,
  1022,
  1069,
  1116 #Karnataka
) #coastline cells
# Read spatial grid data from shapefile

grid_dat <- read.csv("Data/Shp files/Grid_ecol_refined.csv")
grid <- grid |> left_join(grid_dat[, c("id", "distance")], by = "id")

grid_full <- grid
# Filter grid to include only coastline cells

# Function to check if points are inside a polygon

is_in_polygon <- function(ipl, grid) {
  library(sf)
  pnts <- data.frame(
    "x" = ipl$feature_x,
    "y" = ipl$feature_y
  )

  pnts_sf <- sf::st_as_sf(pnts, coords = c("x", "y"), crs = st_crs(grid))
  xx <- st_intersects(pnts_sf, grid)
  yy <- lengths(xx) == 1
  rec_id <- sapply(xx, function(x) if (length(x) == 1) x else NA)
  ret <- data.frame(yy, rec_id)
  return(ret)
}

# Subset trajectory data to include only points close to the coast

traj <- traj[
  which(
    (traj$feature_x > 73.650) &
      (traj$feature_y < 15.76) &
      (traj$feature_y > 14.89)
  ),
]

stm <- Sys.time()
traj_list <- list()
for (i in 1:nrow(traj)) {
  traj_list[[i]] <- traj[i, c("feature_x", "feature_y")]
}

# Create a cluster for parallel processing

clusterExport(cl, c("is_in_polygon", "grid"))
op_pnts <- parLapply(cl, traj_list, function(traj) is_in_polygon(traj, grid))
stopCluster(cl)
op_pnts <- do.call(rbind, op_pnts)
Sys.time() - stm

#save.image("Results/Rdata/estimating_stranding_prob.RData")

traj$strand <- op_pnts$yy
traj$strand_id <- op_pnts$rec_id

part_strand_prob <- traj |>
  group_by(id, strand_id, time_slot, month, year) |>
  summarise(count = sum(strand), .groups = "drop") |>
  mutate(count = ifelse(count > 0, 1, 0))

#stranding prob from each grid for every month
part_strand_prob_from_grid <- part_strand_prob |>
  group_by(id, month, year) |>
  summarise(prob = sum(count), .groups = "drop") |>
  mutate(prob = ifelse(prob > 3, 3, prob), prob = prob / 3) |>
  left_join(grid_full[, c("id", "distance")], by = c("id" = "id")) |>
  group_by(id, distance) |>
  summarise(prob = mean(prob), .groups = "drop")
# Save results to CSV files

write.csv("Results/str_prob_from_grids.csv", row.names = F)

# Plotting stranding probability from each grid cell

left_join(part_strand_prob_from_grid, by = c("id" = "id")) |>
  ggplot() +
  geom_sf(aes(fill = prob)) +
  scale_fill_viridis_c(option = "plasma", na.value = "lightgrey") +
  theme_minimal() +
  labs(fill = "Stranding probability") +
  ggtitle("Stranding probability from each grid cell")

# ggsave("Results/Graphs/str_prob_from_grids.jpg", width = 8, height = 6)

#For all grids
part_strand_prob |>
  group_by(id, month, year) |>
  summarise(prob = sum(count), .groups = "drop") |>
  mutate(prob = ifelse(prob > 3, 3, prob), prob = prob / 3) |>
  group_by(month, year) |>
  summarise(mean_prob = mean(prob), se = sd(prob) / (length(prob))^.5) |>
  write.csv("Results/str_prob_month_year.csv", row.names = F)

#For grids close to the coast (<5 km); these species only occur within this range
part_strand_prob_from_grid <- part_strand_prob |>
  group_by(id, month, year) |>
  summarise(prob = sum(count), .groups = "drop") |>
  mutate(prob = ifelse(prob > 3, 3, prob), prob = prob / 3) |>
  left_join(grid_full[, c("id", "distance")], by = c("id" = "id")) |>
  filter(distance < 10000)

# Save results to CSV files

group_by(month, year) |>
  summarise(mean_prob = mean(prob), se = sd(prob) / (length(prob))^.5) |>
  write.csv("Results/str_prob_month_year_g10km.csv", row.names = F)

# Calculate average stranding probability across months

#part_strand_prob_from_grid = read.csv("Results/str_prob_month_year_g5km.csv")
part_strand_prob_from_grid_months <- # Save results to CSV files
  select(-geometry) |>
  group_by(month) |>
  summarise(
    mean_prob = mean(prob),
    se = sd(prob) / (length(prob))^.5,
    ci5 = quantile(prob, 0.05),
    ,
    ci95 = quantile(prob, 0.95)
  )
fig_strp <-
  part_strand_prob_from_grid_months |>
  ggplot(aes(x = as.factor(month), y = mean_prob)) +
  geom_bar(stat = "identity", fill = "gray") + #, fill = "#9AD9EA", alpha = 1) +
  geom_errorbar(
    aes(ymin = mean_prob - 1.95 * se, ymax = mean_prob + 1.95 * se),
    width = 0.2
  ) +
  #geom_errorbar(aes(ymin = ci5, ymax = ci95), width = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    #axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14)
    #axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))
  ) +
  #labs(y = "Stranding probability\n(mean ± SE)") +
  labs(y = "Stranding probability\n(mean ± 95% CI)") +
  annotate(
    "text",
    x = 1,
    y = Inf,
    label = "a)",
    hjust = 1.1,
    vjust = 1.5,
    size = 6
  )

#ggsave("Results/Graphs/str_prob_months_g5km.jpg", width = 8, height = 5)

ggpubr::ggarrange(fig_strp, fig_pop, ncol = 1) #load fig_pop from the Density_Distance.R file
#ggsave("Results/Graphs/str_prob_pop_distance_CI.jpg", width = 6, height = 8)

#write.csv(part_strand_prob_from, "Results/str_prob_from.csv", row.names = F)
#save.image("Results/Rdata/estimating_stranding_prob.RData")
