# 5_Mapping_mortality.R
# Spatial mapping of latent cetacean mortality risk from stranding records.
#
# The analysis assigns each stranding record to a grid cell, filters records that
# are biologically plausible for the mortality calculation, and then aggregates a
# mortality index across grid cells. The result is a spatial surface showing
# areas of elevated inferred mortality risk.

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

# ----------------------------------------------------------------------------
# 1. Read input data and filter to the focal species
# ----------------------------------------------------------------------------
strand_records <- read.csv("Results/Goa_rev.csv") |>
  arrange(id, desc(Year), desc(Month), desc(Day)) |>
  filter(flag == "FALSE" | flag == "" | is.na(flag)) |>
  filter(Species_co == "Humpback dolphin") |>
  filter(!is.na(lat))

# Spatial layers used for gridding and coastline distance calculations.
goa_grid <- st_read("Data/Shp files/Grid_ecol_refined.shp", quiet = TRUE)
goa_coast <- st_read("Data/Shp files/Goa_coast.shp", quiet = TRUE)

# ----------------------------------------------------------------------------
# 2. Helper functions
# ----------------------------------------------------------------------------
# Function to assign a point to the correct grid cell.
is_in_polygon <- function(point_df, grid_sf) {
  point_sf <- st_as_sf(
    point_df,
    coords = c("lon", "lat"),
    crs = st_crs(grid_sf)
  )

  cell_intersections <- st_intersects(point_sf, grid_sf)
  cell_match <- sapply(cell_intersections, function(x) {
    if (length(x) == 1) x[1] else NA_integer_
  })

  data.frame(
    yy = lengths(cell_intersections) == 1,
    rec_id = cell_match,
    stringsAsFactors = FALSE
  )
}

# This function assigns each record to a grid cell and stores the cell ID.
assign_grid_id <- function(dat, grid_sf) {
  point_assignments <- map_dfr(seq_len(nrow(dat)), function(i) {
    pt <- dat[i, c("lon", "lat")]
    is_in_polygon(pt, grid_sf)
  })

  dat$strand_id <- grid_sf$id[point_assignments$rec_id]
  dat
}

# ----------------------------------------------------------------------------
# 3. Spatial preparation and distance-to-coast check
# ----------------------------------------------------------------------------
# Convert records to sf points to calculate nearest distance to the coast.
point_sf <- st_as_sf(
  data.frame(
    x = strand_records$lon,
    y = strand_records$lat
  ),
  coords = c("x", "y"),
  crs = 4326
)

strand_records <- strand_records |>
  mutate(
    dist_to_coast = as.numeric(st_distance(point_sf, goa_coast))
  )

# ----------------------------------------------------------------------------
# 4. Plot the raw sampling points over the grid and coast
# ----------------------------------------------------------------------------
raw_point_map <- ggplot() +
  geom_sf(data = goa_grid, fill = "white") +
  geom_point(
    data = strand_records,
    mapping = aes(x = lon, y = lat, colour = as.factor(Month)),
    size = 2,
    alpha = 0.7
  ) +
  geom_sf(data = goa_coast, alpha = 0.3) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Month"
  )

raw_point_map

# ----------------------------------------------------------------------------
# 5. Assign each record to a grid cell
# ----------------------------------------------------------------------------
strand_records <- assign_grid_id(strand_records, goa_grid)

# ----------------------------------------------------------------------------
# 6. Apply biological and spatial filters for mortality index calculation
# ----------------------------------------------------------------------------
# The mortality index is calculated only for records that remain plausible under
# the carcass-timing rules and that are sufficiently close to shore.
filtered_records <- strand_records |>
  group_by(id) |>
  filter(adv_time_days <= PMI_max) |>
  filter(dist_to_coast < 5000)

# For each ID, calculate the number of records retained and convert this to a
# mortality contribution of 1/n, matching the original approach.
filtered_records <- filtered_records |>
  group_by(id) |>
  summarise(n = n(), .groups = "drop") |>
  left_join(filtered_records, by = "id") |>
  group_by(id) |>
  mutate(MI = 1 / n)

# ----------------------------------------------------------------------------
# 7. Aggregate mortality into the grid
# ----------------------------------------------------------------------------
# Summed mortality contribution for each grid cell.
cell_mortality <- filtered_records |>
  group_by(strand_id) |>
  summarise(mortality = sum(MI, na.rm = TRUE), .groups = "drop")

mortality_map_data <- goa_grid |>
  left_join(cell_mortality, by = c("id" = "strand_id")) |>
  filter(!is.na(mortality))

# ----------------------------------------------------------------------------
# 8. Create the spatial mortality map
# ----------------------------------------------------------------------------
mortality_map <- ggplot() +
  geom_sf(data = mortality_map_data, aes(fill = mortality)) +
  scale_fill_viridis_c(option = "plasma", na.value = "lightgrey") +
  geom_sf(data = goa_coast, alpha = 0.3) +
  geom_point(
    data = filtered_records |> filter(adv_time_days == 0),
    mapping = aes(x = lon, y = lat),
    colour = "red",
    size = 1,
    alpha = 0.7
  ) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Mortality index"
  ) +
  ggtitle("Spatial pattern of inferred mortality risk")

mortality_map

# ----------------------------------------------------------------------------
# 9. Save outputs
# ----------------------------------------------------------------------------
write.csv(
  mortality_map_data,
  "Results/Goa_rev_IOHD_maxPMI_5km.csv",
  row.names = FALSE
)
