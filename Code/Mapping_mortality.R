library(tidyverse)
library(sf)
library(parallel)

#Function to check if points are inside a polygon
is_in_polygon = function(ipl,grid)
{
  library(sf)
  pnts = data.frame(
    "x" = ipl$lon,
    "y" = ipl$lat)
  
  pnts_sf <- sf::st_as_sf(pnts, coords = c("x", "y"), crs = st_crs(grid))
  xx = st_intersects(pnts_sf, grid)
  yy = lengths(xx) == 1
  rec_id = sapply(xx, function(x) if (length(x) == 1) x else NA)
  ret = data.frame(yy,rec_id)
  return(ret)
}

df_rev = read.csv("Results/2025/Goa_rev.csv")
df_rev = df_rev |> arrange(id,desc(Year),desc(Month),desc(Day)) |> filter(flag == 'FALSE' | flag == "") |> 
  filter(Species_co == "Humpback dolphin") |> 
  #filter(Species_co == "Indo-Pacific finless porpoise") |> 
  filter(!is.na(lat))
grid_full = read_sf("Data/2025/Shp files/Grid_ecol_refined.shp", crs=4326)
goa_coast = read_sf("Data/2025/Shp files/Goa_coast.shp")

points = st_as_sf(x=data.frame(x=df_rev$lon,y=df_rev$lat),coords = c("x", "y"),crs=4326)
df_rev$dist = as.numeric(st_distance(points, goa_coast))

ggplot() +
  geom_sf(data = grid_full, mapping = aes()) +
  geom_point(df_rev, mapping = aes(x=lon, y=lat, color=as.factor(Month)), size=2, alpha=0.7) +
  geom_sf(data = goa_coast, aes(), alpha = 0.3)

stm = Sys.time()
traj_list=list()
for(i in 1:nrow(df_rev))
  traj_list[[i]]=df_rev[i,c("lon","lat")]

cl <- makeCluster(12)
clusterExport(cl, c("is_in_polygon", "grid_full"))
op_pnts = parLapply(cl, traj_list, function(df_rev) is_in_polygon(df_rev, grid_full))
stopCluster(cl)
op_pnts = do.call(rbind, op_pnts)
Sys.time() - stm

df_rev$strand_id = grid_full$id[op_pnts$rec_id]

#important filters
df_rev_min = df_rev |> group_by(id) |> filter(adv_time_days<=PMI_max) |> 
  #filter(PMI_min <= 6) |>
  filter(dist < 5000)
df_rev_min = df_rev_min |> group_by(id) |> count() |> left_join(df_rev_min, by = "id") |> group_by(id) |> mutate(MI = 1/n)

plot_data = 
  grid_full |>
  left_join(df_rev_min |> group_by(strand_id) |> summarise(mortality = sum(MI), .groups = "drop"),
            by = c("id" = "strand_id")) |>
  filter(!is.na(mortality))
plot_data |> 
  ggplot() +
  geom_sf(aes(fill = (mortality))) +
  scale_fill_viridis_c(option = "plasma", na.value = "lightgrey") +
  geom_sf(data = goa_coast, aes(), alpha = 0.3) +
  geom_point(df_rev_min[df_rev_min$adv_time_days == 0,],
             mapping = aes(x=lon, y=lat),
             color = "red", size=1, alpha=0.7) +
  theme_minimal() +
  labs(x = "lon", y = "lat") +
  ggtitle("log(Mortality index)")
  #facet_wrap(~Month)

write.csv(plot_data, "Results/2025/Goa_rev_IOHD_maxPMI_5km.csv", row.names = F)
