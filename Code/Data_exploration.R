library(tidyverse)

dat_str = read.csv("Data/2025/Strandings_flagged_geo_corrected.csv")

dat_str = dat_str |> filter(flag == FALSE | is.na(flag), Year >= 2017)
dat_str |> group_by(Region) |> summarise(n = n()) |> arrange(desc(n))

dat_str = dat_str |> filter(Region == "Goa")
dat_str |> group_by(Species_co) |> summarise(n = n()) |> arrange(desc(n))
#dat_str = dat_str |> filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise", "Small-sized delphinid", "Spinner dolphin")) |> 
  #write.csv("Data/2025/Shp files/Map_data/Strandings_filtered.csv", row.names = F)
dat_str = dat_str |> filter(Species_co %in% c("Humpback dolphin", "Indo-Pacific finless porpoise"))

#Total counts
dat_str |>
  ggplot() + 
  geom_bar(aes(x = as.factor(Month))) +
  facet_wrap(~Species_co, ncol=1)

#Average counts
dat_str |> group_by(Year, Month, Species_co) |> summarise(count = n()) |>
  group_by(Month, Species_co) |> summarise(n = mean(count), n_se = sd(count)/(length(count))^.5, .groups = 'drop') |>
  ggplot() + 
  geom_pointrange(aes(x = as.factor(Month), y = n, ymin = n - n_se, ymax = n + n_se, colour = Species_co), position = position_dodge(width = 0.5)) +
  #facet_wrap(~Species_co, ncol=1) + 
  theme_minimal() +
  scale_colour_manual(
    #values = c("Humpback dolphin" = "#7FCFE0", "Indo-Pacific finless porpoise" = "#B97A57")
    values = c("Humpback dolphin" = "lightsteelblue4", "Indo-Pacific finless porpoise" = "gray")
    ) +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),
        legend.text  = element_text(size = 12), legend.title = element_text(size = 13), legend.position = "bottom") +
  labs(x = "Month", y = "No. of strandings (mean ± SE)", colour = "Species")

ggsave("Results/2025/Graphs/stranding_counts_A1.jpg", width = 8, height = 5)

#counts by age class
dat_str |> filter(Age.Class != "") |> 
  group_by(Year, Month, Species_co, Age.Class) |> summarise(count = n()) |>
  group_by(Month, Species_co, Age.Class) |> summarise(n = mean(count), n_sd = sd(count)) |>
  ggplot() + 
  geom_pointrange(aes(x = as.factor(Month), y = n, ymin = n - n_sd, ymax = n + n_sd)) +
  facet_grid(cols = vars(Species_co), rows = vars(Age.Class))

#proportion by age class
t = dat_str |> filter(Age.Class != "") |> group_by(Month, Species_co, Age.Class) |>
  summarise(count = n())

  dat_str |> filter(Age.Class != "") |> group_by(Month, Species_co) |> summarise(count = n()) |>
  left_join(t, by = c("Month", "Species_co")) |> mutate(prop = count.y/count.x) |>
  group_by(Month, Species_co, Age.Class) |> summarise(prop = mean(prop, na.rm=TRUE), prop_sd = mean(prop, na.rm=TRUE), .groups = "drop") |>
  complete(Month, Species_co, Age.Class, fill = list(prop = NA, prop_sd = NA)) |> 
  ggplot() + 
  geom_bar(aes(x = as.factor(Month), y = prop, fill = Age.Class), stat="identity", position=position_dodge2(2)) +
    facet_grid(cols = vars(Species_co), rows = vars(Age.Class))

#carcass death code counts
#dat_str |> filter(Death_code != "") |> group_by(Species_co) |> count()
dat_str |> filter(Death_code != "") |> group_by(Species_co, Death_code) |> count() |> group_by(Species_co) |> mutate(percentage = n/sum(n)*100)

p1 = 
dat_str |> 
  filter(Death_code %in% c(1,2,3,4)) |> 
  mutate(across(c(PMI_min, PMI_max), as.numeric)) |> 
  mutate(Death_code = as.factor(Death_code)) |> 
  ggplot() + 
  geom_histogram(aes(x = Death_code, fill = Species_co), position = "dodge", stat = "count") + # stat = "count", bins = 10
  theme_minimal() +
  scale_fill_manual(
    #values = c("Humpback dolphin" = "#7FCFE0", "Indo-Pacific finless porpoise" = "#B97A57")
    values = c("Humpback dolphin" = "lightsteelblue4", "Indo-Pacific finless porpoise" = "gray")
    ) +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),
        legend.text  = element_text(size = 12), legend.title = element_text(size = 13)) +
  labs(x = "Decomposition code", y = "Count", fill = "Species")
  #labs(x = "post-morterm interval (lower bound)", y = "Count", fill = "Species")
  #labs(x = "post-morterm interval (upper bound)", y = "Count", fill = "Species")
ggpubr::ggarrange(p1,p2,p3, ncol = 1, common.legend = TRUE, legend = "bottom")
ggsave("Results/2025/Graphs/stranding_decompcode_pmi.jpg", width = 8, height = 10)
  
# Full dataset# Full dataset

dat_str = read.csv("Data/2025/Strandings_flagged_geo_corrected.csv")

dat_str = dat_str |> filter(flag == FALSE | is.na(flag), Year > 2000)
dat_str |> group_by(Region) |> summarise(n = n()) |> arrange(desc(n))
dat_str |> group_by(Species_co) |> summarise(n = n()) |> arrange(desc(n))
dat_str = dat_str |> filter(!(Region %in% c("Sri Lanka", "Pakistan", "Unknown"))) |> 
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

dat_str |>
  #filter(Species_co == "Indo-Pacific finless porpoise") |>
  filter(Region != "Goa") |>
  group_by(Year, species_type, coast) |> reframe(count = n()) |>
  ggplot() + 
  geom_point(aes(x = as.factor(Year), y = count, colour = species_type)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~coast)

#Estimaing congurence in expert vs experimental decemposition and PMI
decomp = read.csv("Data/2025/Decomposition_estimates.csv")
decomp = decomp |> rowwise() |> 
  mutate(m_dc = mean(c_across(starts_with("DC")), na.rm=T)) |> 
  mutate(m_pmil = mean(c_across(starts_with("PMIl")), na.rm=T)) |> 
  mutate(m_pmiu = mean(c_across(starts_with("PMIu")), na.rm=T)) |> 
  ungroup() |> 
  mutate(dc_diff = m_dc - Self_dc,
         pmiu_diff = m_pmiu - Self_PMIu,
         pmil_diff = m_pmil - Self_PMIl) |> 
  mutate(dc_diff_perc = dc_diff/Self_dc * 100,
         pmiu_diff_perc = pmiu_diff/Self_PMIu * 100,
         pmil_diff_perc = pmil_diff/Self_PMIl * 100)
mean(decomp$dc_diff_perc); sd(decomp$dc_diff_perc)#/sqrt(nrow(decomp))
mean(decomp$pmil_diff_perc); sd(decomp$pmil_diff_perc)#/sqrt(nrow(decomp))
mean(decomp$pmiu_diff_perc); sd(decomp$pmiu_diff_perc)#/sqrt(nrow(decomp))

decomp |> 
  rowwise() |> 
  mutate(m = mean(c_across(c(PMIl1, PMIl3, PMIl5, PMIl7, PMIl8)), na.rm = TRUE)) |> select(m)

#Removing the three reviewers with high differences
decomp = decomp |> rowwise() |> 
  mutate(m_dc = mean(c_across(c(DC1, DC3, DC5, DC7, DC8)), na.rm=T)) |> 
  mutate(m_pmil = mean(c_across(c(PMIl1, PMIl3, PMIl5, PMIl7, PMIl8)), na.rm=T)) |> 
  mutate(m_pmiu = mean(c_across(c(PMIu1, PMIu3, PMIu5, PMIu7, PMIu8)), na.rm=T)) |> 
  ungroup() |> 
  mutate(dc_diff = m_dc - Self_dc,
         pmiu_diff = m_pmiu - Self_PMIu,
         pmil_diff = m_pmil - Self_PMIl) |> 
  mutate(dc_diff_perc = dc_diff/Self_dc * 100,
         pmiu_diff_perc = pmiu_diff/Self_PMIu * 100,
         pmil_diff_perc = pmil_diff/Self_PMIl * 100)
mean(decomp$dc_diff_perc); sd(decomp$dc_diff_perc)#/sqrt(nrow(decomp))
mean(decomp$pmil_diff_perc); sd(decomp$pmil_diff_perc)#/sqrt(nrow(decomp))
mean(decomp$pmiu_diff_perc); sd(decomp$pmiu_diff_perc)#/sqrt(nrow(decomp))

#Mapping Fisheries and Tourism from boat surveys
library(sf)
library(purrr)
library(patchwork)
library(tidyverse)

goa_grid = st_read("GIS/Goa_buffer_5km_grid.shp")
goa_coast = st_read("GIS/Goa_coast.shp")

process_gpx_and_s1 = function(gpx_folders, s1_files, goa_grid) {
  results = map2_dfr(gpx_folders, s1_files, function(folder, s1_file) {
    s1 = read.csv(s1_file)
    waypoints_all = list.files(folder, pattern = "Waypoints_.*\\.gpx$", full.names = TRUE) |>
      map_dfr(~ st_read(.x, layer = "waypoints", quiet = TRUE) |>
                mutate(source_file = basename(.x))) |>
      mutate(name = as.numeric(name)) |>
      select(name, geometry) |>
      right_join(s1, by = c("name" = "Waypoint.number"))
    waypoints_all$grid_id <- st_intersects(waypoints_all, goa_grid) %>% as.integer()
    
    # mutate(grid_id = st_intersects(., goa_grid))
    
    waypoints_all |>
      group_by(grid_id) |>
      summarise(
        FBs = sum(Number.of.fishing.boats..small.medium, na.rm = TRUE),
        FBl = sum(Number.of.fishing.boats..large, na.rm = TRUE),
        TB  = sum(Number.of.tourist.boats, na.rm = TRUE)
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

# Call the function
folders = c(
  "C:/Users/imran/OneDrive - Indian Institute of Science/IISc files/IISC/Thesis/Goa_population_estimation/Imran_surveys/GPS/",
  "C:/Users/imran/OneDrive - Indian Institute of Science/IISc files/IISC/Thesis/Goa_population_estimation/Imran_surveys/GPS/Apr 25/",
  "C:/Users/imran/OneDrive - Indian Institute of Science/IISc files/IISC/Thesis/Goa_population_estimation/Imran_surveys/GPS/Nov 25/",
  "C:/Users/imran/OneDrive - Indian Institute of Science/IISc files/IISC/Thesis/Goa_population_estimation/Imran_surveys/GPS/Feb 26/"
)

s1_files <- c("Data/2025/Population/Dolphins_covariates_2312.csv", "Data/2025/Population/Dolphins_covariates_2504.csv", "Data/2025/Population/Dolphins_covariates_2511.csv","Data/2025/Population/Dolphins_covariates_2602.csv")

final_df <- process_gpx_and_s1(folders, s1_files, goa_grid)

dat_map = final_df |> group_by(grid_id) |> summarise(FB_avg = mean(FB, na.rm=T), TB_avg = mean(TB, na.rm=T))

map_TB = 
ggplot() + 
  geom_sf(data = dat_map, aes(fill = TB_avg), alpha = 0.7) +
  geom_sf(data = goa_coast, fill = "NA", colour = "black") +
  scale_fill_gradientn(colours = c("lightsteelblue4", "dodgerblue", "darkblue")) +
  coord_sf() +
  labs(fill = "Average no. of\nTourist Boats") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

map_FN = 
ggplot() + 
  geom_sf(data = dat_map, aes(fill = FB_avg), alpha = 0.7) +
  geom_sf(data = goa_coast, fill = "NA", colour = "black") +
  scale_fill_gradientn(colours = c("lightsteelblue4", "dodgerblue", "darkblue")) +
  coord_sf() +
  labs(fill = "Average no. of\nFishing Gears") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

map_TB + map_FN + plot_layout(ncol = 2)
ggsave("Results/2025/Graphs/Human_activity_maps.jpg", width = 12, height = 6)
