library(tidyverse)
library(Distance)

pop = data.frame()

dol_dat = read.csv("Data/2025/Population/Dolphins_locations_2312.csv")
dol_dat$PDistance_m = dol_dat$Distance_m*sin((abs(dol_dat$Angle)*(pi/180)))
dol_dat$Start = as.Date(dol_dat$Start, "%d-%m-%Y")
dol_dat$Sample.Label <- paste0(dol_dat$Transect.ID, "_R", dol_dat$rep)

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

#Formatting data
dist_data = data.frame(cbind("Region.Label" = "Default", 
                             #"Area" = max(dol_dat$PDistance_m)/1000*2*346, 
                             #"Sample.Label" = dol_dat$rep, 
                             "Sample.Label" = dol_dat$Sample.Label,
                             #"Effort" = 346/2,
                             distance = dol_dat$PDistance_m, 
                             size = dol_dat$N_best,
                             "Boat_speed" = dol_dat$Boat.speed..kmph.,
                             "FB" = dol_dat$Number.of.fishing.boats..small.medium + dol_dat$Number.of.fishing.boats..large,
                             "TB" = dol_dat$Number.of.tourist.boats))#41 is the transect end, effort is per replicate
dist_data[,3:ncol(dist_data)] = sapply(dist_data[,3:ncol(dist_data)], as.numeric)
dist_data[,5:ncol(dist_data)] = sapply(dist_data[,5:ncol(dist_data)], scale)
hist(dist_data$distance, breaks=20)
truncation = 600

# Define region (total area surveyed, km²)
region.table <- data.frame(
  Region.Label = "Default",
  Area = 347/2*truncation/1000*2#346
)

# Define effort per sample (replicate)
sample.table <- data.frame(
  Sample.Label = unique(dol_dat$Sample.Label),
  #Region.Label = sub("_.*", "", unique(dol_dat$Sample.Label)),
  Region.Label = "Default",
  Effort = c(26.3, 20.2, 42.1, 59.2, 55.3-8.6, 42.4+8.6, 63.8, 37.3) #estimated from qgis tracklengths
)

mean(dol_dat$N_best)
conversion <- convert_units("meter", "kilometer", "square kilometer")
mod2 <- ds(dist_data, key="hn",
           adjustment="poly",
           convert_units = conversion,
           truncation = truncation,
           region_table=region.table,
           sample_table=sample.table)
           #formula = ~ Boat_speed + FB + TB)
summary(mod2)
plot(mod2, nc=12, main = "December 2023")
gof_ds(mod2)

pop = rbind(pop, data.frame("Sl" = 1,
                            "Year" = 2023,
                           "Month" = 12,
                           "Season" = "Winter",
                           "Density" = 1.637,#1.7387,
                           "Density_se" = 0.35,#0.2593,
                           "Abundance" = 1.637 * 2 * 130, #coastline is about 130 km including some river edge; dolphins are found up to 2 km from the shore
                           "Abundance_se" = 0.35 * 2 * 130,
                           "lcl" = 1.042 * 2 * 130,
                           "ucl" = 2.571 * 2 * 130
                           ))

#2025 April
dol_dat = read.csv("Data/2025/Population/Dolphins_locations_2504.csv")
dol_dat$PDistance_m = dol_dat$Distance_m*sin((abs(dol_dat$Angle)*(pi/180)))
dol_dat$Start = as.Date(dol_dat$Start, "%d-%m-%Y")
dol_dat$Sample.Label <- paste0(dol_dat$Transect.ID, "_R", dol_dat$rep)

dol_cov = read.csv("Data/2025/Population/Dolphins_covariates_2504.csv")
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

#Formatting data
dist_data = data.frame(cbind("Region.Label" = "Default", 
                             #"Area" = max(dol_dat$PDistance_m)/1000*2*346, 
                             #"Sample.Label" = dol_dat$rep, 
                             "Sample.Label" = dol_dat$Sample.Label,
                             #"Effort" = 346/2,
                             distance = dol_dat$PDistance_m, 
                             size = dol_dat$N_best,
                             "Boat_speed" = dol_dat$Boat.speed..kmph.,
                             "FB" = dol_dat$Number.of.fishing.boats..small.medium + dol_dat$Number.of.fishing.boats..large,
                             "TB" = dol_dat$Number.of.tourist.boats))#41 is the transect end, effort is per replicate
dist_data[,3:ncol(dist_data)] = sapply(dist_data[,3:ncol(dist_data)], as.numeric)
dist_data[,5:ncol(dist_data)] = sapply(dist_data[,5:ncol(dist_data)], scale)
hist(dist_data$distance, breaks=20)
truncation = 500

# Define region (total area surveyed, km²)
region.table <- data.frame(
  Region.Label = "Default",
  Area = 310/2*truncation/1000*2#346
)

# Define effort per sample (replicate)
sample.table <- data.frame(
  Sample.Label = unique(dol_dat$Sample.Label),
  #Region.Label = sub("_.*", "", unique(dol_dat$Sample.Label)),
  Region.Label = "Default",
  Effort = c(28.6, 57.6, 34.4+10.7, 38.5+10.7, 61.1, 23.5, 16.8, 17.1) #estimated from qgis tracklengths
)

mean(dol_dat$N_best)
conversion <- convert_units("meter", "kilometer", "square kilometer")
mod1_2504 <- ds(dist_data, key="hn",
           adjustment="poly",
           convert_units = conversion,
           truncation = truncation,
           region_table=region.table,
           sample_table=sample.table)
           #formula = ~ Boat_speed + TB + FB)
summary(mod1_2504)
plot(mod1_2504, nc=12, main = "April 2025")
gof_ds(mod1_2504)

pop = rbind(pop, data.frame("Sl" = 2,
                            "Year" = 2025,
                            "Month" = 4,
                            "Season" = "Summer",
                            "Density" = 2.0333,#2.5228,
                            "Density_se" = 0.6905,#1.015,
                            "Abundance" = 2.0333 * 2 * 130,
                            "Abundance_se" = 0.6905 * 2 * 130,
                            "lcl" = 0.953 * 2 * 130,
                            "ucl" = 4.337 * 2 * 130
))

#2001 February
dol_dat = read.csv("Data/2025/Population/Dolphins_locations_0202.csv")
dol_dat$PDistance_m = dol_dat$Distance_m*sin((abs(dol_dat$Angle)*(pi/180)))
dol_dat$Start = as.Date(dol_dat$Start, "%d-%m-%Y")
dol_dat$Sample.Label <- paste0(dol_dat$Transect.ID, "_R", dol_dat$rep)

#Formatting data
dist_data = data.frame(cbind("Region.Label" = "Default",
                             #"Region.Label" = dol_dat$Transect.ID,
                             #"Area" = max(dol_dat$PDistance_m)/1000*2*346, 
                             "Sample.Label" = dol_dat$Sample.Label, #dol_dat$rep, 
                             #"Effort" = 346/2,
                             distance = abs(dol_dat$PDistance_m), 
                             size = dol_dat$N_best
                             ))
dist_data[,3:ncol(dist_data)] = sapply(dist_data[,3:ncol(dist_data)], as.numeric)
hist(dist_data$distance, breaks=20)
truncation = 400


# Define region (total area surveyed, km²)
region.table <- data.frame(
  #Region.Label = unique(dist_data$Region.Label),
  Region.Label = "Default",
  #Area = 15,
  Area = 15*20*truncation/1000*2#678
)

# # Define effort per sample (replicate)
# sample.table <- data.frame(
#   Sample.Label = c(1, 2, 3, 4, 5),
#   Region.Label = rep("Default", 5),
#   Effort = rep(678/5,5)#c(28.8,43.2,115.2,86.4,14.4) # see text file in the data folder
# )

sample.table <- data.frame(
  Sample.Label = unique(dol_dat$Sample.Label),
  #Region.Label = sub("_.*", "", unique(dol_dat$Sample.Label)),
  Region.Label = "Default",
  Effort = rep(15, length(unique(dol_dat$Sample.Label))) # each replicate is 15 km
)


mean(dol_dat$N_best)
conversion <- convert_units("meter", "kilometer", "square kilometer")
mod1_0202 <- ds(dist_data, key="hn",
                adjustment="cos",
                convert_units = conversion,
                truncation = truncation,
                region_table=region.table,
                sample_table=sample.table)
summary(mod1_0202)
plot(mod1_0202, nc=12, main = "February 2002")
gof_ds(mod1_0202)

pop = rbind(pop, data.frame("Sl" = 5,
                            "Year" = 2002,
                            "Month" = 2,
                            "Season" = "Spring",
                            "Density" = 3.646,
                            "Density_se" = 0.5939,
                            "Abundance" = 3.646 * 2 * 130,
                            "Abundance_se" = 0.5939 * 2 * 130,
                            "lcl" = 2.635 * 2 * 130,
                            "ucl" = 5.045 * 2 * 130
))

##############################
#Nov 25

dol_dat = read.csv("Data/2025/Population/Dolphins_locations_2511.csv")
dol_dat$PDistance_m = dol_dat$Distance_m*sin((abs(dol_dat$Angle)*(pi/180)))
dol_dat$Start = as.Date(dol_dat$Start, "%d-%m-%Y")
dol_dat$Sample.Label <- paste0(dol_dat$Transect.ID, "_R", dol_dat$rep)

dol_cov = read.csv("Data/2025/Population/Dolphins_covariates_2511.csv")
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

#Formatting data
dist_data = data.frame(cbind("Region.Label" = "Default", 
                             #"Area" = max(dol_dat$PDistance_m)/1000*2*346, 
                             #"Sample.Label" = dol_dat$rep, 
                             "Sample.Label" = dol_dat$Sample.Label,
                             #"Effort" = 346/2,
                             distance = dol_dat$PDistance_m, 
                             size = dol_dat$N_best,
                             "Boat_speed" = dol_dat$Boat.speed..kmph.,
                             "FB" = dol_dat$Number.of.fishing.boats..small.medium + dol_dat$Number.of.fishing.boats..large,
                             "TB" = dol_dat$Number.of.tourist.boats))#41 is the transect end, effort is per replicate
dist_data[,3:ncol(dist_data)] = sapply(dist_data[,3:ncol(dist_data)], as.numeric)
dist_data[,5:ncol(dist_data)] = sapply(dist_data[,5:ncol(dist_data)], scale)
hist(dist_data$distance, breaks=20)
truncation = 500

# Define region (total area surveyed, km²)
region.table <- data.frame(
  Region.Label = "Default",
  Area = 328/2*truncation/1000*2#346
)

# Define effort per sample (replicate)
sample.table <- data.frame(
  Sample.Label = unique(dol_dat$Sample.Label),
  #Region.Label = sub("_.*", "", unique(dol_dat$Sample.Label)),
  Region.Label = "Default",
  Effort = c(26, 56.8, 54.6-10.2, 34.7+10.2, 61.9, 22.7, 14.8, 14.9) #estimated from qgis tracklengths
)

# sample.table <- data.frame(
#   Sample.Label = c(1, 2),
#   Region.Label = "Default",
#   Effort = c(164, 164)
# )

mean(dol_dat$N_best)
conversion <- convert_units("meter", "kilometer", "square kilometer")
mod1_2511 <- ds(dist_data, key="hn",
           adjustment="poly",
           convert_units = conversion,
           truncation = truncation,
           region_table=region.table,
           sample_table=sample.table)
#formula = ~ Boat_speed + FB + TB)
summary(mod1_2511)
plot(mod1_2511, nc=12, main = "November 2025")
gof_ds(mod1_2511)

pop = rbind(pop, data.frame("Sl" = 3,
                            "Year" = 2025,
                            "Month" = 11,
                            "Season" = "Winter",
                            "Density" = 2.0274,
                            "Density_se" = 0.5158,
                            "Abundance" = 2.0274 * 2 * 130, #coastline is about 130 km including some river edge; dolphins are found up to 2 km from the shore
                            "Abundance_se" = 0.5158 * 2 * 130,
                            "lcl" = 1.154 * 2 * 130,
                            "ucl" = 3.561 * 2 * 130
))

##############################
#Feb 26

dol_dat = read.csv("Data/2025/Population/Dolphins_locations_2602.csv")
dol_dat$PDistance_m = dol_dat$Distance_m*sin((abs(dol_dat$Angle)*(pi/180)))
dol_dat$Start = as.Date(dol_dat$Start, "%d-%m-%Y")
dol_dat$Sample.Label <- paste0(dol_dat$Transect.ID, "_R", dol_dat$rep)

dol_cov = read.csv("Data/2025/Population/Dolphins_covariates_2602.csv")
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

#Formatting data
dist_data = data.frame(cbind("Region.Label" = "Default", 
                             #"Area" = max(dol_dat$PDistance_m)/1000*2*346, 
                             #"Sample.Label" = dol_dat$rep, 
                             "Sample.Label" = dol_dat$Sample.Label,
                             #"Effort" = 346/2,
                             distance = dol_dat$PDistance_m, 
                             size = dol_dat$N_best,
                             "Boat_speed" = dol_dat$Boat.speed..kmph.,
                             "FB" = dol_dat$Number.of.fishing.boats..small.medium + dol_dat$Number.of.fishing.boats..large,
                             "TB" = dol_dat$Number.of.tourist.boats))
dist_data[,3:ncol(dist_data)] = sapply(dist_data[,3:ncol(dist_data)], as.numeric)
dist_data[,5:ncol(dist_data)] = sapply(dist_data[,5:ncol(dist_data)], scale)
hist(dist_data$distance, breaks=20)
truncation = 500

# Define region (total area surveyed, km²)
region.table <- data.frame(
  Region.Label = "Default",
  Area = 265/2*truncation/1000*2#this comes from the average transect length for a replicate (of the two)
)

# Define effort per sample (replicate)
sample.table <- data.frame(
  Sample.Label = unique(dol_dat$Sample.Label),
  #Region.Label = sub("_.*", "", unique(dol_dat$Sample.Label)),
  Region.Label = "Default",
  Effort = c(18.5,17.1,23.6,55.8,38,33.7+9,44.4,18.3) #estimated from qgis tracklengths
)

# sample.table <- data.frame(
#   Sample.Label = c(1, 2),
#   Region.Label = "Default",
#   Effort = c(164, 164)
# )

mean(dol_dat$N_best)
conversion <- convert_units("meter", "kilometer", "square kilometer")
mod1_2602 <- ds(dist_data, key="hn",
                adjustment="poly",
                convert_units = conversion,
                truncation = truncation,
                region_table=region.table,
                sample_table=sample.table)
#formula = ~ Boat_speed + FB + TB)
summary(mod1_2602)
plot(mod1_2602, nc=12, main = "February 2026")
gof_ds(mod1_2602)

pop = rbind(pop, data.frame("Sl" = 4,
                            "Year" = 2026,
                            "Month" = 02,
                            "Season" = "Winter",
                            "Density" = 2.55628,
                            "Density_se" = 0.4340,
                            "Abundance" = 2.55628 * 2 * 130, #coastline is about 130 km including some river edge; dolphins are found up to 2 km from the shore
                            "Abundance_se" = 0.4340 * 2 * 130,
                            "lcl" = 1.751 * 2 * 130,
                            "ucl" = 3.731 * 2 * 130
))

#################


fig_pop = 
pop |> 
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Winter"))) |>
  filter(Year>2020) |> 
  complete(Month = 1:12) |> 
  ggplot(aes(x = as.factor(Month), y = Abundance, fill = as.factor(Season))) +
  geom_bar(stat = "identity", fill = "gray", alpha = 1) +#, width = 0.5) +
  #geom_errorbar(aes(ymin = Abundance - Abundance_se, ymax = Abundance + Abundance_se), width = 0.2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14), axis.title.y=element_text(size=14)) +
  #labs(x = "Month", y = "No. of individuals\n(mean ± SE)") +
  labs(x = "Month", y = "No. of individuals\n(mean ± 95% CI)") +
  annotate("text", x = 1, y = Inf, label = "b)", hjust = 1.1, vjust = 1.5, size = 6)
  # scale_x_discrete(labels = c("1" = "December\n2023",
  #                             "2" = "April\n2025",
  #                             "3" = "November\n2025",
  #                             "4" = "February\n2002"))

#ggsave("Results/2025/Graphs/density_distance_final.jpg", width = 6, height = 5)
jpeg("Results/2025/Graphs/density_detection_functions.jpg", width = 1200, height = 1200, units = "px")
