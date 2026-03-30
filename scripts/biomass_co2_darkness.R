library(tidyverse)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(brms)

root_data <- read.csv(file = "data_raw/branch_root_data.csv") %>%
  mutate(treatment = factor(treatment,
                            levels = c("control", "heatwave", "extended"),
                            labels = c("Control", "Heat wave", "Extended season")))
    str(root_data)
    root_data$thin_white_branched <- as.numeric(root_data$thin_white_branched)
    root_data$thin_white_unbranched <- as.numeric(root_data$thin_white_unbranched)
    root_data$thick_white_unbranched <- as.numeric(root_data$thick_white_unbranched)
    root_data$thin_beige_branched <- as.numeric(root_data$thin_beige_branched)
    root_data$thin_beige_unbranched <- as.numeric(root_data$thin_beige_unbranched)
    root_data$thick_beige_unbranched <- as.numeric(root_data$thick_beige_unbranched)
    root_data$thin_brown_branched <- as.numeric(root_data$thin_brown_branched )
    root_data$thin_brown_unbranched <- as.numeric(root_data$thin_brown_unbranched )
    root_data$thick_brown_unbranched <- as.numeric(root_data$thick_brown_unbranched )
    root_data$thin_black_branched  <- as.numeric(root_data$thin_black_branched)
    root_data$thin_black_unbranched  <- as.numeric(root_data$thin_black_unbranched)
    root_data$thick_black_unbranched  <- as.numeric(root_data$thick_black_unbranched)
    
root_morphology <- root_data %>%
  pivot_longer(cols = shrub_biomass:graminoid_biomass,
               names_to = "functional_type_biomass",
               values_to = "pft_biomass") %>%
  filter(pft_biomass != "NA")


co2_flux <- read_csv(file = "data_raw/shifted_co2_data.csv") 


co2_flux_shifted <- co2_flux %>%
  mutate(measurement.week = ifelse(Treatment == "extended",
                                   measurement.week - 3,
                                   measurement.week),
         Treatment = factor(Treatment,
                            levels = c("Control", "Heatwave", "Extended"),
                            labels = c("Control", "Heat wave", "Extended season")))

# calculated co2 flux x darkness for each plant
avg_season_flux <- co2_flux_shifted %>%
  group_by(plant_id, Treatment) %>%
  summarise(total_co2_flux = sum(seasonal_co2_flux, na.rm = TRUE))

  co2_extended_darkness <- avg_season_flux %>%
    filter(plant_id %in% c(1, 8, 6))
  
  co2_extended_notdarkness <- avg_season_flux %>%
    filter(plant_id %in% c(7, 10, 16, 18))

  plant_avg_flux <- co2_flux_shifted %>%
    group_by(plant_id, Treatment) %>%
    summarise(total_co2_flux = sum(seasonal_co2_flux, na.rm = TRUE))

# root biomass x darkness for each plant
root_biomass_extended<- root_morphology %>%
  filter(plant_id %in% c(1, 8, 6, 7, 10, 16, 18)) %>%
  mutate(darkness = "Darkness") %>%
  group_by(plant_id, darkness, treatment) %>%
  summarize(total_root_biomass = median(total_root_biomass))

ratio_data_org_ctrl_hw_sing <- root_morph_separated %>%
  filter(plant_id %in% c(4, 5, 9, 14, 15, 17, 20,	2, 3, 11, 12, 13, 19, 21)) %>%
  mutate(darkness = "Darkness") %>%
  group_by(plant_id, darkness, treatment) %>%
  summarize(total_root_biomass = median(total_root_biomass))

combined_darkness_ratio_all <- bind_rows(ratio_data_org_ext_sing, ratio_data_org_ctrl_hw_sing)   
  