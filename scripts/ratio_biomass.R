library(tidyverse)
library(ggplot2)
library(dunn.test)
library(dplyr)

root_data <- read.csv(file = "data_raw/branch_root_data.csv") %>%
  mutate(treatment = factor(treatment,
                            levels = c("control", "heatwave", "extended"),
                            labels = c("Control", "Heat wave", "Extended season")))

    str(root_data)
    root_data$thin_white_branched <- as.numeric(root_data$thin_white_branched)
    root_data$thin_white_unbranched <- as.numeric(root_data$thin_white_unbranched)
    root_data$thick_white_unbranched <- as.numeric(root_data$thin_white_unbranched)
    root_data$thin_beige_branched <- as.numeric(root_data$thin_beige_branched)
    root_data$thin_beige_unbranched <- as.numeric(root_data$thin_beige_unbranched)
    root_data$thick_beige_unbranched <- as.numeric(root_data$thick_beige_unbranched)
    root_data$thin_brown_branched <- as.numeric(root_data$thin_brown_branched )
    root_data$thin_brown_unbranched <- as.numeric(root_data$thin_brown_unbranched )
    root_data$thick_brown_unbranched <- as.numeric(root_data$thick_brown_unbranched )
    root_data$thin_black_branched  <- as.numeric(root_data$thin_black_branched)
    root_data$thin_black_unbranched  <- as.numeric(root_data$thin_black_unbranched)
    root_data$thick_black_unbranched  <- as.numeric(root_data$thick_black_unbranched)
    root_data$shrub_ratio <- as.numeric(root_data$shrub_ratio)
    root_data$graminoid_ratio <- as.numeric(root_data$graminoid_ratio)

ratio_data <- root_data %>%
  select(plant_id, treatment, shrub_ratio, graminoid_ratio, co2_flux_avg) %>%
  filter(shrub_ratio != "NA" | graminoid_ratio != "NA") %>%
  pivot_longer(cols = shrub_ratio:graminoid_ratio,
               names_to = "functional_type",
               values_to = "ratio") %>%
  filter(ratio != "NA")


# above:below biomass ~ functional type by treatment graph
ggplot(data = ratio_data, aes(x = treatment, y = ratio, fill = treatment)) +
  geom_boxplot(alpha = 0.8) +
  facet_grid(col = vars(functional_type),
             labeller = labeller(functional_type = c("graminoid_ratio" = "Graminoid",
                                                     "shrub_ratio" = "Shrub"))) +
  labs(x = "Treatment",
       y = "Above:Below biomass") +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  theme_classic() +
  theme(legend.position = c(0.12, 0.87),
        legend.title = element_blank(),
        legend.text = element_text(size = 9), 
        strip.background = element_rect(colour = "black", fill = "#EDEDED"),
        strip.text.x = element_text(size = 12, colour = "black"), 
        axis.text = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 13, margin = margin(t = 12)),
        axis.title.y = element_text(size = 13, margin = margin(r = 12)),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ratio_data %>%
  group_by(functional_type) %>%
  summarise(p_value = kruskal.test(ratio ~ treatment)$p.value) 

ratio_data_graminoid <-  ratio_data %>%
  filter(functional_type == "graminoid_ratio") 
dunn.test(ratio_data_graminoid$ratio, ratio_data_graminoid$treatment, method = "bonferroni", kw = TRUE)

ratio_data_shrub <-  ratio_data %>%
  filter(functional_type == "shrub_ratio") 
dunn.test(ratio_data_shrub$ratio, ratio_data_shrub$treatment, method = "bonferroni", kw = TRUE)
