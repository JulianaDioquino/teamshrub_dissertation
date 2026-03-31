library(tidyverse)
library(ggplot2)
library(dplyr)
library(lmerTest)
library(lme4)
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
    
co2_flux <- read_csv(file = "data_raw/shifted_co2_data.csv")

co2_flux_shifted <- co2_flux %>%
  mutate(measurement.week = ifelse(Treatment == "extended",
                                   measurement.week - 3,
                                   measurement.week),
         Treatment = factor(Treatment,
                            levels = c("Control", "Heatwave", "Extended"),
                            labels = c("Control", "Heat wave", "Extended season")))

root_morphology <- root_data %>%
  pivot_longer(cols = shrub_biomass:graminoid_biomass,
               names_to = "functional_type_biomass",
               values_to = "pft_biomass") %>%
  filter(pft_biomass != "NA")


# calculated co2 flux for each plant
avg_season_flux <- co2_flux_shifted %>%
  group_by(plant_id, Treatment) %>%
  summarise(total_co2_flux = sum(seasonal_co2_flux, na.rm = TRUE))


# combining co2 flux with root data frames
flux_root_combine <- left_join(avg_season_flux, root_data, by = 'plant_id') %>%
  filter(Treatment != "NA",
         total_root_biomass != "NA")


# co2 flux ~ root biomass by treatment graph
ggplot(data = flux_root_combine, aes(x = total_root_biomass, y = (total_co2_flux*86400), color = Treatment, fill = Treatment)) +
  geom_smooth(method = lm, alpha = 0.2, size = 1) +
  geom_point(alpha = 0.4) +
  labs(x = "Root biomass (g)",
       y = "Cumulative CO₂ respiration (µmol m⁻²)") +
  scale_color_manual(values = c("Control" = "#6c6563",
                                "Heat wave" = "#b56d5e",
                                "Extended season" = "#bbbc81")) +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  theme_classic() +
  theme(legend.position = c(0.13, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 15, margin = margin(t = 12)),
        axis.title.y = element_text(size = 15, margin = margin(r = 12))) +
  annotate("text",
           x = 0.023,
           y = 0.8e7,
           label = "Treatment: p < 0.001, Root biomass: p < 0.01",
           hjust = 0,
           size = 5.5)


# co2 flux ~ root biomass by treatment analysis
flux_root_combine_anova <- lm(total_co2_flux ~ treatment*total_root_biomass, data = flux_root_combine)
anova(flux_root_combine_anova)
plot(flux_root_combine_anova)
summary(flux_root_combine_anova)






# calculated co2 flux for each plant
avg_season_flux <- co2_flux_shifted %>%
  group_by(plant_id, Treatment) %>%
  summarise(total_co2_flux = sum(seasonal_co2_flux, na.rm = TRUE))


# combining co2 flux with root data frames
flux_pftroot_combine <- left_join(avg_season_flux, root_morphology, by = 'plant_id') %>%
  filter(Treatment != "NA",
         total_root_biomass != "NA")

# co2 flux ~ PFT biomass by treatment graph
ggplot(data = flux_pftroot_combine, aes(x = pft_biomass, y = (total_co2_flux*86400), color = Treatment, fill = Treatment)) +
  geom_smooth(method = lm, alpha = 0.2, size = 1) +
  geom_point(alpha = 0.4) +
  facet_grid(col = vars(functional_type_biomass),
             labeller = labeller(functional_type_biomass = c("graminoid_biomass" = "Graminoid",
                                                             "shrub_biomass" = "Shrub"))) + 
  labs(x = "Root biomass (g)",
       y = "Cumulative CO₂ respiration (µmol m⁻²)") +
  scale_color_manual(values = c("Control" = "#6c6563",
                                "Heat wave" = "#b56d5e",
                                "Extended season" = "#bbbc81")) +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81")) +
  theme_classic() +
  theme(legend.position = c(0.13, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 15, margin = margin(t = 12)),
        axis.title.y = element_text(size = 15, margin = margin(r = 12))) +
  annotate("text",
           x = 0.023,
           y = 0.8e7,
           label = "",
           hjust = 0,
           size = 5.5)

# co2 flux ~ PFT biomass by treatment analysis
flux_pftroot_combine_anova <- lm((total_co2_flux*86400) ~ treatment*pft_biomass*functional_type_biomass, data = flux_pftroot_combine)
  ## treatment has a sig. effect

anova(flux_pftroot_combine_anova)
plot(flux_pftroot_combine_anova)
summary(flux_pftroot_combine_anova)

shrub_data <- flux_pftroot_combine %>%
  filter(functional_type_biomass == "shrub_biomass")

    shrub_data_anova <- lm((total_co2_flux*86400) ~ treatment*pft_biomass, data = shrub_data)
      
      anova(shrub_data_anova)
          # treatment had a significant effect
      
gram_data <- flux_pftroot_combine %>%
  filter(functional_type_biomass == "graminoid_biomass")
    
    gram_data_anova <- lm((total_co2_flux*86400) ~ treatment*pft_biomass, data = gram_data)
    
      anova(gram_data_anova)
        # no significant effect
      

      
      
      
      
      
      
      
      
      
      
      


# bayesian analysis
prior_co2 <- set_prior("normal(5e6, 5e6)", class = "Intercept")

# prior for root biomass
prior_co2 <- c(prior_co2, set_prior("normal(5e7, 5e7)", class = "b", coef = "total_root_biomass"))
# prior for treatment effects
prior_co2 <- c(prior_co2, set_prior("normal(0, 5e6)", class = "b", coef = "TreatmentHeatwave"),
                          set_prior("normal(0, 5e6)", class = "b", coef = "TreatmentExtendedseason"))
# prior for interactions
prior_co2 <- c(prior_co2, set_prior("normal(0, 3e7)", class = "b", coef = "total_root_biomass:TreatmentHeatwave"),
                          set_prior("normal(2e7, 3e7)", class = "b", coef = "total_root_biomass:TreatmentExtendedseason"))
# prior for variance
prior_co2 <- c(prior_co2, set_prior("normal(0.02, 0.05)", class = "sigma"))


bayesian_model_respiration <- brm((total_co2_flux*86400) ~ total_root_biomass*Treatment, # |  trunc(lb=0)  # model formula
                       data = flux_root_combine, # dataset
                       iter = 5000, # number of smapling iteration
                       warmup = 1000, # discarded iterations at the start
                       cores = 3, # a core compute a chain, 3 time faster
                       chains = 3, # number of independant models that we want to converge
                       prior = prior_co2,
                       control = list(adapt_delta = 0.99),
                       family = gaussian(), # distribution
                       #threads = threading(3), # even faster
                       init = 0) # more stable sampling

summary(bayesian_model_respiration, prob = 0.9)
plot(bayesian_model_respiration)
pp_check(bayesian_model_respiration)



# statistical reporting with the fit models
new_data_co2