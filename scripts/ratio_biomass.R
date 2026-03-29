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
  theme(legend.position = c(0.15, 0.87),
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


# ratio biomass ~ functional type by treatment analysis
ratio_data %>%
  group_by(functional_type) %>%
  summarise(p_value = kruskal.test(ratio ~ treatment)$p.value) 

ratio_data_graminoid <-  ratio_data %>%
  filter(functional_type == "graminoid_ratio") 
dunn.test(ratio_data_graminoid$ratio, ratio_data_graminoid$treatment, method = "bonferroni", kw = TRUE)

ratio_data_shrub <-  ratio_data %>%
  filter(functional_type == "shrub_ratio") 
dunn.test(ratio_data_shrub$ratio, ratio_data_shrub$treatment, method = "bonferroni", kw = TRUE)




# shrub above:below bayesian analysis
prior_root_biomass <- set_prior("normal(-0.3, 0.05)" , class = "Intercept")
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.3)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.3)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.01,0.5)" , class = "sigma"))

bayesian_model_shrub <- brm( ratio |  trunc(lb = 0) ~   treatment , # |  trunc(lb=0) 
                            data = ratio_data_shrub,
                            iter = 5000,
                            warmup = 1000,
                            cores = 3,
                            chains = 3,
                            prior = prior_root_biomass ,
                            control = list(adapt_delta = 0.99),
                            family = gaussian(),
                            # threads = threading(3),
                            init = 0)

summary(bayesian_model_shrub, prob = 0.9)
plot(bayesian_model_shrub)
pp_check(bayesian_model_shrub)


# graminoid above:below bayesian analysis
prior_root_biomass <- set_prior("normal(2,2)" , class = "Intercept")
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(2, 2)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(3, 5)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.01,0.5)" , class = "sigma"))

bayesian_model_gram <- brm( ratio |  trunc(lb = 0) ~   treatment , # |  trunc(lb=0) 
                            data = ratio_data_graminoid,
                            iter = 5000,
                            warmup = 1000,
                            cores = 3,
                            chains = 3,
                            prior = prior_root_biomass ,
                            control = list(adapt_delta = 0.99),
                            family = gaussian(),
                            # threads = threading(3),
                            init = 0)

summary(bayesian_model_gram, prob = 0.9)
plot(bayesian_model_gram)
pp_check(bayesian_model_gram)


### SHRUB statistical reporting with the fit models ###
new_data <- data.table(treatment = sort(unique(ratio_data_shrub$treatment)))
preds <- cbind(new_data, fitted(bayesian_model_shrub,
                                newdata = new_data, probs = c(0.05, 0.95)))

preds_full <- data.table(fitted(bayesian_model_shrub,
                                newdata = new_data, summary = F))
colnames(preds_full) <- as.character(new_data$treatment)

preds_delta <- preds_full[,.(`Heat wave` = `Heat wave` - Control ,
                             `Extended season` =`Extended season` - Control)]

# plot set up
preds_full <- melt(preds_full, variable.name = "treatment", value.name = "Estimate")
preds_full$treatment <- as.factor(preds_full$treatment)

preds_delta <- melt(preds_delta, variable.name = "treatment", value.name = "Estimate")

# distribution of the mean plot
ggplot(preds, aes(x = treatment, y = Estimate)) +
  geom_violin(data = preds_full, aes(fill = treatment)) +
  geom_pointrange(aes(ymin = Q5, ymax = Q95), color = "white") + 
  geom_point(data = shrub, aes(y = pft_biomass, x = treatment)) +
  theme_classic() +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))

# how much of the delta control - treatment is positive? 
# if more than 90% i'm confident in treatment effect
preds_delta[,1-mean(Estimate < 0), by = treatment]
preds_delta[,mean(Estimate), by = treatment]

ggplot(preds_delta, aes(x = Estimate, fill = treatment)) + 
  geom_vline(xintercept = 0, lty = 2) +
  geom_density(alpha = 0.5) +
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))+
  facet_wrap(~treatment,scales = "free_y",nrow = 3)

# distribution of the delta control - treatment, JB doesn't like
ggplot(preds_delta,aes(x = treatment, y = Estimate ,fill = treatment ))+
  geom_hline(yintercept = 0,lty = 2)+
  geom_violin(alpha = 0.95)+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))



### GRAMINOID statistical reporting with the fit models ###
new_data <- data.table(treatment = sort(unique(ratio_data_graminoid$treatment)))
preds <- cbind(new_data, fitted(bayesian_model_gram,
                                newdata = new_data, probs = c(0.05, 0.95)))

preds_full <- data.table(fitted(bayesian_model_gram,
                                newdata = new_data, summary = F))
colnames(preds_full) <- as.character(new_data$treatment)

preds_delta <- preds_full[,.(`Heat wave` = `Heat wave` - Control ,
                             `Extended season` =`Extended season` - Control)]

# plot set up
preds_full <- melt(preds_full, variable.name = "treatment", value.name = "Estimate")
preds_full$treatment <- as.factor(preds_full$treatment)

preds_delta <- melt(preds_delta, variable.name = "treatment", value.name = "Estimate")

# distribution of the mean plot
ggplot(preds, aes(x = treatment, y = Estimate)) +
  geom_violin(data = preds_full, aes(fill = treatment)) +
  geom_pointrange(aes(ymin = Q5, ymax = Q95), color = "white") + 
  geom_point(data = graminoid, aes(y = pft_biomass, x = treatment)) +
  theme_classic() +
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))

# how much of the delta control - treatment is positive? 
# if more than 90% i'm confident in treatment effect
preds_delta[,1-mean(Estimate < 0), by = treatment]
preds_delta[,mean(Estimate), by = treatment]

ggplot(preds_delta, aes(x = Estimate, fill = treatment)) + 
  geom_vline(xintercept = 0, lty = 2) +
  geom_density(alpha = 0.5) +
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))+
  facet_wrap(~treatment,scales = "free_y",nrow = 3)

# distribution of the delta control - treatment, JB doesn't like
ggplot(preds_delta,aes(x = treatment, y = Estimate ,fill = treatment ))+
  geom_hline(yintercept = 0,lty = 2)+
  geom_violin(alpha = 0.95)+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))

