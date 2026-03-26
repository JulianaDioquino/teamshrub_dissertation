library(tidyverse)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(data.table)

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


root_morphology <- data.table(root_morphology)
root_morphology[,pft_biomass := as.numeric(pft_biomass)]
root_morphology <- root_morphology[!is.na(pft_biomass),]
root_morphology <- root_morphology[!duplicated(root_morphology),]

# root biomass ~ functional type by treatment graph
ggplot(data = root_morphology, aes(x = treatment, y = pft_biomass, fill = treatment)) +
  geom_boxplot(alpha = 0.8,outliers = F) +
  facet_grid(col = vars(functional_type_biomass),
             labeller = labeller(functional_type_biomass = c("graminoid_biomass" = "Graminoid",
                                                             "shrub_biomass" = "Shrub"))) + 
  labs(x = "Treatment",
       y = "Root biomass (g)") +
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




# root biomass ~ functional type by treatment analysis
root_morphology %>%
  group_by(functional_type_biomass) %>%
  summarise(p_value = kruskal.test(pft_biomass ~ treatment)$p.value)

graminoid <-  root_morphology %>%
  filter(functional_type_biomass == "graminoid_biomass") 
dunn.test(graminoid$pft_biomass, graminoid$treatment, method = "bonferroni")

shrub <-  root_morphology %>%
  filter(functional_type_biomass == "shrub_biomass") 
dunn.test(shrub$pft_biomass, shrub$treatment, method = "bonferroni")

library(brms)

prior_root_biomass <- set_prior("normal(-0.03,0.05)" , class = "Intercept")
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.02)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.02)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.01,0.3)" , class = "sigma"))

hist(graminoid$pft_biomass)
hist(shrub$pft_biomass)

bayesian_model <- brm( pft_biomass |  trunc(lb = 0) ~   treatment , # |  trunc(lb=0) 
                       data = graminoid,
                       iter = 5000,
                       warmup = 1000,
                       cores = 3,
                       chains = 3,
                       prior = prior_root_biomass ,
                       family = gaussian(),
                       threads = threading(3),
                       init = 0)



prior_root_biomass <- set_prior("normal(-0.05,0.1)" , class = "Intercept")
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.04)" , class = "b" , coef  = "treatmentHeatwave"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0,0.07)" , class = "b" , coef  = "treatmentExtendedseason"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(0.02,0.3)" , class = "sigma"))
prior_root_biomass <- c(prior_root_biomass,set_prior("normal(5,3)" , class = "alpha"))


hist(graminoid$pft_biomass)

bayesian_model <- brm( pft_biomass |  trunc(lb=0)   ~  treatment , # |  trunc(lb=0)  # model formula
                       data = shrub, # dataset
                       iter = 5000, # number of smapling iteration
                       warmup = 1000, # discarded iterations at the start
                       cores = 3, # a core compute a chain, 3 time faster
                       chains = 3, # number of independant models that we want to converge
                       prior = prior_root_biomass ,
                       family = gaussian(), # distribution
                       #threads = threading(3), # even faster
                       init = 0) # more stable sampling

summary(bayesian_model,prob = 0.9)

new_data <-  data.table(treatment = sort(unique(shrub$treatment)))

preds <- cbind(new_data,fitted(bayesian_model,newdata = new_data,probs = c(0.05,0.95)))

preds_full <- data.table(fitted(bayesian_model,newdata = new_data,summary = F))
colnames(preds_full) <- as.character(new_data$treatment)

preds_delta <- preds_full[,.(`Heat wave` = `Heat wave` - Control ,
                             `Extended season` =`Extended season` - Control)]

preds_full <- melt(preds_full,variable.name = "treatment",value.name = "Estimate")
preds_full$treatment <- as.factor(preds_full$treatment)

preds_delta <- melt(preds_delta,variable.name = "treatment",value.name = "Estimate")

ggplot(preds,aes(x = treatment,y = Estimate))+
  geom_violin(data = preds_full,aes(fill = treatment))+
  geom_pointrange(aes(ymin =  Q5, ymax = Q95))+
  geom_point(data = root_data,aes(y = ,x = treatment))+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                             "Heat wave" = "#b56d5e",
                              "Extended season" = "#bbbc81")) 

preds_delta[,1-mean(Estimate>0),by = treatment]
preds_delta[,mean(Estimate),by = treatment]

ggplot(preds_delta,aes(x = Estimate,fill = treatment ))+
  geom_vline(xintercept = 0,lty = 2)+
  geom_density(alpha = 0.5)+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))


ggplot(preds_delta,aes(x = treatment, y = Estimate ,fill = treatment ))+
  geom_hline(yintercept = 0,lty = 2)+
  geom_violin(alpha = 0.95)+
  theme_classic()+  
  scale_fill_manual(values = c("Control" = "#6c6563",
                               "Heat wave" = "#b56d5e",
                               "Extended season" = "#bbbc81"))

plot(bayesian_model)
pp_check(bayesian_model)

# https://distribution-explorer.github.io/

