library(tidyverse)
library(glmmTMB)

#### Read data ####
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

# Center and scale covariates
lyrs <- read_csv("output/pigs_raster_mean_sds.csv")

pigs <- pigs %>%
  mutate(hardwoods=(hardwoods-lyrs$mean[1]) / lyrs$sd[1],
         shrubs=(shrubs-lyrs$mean[2]) / lyrs$sd[2],
         gramanoids=(gramanoids-lyrs$mean[3]) / lyrs$sd[3],
         foodcrops=(foodcrops-lyrs$mean[4]) / lyrs$sd[4],
         developed=(developed-lyrs$mean[5]) / lyrs$sd[5],
         dist_water=(dist_water-lyrs$mean[6]) / lyrs$sd[6]) %>%
  rename(weights=weight)

## Set up to loop by individual
un.id <- unique(pigs$key)

cor(pigs[,c("hardwoods","shrubs","gramanoids","developed","dist_water")])

## Run the mixed-effects RSF (~40 min)
system.time(
  rsf_model <- glmmTMB(
    formula = case ~ hardwoods + gramanoids + developed + dist_water + I(dist_water^2) + (1 + hardwoods + gramanoids + developed + dist_water | key),  # Random intercept + slopes
    data = pigs,
    family = binomial(link = "logit"),  # Logistic regression
    weights = weights,
  )
)

summary(rsf_model)

saveRDS(rsf_model, "output/pigs_mixed_model_rsf.RDS")

resid_all <- residuals(rsf_model, type = "pearson")
pigs$residuals <- resid_all

ggplot(pigs, aes(x = key, y = residuals)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Residuals by Individual", y = "Pearson Residual") +
  theme(axis.text.x=element_text(angle=90))

ggplot(pigs, aes(x = key, y = residuals)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal()

beta <- fixef(rsf_model)$cond
beta <- as.data.frame(beta)
beta <- rownames_to_column(beta, var="covariate")

beta <- beta %>% filter(covariate!='(Intercept)') %>%
  mutate(covariate=if_else(covariate=="I(dist_water^2)", "Water2", covariate))

write_csv(beta, "output/pigs_mixed_effects_betas.csv")
