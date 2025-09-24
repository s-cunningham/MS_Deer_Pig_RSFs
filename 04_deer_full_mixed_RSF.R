
library(tidyverse)
library(glmmTMB)

#### Read data ####
# deer <- read_csv("output/deer_used_avail_covariates_ind-wts_180m.csv")
deer <- read_csv("output/deer_used_avail_covariates.csv")

# Center and scale covariates
lyrs <- read_csv("output/deer_raster_mean_sds.csv")

deer <- deer |>
  mutate(allhardwoods=(allhardwoods-lyrs$mean[1]) / lyrs$sd[1],
         shrubs=(shrubs-lyrs$mean[2]) / lyrs$sd[2],
         gramanoids=(gramanoids-lyrs$mean[3]) / lyrs$sd[3],
         foodcrops=(foodcrops-lyrs$mean[4]) / lyrs$sd[4],
         developed=(developed-lyrs$mean[5]) / lyrs$sd[5],
         water_dist=(water_dist-lyrs$mean[6]) / lyrs$sd[6]) 

## Set up to loop by individual
un.id <- unique(deer$id)

## Run the full mixed-effects RSF
# Set up so we can fix the intercept variance
rsf.tmp <- glmmTMB(case ~ allhardwoods + gramanoids + foodcrops + shrubs + developed + water_dist + I(water_dist^2) + 
                     (1 | id) + # Random intercept
                     (0 + allhardwoods | id) + # Random slope
                     (0 + gramanoids | id) + # Random slope
                     (0 + foodcrops | id) + # Random slope
                     (0 + shrubs | id) + # Random slope
                     (0 + developed | id) + # Random slope
                     (0 + water_dist | id), # Random slope
               data=deer, family=binomial(link = "logit"), doFit = FALSE, weights=weight)

# Fix standard deviation of first random term, equal to variance of 10^6
rsf.tmp$parameters$theta[1] = log(1e3)

# We need to tell `glmmTMB` not to change the first entry of the vector of variances, and give all other variances another indicator to make sure they can be freely estimated:
rsf.tmp$mapArg = list(theta=factor(c(NA,1:6)))

# Fit the model
system.time(
rsf <- fitTMB(rsf.tmp)
)

summary(rsf)

saveRDS(rsf, "output/deer_mixed_model_rsf.RDS")
rsf <- readRDS("output/deer_mixed_model_rsf.RDS")

resid_all <- residuals(rsf, type = "pearson")
deer$residuals <- resid_all

ggplot(deer, aes(x = key, y = residuals)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Residuals by Individual", y = "Pearson Residual") +
  theme(axis.text.x=element_text(angle=90))

# Get fixed effect standard errors
se <- summary(rsf)$coefficients$cond[,"Std. Error"]
se <- se[-1]
se <- t(t(se))
se <- as.data.frame(se)
se <- rownames_to_column(se, var="covariate")
write_csv(se, "output/deer_mixed_effects_se.csv")

# Get fixed effect coefficients
beta <- fixef(rsf)$cond

# drop intercept and save named vector
beta <- beta[-1]
saveRDS(beta, "output/deer_mixed_effects_beta.RDS")

beta <- as.data.frame(beta)
beta <- rownames_to_column(beta, var="covariate")

beta <- beta |> filter(covariate!='(Intercept)') |>
  mutate(covariate=if_else(covariate=="I(dist_water^2)", "Water2", covariate))

write_csv(beta, "output/deer_mixed_effects_betas.csv")

# Covariance matrix of fixed effects
vcov_beta <- vcov(rsf)$cond
saveRDS(vcov_beta, "output/deer_mixed_effects_vcov.RDS")

