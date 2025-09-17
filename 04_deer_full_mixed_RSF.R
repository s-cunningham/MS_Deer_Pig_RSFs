
library(tidyverse)
library(glmmTMB)

#### Read data ####
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
un.id <- unique(deer$key)

#
isw <- deer |> group_by(key, case) |> count() |> pivot_wider(names_from="case", values_from="n") |>
  mutate(ratio=`1`/`0`) |> mutate(weight = 5000*ratio) |>
  select(key, weight)

deer <- left_join(deer, isw, by="key")
deer <- deer |> mutate(weight.x = if_else(case==0, weight.y, weight.x)) |>
  select(X:water_dist) |>
  rename(weight=weight.x)

## Run the full mixed-effects RSF
system.time(
  rsf <- glmmTMB(case ~ allhardwoods + gramanoids + foodcrops + shrubs + developed + water_dist + I(water_dist^2) +
                   (1 + allhardwoods + gramanoids + foodcrops + shrubs + developed + water_dist | key), 
                 data=deer, family=binomial(link = "logit"), weight=weight)
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

