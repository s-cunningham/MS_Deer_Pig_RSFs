
library(tidyverse)
library(glmmTMB)

#### Read data ####
deer <- read_csv("output/deer_used_avail_covariates.csv")

# Center and scale covariates
lyrs <- read_csv("output/deer_raster_mean_sds.csv")

deer <- deer %>%
  mutate(allhardwoods=(allhardwoods-lyrs$mean[1]) / lyrs$sd[1],
         shrubs=(shrubs-lyrs$mean[2]) / lyrs$sd[2],
         gramanoids=(gramanoids-lyrs$mean[3]) / lyrs$sd[3],
         foodcrops=(foodcrops-lyrs$mean[4]) / lyrs$sd[4],
         developed=(developed-lyrs$mean[5]) / lyrs$sd[5],
         water_dist=(water_dist-lyrs$mean[6]) / lyrs$sd[6]) 

## Set up to loop by individual
un.id <- unique(deer$key)


## Run the full mixed-effects RSF
system.time(
  rsf <- glmmTMB(case ~ allhardwoods + gramanoids + foodcrops + shrubs + developed + water_dist + I(water_dist^2) +
                 (allhardwoods + gramanoids + foodcrops + shrubs + developed + water_dist | key), 
               data=deer, family=binomial(link = "logit"), weight=weight)
)


summary(rsf)

saveRDS(rsf, "output/deer_mixed_model_rsf.RDS")


beta <- fixef(rsf)$cond
beta <- as.data.frame(beta)
beta <- rownames_to_column(beta, var="covariate")

beta <- beta %>% filter(covariate!='(Intercept)') %>%
  mutate(covariate=if_else(covariate=="I(dist_water^2)", "Water2", covariate))

write_csv(beta, "output/deer_mixed_effects_betas.csv")

