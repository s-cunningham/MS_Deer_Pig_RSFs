library(tidyverse)
library(glmmTMB)

#### Read data ####
pigs <- read_csv("output/pigs_used_avail_covariates_indwt_30m.csv")

## Per-individual scaling
# Check ratio of used:avail
ratios <- pigs |> group_by(key, case) |> count() |>
  pivot_wider(names_from=case, values_from=n) |>
  rename(used=`1`, avail=`0`) |>
  mutate(n_avail_used=avail/used) 

# total weight per ID
W <- 5000

pigs <- pigs |>
  left_join(ratios,by="key") |>
  mutate(weight = if_else(case==1, 1, W/avail))

# Center and scale covariates
lyrs <- read_csv("output/pigs_raster_mean_sds.csv")

pigs <- pigs |>
  mutate(shrubs=(shrubs-unlist(unname(lyrs[lyrs$layer=="shrubs", 2]))) / unlist(unname(lyrs[lyrs$layer=="shrubs", 3])),
         gramanoids=(gramanoids-unlist(unname(lyrs[lyrs$layer=="gramanoids", 2]))) / unlist(unname(lyrs[lyrs$layer=="gramanoids", 3])),
         developed=(developed-unlist(unname(lyrs[lyrs$layer=="developed", 2]))) / unlist(unname(lyrs[lyrs$layer=="developed", 3])),
         hardwoods=(hardwoods-unlist(unname(lyrs[lyrs$layer=="hardwoods", 2]))) / unlist(unname(lyrs[lyrs$layer=="hardwoods", 3])),
         dist_water=(dist_water-unlist(unname(lyrs[lyrs$layer=="dist_water", 2]))) / unlist(unname(lyrs[lyrs$layer=="dist_water", 3]))) 

## Set up to loop by individual
un.id <- unique(pigs$key)

cor(pigs[,c("hardwoods","shrubs","gramanoids","developed","dist_water")])

## Run the full mixed-effects RSF
# Set up so we can fix the intercept variance
rsf.tmp <- glmmTMB(case ~ hardwoods + gramanoids + shrubs + developed + dist_water + I(dist_water^2) + 
                     (1|key) + # Random intercept
                     (0 + hardwoods | key) + # Random slope
                     (0 + gramanoids | key) + # Random slope
                     (0 + shrubs | key) + # Random slope
                     (0 + developed | key) + # Random slope
                     (0 + dist_water | key) + # Random slope
                     (0 + I(dist_water^2) | key), 
                   data=pigs, family=binomial(link = "logit"), doFit = FALSE, weights=weight)

# Fix standard deviation of first random term, equal to variance of 10^6
rsf.tmp$parameters$theta[1] = log(1e3)

# We need to tell `glmmTMB` not to change the first entry of the vector of variances, and give all other variances another indicator to make sure they can be freely estimated:
rsf.tmp$mapArg = list(theta=factor(c(NA,1:6)))

# Fit the model
system.time(
  rsf <- fitTMB(rsf.tmp)
)

summary(rsf)

saveRDS(rsf, "output/pigs_mixed_model_rsf_indwt_30m.RDS")
# rsf <- readRDS("output/pigs_mixed_model_rsf.RDS")

resid_all <- residuals(rsf, type = "pearson")
pigs$residuals <- resid_all

resid_summary <- pigs |>
  group_by(key) |>
  reframe(
    mean_resid = mean(residuals, na.rm = TRUE),
    sd_resid   = sd(residuals, na.rm = TRUE),
    max_resid  = max(residuals, na.rm = TRUE),
    min_resid  = min(residuals, na.rm = TRUE),
    n_large    = sum(abs(residuals) > 2, na.rm = TRUE),
    n_total    = n()
  ) |>
  arrange(desc(n_large))

# Get fixed effect standard errors
se <- summary(rsf)$coefficients$cond[,"Std. Error"]
se <- se[-1]
se <- t(t(se))
se <- as.data.frame(se)
se <- rownames_to_column(se, var="covariate")
write_csv(se, "output/pigs_mixed_effects_se_indwt_30m.csv")

# Get fixed effect coefficients
beta <- fixef(rsf)$cond

# drop intercept and save named vector
beta <- beta[-1]
saveRDS(beta, "output/pigs_mixed_effects_beta_indwt_30m.RDS")

beta <- as.data.frame(beta)
beta <- rownames_to_column(beta, var="covariate")

beta <- beta |> filter(covariate!='(Intercept)') |>
  mutate(covariate=if_else(covariate=="I(dist_water^2)", "Water2", covariate))

write_csv(beta, "output/pigs_mixed_effects_betas_indwt_30m.csv")

# Covariance matrix of fixed effects
vcov_beta <- vcov(rsf)$cond
saveRDS(vcov_beta, "output/pigs_mixed_effects_vcov_indwt_30m.RDS")

# look at coefficient SDs
summary(rsf)$varcor