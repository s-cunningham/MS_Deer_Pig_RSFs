# ----------------------------------------------------------------------------------------------------------------------------
# Name: Resource Selection functions for Deer
# Author: Stephanie Cunningham
# Objective: Run the resource selection function models for deer via elastic net and estimate coefficients from a GLM
# Input: Used and available points, with associated covariates extracted from rasters
# Output: Coefficients (beta and SE of beta) from RSF models
# Required packages: tidyverse, arm
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(arm)

#### Read data ####
deer <- read_csv("output/deer_used_avail_covariates.csv")

# Center and scale covariates
deer[,7:19] <- scale(deer[,7:19])

## Set up to loop by individual
un.id <- unique(deer$key)

# Check ratio of used:available for each deer
ua <- deer %>% 
  dplyr::select(key, case) %>% 
  group_by(key, case) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(names_from="case", values_from="n") %>%
  mutate(used_avail=`1`/`0`,
         npts=`1`+`0`) %>%
  dplyr::select(key, npts, used_avail)


#### Run the RSF ####
# Create list to save models
deer_rsf <- list()

#### Run RSF with elastic net (which might sometimes end up as a lasso or ridge) ####
# Loop to run LASSO over each individual HR
set.seed(1)
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- deer %>% filter(key==un.id[i])

  # Run the LASSO with the optimal lambda and alpha
  rsf <- bayesglm(case ~ allhardwoods + gramanoids + shrubs + foodcrops + developed + water_dist + I(water_dist^2), data=temp, family=binomial(link = "logit"), weight=weight)

  # add model object to list
  deer_rsf[[i]] <- rsf
  
}

# add names to list

saveRDS(deer_rsf, "output/deer_glms.RDS")
# reduced_glms <- readRDS("output/deer_reduced_glms.RDS")

#### Summarize coefficients ####
# extract mean coeficients and combine into a single table
r_glms <- lapply(deer_rsf, coef)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- lapply(r_glms, t)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- do.call(bind_rows, r_glms)

r_glms <- r_glms %>% 
  as_tibble() %>%
  rename(intercept=`(Intercept)`, water_dist2=`I(water_dist^2)`) %>%  
  # drop intercept
  dplyr::select(-intercept) %>%
  # Reorder
  dplyr::select(allhardwoods:water_dist2) %>%
  # fill with 0
  mutate(across(allhardwoods:water_dist2, \(x) coalesce(x, 0)))

# Add IDs 
r_glms$id <- un.id

# Calculate mean across all individuals
glm_betas <- r_glms %>% 
  reframe(across(allhardwoods:water_dist2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_betas, "output/deer_glm_betas.csv")


# Extract SE for confidence intervals
se <- lapply(deer_rsf, se.coef)

# combine into a tibble
se <- do.call(bind_rows, se)

# rename intercept and add ID column
se <- se %>%
  rename(intercept=`(Intercept)`, water_dist2=`I(water_dist^2)`) %>%  
  # drop intercept
  dplyr::select(-intercept) %>%
  # Reorder
  dplyr::select(allhardwoods:water_dist2) %>%
  # fill with 0
  mutate(across(allhardwoods:water_dist2, \(x) coalesce(x, 0)))

# Add IDs
se$id <- un.id

# Calculate mean across all individuals
glm_se <- se %>% 
  reframe(across(allhardwoods:water_dist2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_se, "output/deer_glm_se.csv")
