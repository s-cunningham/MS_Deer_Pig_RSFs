# ----------------------------------------------------------------------------------------------------------------------------
# Name: Resource Selection functions for Pigs
# Author: Stephanie Cunningham
# Objective: Run the resource selection function models for wild pigs via elastic net and estimate coefficients from a GLM
# Input: Used and available points, with associated covariates extracted from rasters
# Output: Coefficients (beta and SE of beta) from RSF models
# Required packages: tidyverse, glmnet, glmnetUtils, recipes, arm
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(arm)

#### Read data ####
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

# Center and scale covariates
pigs[,7:22] <- scale(pigs[,7:22])

## Set up to loop by individual
un.id <- unique(pigs$key)

# Create list to save models
pigs_rsf <- list()

#### Run RSF with elastic net (which might sometimes end up as a lasso or ridge) ####
# Loop to run LASSO over each individual HR
set.seed(1)
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- pigs %>% filter(key==un.id[i]) 
  
  # Set up data for LASSO in glmnet
  # x <- model.matrix(case ~ bottomland + deciduous + gramanoids + shrubs + developed + water, temp)[, -1]
  # y <- temp$case
  # wts <- temp$weight
  
  rsf <- bayesglm(case ~ hardwoods + shrubs + gramanoids + developed + dist_water + I(dist_water^2), 
                        data=temp, family=binomial(link = "logit"), weight=weight)
  
  pigs_rsf[[i]] <- rsf

}

saveRDS(pigs_rsf, "output/pigs_glms.RDS")
# reduced_glms <- readRDS("output/pigs_reduced_glms.RDS")

#### Summarizing coefficients ####
## extract mean coeficients and combine into a single table
r_glms <- lapply(pigs_rsf, coef)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- lapply(r_glms, t)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- do.call(bind_rows, r_glms)

r_glms <- r_glms %>% 
  as_tibble() %>%
  rename(intercept=`(Intercept)`, water2=`I(dist_water^2)`) %>%  
  # drop intercept
  dplyr::select(-intercept) %>%
  # Reorder
  dplyr::select(hardwoods:water2) %>%
  # fill with 0
  mutate(across(hardwoods:water2, \(x) coalesce(x, 0)))
  
# Create list
beta_list <- lapply(1:nrow(r_glms), function(i) {
  b <- as.numeric(r_glms[i, ])
  names(b) <- colnames(r_glms)
  b
})

# Add IDs 
r_glms$id <- un.id

# Calculate mean across all individuals
glm_betas <- r_glms %>% 
  reframe(across(hardwoods:water2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_betas, "output/pigs_glm_betas.csv")

## Find between / within variance
# get covariance matrices
vcov_mat <- lapply(pigs_rsf, vcov)

# Remove intercept row and column from each covariance matrix
vcov_mat <- lapply(vcov_mat, function(v) {
  idx <- which(rownames(v) != "(Intercept)")
  v[idx, idx, drop = FALSE]
})

# 1. Stack beta vectors into a matrix
B <- do.call(rbind, beta_list)  # N x K
B_bar <- colMeans(B)

# 2. Within-individual variance
within_var <- Reduce(`+`, vcov_mat) / length(beta_list)

# 3. Between-individual variance
B_centered <- sweep(B, 2, B_bar)
between_var <- crossprod(B_centered) / (nrow(B) - 1)

# 4. Total variance of average coefficients
total_var <- within_var / length(beta_list) + between_var

# Final output
beta_avg <- B_bar
vcov_avg <- total_var

final <- list(beta_avg, vcov_avg)

# save cov matrix
saveRDS(final, "output/pigs_vcov.rds")

## Extract SE for confidence intervals
se <- lapply(pigs_rsf, se.coef)

# combine into a tibble
se <- do.call(bind_rows, se)

# rename intercept and add ID column
se <- se %>%
  rename(intercept=`(Intercept)`, water2=`I(dist_water^2)`) %>%  
  # drop intercept
  dplyr::select(-intercept) %>%
  # Reorder
  dplyr::select(hardwoods:water2) %>%
  # fill with 0
  mutate(across(hardwoods:water2, \(x) coalesce(x, 0)))

# Add IDs
se$id <- un.id

# Calculate mean across all individuals
glm_se <- se %>% 
  reframe(across(hardwoods:water2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_se, "output/pigs_glm_se.csv")


