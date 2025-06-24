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

  # Fit model and print the iteration for IDs that have warnings
  rsf <- tryCatch(bayesglm(case ~ allhardwoods + gramanoids + foodcrops + shrubs + developed + water_dist + I(water_dist^2), 
                           data=temp, family=binomial(link = "logit"), weight=weight), warning=function(w) print(i))
  
  # add model object to list
  deer_rsf[[i]] <- rsf
  
}

# Drop individuals that have complete separation
deer_rsf <- deer_rsf[-c(5,42)]
un.id <- un.id[-c(5,42)]

#### Summarize coefficients ####
# get covariance matrices
vcov_mat <- lapply(deer_rsf, vcov)

# Remove intercept row and column from each covariance matrix
vcov_mat <- lapply(vcov_mat, function(v) {
  idx <- which(rownames(v) != "(Intercept)")
  v[idx, idx, drop = FALSE]
})

flag <- lapply(vcov_mat, function(x) sum(diag(x) > 100)) %>% unlist()
flag <- ifelse(flag==2, 1, flag)
which(flag==1)

# Drop individuals that have complete separation
deer_rsf <- deer_rsf[-c(4,67,72)]
un.id <- un.id[-c(4,67,72)]

saveRDS(deer_rsf, "output/deer_glms.RDS")
# reduced_glms <- readRDS("output/deer_reduced_glms.RDS")

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
  reframe(across(allhardwoods:water_dist2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_betas, "output/deer_glm_betas.csv")

## Find between / within variance
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
saveRDS(final, "output/deer_vcov.rds")

## Extract SE for confidence intervals
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
