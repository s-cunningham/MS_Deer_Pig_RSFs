# ----------------------------------------------------------------------------------------------------------------------------
# Name: 
# Author: Stephanie Cunningham
# Objective: Plot functional response
# Input: RSF predictions (30x30)
# Output: 
# Required packages: tidyverse
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)

## Read in deer data
deer <- read_csv("output/deer_used_avail_covariates.csv")

# Center and scale covariates
deer[,7:19] <- scale(deer[,7:19])

## Set up to loop by individual
un.id <- unique(deer$key)



#### Summarize coefficients ####
# read in models
deer_rsf <- readRDS("output/deer_glms.RDS")

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


plot(r_glms$allhardwoods)
