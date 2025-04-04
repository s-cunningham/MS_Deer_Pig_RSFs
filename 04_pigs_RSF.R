# ----------------------------------------------------------------------------------------------------------------------------
# Name: Resource Selection functions for Pigs
# Author: Stephanie Cunningham
# Objective: Run the resource selection function models for wild pigs via elastic net and estimate coefficients from a GLM
# Input: Used and available points, with associated covariates extracted from rasters
# Output: Coefficients (beta and SE of beta) from RSF models
# Required packages: tidyverse, glmnet, glmnetUtils, recipes, arm
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(glmnetUtils)

#### Helper functions for cva.glmnet ####
# Copied from https://stackoverflow.com/questions/54803990/extract-the-best-parameters-from-cva-glmnet-object
# Get alpha.
get_alpha <- function(fit) {
  alpha <- fit$alpha
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  alpha[which.min(error)]
}

# Get all parameters.
get_model_params <- function(fit) {
  alpha <- fit$alpha
  lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
  error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
             lambdaSE = lambdaSE[best], eror = error[best])
}

#### Read data ####
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

## Set up to loop by individual
un.id <- unique(pigs$key)

# Create list to save models
pigs_rsf <- list()

# list to save coefficients
cfs <- list()

# # create a list to save VIF from GLMs
# glm_vifs <- list()

# Save GLMS with reduced variable sets
reduced_glms <- list()

#### Run RSF with elastic net (which might sometimes end up as a lasso or ridge) ####
# Loop to run LASSO over each individual HR
set.seed(1)
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- pigs %>% filter(key==un.id[i])
  
  # Set up data for LASSO in glmnet
  x <- model.matrix(case ~ bottomland + deciduous + evergreen + mixed + gramanoids + shrubs + foodcrops + water + I(water^2), temp)[, -1]
  y <- temp$case
  wts <- temp$weight
  
  # Run cross validation for alpha and lambda
  cv.out <- cva.glmnet(x, y, data=temp, weights=wts, family="binomial", type.measure="deviance")
  
  # Extract "best" paramters
  cv.out.p <- get_model_params(cv.out)
  
  # Print for status update
  print(paste0("Alpha for ", un.id[i], " is ", cv.out.p$alpha[1]))
  
  # Run the LASSO with the optimal lambda and alpha
  rsf <- glmnet(x, y, family="binomial", alpha=cv.out.p$alpha[1], weights=wts, lambda=cv.out.p$lambdaMin[1])
  
  # add model object to list
  pigs_rsf[[i]] <- rsf
  
  # Extract coefficients
  l.coef <- predict(rsf, type="coefficients", s=cv.out.p$lambdaMin[1])
  
  # Save LASSO coefficients
  cfs[[i]] <- as.matrix(t(l.coef))
  
  # Extract the non-zero coefficients
  coefs <- as.matrix(l.coef) %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(rowname!="(Intercept)" & s1!=0) %>%
    select(rowname) %>% 
    as.vector() %>% 
    unname() %>%
    unlist() %>%
    paste(collapse=" + ")
  
  if (str_length(coefs) > 0) {
    # Create a formula with the non-zero coefficients  
    mod_form <- formula(paste("case ~", coefs))
    
    # run regression model with selected covariates
    test <- glm(mod_form, data=temp, family=binomial(link = "logit"), weight=weight)
    
    # Save GLM with the reduced covariate set
    reduced_glms[[i]] <- test
    
    # glm_vifs[[i]] <- car::vif(test)
  }
}

saveRDS(reduced_glms, "output/pigs_reduced_glms.RDS")
# reduced_glms <- readRDS("output/pigs_reduced_glms.RDS")

#### reduced GLMs ####
# extract mean coeficients and combine into a single table
r_glms <- lapply(reduced_glms, coef)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- lapply(r_glms, t)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- do.call(bind_rows, r_glms)

r_glms <- r_glms %>% 
  as_tibble() %>%
  rename(intercept=`(Intercept)`, water2=`I(water^2)`) %>%
  # drop intercept
  select(-intercept) %>%
  # Reorder
  select(deciduous:gramanoids,foodcrops,water,water2) %>%
  # fill with 0
  mutate(across(deciduous:water2, \(x) coalesce(x, 0))) %>%
  mutate(water2 = if_else(water==0, 0, water2))
  
# Add IDs 
r_glms$id <- un.id

# Drop those with unreasonable SEs (probably because the AKDE was too big)
r_glms <- r_glms %>%
filter(id!="30251_Eastern_2" & id!="35493_Noxubee_6" & id!="35490_Noxubee_2" & 
         id!="19212_Delta_2" & id!="19212_Delta_3" & id!="19219_Delta_2") %>%
select(-id)

# Calculate mean across all individuals
glm_betas <- r_glms %>% 
  reframe(across(deciduous:water2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_betas, "output/pigs_glm_betas.csv")


# Extract SE for confidence intervals
se <- lapply(reduced_glms, arm::se.coef)

# combine into a tibble
se <- do.call(bind_rows, se)

# rename intercept and add ID column
se <- se %>%
  rename(intercept=`(Intercept)`, water2=`I(water^2)`) %>%
  # drop intercept
  select(-intercept) %>%
  # Reorder
  select(deciduous:gramanoids,foodcrops,water,water2) %>%
  # fill with 0
  mutate(across(deciduous:water2, \(x) coalesce(x, 0))) %>%
  mutate(water2 = if_else(water==0, 0, water2))

# Add IDs
se$id <- un.id

# Drop those with unreasonable SEs (probably because the AKDE was too big)
se <- se %>%
  filter(id!="30251_Eastern_2" & id!="35493_Noxubee_6" & id!="35490_Noxubee_2" & 
           id!="19212_Delta_2" & id!="19212_Delta_3" & id!="19219_Delta_2") %>%
  select(-id)

# Calculate mean across all individuals
glm_se <- se %>% 
  reframe(across(deciduous:water2, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_se, "output/pigs_glm_se.csv")

