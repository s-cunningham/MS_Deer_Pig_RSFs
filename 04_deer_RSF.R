# ----------------------------------------------------------------------------------------------------------------------------
# Name: Resource Selection functions for Deer
# Author: Stephanie Cunningham
# Objective: Run the resource selection function models for deer via elastic net and estimate coefficients from a GLM
# Input: Used and available points, with associated covariates extracted from rasters
# Output: Coefficients (beta and SE of beta) from RSF models
# Required packages: tidyverse, glmnet, glmnetUtils, recipes, arm
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(glmnetUtils)
library(recipes)

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
deer <- read_csv("output/deer_used_avail_covariates.csv")

# Center and scale covariates
deer[,9:19] <- scale(deer[,9:19])

## Set up to loop by individual
un.id <- unique(deer$key)

# Check ratio of used:available for each deer
ua <- deer %>% 
  select(key, case) %>% 
  group_by(key, case) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(names_from="case", values_from="n") %>%
  mutate(used_avail=`1`/`0`,
         npts=`1`+`0`) %>%
  select(key, npts, used_avail)


#### Run the RSF ####
# Create list to save models
deer_rsf <- list()

# list to save coefficients
cfs <- list()

# # create a list to save VIF from GLMs
glm_vifs <- list()

# Save GLMS with reduced variable sets
reduced_glms <- list()
glm_names <- character()

# Create vector to indicate when warnings are thrown
fail <- numeric(length(un.id))
lasso_fail <- numeric(length(un.id))

#### Run RSF with elastic net (which might sometimes end up as a lasso or ridge) ####
# Loop to run LASSO over each individual HR
set.seed(1)
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- deer %>% filter(key==un.id[i])
  
  # Set up data for LASSO in glmnet
  x <- model.matrix(case ~ deciduous + evergreen + bottomland + gramanoids + shrubs + foodcrops + developed, temp)[, -1]
  y <- temp$case
  wts <- temp$weight
  
  # Run cross validation for alpha and lambda
  tryCatch(cv.out <- cva.glmnet(x, y, data=temp, weights=wts, family="binomial", type.measure="deviance"),
           warning = function(w) {
             print(w)
             lasso_fail[i] <<- 1
           })
  
  # Extract "best" paramters
  cv.out.p <- get_model_params(cv.out)
  
  # Print for status update
  print(paste0("Alpha for ", un.id[i], " is ", cv.out.p$alpha[1]))
  
  # Run the LASSO with the optimal lambda and alpha
  rsf <- glmnet(x, y, family="binomial", alpha=cv.out.p$alpha[1], weights=wts, lambda=cv.out.p$lambdaMin[1]) #l

  # add model object to list
  deer_rsf[[i]] <- rsf
  
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
    tryCatch(test <- glm(mod_form, data=temp, family=binomial(link = "logit"), weight=weight),
             warning = function(w) {
               print(w)
               fail[i] <<- 1
             })
    
    # test2 <- arm::bayesglm(mod_form, data=temp, family=binomial(link = "logit"), weight=weight)
    
    # Save GLM with the reduced covariate set
    reduced_glms[[i]] <- test
    
    glm_names <- c(glm_names, un.id[i])
    glm_vifs[[i]] <- car::vif(test)
  }
}

# add names to list

saveRDS(reduced_glms, "output/deer_reduced_glms.RDS")
# reduced_glms <- readRDS("output/deer_reduced_glms.RDS")

glm_summary <- lapply(reduced_glms, summary)

#### reduced GLMs ####
## Extract SE for confidence intervals
se <- lapply(reduced_glms, arm::se.coef)   # is ordering right??

# combine into a tibble
se <- do.call(bind_rows, se)

# rename intercept and add ID column
se <- se %>%
  rename(intercept=`(Intercept)`) %>%  
  # drop intercept
  select(-intercept) %>%
  # fill with 0
  mutate(across(deciduous:developed, \(x) coalesce(x, 0)))

# Add IDs 
se$id <- glm_names

# Drop those with bad models
se <- se[which(fail==0),]

# Calculate mean across all individuals
glm_se <- se %>% 
  reframe(across(deciduous:developed, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_se, "output/deer_glm_se.csv")



## extract mean beta coeficients and combine into a single table
r_glms <- lapply(reduced_glms, coef)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- lapply(r_glms, t)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- do.call(bind_rows, r_glms)

# Add IDs 
r_glms$id <- glm_names

r_glms <- r_glms %>% 
  as_tibble() %>%
  rename(intercept=`(Intercept)`) %>%
  # drop intercept
  select(-intercept) %>%
  # fill with 0
  mutate(across(deciduous:developed, \(x) coalesce(x, 0))) 

# Drop those with bad models
r_glms <- r_glms[which(fail==0),]

# Calculate mean across all individuals
glm_betas <- r_glms %>% 
  reframe(across(deciduous:developed, \(x) mean(x, na.rm = TRUE)))

write_csv(glm_betas, "output/deer_glm_betas.csv")


