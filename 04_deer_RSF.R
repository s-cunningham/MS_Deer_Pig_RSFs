## Run the resource selection function for deer

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
deer <- read_csv("output/deer_used_avail_covariates.csv")

## Set up to loop by individual
deer <- read_csv("output/deer_used_avail_covariates.csv")

## Set up to loop by individual
un.id <- unique(deer$key)

# Create list to save models
deer_rsf <- list()

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
  temp <- deer %>% filter(key==un.id[i])
  
  # Set up data for LASSO in glmnet
  x <- model.matrix(case ~ deciduous + evergreen + gramanoids + shrubs + foodcrops + water + I(water^2), temp)[, -1]
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
    test <- glm(mod_form, data=temp, family=binomial(link = "logit"), weight=weight)
    
    # Save GLM with the reduced covariate set
    reduced_glms[[i]] <- test
    
    # glm_vifs[[i]] <- car::vif(test)
  }
}

#### Summarize coefficients ####
# combine into a tibble
cfs <- lapply(cfs, as.data.frame)
cfs <- do.call(bind_rows, cfs)

# rename columns and add ID
cfs <- cfs %>% 
  as_tibble() %>%
  rename(water2=`I(water^2)`, intercept=`(Intercept)`) %>%
  mutate(id=un.id) %>%
  select(id, intercept:water2)

# calculate column means
betas <- cfs %>% 
  reframe(across(deciduous:water2, \(x) mean(x, na.rm = TRUE)))
exp(betas)

write_csv(betas, "output/deer_lasso_betas.csv")


# calculate betas and sd for confidence interval
betas <- cfs %>% 
  reframe(across(deciduous:water2, list(mean=mean, sd=sd)))


betas <- betas %>% 
  # pivot to be longer, instead of single row for all info
  pivot_longer(1:14, names_to="var", values_to="value") %>%
  # Split column so that there is one for the covariate layer and one for the metric
  separate(var, into=c("covar", "metric"), sep="_") %>%
  # sort by metric (i.e. group means & SDs)
  arrange(metric) %>%
  # Expand so that means and SD are in their own columns
  pivot_wider(id_cols=1, names_from="metric", values_from="value") %>%
  # Calculate CIs
  mutate(lci = mean - 1.96*(sd/sqrt(120)),
         uci = mean + 1.96*(sd/sqrt(120)))

write_csv(betas, "output/deer_lasso_betas_uncert.csv")


