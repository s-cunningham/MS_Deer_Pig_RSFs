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

#### Read data and set up for loop ####
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
lam.grid <- 10^seq(20, -2, length=200)
# Loop to run LASSO over each individual HR
set.seed(1)
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- pigs %>% filter(key==un.id[i])

  # Set up data for LASSO in glmnet
  x <- model.matrix(case ~ bottomland + water + I(water^2), temp)[, -1]
  y <- temp$case
  wts <- temp$weight
  
  # Run cross validation for alpha and lambda
  cv.out <- cva.glmnet(x, y, data=temp, weights=wts, family="binomial", type.measure="deviance")
  
  # Extract "best" paramters
  cv.out.p <- get_model_params(cv.out)

  # Print for status update
  print(paste0("Alpha for ", un.id[i], " is ", cv.out.p$alpha[1]))
  
  # cv.out <- cv.glmnet(x, y, data=temp, alpha=0, weights=wts, family="binomial", type.measure="deviance", lambda=lam.grid)
  
  # Run the LASSO with the optimal lambda and alpha
  rsf <- glmnet(x, y, family="binomial", alpha=cv.out.p$alpha[1], weights=wts, lambda=cv.out.p$lambdaMin[1])
  # rsf <- glmnet(x, y, family="binomial", alpha=0, weights=wts, lambda=cv.out$lambda.min)
  
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
  reframe(across(bottomland:water2, \(x) mean(x, na.rm = TRUE)))

# WRite to csv
write_csv(betas, "output/pigs_lasso_betas.csv")


plot(cfs$decidmixed, pch=16, col="blue")
abline(h=0)



# # get coefficients from all the models
# cfs <- lapply(pigs_rsf, coef)
# 
# # combine into a tibble
# cfs <- do.call(bind_rows, cfs)
# 
# # rename intercept and add ID column
# cfs <- cfs %>%
#   rename(intercept=`(Intercept)`) %>%
#   mutate(id = un.id) %>%
#   select(id, intercept:`I(water^2)`)
# 
# # calculate column means
# betas <- cfs %>% 
#   reframe(across(bottomlandhw:evergreen, \(x) mean(x, na.rm = TRUE)))
# 
# # Extract SE for confidence intervals
# se <- lapply(pigs_rsf, arm::se.coef)
# 
# # combine into a tibble
# se <- do.call(bind_rows, se)
# 
# # rename intercept and add ID column
# se <- se %>%
#   rename(intercept=`(Intercept)`) %>%
#   mutate(id = un.id) %>%
#   select(id, intercept:evergreen)


#### Functional response ####

# Read raster with land cover. Is already same projection as layers
lulc <- rast("data/landscape_data/cdl23_classes.tif")

# create vector of class names
lc <- c("palatablecrops", "barren", "gramanoids", "decidmixed", "evergreen", 
        "developed", "bottomlandhw", "water", "othercrops", "shrublands")

# Set up NLCD classes and class values
cdl_values <- c(1,2,3,4,5,6,7,8,9,10)
lc <- c("foodcrops", "barren", "gramanoids", "decidmixed", "evergreen", 
        "developed", "bottomlandhw", "water", "othercrops", "shrublands")

# Add class names and numbers to the raster
levels(lulc) <- list(data.frame(ID = cdl_values,
                                landcov = lc))

# extract classes at each point
pigs_lc <- extract(lulc, pigs_v) 

# join data so we have a tibble with points and land cover
hr_lc <- pigs %>% 
  bind_cols(pigs_lc) %>%
  select(key, id, case, landcov)

# Summarize by key...want the number of points in available for each landcover class
hr_lc <- hr_lc %>%
  # just want the available points
  filter(case==0) %>%
  # group by key
  group_by(key, landcov) %>%
  # summarize factor levels
  count() 

# Calculate how many points (cells) in each home range
ncells <- pigs %>% 
  filter(case==0) %>%
  group_by(key) %>%
  count() %>%
  rename(ntotal=n)

# join home range cells to hr_lc
hr_lc <- hr_lc %>%
  left_join(ncells, by="key") %>%
  ungroup() %>%
  # filter so we only have the rows that were in the model
  filter(landcov %in% m.covars) %>%
  # calculate % of home range
  mutate(pct=n/ntotal)


# pivot to long format with a column for landcov and a column for beta (the value)
cfs <- cfs %>% 
  # drop the intercept
  select(-intercept) %>%
  # pivot
  pivot_longer(2:7, names_to="landcov", values_to="beta") %>%
  # rename to match other data
  rename(key=id)

# now join betas to home range info
hr_lc <- left_join(hr_lc, cfs, by=c("key", "landcov"))














