library(tidyverse)
library(glmnet)
library(glmnetUtils)

## Read data 
deer <- read_csv("output/deer_used_avail_covariates.csv")

## Set up to loop by individual
un.id <- unique(deer$key)

# Create list to save models
deer_rsf <- list()

# List to save best lambdas
lambda_vals <- list()

# list to save coefficients
cfs <- list()

# # create a list to save VIF from GLMs
# glm_vifs <- list()

# Save GLMS with reduced variable sets
reduced_glms <- list()

# Set up lasso lambda grid
set.seed(1)
lam.grid <- 10^seq(20, -2, length=200)

# alpha.grid <- seq(0,1,0.02)

# Loop to run LASSO over each individual HR
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- deer %>% filter(key==un.id[i])
  
  # Set up data for LASSO in glmnet
  x <- model.matrix(case ~ deciduous + evergreen + gramanoids + shrubs + foodcrops + water + I(water^2), temp)[, -1]
  y <- temp$case
  wts <- temp$weight

  cv.out <- cv.glmnet(x, y, family="binomial", alpha=1, weights=wts, type.measure="deviance")
  
  bestlam <- cv.out$lambda.min
  lambda_vals[[i]] <- bestlam
  
  # Run the LASSO with the optimal lambda and alpha
  rsf <- glmnet(x, y, family="binomial", alpha=1, weights=wts, lambda=bestlam)
  
  # add model object to list
  deer_rsf[[i]] <- rsf
  
  # Extract coefficients
  l.coef <- predict(rsf, type="coefficients", s=bestlam)
  
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
  # fill with 0
  mutate(across(deciduous:water2, \(x) coalesce(x, 0)))

# how many do we have?
nrow(r_glms)

# add rows of 0s to IDs that dropped all coeficients
zero_rows <- tibble(deciduous=rep(0,3),
                      evergreen=rep(0,3),
                      gramanoids=rep(0,3), 
                      shrubs=rep(0,3),
                      foodcrops=rep(0,3),
                      water=rep(0,3),
                      water2=rep(0,3))
r_glms <- bind_rows(r_glms, zero_rows)

glm_betas <- r_glms %>% 
  reframe(across(deciduous:water2, \(x) mean(x, na.rm = TRUE)))






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
deer_lc <- extract(lulc, deer_v) 

# join data so we have a tibble with points and land cover
hr_lc <- deer %>% 
  bind_cols(deer_lc) %>%
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
ncells <- deer %>% 
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


btmhw <- hr_lc %>% filter(landcov=="bottomlandhw")  

ggplot(btmhw, aes(x=pct, y=beta)) +
  coord_cartesian(xlim=c(0,1)) +
  geom_hline(yintercept=0) +
  geom_point() + 
  geom_smooth(method='lm', formula= y~x)

s1 <- hr_lc %>% filter(landcov=="evergreen")  

ggplot(s1, aes(x=pct, y=beta)) +
  geom_hline(yintercept=0) +
  geom_point() + 
  geom_smooth(method='lm', formula= y~x)




