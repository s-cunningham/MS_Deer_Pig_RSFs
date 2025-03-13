library(tidyverse)
library(terra)
library(lme4)
library(glmnet)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
deer <- read_csv("output/deer_used_avail_locations.csv") %>%
  # Add column for weight
  mutate(weight=if_else(case==1, 1, 5000))

# Rasters
rast_list <- c("data/landscape_data/bottomlandHW_180m_sum.tif",
               "data/landscape_data/deciduous_180m_sum.tif",
               "data/landscape_data/mixed_180m_sum.tif",
               "data/landscape_data/shrublands_180m_sum.tif",
               "data/landscape_data/othercrops_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif", 
               "data/landscape_data/evergreen_180m_sum.tif",
               "data/landscape_data/barren_180m_sum.tif", 
               "data/landscape_data/developed_180m_sum.tif",
               "data/landscape_data/palatable_crops_180m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Center and scale continuous rasters
layers <- scale(layers)

# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
deer_v <- vect(deer, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Rename layers
names(layers) <- c("bottomlandhw", "decidmixed", "othercrops", "gramanoids", "evergreen", "barren", "developed", "foodcrops", "water")

#### Extract covariates and used and available locations ####

# Extract distance values at used and available locations
dat_deer <- extract(layers, deer_v)

# Join extracted data back to location data frame
deer <- bind_cols(deer, dat_deer)

# correlation matrix
cor(deer[,c(8:16)])

#### Run RSF ####

# create key column
deer <- deer %>%
  unite("key", c("id", "burst"), sep="_", remove=FALSE)

deer$burst <- as.character(deer$burst)

## Set up to loop by individual
un.id <- unique(deer$id)

# Create list to save models
deer_rsf <- list()

# Set up lasso lambda grid
set.seed(1)
grid <- 10^seq(10, -2, length=100)

# Loop to run LASSO over each individual HR
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- deer %>% filter(id==un.id[i])
  
  # # Run rsf (but maybe change this to LASSO)
  # rsf <- glm(case ~ bottomlandhw + decidmixed + gramanoids + foodcrops + evergreen, 
  #              data=temp, family=binomial(link = "logit"), weight=weight)
  # summary(rsf)
  # print(car::vif(rsf))
 
  # Set up data for LASSO in glmnet

  x <- model.matrix(case ~ bottomlandhw + decidmixed + gramanoids + foodcrops + evergreen, temp)[, -1]
  y <- temp$case
  wts <- temp$weight

  lasso.mod <- glmnet(x, y, family="binomial", alpha=1, weights=wts, lambda=grid)
  plot(lasso.mod)

  cv.out <- cv.glmnet(x, y, family="binomial", alpha=1, weights=wts)
  bestlam <- cv.out$lambda.min

  rsf <- glmnet(x, y, family="binomial", alpha=1, weights=wts, lambda=bestlam)
  coef(rsf, s = bestlam)

  lasso.coef <- predict(rsf, type="coefficients", s=bestlam)
  lasso.coef[lasso.coef != 0]
  # ### LASSO
  
  
  # add model object to list
  deer_rsf[[i]] <- rsf
  
}

m.covars <- c("bottomlandhw", "decidmixed", "gramanoids", "foodcrops", "developed", "evergreen")

#### Summarize coefficients ####
# get coefficients from all the models
cfs <- lapply(deer_rsf, coef)

# combine into a tibble
cfs <- do.call(bind_rows, cfs)

# rename intercept and add ID column
cfs <- cfs %>%
  rename(intercept=`(Intercept)`) %>%
  mutate(id = un.id) %>%
  select(id, intercept:evergreen)

# calculate column means
betas <- cfs %>% 
    reframe(across(bottomlandhw:evergreen, \(x) mean(x, na.rm = TRUE)))

# Extract SE for confidence intervals
se <- lapply(deer_rsf, arm::se.coef)

# combine into a tibble
se <- do.call(bind_rows, se)

# rename intercept and add ID column
se <- se %>%
  rename(intercept=`(Intercept)`) %>%
  mutate(id = un.id) %>%
  select(id, intercept:evergreen)


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



