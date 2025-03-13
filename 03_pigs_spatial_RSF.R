library(tidyverse)
library(terra)
library(lme4)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
pigs <- read_csv("output/pigs_used_avail_locations.csv") %>%
            # Add column for weight
            mutate(weight=if_else(case==1, 1, 5000))

# Rasters
rast_list <- c("data/landscape_data/bottomlandHW_210m_sum.tif",
               "data/landscape_data/decidmixed_210m_sum.tif",
               "data/landscape_data/othercrops_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/evergreen_210m_sum.tif",
               "data/landscape_data/barren_210m_sum.tif", 
               "data/landscape_data/developed_210m_sum.tif",
               "data/landscape_data/palatable_crops_210m_sum.tif") 
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

#read tree camopy cover
tcc <- rast("data/landscape_data/nlcd_tcc2021_ms50km_range210m.tif")
tcc <- project(tcc, crs(layers))
tcc <- resample(tcc, layers)
ext(tcc) <- ext(layers)

layers <- c(layers, tcc)

# Center and scale continuous rasters
layers <- scale(layers)

# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
pigs_v <- vect(pigs, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Rename layers
names(layers) <- c("bottomlandhw", "decidmixed", "othercrops", "gramanoids", "evergreen", "barren", "developed", "foodcrops", "water", "canopy")

#### Extract covariates and used and available locations ####

# Extract distance values at used and available locations
dat_pigs <- extract(layers, pigs_v)

# Join extracted data back to location data frame
pigs <- bind_cols(pigs, dat_pigs)

ggplot(pigs,aes(x=water, y=case)) +
  geom_point(alpha=0.3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) 

ggplot(pigs, aes(x=water, y=foodcrops, color=id)) +
  geom_point(alpha=0.1) +
  facet_wrap(vars(case)) +
  theme(legend.position="none")


#### Run RSF ####
cor(pigs[,c(9:18)])

# create key column
pigs <- pigs %>%
          unite("key", c("id", "burst"), sep="_", remove=FALSE)

## Set up to loop by individual
un.id <- unique(pigs$key)

# Create list to save models
pigs_rsf <- list()

# Create list to save VIF

# Run loop
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- pigs %>% filter(key==un.id[i])

  # Run rsf (but maybe change this to LASSO)
  rsf <- glm(case ~ canopy + I(canopy^2) + water + I(water^2), data=temp, family=binomial(link = "logit"), weight=weight)
  summary(rsf)
  print(car::vif(rsf))
 
  # add model object to list
  pigs_rsf[[i]] <- rsf
  
}

m.covars <- c("bottomlandhw", "decidmixed", "gramanoids", "canopy", "canopy^2", "water", "water^2")

#### Summarize coefficients ####
# get coefficients from all the models
cfs <- lapply(pigs_rsf, coef)

# combine into a tibble
cfs <- do.call(bind_rows, cfs)

# rename intercept and add ID column
cfs <- cfs %>%
  rename(intercept=`(Intercept)`) %>%
  mutate(id = un.id) %>%
  select(id, intercept:`I(water^2)`)

# calculate column means
betas <- cfs %>% 
  reframe(across(bottomlandhw:evergreen, \(x) mean(x, na.rm = TRUE)))

# Extract SE for confidence intervals
se <- lapply(pigs_rsf, arm::se.coef)

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








