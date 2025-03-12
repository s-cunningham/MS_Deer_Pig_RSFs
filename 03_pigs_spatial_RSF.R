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
cor(pigs[,c(8:17)])

# create key column
pigs <- pigs %>%
          unite("key", c("id", "burst"), sep="_", remove=FALSE)

## Set up to loop by individual
un.id <- unique(pigs$key)

# Create list to save models
pigs_rsf <- list()

# Run loop
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- pigs %>% filter(key==un.id[i])

  # Run rsf (but maybe change this to LASSO)
  rsf <- glm(case ~ bottomlandhw + decidmixed + gramanoids + barren + water + I(water^2), data=temp, family=binomial(link = "logit"), weight=weight)
  summary(rsf)
  print(car::vif(rsf))
 
  
  # add model object to list
  pigs_rsf[[i]] <- rsf
  
}

coef(summary(pigs_rsf[[1]]))


cfs <- lapply(pigs_rsf, coef)

ggplot(temp, aes(x=foodcrops, y=case)) +
  geom_point()










