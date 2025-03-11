library(tidyverse)
library(terra)
library(lme4)

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
               "data/landscape_data/palatable_crops_210m_sum.tif",
               "data/landscape_data/water_210m_sum.tif") # replace with water (distance)
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# Center and scale continuous rasters
layers <- scale(layers)

# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
pigs_v <- vect(pigs, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Rename layers
names(layers) <- c("bottomlandhw", "decidmixed", "othercrops", "gramanoids", "evergreen", "barren", "developed", "foodcrops", "water")

#### Extract covariates and used and available locations ####

# Extract distance values at used and available locations
dat_pigs <- extract(layers, pigs_v)

# Join extracted data back to location data frame
pigs <- bind_cols(pigs, dat_pigs)

#### Run RSF ####
cor(pigs[,c(8:16)])

# create key column
pigs <- pigs %>%
          unite("key", c("id", "burst"), sep="_", remove=FALSE)

## Set up to loop by individual
un.id <- unique(pigs$key)

# Create list to save models
pigs_rsf <- list()

library(glmnet)

# Run loop
for (i in 1:length(un.id)) {
  
  # Filter by ID
  temp <- pigs %>% filter(key==un.id[i])
  

  # Run rsf (but maybe change this to LASSO)
  rsf <- glm(case ~ bottomlandhw + decidmixed + foodcrops, data=temp, family=binomial(link = "logit"), weight=weight)
  summary(rsf)
  
  set.seed(1)
  grid <- 10^seq(10, -2, length=100)
  
  x <- model.matrix(case ~ bottomlandhw + decidmixed + water + foodcrops + othercrops + barren + gramanoids, temp)[, -1]
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
  
  # add model object to list
  pigs_rsf[[i]] <- rsf
  
}


ggplot(temp, aes(x=foodcrops, y=case)) +
  geom_point()










