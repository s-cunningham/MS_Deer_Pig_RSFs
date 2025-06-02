# ----------------------------------------------------------------------------------------------------------------------------
# Name: 
# Author: Stephanie Cunningham
# Objective: Plot functional response
# Input: RSF predictions (30x30)
# Output: 
# Required packages: tidyverse
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(exactextractr)

## Read in deer data
deer <- read_csv("output/deer_used_avail_covariates.csv")

# Center and scale covariates
deer[,7:19] <- scale(deer[,7:19])

## Set up to loop by individual
un.id <- unique(deer$key)

#### AKDE polygons and CDL #### 
## read cropland data layer
cdl <- rast("data/landscape_data/CDL2023_MS50km_30m_nlcdproj.tif")
names(cdl) <- "value"

# Reclassify so that there is only 1 id for each class
m <- cbind(c(141,190,1,3,5,26,241,225,12,152,4,23,24,26,27,28,29,36,37,176,205,236,121,122,123,124),
           c(1,1,2,2,2,2,2,2,2,3,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5))
cdl <- classify(cdl, m)
remv <- unique(cdl)
remv <- remv[remv > 5]
remv <- c(0, remv)
fillNAs <- rep(NA, length(remv))
m <- cbind(remv,
           fillNAs)
cdl <- classify(cdl, m)
unique(cdl)

# Set up NLCD classes and class values
cdl_values <- c(1,2,3,4,5)
cdl_class <- c("Hardwoods","Crops","Shrubs","Gramanoids","Developed")

# Add class names and numbers to the raster
levels(cdl) <- list(data.frame(value = cdl_values,
                                landcov = cdl_class))


# Read in polygons of AKDES
ud <- vect("output/deer_AKDE_est.shp") 
# check projections
crs(ud)==crs(cdl)
# Reproject
ud <- project(ud, cdl)
# Convert SpatVector to sf
ud <- sf::st_as_sf(ud) %>%
  unite("key", c("id", "burst"), sep="_", remove=FALSE)

## Extract values from NLCD raster based on buffer using exactextractr
landcov_fracs <- exact_extract(cdl, ud, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(key, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'key', progress = FALSE)

# Clean up data frame a little
landcov_fracs <- landcov_fracs %>% 
  # drop NAs
  filter(!is.na(value)) %>%
  # replace value column with covariate names
  mutate(value=as.character(value),
         value=case_when(value==1 ~ "Hardwoods",
                         value==2 ~ "Crops",
                         value==3 ~ "Shrubs",
                         value==4 ~ "Gramanoids",
                         value==5 ~ "Developed"))

#### Coefficients ####
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
  dplyr::select(allhardwoods:developed) %>%
  # fill with 0
  mutate(across(allhardwoods:developed, \(x) coalesce(x, 0))) %>%
  # rename to match landcover data
  rename(Hardwoods=allhardwoods, Gramanoids=gramanoids, Shrubs=shrubs, Crops=foodcrops, Developed=developed)

# Add IDs 
r_glms$key <- un.id

# Pivot longer
r_glms <- r_glms %>%
  pivot_longer(1:5, names_to="value", values_to="beta")

## Join pct cover and glm betas
dat <- left_join(r_glms, landcov_fracs, by=c("key", "value"))

#### Plot functional response
ggplot(dat) +
  geom_hline(yintercept=0) +
  geom_point(aes(x=freq, y=beta)) +
  facet_wrap(vars(value), scales="free", ncol=5) +
  theme_classic()

