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
                         value==4 ~ "Graminoids",
                         value==5 ~ "Developed"))

# Average distance to water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, cdl)
ext(water) <- ext(cdl)

# Reclassify missing data to 0
m <- rbind(c(0, NA))
water <- classify(water, m)

# Average distance to water in home range
mean_dist <- exact_extract(water, ud, 'mean')

# Add to data frame with individual ID
mean_dist <- data.frame(key=ud$key, Water=mean_dist, Water2=mean_dist^2)

# pivot so it can have rows added to landcover_fracs
mean_dist <- mean_dist %>%
  as_tibble() %>%
  pivot_longer(cols=2:3, names_to="value", values_to="freq")

# Bind rows to landcov_fracs
landcov_fracs <- bind_rows(landcov_fracs, mean_dist)

#### Coefficients ####
# read in models
deer_rsf <- readRDS("output/deer_glms.RDS")
deerSD <- read_csv("output/deer_used_avail_covariates.csv") %>%
  # Select just the variables from the model
  select("shrubs", "gramanoids", "foodcrops", "developed", "allhardwoods", "water_dist") %>%
  # Create squared water
  mutate(Water2=water_dist^2) %>%
  # rename to match landcover data
  rename(Hardwoods=allhardwoods, Graminoids=gramanoids, Shrubs=shrubs, Crops=foodcrops, Developed=developed,
         Water=water_dist) %>%
  # reorder to match glms
  select(Hardwoods, Graminoids, Shrubs, Crops, Developed, Water, Water2)

# Center and scale
deerSD <- scale(deerSD)

deerSD <- attr(deerSD, "scaled:scale")

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
  dplyr::select(allhardwoods:water_dist2) %>%
  # fill with 0
  mutate(across(allhardwoods:water_dist2, \(x) coalesce(x, 0))) %>%
  # rename to match landcover data
  rename(Hardwoods=allhardwoods, Graminoids=gramanoids, Shrubs=shrubs, Crops=foodcrops, Developed=developed,
         Water=water_dist, Water2=water_dist2)

# Unscale
r_glms <- r_glms %>%
  mutate(Hardwoods=Hardwoods / deerSD[1],
         Graminoids=Graminoids / deerSD[2],
         Shrubs=Shrubs / deerSD[3],
         Crops=Crops / deerSD[4],
         Developed=Developed / deerSD[5],
         Water=Water / deerSD[6],
         Water2=Water2 / deerSD[7])

# Add IDs 
r_glms$key <- un.id

# Pivot longer
r_glms <- r_glms %>%
  pivot_longer(1:7, names_to="value", values_to="beta")

## Join pct cover and glm betas
deer_fr <- left_join(r_glms, landcov_fracs, by=c("key", "value"))

# Fill in 0s
deer_fr <- deer_fr %>% mutate(freq=coalesce(freq, 0))

# Save CSV
write_csv(deer_fr, "output/deer_functional.csv")

#### Plot functional response
ggplot(deer_fr) +
  geom_hline(yintercept=0, linewidth=0.5, color="gray60") +
  geom_point(aes(x=freq, y=beta), alpha=0.5) +
  facet_wrap(vars(value), scales="free", ncol=7, strip.position="bottom") +
  labs(y = expression(beta)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0, vjust = 0.5, size=16),
        axis.text=element_text(size=10),
        strip.placement="outside",
        strip.background=element_rect(color=NA, fill=NA))

