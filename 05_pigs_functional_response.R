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

## Read in pigs data
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

# Center and scale covariates
lyrs <- read_csv("output/pigs_raster_mean_sds.csv")

pigs <- pigs %>%
  mutate(hardwoods=(hardwoods-lyrs$mean[4]) / lyrs$sd[4],
         shrubs=(shrubs-lyrs$mean[1]) / lyrs$sd[1],
         gramanoids=(gramanoids-lyrs$mean[2]) / lyrs$sd[2],
         developed=(developed-lyrs$mean[3]) / lyrs$sd[3],
         dist_water=(dist_water-lyrs$mean[5]) / lyrs$sd[5]) 

## Set up to loop by individual
un.id <- unique(pigs$key)

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
ud <- vect("output/pigs_AKDE_est.shp") 
# check projections
crs(ud)==crs(cdl)
# Reproject
ud <- project(ud, cdl)
# Convert SpatVector to sf
ud <- sf::st_as_sf(ud) #%>%
  # unite("key", c("id", "burst"), sep="_", remove=FALSE)

## Extract values from NLCD raster based on buffer using exactextractr
landcov_fracs <- exact_extract(cdl, ud, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(id, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'id', progress = FALSE)

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
                         value==5 ~ "Developed")) %>%
  rename(key=id)

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
mean_dist <- data.frame(key=ud$id, Water=mean_dist, Water2=mean_dist^2)

# pivot so it can have rows added to landcover_fracs
mean_dist <- mean_dist %>%
  as_tibble() %>%
  pivot_longer(cols=2:3, names_to="value", values_to="freq")

# Bind rows to landcov_fracs
landcov_fracs <- bind_rows(landcov_fracs, mean_dist)

#### Coefficients ####
# read in models
pigs_rsf <- readRDS("output/pigs_glms.RDS")

# Center and scale
# extract mean coeficients and combine into a single table
r_glms <- lapply(pigs_rsf, coef)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- lapply(r_glms, t)
r_glms <- lapply(r_glms, as.data.frame)
r_glms <- do.call(bind_rows, r_glms)

r_glms <- r_glms %>% 
  as_tibble() %>%
  rename(intercept=`(Intercept)`, water_dist2=`I(dist_water^2)`) %>%  
  # drop intercept
  dplyr::select(-intercept) %>%
  # Reorder
  dplyr::select(hardwoods:water_dist2) %>%
  # fill with 0
  mutate(across(hardwoods:water_dist2, \(x) coalesce(x, 0))) %>%
  # rename to match landcover data
  rename(Hardwoods=hardwoods, Graminoids=gramanoids, Shrubs=shrubs, Developed=developed,
         Water=dist_water, Water2=water_dist2)

# Unscale
r_glms <- r_glms %>%
  mutate(Hardwoods=Hardwoods / lyrs$sd[4],
         Shrubs=Shrubs / lyrs$sd[1],
         Graminoids=Graminoids / lyrs$sd[2],
         Developed=Developed / lyrs$sd[3],
         Water=Water / lyrs$sd[5],
         Water2=Water2 / lyrs$sd[6])

# Add IDs 
r_glms$key <- un.id[-c(13,18,28)]

# Pivot longer
r_glms <- r_glms %>%
  pivot_longer(1:6, names_to="value", values_to="beta")

## Join pct cover and glm betas
pigs_fr <- left_join(r_glms, landcov_fracs, by=c("key", "value"))

# Fill in 0s
pigs_fr <- pigs_fr %>% mutate(freq=coalesce(freq, 0))

# Save CSV
write_csv(pigs_fr, "output/pigs_functional.csv")


#### Plot functional response
ggplot(pigs_fr) +
  geom_hline(yintercept=0) +
  geom_point(aes(x=freq, y=beta)) +
  xlab("% landcover in home range") +
  facet_wrap(vars(value), scales="free", ncol=6, strip.position="bottom") +
  labs(y = expression(beta)) +
  theme_classic()+
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.y=element_text(angle = 0, vjust = 0.5, size=16),
        axis.text=element_text(size=10),
        strip.placement="outside",
        strip.background=element_rect(color=NA, fill=NA))


