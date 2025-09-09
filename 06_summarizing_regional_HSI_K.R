library(tidyverse)
library(terra)
library(sf)
library(exactextractr)

# Read region polygon
regions <- st_read("data/landscape_data/edited_ms_soil_regions.shp")

# Read suitability level raster
pigs_area <- rast("data/landscape_data/pigs_suitability_levels.tif")
deer_area <- rast("data/landscape_data/deer_suitability_levels.tif")

# Read RSF predictions (habitat suitability)
pigs_rsf <- rast("results/predictions/pigs_rsf_predicted_90m.tif")
deer_rsf <- rast("results/predictions/deer_rsf_predicted_90m.tif")

# reproject regions to raster CRS
regions <- st_transform(regions, crs=st_crs(deer_rsf))

# Extract RSF value
p_avgRSF <- exact_extract(pigs_rsf, regions, 'mean', append_cols = 'name', progress = FALSE)
d_avgRSF <- exact_extract(deer_rsf, regions, 'mean', append_cols = 'name', progress = FALSE)

p_totRSF <- exact_extract(pigs_rsf, regions, 'sum', append_cols = 'name', progress = FALSE)
d_totRSF <- exact_extract(deer_rsf, regions, 'sum', append_cols = 'name', progress = FALSE)

p_totRSF$K <- p_totRSF$sum * 0.219
d_totRSF$K <- d_totRSF$sum * 0.116

# round to 2 digits
d_avgRSF$mean <- round(d_avgRSF$mean, digits=2) 
p_avgRSF$mean <- round(p_avgRSF$mean, digits=2) 

# calculate fraction of area suitable for pigs
regions <- st_transform(regions, crs=st_crs(deer_area))
suitability_deer <- exact_extract(deer_area, regions, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)


expanse(deer_area, unit="km")
