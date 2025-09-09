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

# Calculate area of regions
reg_area <- data.frame(regions$name, area_km2=st_area(regions)/1000/1000)
reg_area <- reg_area |>
  mutate(area_km2 =as.numeric(area_km2)) |>
  rename(name=regions.name)

regions <- st_transform(regions, crs=st_crs(deer_area))

# Calcuate percentages
## deer
suitability_deer <- exact_extract(deer_area, regions, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)

suitability_deer <- suitability_deer |>
  ungroup() |>
  # drop NAs 
  filter(!is.na(value)) |>
  # pivot
  pivot_wider(names_from = "value", values_from = "freq") |>
  rename(marginal=`1`, moderate=`2`, core=`3`) 
  
suitability_deer <- left_join(suitability_deer, reg_area, by="name")

suitability_deer <- suitability_deer |>
  mutate(core_km2 = core*area_km2,
         moderate_km2 = moderate*area_km2,
         marginal_km2 = marginal*area_km2)

## pigs
suitability_pigs <- exact_extract(pigs_area, regions, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)

suitability_pigs <- suitability_pigs |>
  ungroup() |>
  # drop NAs 
  filter(!is.na(value)) |>
  # pivot
  pivot_wider(names_from = "value", values_from = "freq") |>
  rename(marginal=`1`, moderate=`2`, core=`3`) 

suitability_pigs <- left_join(suitability_pigs, reg_area, by="name")

suitability_pigs <- suitability_pigs |>
  mutate(core_km2 = core*area_km2,
         moderate_km2 = moderate*area_km2,
         marginal_km2 = marginal*area_km2)


# area-adjusted total suitability?


d_totRSF <- left_join(d_totRSF, reg_area, by="name")
d_totRSF <- d_totRSF |>
  mutate(a_adj_totalHSI = sum/area_km2)
