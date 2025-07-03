library(tidyverse)
library(terra)
library(tidyterra)

# Load prediction maps (0-1, 90m)
deer <- rast("results/predictions/deer_rsf_predicted_90m.tif")
pigs <- rast("results/predictions/pigs_rsf_predicted_90m.tif")

# Load patch rasters
d_core <- rast("results/exported_patches/deer_core.tif")
p_core <- rast("results/exported_patches/pigs_core.tif")

d_marg <- rast("results/exported_patches/deer_marginal.tif")
p_marg <- rast("results/exported_patches/pigs_marginal.tif")

d_very <- rast("results/exported_patches/deer_highlymarginal.tif")
p_very <- rast("results/exported_patches/pigs_highlymarginal.tif")


# Mask rasters
dc_mask <- mask(deer, d_core)
pc_mask <- mask(pigs, p_core)

dm_mask <- mask(deer, d_marg)
pm_mask <- mask(pigs, p_marg)

dhm_mask <- mask(deer, d_very)
phm_mask <- mask(pigs, p_very)


## Summarize values
global(dc_mask, "mean", na.rm=TRUE)
global(pc_mask, "mean", na.rm=TRUE)

global(dm_mask, "mean", na.rm=TRUE)
global(pm_mask, "mean", na.rm=TRUE)

global(dhm_mask, "mean", na.rm=TRUE)
global(phm_mask, "mean", na.rm=TRUE)
