
library(terra)

# Read in raster
deer_rsf <- rast("results/predictions/deer_rsf_predicted_90m.tif")

# Replace all numerical values with 1
deer_rsf[deer_rsf >= 0] <- 1

# Replace all NAs with 0
r_no_na <- subst(deer_rsf, NA, 0)

# write new binary raster
writeRaster(r_no_na, "data/spatial_data/useable_land.tif")
