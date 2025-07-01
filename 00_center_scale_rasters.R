library(terra)

#### Deer ####
# Read in rasters
rast_list <- c("data/landscape_data/shrublands_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif",
               "data/landscape_data/developed_180m_sum.tif",
               "data/landscape_data/palatable_crops_180m_sum.tif",
               "data/landscape_data/allhardwoods_180m_sum.tif")
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to %
layers <- layers / 113

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

w2 <- water^2

layers <- c(layers, w2)

# Rename layers
names(layers) <- c("shrubs", "gramanoids", "developed", "foodcrops", "hardwoods", "water_dist", "water2")

# Get mean and sd from entire covariate rasters
mean_vals <- global(layers, "mean", na.rm = TRUE)
sd_vals   <- global(layers, "sd", na.rm = TRUE)

rast_cs <- dplyr::bind_cols(mean_vals, sd_vals)
rast_cs <- rownames_to_column(rast_cs, var="layer")
write_csv(rast_cs, "output/deer_raster_mean_sds.csv")

# Apply to covariate rasters
layers <- (layers - mean_vals$mean) / sd_vals$sd

# Save deer center & scaled
# Define the output file names. 
lnames <- c("shrubs", "gramanoids", "developed","foodcrops", "hardwoods", "water_dist")
output_files <- paste0("data/landscape_data/scaled_rasters/deer_scaled_", lnames, ".tif")

# # Write each layer to a separate file
writeRaster(layers, output_files, overwrite=TRUE)

#### Pigs ####
# Read in rasters
rast_list <- c("data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/developed_210m_sum.tif",
               "data/landscape_data/allhardwoods_210m_sum.tif") 
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

#
w2 <- water^2
layers <- c(layers, w2)

# Rename layers
names(layers) <- c("shrubs", "gramanoids", "developed", "hardwoods", "dist_water", "water2")

# Get mean and sd from entire covariate rasters
mean_vals <- global(layers, "mean", na.rm = TRUE)
sd_vals   <- global(layers, "sd", na.rm = TRUE)

rast_cs <- dplyr::bind_cols(mean_vals, sd_vals)
rast_cs <- rownames_to_column(rast_cs, var="layer")
write_csv(rast_cs, "output/pigs_raster_mean_sds.csv")


# Apply to covariate rasters
layers <- (layers - mean_vals$mean) / sd_vals$sd

# Save deer center & scaled
# Define the output file names. 
lnames <- c("shrubs", "gramanoids", "developed", "hardwoods", "dist_water")
output_files <- paste0("data/landscape_data/scaled_rasters/pigs_scaled_", lnames, ".tif")

# # Write each layer to a separate file
writeRaster(layers, output_files, overwrite=TRUE)

