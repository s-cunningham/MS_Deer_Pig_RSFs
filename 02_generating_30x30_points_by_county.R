library(tidyverse)
library(terra)

# Read in a 30x30 raster as a template
crops <- rast("data/landscape_data/palatable_crops_180m_sum.tif")

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
crops <- classify(crops, m)

# Convert to percent
crops <- crops / 113

# Read shapefile of MS counties and reproject
counties <- vect("data/landscape_data/county_nrcs_a_ms.shp")
counties <- project(counties, crops)

# Create a list where each element is a county
split_counties <- split(counties, "COUNTYNAME") 

# Loop over each county and convert cells to points
for (i in 1:length(split_counties)) {
  
  # Subset just a single county
  cty <- split_counties[[i]]
  
  # Extract county name and create a filename
  n <- unique(cty$COUNTYNAME)
  filename <- paste0("data/landscape_data/county_points/", n, "_points.shp")
  
  # Buffer the county by 5 km
  # ctyB <- buffer(cty, width=5000)
  
  # Mask and crop raster to county
  temp <- mask(crops, cty)
  temp <- crop(temp, cty)
  
  # Convert to points
  pts <- as.points(temp, values=FALSE, na.rm=TRUE, na.all=FALSE)
  
  # Write shapefile
  writeVector(pts, filename, overwrite=TRUE)
  
}





