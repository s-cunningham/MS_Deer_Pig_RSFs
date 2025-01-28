library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(ggspatial)

# Read in MS shapefile
ms <- vect("data/landscape_data/mississippi_ACEA.shp")

# List patch .asc files
patches <- list.files(path="C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/patches/",
                      pattern=".ASC", full.names=TRUE)  # note that the pattern argument is case-sensitive, and RAMAS writes ASCII files as .ASC instead of .asc like R
# Make raster "stack"
patches <- rast(patches)

# Define projection
crs(patches) <- crs(ms)

# Shorten the names a bit (taking of "_patches")
names(patches) <- gsub("_patches", "", names(patches))

# Reclassify so that 0 becomes NA (backround is 0, patches are numbered starting with 1), and all patches receive same class value (1)
# If you want to show different patches, don't do the second part
# Replace 0 with NA
patches <- ifel(patches == 0, NA_integer_, patches)

# Reclassify all patches as 1
patches <- ifel(patches>=1, 1, patches)

# loop over patches and write as .tif
for (i in 1:dim(patches)[3]) {
  # Subset to each layer
  patch <- patches[[i]]
  
  # create filename
  filename <- paste0("results/habitat_patches/", names(patch), ".tif")
  
  # write raster
  writeRaster(patch, filename)
}

## Calculate area (in hectares)
p_area <- expanse(patches, unit="ha", transform=TRUE)

# convert ha to km2, then add the layer names to the layer column
p_area <- p_area %>% mutate(areakm2 = area/100) %>%
              mutate(layer=names(patches))

# Convert rasters to polygones
polylist <- lapply(as.list(patches), as.polygons)





