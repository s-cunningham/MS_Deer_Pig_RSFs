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

# Loop over patches and write as .tif
for (i in 1:dim(patches)[3]) {
  # Subset to each layer
  patch <- patches[[i]]
  
  # create filename
  filename <- paste0("results/habitat_patches/", names(patch), ".tif")
  
  # write raster
  writeRaster(patch, filename)
}

## Calculate area (in hectares)
p_area <- expanse(patches, unit="km", transform=TRUE)

# Add the layer names to the layer column
p_area <- p_area %>% mutate(layer=names(patches))

# Convert rasters to polygones
polylist <- lapply(as.list(patches), as.polygons)

# Make a named list so layers can by extracted by name
names(polylist) <- names(patches)

## Calculate area of overlap between Deer core patches and pigs
# Need to add terra:: before the union function, it must be masked by something else

# Core deer & core pig habitat overlap 
ovp1 <- terra::union(polylist$deer_c_mean, polylist$pig_c_mean)
expanse(ovp1, unit="km")[[3]]

# Core deer & core pig habitat overlap - upper 95% CI
ovp2 <- terra::union(polylist$deer_c_uci, polylist$pig_c_uci)
expanse(ovp2, unit="km")[[3]]

# Core deer & core pig habitat overlap - lower 95% CI
ovp3 <- terra::union(polylist$deer_c_lci, polylist$pig_c_lci)
expanse(ovp3, unit="km")[[3]]

# Core deer & marginal pig habitat overlap
ovp4 <- terra::union(polylist$deer_c_mean, polylist$pig_m_mean)
expanse(ovp4, unit="km")[[3]]

ovp5 <- terra::union(polylist$deer_c_lci, polylist$pig_m_lci)
expanse(ovp5, unit="km")[[3]]

ovp6 <- terra::union(polylist$deer_c_uci, polylist$pig_m_uci)
expanse(ovp6, unit="km")[[3]]





#### Create patch plot for manuscript ####
## For each species, combine into multiclass raster (mean predictions)
# Will use the patch values to reclassify
patch_vals <- c(3,2,1)
patch_class <- c("Core", "Marginal", "Highly Marginal")

## Mosaic rasters
deer_mn <- mosaic(patches["deer_hm_mean"], patches["deer_m_mean"], fun = "sum")
deer_mn <- mosaic(deer_mn, patches["deer_c_mean"], fun="sum")

pigs_mn <- mosaic(patches["pig_hm_mean"], patches["pig_m_mean"], fun = "sum")
pigs_mn <- mosaic(pigs_mn, patches["pig_c_mean"], fun="sum")

## Reclassify
# Add patch types to raster
levels(deer_mn) <- list(data.frame(ID = patch_vals,
                                patch = patch_class))

levels(pigs_mn) <- list(data.frame(ID = patch_vals,
                                   patch = patch_class))


## Plot
deer_patches <- ggplot() +
  geom_spatraster(data=deer_mn) +
  scale_fill_viridis_d(option="C", direction=1, na.value="transparent", name="Patch Type:", na.translate=FALSE) +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.85, "npc"),
    style = north_arrow_minimal()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.07, "npc"),
    text_cex = .8) 

pig_patches <- ggplot() +
  geom_spatraster(data=pigs_mn) +
  scale_fill_viridis_d(option="C", direction=1, name="Patch Type:", na.value="transparent", na.translate=FALSE) +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() 

deer_patches + pig_patches +
  plot_annotation(tag_levels = 'A', tag_prefix="(", tag_suffix=")") +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', 
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))

ggsave(file="figs/patches.svg")
# Saving 10.6 x 9.06 in image
