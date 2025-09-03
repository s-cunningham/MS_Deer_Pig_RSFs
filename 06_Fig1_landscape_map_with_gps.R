library(tidyverse)
library(terra)
library(tidyterra)
library(ggspatial)

## Read region map
regions <- vect("data/landscape_data/ms_soil_ecoregions.shp")

# add region names
regions <- regions |>
  mutate(name=c("Delta", "Interior Flatwoods", "Upper Coastal Plain", "Blackland Prairie", "Coastal Flatwoods",
                "Lower Coastal Plain", "Thin Loess", "Thick Loess"))

# Check the names are correct
ggplot() +
  geom_spatvector(data=regions, aes(fill=name))

## Read GPS locations
# Pigs
pigs <- read_csv("data/location_data/pigs_filtered.csv")
pigs <- vect(pigs, geom=c("X", "Y"), crs="epsg:32616")

# Deer
deer <- read_csv("data/location_data/deer_filtered.csv")
deer <- vect(deer, geom=c("X", "Y"), crs="epsg:32616")

## Create 1km buffers
# Pigs
pigs <- buffer(pigs, width=1000)
pigs <- aggregate(pigs) # Dissolve buffers

# Deer
deer <- buffer(deer, width=1000)
deer <- aggregate(deer) # Dissolve buffers

## Read spatial data
nlcd <- rast("data/landscape_data/NLCD2023_MS50kmbuff_90m.tif")
ms <- vect("data/landscape_data/mississippi_ACEA.shp")

## Create empty raster
temp_rast <- rast(ext(ms), resolution=1000) # create template raster 1 km x 1 km

# Mask to Mississippi
nlcd <- mask(nlcd, ms)

# Add color table for plotting
FedData::nlcd_colors()
coltab(nlcd) <- data.frame(value=c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95), 
                           col=c("#5475A8","#E8D1D1","#E29E8C","#ff0000","#B50000",
                                                   "#D2CDC0","#85C77E","#38814E","#D4E7B0","#DCCA8F",
                                                   "#FDE9AA","#FBF65D","#CA9146","#C8E6F8","#64B3D5"))

# Set up NLCD classes and class values
nlcd_values <- c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95)
nlcd_class <- c("Open Water", "Developed, Open Space", "Developed, Low Intensity", 
                "Developed, Medium Intensity", "Developed, High Intensity", "Barren", 
                "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                "Grassland/Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                "Emergent Herbaceous Wetlands")

# Add class names and numbers to the raster
levels(nlcd) <- list(data.frame(ID = nlcd_values,
                                landcov = nlcd_class))

# Resample to 1 x 1 km for plot
# nlcd <- resample(nlcd, temp_rast, method="near")
nlcd <- mask(nlcd, ms)
nlcd <- crop(nlcd, ext(ms))

# Reproject polygons
deer <- project(deer, nlcd)
pigs <- project(pigs, nlcd)
regions <- project(regions, nlcd)

map <- ggplot() +
  geom_spatraster(data=nlcd) +
  scale_fill_coltab(data=nlcd, na.translate = FALSE,
                    na.value="transparent", name="NLCD") +
  geom_spatvector(data=regions, col="black", fill=NA, linewidth=0.4) +
  geom_spatvector(data=pigs, col="black", fill="black", alpha=0.8, linewidth=0.7) +
  geom_spatvector(data=deer, col="black", fill="white", alpha=0.8, linewidth=0.7) +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.85, "npc"),
    style = north_arrow_minimal()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.15, "npc"),
    text_cex = .8) +
  theme_bw() +
  theme(panel.border=element_rect(fill=NA, color=NA),
        legend.margin=margin(t = 12, unit='cm')) 

# Create inset
us <- vect("data/landscape_data/north_am_states.shp")
us <- project(us, deer)

inset <- ggplot() +
  ylim(500000,1800000) + xlim(-500000,1500000) +
  geom_spatvector(data=us) +
  geom_spatvector(data=ms, color="gray30", fill="gray30") +
  theme_bw() +
  theme(panel.border=element_rect(fill=NA, color="black"))

## Create 2nd inset (regions)
# get centroid of each polygon
r_centers <- centroids(regions)


ggplot() +
  geom_spatvector(data=regions, aes(fill=name), color=NA) + #linewidth=0.8,
  scale_fill_viridis_d(option="A", name="Region") +
  theme_void()



## Combine with patchwork
library(patchwork)
layout <- c(
  area(t = 1, l = 1, b = 15, r = 4),  # map
  area(t = 0, l = 3, b = 4, r = 4)  # inset?
)

map + inset +
  plot_layout(design = layout)

# ggsave("figs/figure1_map.svg")
# Saving 8.47 x 9.06 in image
