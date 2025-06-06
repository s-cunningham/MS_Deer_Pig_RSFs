library(tidyverse)
library(terra)
library(tidyterra)
library(ggspatial)
library(patchwork)

theme_set(theme_void())

# Read MS shapefile
ms <- vect("data/landscape_data/mississippi_ACEA.shp")

# Read in rasters
deer <- rast("results/predictions/deer_rsf_predicted_30m.tif")
pigs <- rast("results/predictions/pigs_rsf_predicted_30m.tif")

# Binary classification for selection
m <- c(minmax(deer)[1], 0, 0,
       0, minmax(deer)[2], 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
deer <- classify(deer, rclmat, include.lowest=TRUE)

m <- c(minmax(pigs)[1], 0, 0,
       0, minmax(pigs)[2], 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
pigs <- classify(pigs, rclmat, include.lowest=TRUE)

# Convert to categorical
deer <- as.factor(deer)
pigs <- as.factor(pigs)

# Reproject vector
ms <- project(ms, pigs)


d_plot <- ggplot() +
  geom_spatraster(data=deer) +
  geom_spatvector(data=ms, fill=NA, color="black") +
  scale_fill_grass_d(palette = "viridis", na.translate=FALSE) +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.02, "npc"),
    style = north_arrow_orienteering())

p_plot <- ggplot() +
  geom_spatraster(data=pigs) +
  geom_spatvector(data=ms, fill=NA, color="black") +
  scale_fill_grass_d(palette = "viridis", na.translate=FALSE) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.07, "npc"),
    text_cex = .8) 


d_plot + p_plot + plot_layout(guides = "collect") & theme(legend.position="bottom")

ggsave("figs/new_figure4.svg")
