library(tidyverse)
library(terra)
library(tidyterra)
library(ggspatial)
library(patchwork)

theme_set(theme_void())

# Read in SE rasters (clamped at 10)
deer <- rast("results/predictions/deer_rsf_se_30m.tif")
pigs <- rast("results/predictions/pigs_rsf_se_30m.tif")

# Read in MS shapefile (to drop islands)
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
ms <- project(ms, pigs)

# Load template raster (cells must be *exactly* the same width/length for RAMAS)
temp_rast <- rast("data/landscape_data/CDL2023_90mACEA_mask.tif")

# Resample
deer <- resample(deer, temp_rast)
deer <- crop(deer, ms)

pigs <- resample(pigs, temp_rast)
pigs <- crop(pigs, ms)


both <- c(deer, pigs)
names(both) <- c('deer', 'pigs')

ggplot() +
  geom_spatraster(data=both) +
  geom_spatvector(data=ms, fill=NA, color="black", linewidth=0.5) +
  facet_wrap(~lyr) +
  scale_fill_grass_c(palette="viridis", direction=1, name="Standard Error") +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1)) +
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
    text_cex = .8) +
  theme(legend.position="bottom",
        legend.title=element_text(vjust=0.9, size=12),
        legend.text=element_text(size=11),
        strip.text=element_blank())


# save figure
ggsave("figs/mapping_standard_error.svg")
# Saving 8.91 x 8.05 in image




# Flag pixels with covariates outside training range
extrap_mask <- sum(abs(cov_stack) > 3, na.rm = TRUE) > 0
plot(extrap_mask, main = "Covariate Extrapolation Areas")