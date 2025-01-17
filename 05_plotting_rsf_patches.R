library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(ggspatial)

# Misissippi polygon
ms <- vect("data/landscape_data/mississippi_ACEA.shp")

# RSF prediction rasters
pigsRSF <- rast("data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.tif")
deerRSF <- rast("data/predictions/deer_glmm_rsf_FINALFINALFINALFINAL.tif")

# Mask to Mississippi
pigsRSF <- mask(pigsRSF, ms)
deerRSF <- mask(deerRSF, ms)

pigs_rsf <- ggplot() +
  geom_spatraster(data=pigsRSF) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability") +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void()

deer_rsf <- ggplot() +
  geom_spatraster(data=deerRSF) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability") +
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

deer_rsf + pigs_rsf +
  plot_annotation(tag_levels = 'A', tag_prefix="(", tag_suffix=")") +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', 
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        # legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"))

ggsave(file="figs/both_rsf_maps.svg")
# Saving 7.43 x 6.7 in image

## Plotting habitat patches
# Core habitat patches
p_core <- rast("data/habitat_patches/pigs_core_patches.tif")
d_core <- rast("data/habitat_patches/deer_core_patches.tif")

# Marginal habitat patches
p_marg <- rast("data/habitat_patches/pigs_marginal_patches.tif")
d_marg <- rast("data/habitat_patches/deer_marginal_patches.tif")

# highly marginal habitat patches
p_hmar <- rast("data/habitat_patches/pigs_highly_marginal_patches.tif")
d_hmar <- rast("data/habitat_patches/deer_highly_marginal_patches.tif")

## combine patch rasters into single, 3-class raster
deer <- mosaic(d_hmar, d_marg, fun = "sum")
deer <- mosaic(deer, d_core, fun="sum")
# Reclassify
patch_vals <- c(3,2,1)
patch_class <- c("Core", "Marginal", "Highly Marginal")

# Add patch types to raster
levels(deer) <- list(data.frame(ID = patch_vals,
                                patch = patch_class ))

pigs <- mosaic(p_hmar, p_marg, fun = "sum")
pigs <- mosaic(pigs, p_core, fun="sum")
# Add class names and numbers to the raster
levels(pigs) <- list(data.frame(ID = patch_vals,
                                patch = patch_class ))

## Plot
deer_patches <- ggplot() +
  geom_spatraster(data=deer) +
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
  geom_spatraster(data=pigs) +
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
# Saving 7.43 x 6.7 in image

## subtract core and marginal area from RSF
## Deer
m <- rbind(c(0.73,1,0))
deerNoCore <- classify(deerRSF, m)
plot(deerNoCore)

writeRaster(deerNoCore, "data/predictions/deer_glmm_rsf_FINAL_noCore.asc", NAflag=-9999, overwrite=TRUE)

m <- rbind(c(0.365,1,0))
deerNoMarg <- classify(deerRSF, m)
plot(deerNoMarg)
writeRaster(deerNoMarg, "data/predictions/deer_glmm_rsf_FINAL_noCoreMarg.asc", NAflag=-9999, overwrite=TRUE)

## Pigs
m <- rbind(c(0.53,1,0))
pigsNoCore <- classify(pigsRSF, m)
plot(pigsNoCore)

writeRaster(pigsNoCore, "data/predictions/pigs_glmm_rsf_FINAL_noCore.asc", NAflag=-9999, overwrite=TRUE)

m <- rbind(c(0.265,1,0))
pigsNoMarg <- classify(pigsRSF, m)
plot(pigsNoMarg)
writeRaster(pigsNoMarg, "data/predictions/pigs_glmm_rsf_FINAL_noCoreMarg.asc", NAflag=-9999, overwrite=TRUE)

