library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(ggspatial)

# Misissippi polygon
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
lakes <- vect("data/landscape_data/ne_ms_lakes.shp")
lakes <- project(lakes, ms)

# RSF prediction rasters
pigsRSF <- rast("data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.tif")
deerRSF <- rast("data/predictions/deer_glmm_rsf_FINALFINALFINALFINAL.tif")

# Mask to Mississippi
pigsRSF <- mask(pigsRSF, ms)
deerRSF <- mask(deerRSF, ms)

# Mask lakes
pigsRSF <- mask(pigsRSF, lakes, inverse=TRUE)
deerRSF <- mask(deerRSF, lakes, inverse=TRUE)

pigs_rsf <- ggplot() +
  geom_spatraster(data=pigsRSF) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability") +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  theme(legend.position="none")

deer_rsf <- ggplot() +
  geom_spatraster(data=deerRSF) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability") +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  theme(legend.position = 'bottom',
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        # legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm")) +
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
  plot_annotation(tag_levels = 'A', tag_prefix="(", tag_suffix=")") 

# ggsave(file="figs/both_rsf_maps.svg")
# Saving 7.43 x 6.7 in image
# Saving 10.1 x 9.06 in image (2025-01-24)
# Saving 10.6 x 9.06 in image (2025-01-28)

## Load confidence intervals
deer_lci <- rast("data/predictions/deer_glmm_rsf_LowerCI.asc")
deer_uci <- rast("data/predictions/deer_glmm_rsf_UpperCI.asc")
pigs_lci <- rast("data/predictions/pigs_glmm_rsf_LowerCI.asc")
pigs_uci <- rast("data/predictions/pigs_glmm_rsf_UpperCI.asc")

# Mask to Mississippi
deer_lci <- mask(deer_lci, ms)
deer_uci <- mask(deer_uci, ms)
pigs_lci <- mask(pigs_lci, ms)
pigs_uci <- mask(pigs_uci, ms)

# Mask lakes
deer_lci <- mask(deer_lci, lakes, inverse=TRUE)
deer_uci <- mask(deer_uci, lakes, inverse=TRUE)
pigs_lci <- mask(pigs_lci, lakes, inverse=TRUE)
pigs_uci <- mask(pigs_uci, lakes, inverse=TRUE)

deerLCI <- ggplot() +
  geom_spatraster(data=deer_lci) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability", limits=c(0, 0.9997315)) +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  theme(legend.position = 'bottom',
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        # legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm")) +
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

deerUCI <- ggplot() +
  geom_spatraster(data=deer_uci) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability", limits=c(0, 0.9997315)) +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  theme(legend.position = 'none',
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))
        # legend.key.height = unit(1, "cm"),
        # legend.key.width = unit(2, "cm"))

pigsLCI <- ggplot() +
  geom_spatraster(data=pigs_lci) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability", limits=c(0, 0.863947)) +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  theme(legend.position = 'bottom',
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        # legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2, "cm"))

pigsUCI <- ggplot() +
  geom_spatraster(data=pigs_uci) +
  scale_fill_viridis_c(option="C", direction=-1, na.value="transparent", name="Suitability", limits=c(0, 0.863947)) +
  geom_spatvector(data=ms, color="black", fill=NA, linewidth=0.6) +
  theme_void() +
  theme(legend.position = 'none',
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))
        # legend.key.height = unit(1, "cm"),
        # legend.key.width = unit(2, "cm"))



(deerLCI + deerUCI) / (pigsLCI + pigsUCI) + #plot_layout(guides="collect") + 
  plot_annotation(tag_levels = 'a', tag_prefix="(", tag_suffix=")") #&
  # theme(legend.position="bottom")
ggsave("figs/rsf_ci_plot.svg", width=5.7, height=9)




## Removing Core and marginal from HSI rasters
## subtract core and marginal area from RSF
## Deer Core
m <- rbind(c(0.70,1,0))
deerNoCore <- classify(deerRSF, m)
plot(deerNoCore)

writeRaster(deerNoCore, "data/predictions/deer_glmm_rsf_FINAL_noCore_mean.asc", NAflag=-9999, overwrite=TRUE)

deerNoCoreLCI <- classify(deer_lci, m)
deerNoCoreUCI <- classify(deer_uci, m)



## Deer Marginal
m <- rbind(c(0.35,1,0))
deerNoMarg <- classify(deerRSF, m)
plot(deerNoMarg)
writeRaster(deerNoMarg, "data/predictions/deer_glmm_rsf_FINAL_noCoreMarg.asc", NAflag=-9999, overwrite=TRUE)

m <- rbind(c(0.35,1,0))
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

