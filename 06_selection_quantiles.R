library(tidyverse)
library(terra)
library(tidyterra)


# Load template raster (cells must be *exactly* the same width/length for RAMAS)
temp_rast <- rast("data/landscape_data/CDL2023_90mACEA_mask.tif")

# Load 30 m rasters
pigs30 <- rast("results/predictions/pigs_rsf_predicted_30m.tif")
deer30 <- rast("results/predictions/deer_rsf_predicted_30m.tif")

# Resample
deer90 <- resample(deer30, temp_rast)
deer90 <- crop(deer90, deer30)

pigs90 <- resample(pigs30, temp_rast)
pigs90 <- crop(pigs90, pigs30)

# Save only cells > 0
vd90 <- values(deer90, mat = TRUE)
vd90 <- vd90 %>% na.omit()
vd90 <- vd90 %>% as_tibble() %>%
  filter(RSF > 0)

quantile(vd90$RSF, probs=c(0.5))


vp90 <- values(pigs90, mat = TRUE)
vp90 <- vp90 %>% na.omit()
vp90 <- vp90 %>% as_tibble() %>%
  filter(RSF > 0)

quantile(vp90$RSF, probs=c(0.5))


# linear stretch
# deer
(1.429825 - minmax(deer90)[1])/(minmax(deer90)[2]-minmax(deer90)[1])

# pigs
(1.877943 - minmax(pigs90)[1])/(minmax(pigs90)[2]-minmax(pigs90)[1])


## full
pigs90ls <- (pigs90 - minmax(pigs90)[1])/(minmax(pigs90)[2]-minmax(pigs90)[1])
deer90ls <- (deer90 - minmax(deer90)[1])/(minmax(deer90)[2]-minmax(deer90)[1])

writeRaster(pigs90ls, "results/predictions/pigs_rsf_predicted.asc", NAflag=-9999, overwrite=TRUE)
##
deer90 <- rast("results/predictions/deer_rsf_predicted_90m.tif")
global(deer90ls, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)

pigs90 <- rast("results/predictions/pigs_rsf_predicted_90m.tif")
global(pigs90ls, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)
