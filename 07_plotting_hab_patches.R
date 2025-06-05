library(tidyverse)
library(terra)
library(tidyterra)

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

# area/cell frequency
fdeer <- freq(deer)
fdeer$area <- fdeer$count * 0.0009

fpigs <- freq(pigs)
fpigs$area <- fpigs$count * 0.0009





