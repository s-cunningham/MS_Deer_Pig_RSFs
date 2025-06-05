library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(ggspatial)

# set ggplot theme
theme_set(theme_classic())

## Read predicted maps (90 x 90)
deer <- rast("results/predictions/deer_rsf_bins_30m.tif")
pigs <- rast("results/predictions/pigs_rsf_bins_30m.tif")

# Convert to data frame
deer <- as.data.frame(deer, cells=TRUE)
pigs <- as.data.frame(pigs, cells=TRUE)

# Rename the suitability column
names(deer)[2] <- "deer" 
names(pigs)[2] <- "pigs"

# Combine
dat <- left_join(deer, pigs, by="cell")
dat$deer <- factor(dat$deer, levels=1:10, labels=c("binA","binB","binC","binD","binE","binF","binG","binH","binI","binJ"))
dat$pigs <- factor(dat$pigs, levels=1:10, labels=c("binA","binB","binC","binD","binE","binF","binG","binH","binI","binJ"))

# Create "confusion matrix"
cm <- table(Deer = as.vector(dat$deer), Pigs = as.vector(dat$pigs))

# Get upper triangle
cm[upper.tri(cm)]
sum(cm[upper.tri(cm)])/nrow(dat)

# Diagonal
sum(diag(cm))/nrow(dat)

# get lower triange
cm[lower.tri(cm)]
sum(cm[lower.tri(cm)])/nrow(dat)




