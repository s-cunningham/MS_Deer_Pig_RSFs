library(tidyverse)
library(terra)
library(tidyterra)
library(ggspatial)
library(patchwork)

theme_set(theme_void())

# Read in MS shapefile (to drop islands)
ms <- vect("data/landscape_data/mississippi_ACEA.shp")

# Read in RSF raster
mean_d <- rast("results/predictions/deer_rsf_predicted_30m.tif")
mean_p <- rast("results/predictions/pigs_rsf_predicted_30m.tif")

# Read in LCI
lower_d <- rast("results/predictions/deer_rsf_LCI_30m.tif")
lower_p <- rast("results/predictions/pigs_rsf_LCI_30m.tif")

# Read in UCI
upper_d <- rast("results/predictions/deer_rsf_UCI_30m.tif")
upper_p <- rast("results/predictions/pigs_rsf_UCI_30m.tif")

# project shapefile
ms <- project(ms, mean_d)

# Calculate different between upper and lower to show areas with high variability
d_diff <- upper_d - lower_d
p_diff <- upper_p - lower_p

diffs <- c(d_diff, p_diff)

ggplot() +
  geom_spatraster(data=d_diff) +
  scale_fill_gradient2(low="white", high="red", na.value="transparent")

# Step 1: Aggregate to reduce resolution for plotting
agg_factor <- 30  # adjust this to balance quality and performance
mean_d <- aggregate(mean_d, fact = agg_factor, fun = mean)
lower_d <- aggregate(lower_d, fact = agg_factor, fun = mean)
upper_d <- aggregate(upper_d, fact = agg_factor, fun = mean)
mean_p <- aggregate(mean_p, fact = agg_factor, fun = mean)
lower_p <- aggregate(lower_p, fact = agg_factor, fun = mean)
upper_p <- aggregate(upper_p, fact = agg_factor, fun = mean)

# Combined into single raster stack
all_d <- c(lower_d, mean_d, upper_d)
names(all_d) <- c("LCI", "RSF", "UCI")

all_p <- c(lower_p, mean_p, upper_p)
names(all_p) <- c("LCI", "RSF", "UCI")

global(all_d, quantile, probs=seq(0.1, 1, by=0.1), na.rm=TRUE)
global(all_p, quantile, probs=seq(0.1, 1, by=0.1), na.rm=TRUE)

# Plot all rasters
ggplot() +
  geom_spatraster(data=all_d) +
  scale_fill_grass_c(palette="viridis", direction=1) +
  facet_wrap(~lyr) +
  theme_void()

# Convert pixels to data frame
df_d <- as.data.frame(all_d, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = c(LCI, RSF, UCI), names_to = "type", values_to = "value")

df_p <- as.data.frame(all_p, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = c(LCI, RSF, UCI), names_to = "type", values_to = "value")

# Compute global quantile breaks (based on downsampled data only)
quantiles_d <- quantile(df_d$value, probs = seq(0, 1, 0.1), na.rm = TRUE)
quantiles_p <- quantile(df_p$value, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Bin values by quantiles
df_d$bin <- cut(df_d$value, breaks=quantiles_d, labels=1:10, include.lowest=TRUE)
df_d$bin <- as.numeric(as.character(df_d$bin))

df_p$bin <- cut(df_p$value, breaks=quantiles_p, labels=1:10, include.lowest=TRUE)
df_p$bin <- as.numeric(as.character(df_p$bin))

# Pivot wider and drop values
df_d <- df_d %>% 
  select(-value) %>%
  pivot_wider(names_from="type", values_from="bin")

df_p <- df_p %>% 
  select(-value) %>%
  pivot_wider(names_from="type", values_from="bin")

# Convert the data frame to a multi-layer SpatRaster
ci_bins_d <- rast(df_d, type = "xyz",
                crs = crs(mean_d), # Example CRS
                extent = ext(mean_d))

ci_bins_p <- rast(df_p, type = "xyz",
                  crs = crs(mean_p), # Example CRS
                  extent = ext(mean_p))

# Convert to factor
ci_bins_d <- as.factor(ci_bins_d)
ci_bins_p <- as.factor(ci_bins_p)

# Plot all rasters
deer_plt <- ggplot() +
  geom_spatraster(data=ci_bins_d) +
  scale_fill_grass_d(palette="viridis", direction=1, name="Standard Error") +
  geom_spatvector(data=ms, fill=NA, color="black", linewidth=0.2) +
  facet_wrap(~lyr) +
  theme_void() +
  guides(fill=guide_legend(nrow=1,byrow=TRUE,title="Quantile")) +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        strip.text=element_text(face="bold", size=12)) +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.85, "npc"),
    height= unit(0.08, "npc"),
    style = north_arrow_minimal()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.1, "npc"),
    text_cex = .7)

pig_plt <- ggplot() +
  geom_spatraster(data=ci_bins_p) +
  scale_fill_grass_d(palette="viridis", direction=1, name="Standard Error") +
  geom_spatvector(data=ms, fill=NA, color="black", linewidth=0.2) +
  facet_wrap(~lyr) +
  theme_void() +
  guides(fill=guide_legend(nrow=1,byrow=TRUE,title="Quantile")) +
  theme(legend.position="bottom",
        legend.text=element_blank(),
        strip.text=element_blank()) 


## plot together

deer_plt / pig_plt + plot_layout(guides="collect") & theme(legend.position="bottom")

ggsave(filename="figs/conf_int_maps.svg")
# Saving 8.59 x 8.15 in image


#### Standard Error
deer_se <- rast("results/predictions/deer_rsf_se_30m.tif")
pigs_se <- rast("results/predictions/pigs_rsf_se_30m.tif")

# Rename 
names(deer_se) <- "Deer"
names(pigs_se) <- "Pigs"

# Combine rasters
se_rast <- c(deer_se, pigs_se)

# project shapefile
ms <- project(ms, pigs_se)

deer_se_c <- clamp(deer_se, upper=4, values=TRUE)
pigs_se_c <- clamp(pigs_se, upper=4, values=TRUE)

# Rename 
names(deer_se_c) <- "Deer"
names(pigs_se_c) <- "Pigs"

# Combine clamped rasters
se_rast_c <- c(deer_se_c, pigs_se_c)


ggplot() +
  geom_spatraster(data=se_rast_c) +
  scale_fill_grass_c(palette="viridis", direction=1, name="Standard Error") +
  geom_spatvector(data=ms, fill=NA, color="black", linewidth=0.2) +
  facet_wrap(~lyr) +
  theme(legend.position="bottom",
        strip.text=element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12)) + 
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.85, "npc"),
    height= unit(0.08, "npc"),
    style = north_arrow_minimal()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.1, "npc"),
    text_cex = 0.9)
ggsave("figs/clamped_SE_plot.svg")
# Saving 9.81 x 9 in image

ggplot() +
  geom_spatraster(data=se_rast) +
  scale_fill_grass_c(palette="viridis", direction=1, name="Standard Error") +
  geom_spatvector(data=ms, fill=NA, color="black", linewidth=0.2) +
  facet_wrap(~lyr) +
  theme(legend.position="bottom",
        strip.text=element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12)) + 
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.85, "npc"),
    height= unit(0.08, "npc"),
    style = north_arrow_minimal()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.1, "npc"),
    text_cex = 0.9)
ggsave("figs/SE_plot.svg")

