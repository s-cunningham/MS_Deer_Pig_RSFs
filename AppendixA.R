library(dplyr)
library(ggplot2)
library(terra)
library(tidyterra)

## read in RSF rasters

deer <- rast("data/deer_rsf_predicted.asc")
crs(deer) <- "EPSG:32616"

deer_mar <- ifel(deer>0.34, deer, NA)
deer_mod <- ifel(deer>0.68, deer, NA)
deer_core <- ifel(deer>0.92, deer, NA)


global(deer_mar, fun="sum", na.rm=TRUE)
global(deer_mod, fun="sum", na.rm=TRUE)
global(deer_core, fun="sum", na.rm=TRUE)





## Pigs
pigs <- rast("data/pigs_rsf_predicted.asc")

pig_mar <- ifel(pigs>=0.155, pigs, NA)
pig_mod <- ifel(pigs>0.31, pigs, NA)
pig_core <- ifel(pigs>0.36, pigs, NA)

global(pig_mar, fun="sum", na.rm=TRUE)
global(pig_mar, fun="mean", na.rm=TRUE)

global(pig_mod, fun="sum", na.rm=TRUE)
global(pig_core, fun="sum", na.rm=TRUE)

sum(values(deer_mar), na.rm=TRUE)

## plot
ggplot() +
  geom_spatraster(data=deer)



# Set crop extent  ext(xmin, xmax, ymin, ymax)
crop_ext <- ext(650000,652500,1020000,1022500)

deer2 <- crop(deer, crop_ext)
pigs2 <- crop(pigs, crop_ext)
plot(pigs2, pax=list(cex.axis=2))

ggplot() +
  geom_spatraster(data=pigs2) +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(title="Suitability",
                               title.position="top",
                                barwidth = unit(20, "cm"), 
                                barheight = unit(1, "cm"))) +
  theme_void() +
  theme(legend.position="bottom",
        legend.text=element_text(size=14),
        legend.title=element_text(size=18, face="bold"))


ggsave("figs/example_hsi_supplement.pdf") #Saving 10.8 x 11.5 in image


global(pigs2, fun="sum", na.rm=TRUE)




# Check new quantiles
global(deer, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)


global(pigs, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)



deer_df <- as.data.frame(deer,na.rm = TRUE)
deer_df$species <- "White-tailed Deer"
names(deer_df)[1] <- "suitability"

pigs_df <- as.data.frame(pigs, na.rm = TRUE)
pigs_df$species <- "Wild Pigs"
names(pigs_df)[1] <- "suitability"

df <- bind_rows(deer_df, pigs_df)

thresholds <- data.frame(species=rep(c("White-tailed Deer", "Wild Pigs"), each=3), 
                         levels=rep(c("Core", "Moderate", "Marginal"), 2),
                         min=c(0.92, 0.68, 0.34, 0.36, 0.31, 0.155),
                         max=c(1, 0.92, 0.68, 1, 0.36, 0.31))

thresholds$levels <- factor(thresholds$levels, levels=c("Core", "Moderate", "Marginal"))



ggplot(df) +
  geom_rect(data=thresholds, aes(xmin=min, xmax=max, ymin=0, ymax=5E6+50000, fill=levels), alpha=0.8) +
  geom_histogram(aes(x=suitability), fill="white", color="black", alpha=0.6) +
  # geom_vline(data=thresholds, aes(xintercept=min, color=levels), linewidth=1) +
  guides(
    fill=guide_legend(title="Suitability Levels", position="inside")#,
    # color=guide_legend(title="Suitability Levels", position="inside")
  ) +
  scale_fill_manual(values=c("#fde725","#21918c","#440154")) +
  # scale_color_manual(values=c("#fde725","#21918c","#440154")) +
  facet_wrap(vars(species), ncol=1) +
  ylab("Number of raster cells") + xlab("Suitability") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        strip.text=element_text(size=14),
        panel.grid=element_blank(),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.position=c(0.01, 0.98),
        legend.justification=c(0,1))

ggsave("figs/suitability_histograms.png", dpi=350) #Saving 7.81 x 9.44 in image
