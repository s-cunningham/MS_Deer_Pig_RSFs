library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(ggspatial)

# Function
beta_seq <- function(y) {
  if (y < 0) {
    x <- seq(0, y, by = -0.001)
  } else {
    x <- seq(0, y, by = 0.001)
  }
  return(x)
}

# set ggplot theme
theme_set(theme_classic())

## Read predicted maps (90 x 90)
deer <- rast("results/predictions/deer_rsf_bins_30m.tif")
pigs <- rast("results/predictions/pigs_rsf_bins_30m.tif")

# Read in MS polygon
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
ms <- project(ms, deer)


d_map <- ggplot() +
  geom_spatraster(data=deer) +
  scale_x_continuous(expand=expansion(0)) + 
  scale_y_continuous(expand=expansion(0)) +
  scale_fill_grass_d(palette="viridis", direction=1, na.translate=FALSE) +
  guides(fill=guide_legend(nrow=1,direction="horizontal")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_blank()) +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.01, "npc"),
    pad_y = unit(0.80, "npc"),
    style = north_arrow_minimal()) 

p_map <- ggplot() +
  geom_spatraster(data=pigs) +
  scale_x_continuous(expand=expansion(0)) + 
  scale_y_continuous(expand=expansion(0)) +
  scale_fill_grass_d(palette="viridis", direction=1, na.translate=FALSE) +
  guides(fill=guide_legend(nrow=1,direction="horizontal")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_blank()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.5,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.07, "npc"),
    text_cex = .8) 

maps <- d_map + p_map + plot_layout(guides = 'collect') & theme(legend.position = 'bottom') 


## Read regression coefficients
beta_d <- read_csv("output/deer_mixed_effects_betas.csv") %>%
  # switch format
  pivot_wider(names_from="covariate", values_from="beta") %>%
  # rename all hardwoods column to match pigs
  rename(hardwoods=allhardwoods, water=water_dist, water2=`I(water_dist^2)`) %>%
  # add column for species
  mutate(species="deer")
beta_p <- read_csv("output/pigs_mixed_effects_betas.csv") %>%
  # switch format
  pivot_wider(names_from="covariate", values_from="beta") %>%
  # reorder columns to match deer
  select(hardwoods, gramanoids, shrubs, developed:Water2) %>%
  # rename water
  rename(water=dist_water, water2=Water2) %>%
  # add column for species
  mutate(species="pigs")

# Read beta standard errors
se_d <- read_csv("output/deer_mixed_effects_se.csv") %>%
  # switch format
  pivot_wider(names_from="covariate", values_from="V1") %>%
  # rename all hardwoods column to match pigs
  rename(hardwoods=allhardwoods, water=water_dist, water2=`I(water_dist^2)`) %>%
  # add column for species
  mutate(species="deer")
se_p <- read_csv("output/pigs_mixed_effects_se.csv") %>%
  # switch format
  pivot_wider(names_from="covariate", values_from="V1") %>%
  # rename water
  rename(water=dist_water, water2=`I(dist_water^2)`) %>%
  # reorder columns to match deer
  select(hardwoods, gramanoids, shrubs, developed:water2) %>%
  # add column for species
  mutate(species="pigs")

## Combine species and coefficients
# beta coefficients
beta <- bind_rows(beta_d, beta_p)  

# Pivot so that covariates are in a in a column, and beta coefficients in another column
beta <- beta %>% pivot_longer(1:(ncol(.)-1), names_to="covariate", values_to="beta") %>%
  mutate(beta=round(beta, digits=3)) %>%
  # drop row for pigs + foodcrops
  filter(!is.na(beta)) %>%
  # create single column for species & covariate
  unite("sp_cov", c("species", "covariate"), sep="_", remove=FALSE)

# standard errors
se <- bind_rows(se_d, se_p)

# Pivot so that SEs are in a in a column, and beta coefficients in another column
se <- se %>% pivot_longer(1:(ncol(.)-1), names_to="covariate", values_to="se")

se <- left_join(se, beta, by=c("species", "covariate")) %>%
  # drop row for pigs + foodcrops
  filter(!is.na(se)) %>%
  # create single column for species & covariate
  unite("covar", c("species", "covariate"), sep="_", remove=FALSE) %>%
  separate(covar, into=c("species", "covariate"), sep="_", remove=FALSE) %>%
  # refactor covariates for plotting (create numeric indicator for covariate) 
  mutate(x_covar = fct_recode(covar, p1="deer_hardwoods",p2="deer_gramanoids",p3="deer_shrubs",p4="deer_foodcrops",p5="deer_developed",
                              p6="deer_water",p7="deer_water2",p8="pigs_hardwoods",p9="pigs_gramanoids",p10="pigs_shrubs",
                              p11="pigs_developed",p12="pigs_water",p13="pigs_water2"),
         x_covar = as.character(x_covar)) %>%
  # convert x_covar to number
  mutate(x_covar=parse_number(x_covar)) %>%
  mutate(se_min=if_else(beta>0, beta-se, beta-(-1*se)),
         se_max=if_else(beta>0, beta+se, beta+(-1*se))) %>%
  select(covar,species, covariate, x_covar,se, se_min, beta, se_max)

# Interpolate values from 0 to beta
vals <- lapply(beta$beta, beta_seq)

# Create corresponding number of x values
y <- rep(beta$sp_cov, lengths(vals))
vals <- unlist(vals)

# Combine into data frame
coefs <- data.frame(covar=y, beta=vals)

# Convert to tibble and split species and covariate
coefs <- coefs %>% as_tibble() %>%
  separate(covar, into=c("species", "covariate"), sep="_", remove=FALSE) %>%
  # refactor covariates for plotting (create numeric indicator for covariate) 
  mutate(x_covar = fct_recode(covar, p1="deer_hardwoods",p2="deer_gramanoids",p3="deer_shrubs",p4="deer_foodcrops",p5="deer_developed",
                              p6="deer_water",p7="deer_water2",p8="pigs_hardwoods",p9="pigs_gramanoids",p10="pigs_shrubs",
                              p11="pigs_developed",p12="pigs_water",p13="pigs_water2"),
         x_covar = as.character(x_covar)) %>%
  # convert x_covar to number
  mutate(x_covar=parse_number(x_covar)) %>%
  # create x values
  mutate(x = x_covar - 0.4,
         xend = x_covar + 0.4)

covar_labels <- unique(coefs$covar)
covar_labels <- gsub("pigs_", "", covar_labels)
covar_labels <- gsub("deer_", "", covar_labels)
covar_labels <- c("Hardwoods", "Graminoids", "Shrubs", "Crops", "Developed",
                  "Dist. Water", expression(Dist.~Water^2), "Hardwoods", "Graminoids",
                     "Shrubs", "Developed", "Dist. Water", expression(Dist.~Water^2))

# Create plot showing beta coefficients
coef_plt <- ggplot(se) +
  geom_hline(yintercept=0, color="gray") +
  geom_errorbar(aes(x=x_covar, ymin=se_min, ymax=se_max), width=0.1, color="black") +
  geom_point(aes(x=x_covar, y=beta, fill=beta), shape=22, size=4) +
  scale_fill_viridis_c(option="D") +
  scale_x_continuous(breaks=1:13, labels=covar_labels) +
  facet_wrap(vars(species), scales="free") +
  labs(y = "Selection Coefficient") +
  theme(#panel.border=element_rect(fill=NA, color="black", linewidth=1),
    legend.position="none",
    strip.text=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=12),
    axis.text.x=element_text(angle=25, vjust=0.8, hjust=0.8, size=11),
    axis.text.y=element_text(size=11))

# Set up layout
layout <- "
AAAA
AAAA
AAAA
BBBB
"
# Combine figures
maps / free(coef_plt) + plot_layout(design = layout) #+ plot_annotation(tag_levels = 'a')

# save figure
ggsave("figs/Fig2_maps_betas.svg")
# Saving 7.71 x 9 in image





