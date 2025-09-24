## This script reads in text files with RAMAS results (Area and Carrying Capacity), then summarizes total K and area

library(tidyverse)
library(ggtext)

source("00_functions.R")

### Carrying capacity ####
K <- list.files(path="results/simulations/", pattern="_K.txt$", full.names=TRUE)
K <- lapply(K, ramas_K)
K <- do.call("bind_rows", K)

K <- K |>
  # drop patches that have K of < 5
  filter(K >= 5) |>
  # filter(value=="mean") |>
  # summarize by species, density 
  group_by(species, density, suitability, value) |>
  reframe(K=sum(K), avgHS=mean(AvgHS), totHS=sum(TotalHS)) |> # might need to recalculate avgHS
  # replace actual density number with low and high
  mutate(density=if_else(density==14 | density==8, "low", "high")) 

# Relabel species
K <- K |>
  mutate(species=if_else(species=="deer", "White-tailed Deer", "Wild Pigs"))

## Look at patch types
Kpatch <- K |>
  # Drop initial abundance and average suitability value
  select(-avgHS, -totHS) |>
  # pivot to wide format
  pivot_wider(names_from="value", values_from="K") |>
  # reorder columns
  select(species, density, suitability, LCI, mean, UCI)

# Fill in 0 for the deer one that couldn't find patches (no core)
Kpatch$LCI[Kpatch$species=="White-tailed Deer" & Kpatch$suitability=="core"] <- 0

# Reorder patch types so that they're in order of size, not alphabetical
Kpatch$suitability <- factor(Kpatch$suitability, levels=c("core", "moderate", "marginal"), labels=c("Core", "Moderate", "Marginal"))
# Relabel density to be more informative for plot
Kpatch$density <- factor(Kpatch$density, levels=c("low", "high"), labels=c("14.3 deer | 8 pigs", "22 deer | 27 pigs"))

## Add dashed line for MDWFP estimate
Kpatch <- Kpatch |>
  mutate(mdwfp_est = if_else(species=="Wild Pigs", 1000000, 1610000))


ggplot(Kpatch) +
  coord_cartesian(ylim=c(0, 2500000)) +
  geom_hline(aes(yintercept=mdwfp_est), linetype=2) +
  geom_linerange(aes(x=suitability, ymin=UCI, ymax=LCI, color=density),position=position_dodge(.2), linewidth=1) +
  geom_point(aes(x=suitability, y=mean, color=density), position=position_dodge(.2), size=2) +
  scale_color_manual(values=c("#ed7953","#9c179e"), name='Density (animals/km<sup>2</sup>)') +
  scale_y_continuous(breaks=seq(0, 2500000, by=500000), labels=c(0,0.5,1,1.5,2,2.5)) +
  xlab("Suitability Level") + ylab("Carrying Capacity (in millions)") +
  guides(
    colour = guide_legend(position = "inside")
  ) + 
  facet_grid(~species) +
  theme_classic() +
  theme(legend.position.inside=c(0,1),
        legend.justification=c(0,1),
        legend.title = element_markdown(size=12),
        legend.text=element_text(size=11),
        legend.background = element_rect(fill=NA),
        strip.text=element_blank(),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))


ggsave("figs/patch_carrying_capacity.svg")
# Saving 7.18 x 4.26 in image

## Look at suitability in mean patches
K <- K |>
  filter(value=="mean")

### Area ####
A <- list.files(path="results/simulations/", pattern="_area.txt$", full.names=TRUE)
A <- lapply(A, ramas_area)
A <- do.call("bind_rows", A)

# Summarize small patches
A <- A |>
  group_by(species, density, level, value) |>
  reframe(Ncells=sum(Ncells), Area_km2=sum(Area_km2))

A |> filter(value=="mean")
