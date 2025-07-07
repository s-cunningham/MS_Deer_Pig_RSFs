## This script reads in text files with RAMAS results (Area and Carrying Capacity), then summarizes total K and area

library(tidyverse)
library(ggtext)

source("00_functions.R")

## abundance
N <- list.files(path="results/simulations/", pattern="_pop1.txt$", full.names=TRUE)
N <- lapply(N, ramas_abun)
N <- do.call("bind_rows", N)

N <- N |> 
  # Filter to keep only rows with 0 or no digits
  filter(str_detect(bin, "0") | !str_detect(bin, "\\d")) |>
  # drop marginal
  filter(bin!="marginal") |>
  # make 'above0' layer 'marginal'
  mutate(bin=if_else(bin=="0", "marginal", bin)) |>
  # drop mean column (didn't do CIs for pop simulation)
  select(-value) |>
  # convert density to categorical
  mutate(density=if_else(density==14 | density==8, "low", "high")) |>
  # keep just the high densityies
  filter(density=="high")

# Reorder so patches are not in alphabetical order
N$bin <- factor(N$bin, levels=c("core", "marginal", "highlymarginal"), labels=c("Core", "Marginal", "Highly\nMarginal"))

ggplot(N) +
  geom_ribbon(aes(x=time, ymin=lowerSD, ymax=upperSD, group=bin, color=bin, fill=bin), alpha=0.3) +
  scale_y_continuous(breaks=seq(0, 30000000, by=5000000), labels=seq(0,30, by=5)) +
  geom_line(aes(x=time, y=average, group=bin, color=bin)) +
  scale_color_manual(values=c("#fde725","#21918c","#440154"), name="Patch Type") + 
  scale_fill_manual(values=c("#fde725","#21918c","#440154"), name="Patch Type") + 
  facet_grid(.~species) +
  guides(
    colour = guide_legend(position = "inside"),
    fill = guide_legend(position = "inside")
  ) + 
  ylab("Abundance (in millions)") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA),
        legend.position.inside=c(0,1),
        legend.justification=c(0,1),
        legend.background=element_rect(fill=NA))

# Calculate lambda and geometric mean for each scenario
pop_growth <- N |>
  group_by(species, density, bin) |>
  reframe(lambda=lambda_calc(average)) |>
  group_by(species, density, bin) |>
  reframe(lambda_mgm=gm_mean(lambda)) 

N <- N |>
  mutate(species=if_else(species=="deer", "White-tailed Deer", "Wild Pigs"))

## Carrying capacity
K <- list.files(path="results/simulations/", pattern="_K.txt$", full.names=TRUE)
K <- lapply(K, ramas_K)
K <- do.call("bind_rows", K)

K <- K |>
  # drop patches that have K of < 5
  # filter(K >= 5) |>
  # filter(value=="mean") |>
  # summarize by species, density and bin
  group_by(species, density, bin, value) |>
  reframe(K=sum(K), N0=sum(N0), avgHS=mean(AvgHS)) |> # might need to recalculate avgHS
  # replace actual density number with low and high
  mutate(density=if_else(density==14 | density==8, "low", "high")) 

# Relabel species
K <- K |>
  mutate(species=if_else(species=="deer", "White-tailed Deer", "Wild Pigs"))

## Look at patch types
Kpatch <- K |> 
  # Filter to keep only rows with 0 or no digits
  filter(str_detect(bin, "0") | !str_detect(bin, "\\d")) |>
  # Drop initial abundance and average suitability value
  select(-N0, -avgHS) |>
  # drop marginal
  filter(bin!="marginal") |>
  # make 'above0' layer 'marginal'
  mutate(bin=if_else(bin=="0", "marginal", bin)) |>
  # pivot to wide format
  pivot_wider(names_from="value", values_from="K") |>
  # reorder columns
  select(species, density, bin, UCI, mean, LCI)

# Reorder patch types so that they're in order of size, not alphabetical
Kpatch$bin <- factor(Kpatch$bin, levels=c("core", "marginal", "highlymarginal"), labels=c("Core", "Marginal", "Highly\nMarginal"))
# Relabel density to be more informative for plot
Kpatch$density <- factor(Kpatch$density, levels=c("high", "low"), labels=c("22 deer | 27 pigs", "14.3 deer | 8 pigs"))

ggplot(Kpatch) +
  geom_linerange(aes(x=bin, ymin=UCI, ymax=LCI, color=density),position=position_dodge(.2), linewidth=1) +
  geom_point(aes(x=bin, y=mean, color=density), position=position_dodge(.2), size=2) +
  scale_color_manual(values=c("#21918c","#440154"), name='Density (animals/km<sup>2</sup>)') +
  scale_y_continuous(breaks=seq(0, 3000000, by=500000), labels=c(0,0.5,1,1.5,2,2.5,3)) +
  xlab("Patch Type") + ylab("Carrying Capacity (in millions)") +
  guides(
    colour = guide_legend(position = "inside")
  ) + 
  facet_grid(~species) +
  theme_classic() +
  theme(legend.position.inside=c(0,1),
        legend.justification=c(0,1),
        legend.title = element_markdown(),
        legend.text=element_text(size=10),
        legend.background = element_rect(fill=NA),
        strip.text=element_blank(),
        axis.title=element_text(size=11),
        axis.text=element_text(size=10),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))
ggsave("figs/patch_carrying_capacity.svg")
# Saving 7.18 x 4.26 in image


Kpatch <- Kpatch %>% arrange(species, density, bin)


Nt <- N %>%
  filter(species=="White-tailed Deer" & density=="high" & bin=="Highly\nMarginal")

y <- Nt$average[2:nrow(Nt)] - Nt$average[1:(nrow(Nt)-1)]

x <- Nt$average[1:(nrow(Nt)-1)]

plot(y=y, x=x)

