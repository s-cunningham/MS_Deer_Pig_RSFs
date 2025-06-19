## This script reads in text files with RAMAS results (Area and Carrying Capacity), then summarizes total K and area

library(tidyverse)


source("00_functions.R")

# Read file with threshold values
dat <- read_csv("data/logrsf_cutoffs.csv") %>%
  mutate(bin=gsub(" ", "",bin))

## abundance
N <- list.files(path="results/simulations/", pattern="_abundance.txt$", full.names=TRUE)
N <- lapply(N, ramas_abun)
N <- do.call("bind_rows", N)

# Drop density dependence
N <- N %>% filter(densDep=="exp" & bin!="10") %>% select(-densDep)

ggplot(N) +
  geom_line(aes(x=time, y=average, group=factor(bin), color=factor(bin))) +
  facet_wrap(vars(species, density))

# Calculate lambda and geometric mean for each scenario
N <- N %>%
  group_by(species, density, bin) %>%
  reframe(lambda=lambda_calc(average)) %>%
  group_by(species, density, bin) %>%
  reframe(lambda_mgm=gm_mean(lambda)) %>%
  mutate(density=if_else(density==14 | density==8, "low", "high")) 

N <- left_join(N, dat, by=c("species", "bin"))

N <- N %>%
  mutate(species=if_else(species=="deer", "White-tailed Deer", "Wild Pigs"))

zero_line <- data.frame(species=c("White-tailed Deer", "Wild Pigs"),
                        cutoff=c(0.884, 0.808))

ggplot(N) +
  geom_hline(yintercept=1, linewidth=0.5) +
  geom_vline(data=zero_line, aes(xintercept=cutoff), linewidth=0.2) +
  geom_point(aes(x=cutoff, y=lambda_mgm, group=density, color=density)) +
  facet_grid(species~density) +
  theme_bw()




## Patch area
area <- list.files(path="results/simulations/", pattern="_area.txt$", full.names=TRUE)
area <- lapply(area, ramas_area) %>% suppressWarnings()
area <- do.call("bind_rows", area)

area <- area %>%
  group_by(species, density, bin) %>%
  reframe(Areakm2=sum(Area_km2)) %>%
  # replace actual density number with low and high
  mutate(density=if_else(density==14 | density==8, "low", "high"))
  
area <- left_join(area, dat, by=c("species", "bin"))

ggplot(area) +
  geom_point(aes(x=cutoff, y=Areakm2, group=density, color=density)) +
  facet_wrap(vars(species)) +
  theme_bw()

## Carrying capacity
K <- list.files(path="results/simulations/", pattern="_K.txt$", full.names=TRUE)
K <- lapply(K, ramas_K)
K <- do.call("bind_rows", K)

K <- K %>%
  # drop patches that have K of < 5
  filter(K >= 5) %>%
  # summarize by species, density and bin
  group_by(species, density, bin) %>%
  reframe(K=sum(K), N0=sum(N0), avgHS=mean(AvgHS)) %>% # might need to recalculate avgHS
  # replace actual density number with low and high
  mutate(density=if_else(density==14 | density==8, "low", "high")) 

# Join to cutoff value
K <- left_join(K, dat, by=c("species", "bin"))

K <- K %>%
  mutate(species=if_else(species=="deer", "White-tailed Deer", "Wild Pigs"))

zero_line <- data.frame(species=c("White-tailed Deer", "Wild Pigs"),
                        cutoff=c(0.884, 0.808))

ggplot() +
  geom_vline(data=zero_line, aes(xintercept=cutoff), linewidth=0.2) +
  geom_point(data=K, aes(x=cutoff, y=K, group=density, shape=density, color=density)) +
  scale_shape_manual(values=c(21, 23)) +
  facet_wrap(vars(species)) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.5))

