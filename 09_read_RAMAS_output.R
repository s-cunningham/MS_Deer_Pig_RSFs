## This script reads in text files with RAMAS results (Area and Carrying Capacity), then summarizes total K and area

library(tidyverse)


source("00_functions.R")

## abundance
N <- list.files(path="results/simulations/", pattern="_abundance.txt$", full.names=TRUE)
N <- lapply(N, ramas_abun)
N <- do.call("bind_rows", N)

# Calculate lambda and geometric mean for each scenario
N <- N %>%
  group_by(species, density, bin) %>%
  reframe(lambda=lambda_calc(average)) %>%
  group_by(species, density, bin) %>%
  reframe(lamda_mgm=gm_mean(lambda))


ggplot(N) +
  geom_line(aes(x=time, y=average, group=factor(bin), color=factor(bin))) +
  facet_wrap(vars(species, density))


## Patch area
area <- list.files(path="results/simulations/", pattern="_area.txt$", full.names=TRUE)
area <- lapply(area, ramas_area)
area <- do.call("bind_rows", area)

area <- area %>%
  group_by(species, density, bin) %>%
  reframe(Areakm2=sum(Area_km2))
  
ggplot(area) +
  geom_point(aes(x=factor(bin), y=Areakm2, group=density, color=density)) +
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

ggplot(K) +
  geom_point(aes(x=factor(bin), y=K, group=density, color=density)) +
  facet_wrap(vars(species)) +
  theme_bw()
