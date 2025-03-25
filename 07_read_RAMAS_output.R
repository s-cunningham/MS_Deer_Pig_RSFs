## This script reads in text files with RAMAS results (Area and Carrying Capacity), then summarizes total K and area

library(tidyverse)

#### Functions to read .txt files from RAMAS ####
area_RAMAS <- function(x) {
  
  # Read in .txt file
  dat <- read_table(x) %>%
    select(Ptch, ncells, Areakm2, CoreA.) %>%
    rename(PatchNo=Ptch)
  
  # Extract pieces from file name
  spp <- gsub("[257]+", "", str_split(x, "[:punct:]")[[1]][4])
  pct <- gsub("[a-z]+", "", str_split(x, "[:punct:]")[[1]][4])
  patch <- str_split(x, "[:punct:]")[[1]][5] 
  density <- str_split(x, "[:punct:]")[[1]][6] 
  
  # Add to data frame
  dat$spp <- spp
  dat$density <- density
  dat$pct <- pct
  dat$patch <- patch

  # Rearrange columns
  dat <- dat %>%
    select(spp:patch, PatchNo:CoreA.)
  
  return(dat)
}
K_RAMAS <- function(x) {
  
  # Read in .txt file
  dat <- read_table(x) %>%
    select(Ptch:N0) %>%
    rename(PatchNo=Ptch)
  
  # Extract pieces from file name
  spp <- gsub("[257]+", "", str_split(x, "[:punct:]")[[1]][4])
  pct <- gsub("[a-z]+", "", str_split(x, "[:punct:]")[[1]][4])
  patch <- str_split(x, "[:punct:]")[[1]][5] 
  density <- str_split(x, "[:punct:]")[[1]][6] 
  
  # Add to data frame
  dat$spp <- spp
  dat$density <- density
  dat$pct <- pct
  dat$patch <- patch
  
  # Rearrange columns
  dat <- dat %>%
    select(spp:patch, PatchNo:N0)
  
  return(dat)
}

#### Read in data ####

## Read in area
# List all the files in the ramas_output folder that end in _area.txt
files <- list.files(path="data/ramas_output/", pattern="_area.txt", full.names=TRUE)
# Read all files in the list
files <- lapply(files, area_RAMAS)
# Reduce the list to a data frame (tibble)
area <- do.call("bind_rows", files)
# convert density and percentage columns to numeric
area <- area %>% mutate(across(2:3, as.numeric))

## Read in carrying capacity
# List all the files in the ramas_output folder that end in _K.txt
files <- list.files(path="data/ramas_output/", pattern="_K.txt", full.names=TRUE)
# Read all files in the list
files <- lapply(files, K_RAMAS)
# Reduce the list to a data frame (tibble)
k <- do.call("bind_rows", files)
# convert density and percentage columns to numeric
k <- k %>% mutate(across(2:3, as.numeric))

## Join K and area
dat <- left_join(k, area, by=c("spp", "density", "pct", "patch", "PatchNo"))

# Drop Area and recalculate
dat <- dat %>%
  # drop area columns
  select(-Areakm2, -CoreA.) %>%
  # Calculate area of 90 x 90 m cells, convert to km2
  mutate(Areakm2 = (ncells * 90 * 90)/1000/1000) 

# Drop patches that have K < 6 (might be variable between densities)
# dat <- dat %>% 
#   filter(K >= 6)

# Summarize K and area of reduced patch size
dat %>% group_by(density, pct, patch) %>%
  reframe(K=sum(K), Area=sum(Areakm2)) %>%
  arrange(patch)


