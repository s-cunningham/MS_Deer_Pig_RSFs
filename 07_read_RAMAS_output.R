library(tidyverse)


rd_txt <- function(x) {
  
  # Read in .txt file
  dat <- read_table(x) %>%
    select(Ptch, ncells, Areakm2, CoreA.)
  
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
    select(spp:patch, Ptch:CoreA.)
  
  return(dat)
}

## 
files <- list.files(path="data/ramas_output/", pattern="_area.txt", full.names=TRUE)


test <- read_table(files[[1]])


str_split(files[[15]], "[:punct:]")[[1]][4]
spp <- gsub("[257]+", "", str_split(files[[15]], "[:punct:]")[[1]][4])
pct <- gsub("[a-z]+", "", str_split(files[[15]], "[:punct:]")[[1]][4])

