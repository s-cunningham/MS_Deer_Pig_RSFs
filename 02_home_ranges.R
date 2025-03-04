library(tidyverse)
library(amt)
library(sf)
library(adehabitatLT)

## Create template raster for 


#### Deer ####
set.seed(1)
# Read in deer location data
deer <- read_csv("data/location_data/deer_filtered.csv") %>%
          mutate(month=as.numeric(format(timestamp, "%m")),
                 year=as.numeric(format(timestamp, "%Y"))) %>%
          dplyr::select(DeerID:study,month,year,timestamp,X,Y) %>%
          unite("DeerID", c(1,5), sep="_") # separate by year

# Create track of all individuals
deer_tr <- make_track(deer, X, Y, timestamp, id=DeerID, crs=32616)

# Check sampling rate of all individuals
trk_info <- deer_tr %>% summarize_sampling_rate_many(c("id"), time_unit="min")

# List individual IDs
un.id <- unique(deer$DeerID)

## Loop over each ID and 
# 1. Resample to every 4 hours
# 2. Build ctmm home range
# 3. Extract available points from home range

## Loop over individuals
avail <- data.frame()
deer_res <- data.frame()
size <- data.frame()
dist <- numeric()
for (i in 1:lenght(un.id)) {
  
  # Filter to ID
  dat <- deer_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  
  # Convert to ltraj
  df <- as.data.frame(dat)
  df <- as.ltraj(xy=df[,c("x_","y_")], date=df$t_, id=df$id)
  
  test <- lavielle(df, Lmin=10, which="R2n", Kmax=5, plotit=TRUE)
  fp <- findpath(test, 5)
  
  ## Regularize trajectories, round values
  refda <- min(df[[1]]$date)
  df_NA <- setNA(df, refda, 240, units="min") 
  df <- sett0(df_NA, refda, 240, units="min")
  
  # Save distance between locations
  dist <- c(dist, df[[1]]$dist)
  
}



#### Pigs ####
# Read in deer location data
pigs <- read_csv("data/location_data/pigs_filtered.csv") %>%
  mutate(month=as.numeric(format(timestamp, "%m")),
         year=as.numeric(format(timestamp, "%Y"))) %>%
  dplyr::select(PigID:study,month,year,timestamp,X,Y) %>%
  unite("PigID", c(1,5), sep="_") %>% # separate by year
  filter(PigID != "35493_Noxubee_2022")
