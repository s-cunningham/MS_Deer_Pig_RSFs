library(tidyverse)
library(amt)

set.seed(1)

# Read in deer location data
deer <- read_csv("data/location_data/deer_filtered.csv") %>%
          mutate(month=as.numeric(format(timestamp, "%m")),
                 year=as.numeric(format(timestamp, "%Y"))) %>%
          select(DeerID:study,month,year,timestamp,X,Y) %>%
          unite("DeerID", c(1,5), sep="_") # separate by year

# Create track of all individuals
deer_tr <- make_track(deer, X, Y, timestamp, id=DeerID, crs=32616)

# Check sampling rate of all individuals
trk_info <- deer_tr %>% summarize_sampling_rate_many(c("id"), time_unit="min")

range(trk_info$mean)

# Loop over each individual, building home range and extracting points
un.id <- unique(deer$DeerID)
avail <- data.frame()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- deer_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(15))
  
  # Make template raster
  trast <- make_trast(dat, res=10)

  # Create HR
  hr_deer <- hr_kde(dat, trast=trast, level=0.95)

  # Generate random points
  r1 <- random_points(hr_deer, n = 5*nrow(dat)) 
  # plot(r1)

  # Add ID to each location
  r1 <- r1 %>% mutate(DeerID=un.id[i])
  
  # Save points
  avail <- bind_rows(avail, r1)
  
}

deer <- deer %>% select(DeerID, X, Y) %>%
          mutate(type=1)

avail <- avail %>% as_tibble() %>% 
            select(DeerID,  x_, y_, case_) %>%
            rename(X=x_, Y=y_, type=case_) %>%
            mutate(type=as.numeric(type))

deer <- bind_rows(deer, avail)

write_csv(deer, "data/location_data/deer_used_avail.csv")
