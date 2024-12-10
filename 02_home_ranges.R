library(tidyverse)
library(amt)

#### Deer ####
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
deer_res <- data.frame()
size <- data.frame()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- deer_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  
  # Save resampled data
  deer_res <- bind_rows(deer_res, dat)
  
  # Make template raster
  trast <- make_trast(dat, res=10)

  # Create HR
  hr_deer <- hr_akde(dat, trast=trast, level=0.95)
 
  # Calculate HR area
  area <- hr_deer %>% hr_area() %>% mutate(areakm2=area/(1000^2)) %>%
              select(-level, -area) %>%
              pivot_wider(names_from = "what", values_from="areakm2")
  
  # Generate random points
  r1 <- random_points(hr_deer, n = 5*nrow(dat)) 
  # plot(r1)

  # Add ID to each location
  r1 <- r1 %>% mutate(DeerID=un.id[i])

  # Save random points
  avail <- bind_rows(avail, r1)
  
  # Save HR area
  area$DeerID <- un.id[i]
  size <- bind_rows(size, area)
  
}

deer_res <- deer_res %>% rename(DeerID=id, X=x_, Y=y_) %>% 
          select(DeerID, X, Y) %>%
          mutate(type=1) %>% as_tibble()

avail <- avail %>% as_tibble() %>% 
            select(DeerID,  x_, y_, case_) %>%
            rename(X=x_, Y=y_, type=case_) %>%
            mutate(type=as.numeric(type))

all <- bind_rows(deer_res, avail)

write_csv(all, "data/location_data/deer_used_avail.csv")
write_csv(size, "data/location_data/deer_est_akde_homeranges.csv")

#### Pigs ####
set.seed(1)
# Read in deer location data
pigs <- read_csv("data/location_data/pigs_filtered.csv") %>%
  mutate(month=as.numeric(format(timestamp, "%m")),
         year=as.numeric(format(timestamp, "%Y"))) %>%
  select(PigID:study,month,year,timestamp,X,Y) %>%
  unite("PigID", c(1,5), sep="_") %>% # separate by year
  filter(PigID != "35493_Noxubee_2022")

# Create track of all individuals
pigs_tr <- make_track(pigs, X, Y, timestamp, id=PigID, crs=32616)

# Check sampling rate of all individuals
trk_info <- pigs_tr %>% summarize_sampling_rate_many(c("id"), time_unit="min")

range(trk_info$mean)

# Loop over each individual, building home range and extracting points
un.id <- unique(pigs$PigID)
avail <- data.frame()
pigs_res <- data.frame()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- pigs_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  
  # Save resampled data
  pigs_res <- bind_rows(pigs_res, dat)
  
  # Make template raster
  trast <- make_trast(dat, res=10)
  
  # Create HR
  hr_pigs <- hr_akde(dat, trast=trast, level=0.95)
  
  # Calculate HR area
  area <- hr_pigs %>% hr_area() %>% mutate(areakm2=area/(1000^2)) %>%
    select(-level, -area) %>%
    pivot_wider(names_from = "what", values_from="areakm2")
  
  # Generate random points
  r1 <- random_points(hr_pigs, n = 5*nrow(dat)) 
  # plot(r1)
  
  # Add ID to each location
  r1 <- r1 %>% mutate(PigID=un.id[i])
  
  # Save random points
  avail <- bind_rows(avail, r1)
  
  # Save HR area
  area$PigID <- un.id[i]
  size <- bind_rows(size, area)
}

pigs_res <- pigs_res %>% rename(PigID=id, X=x_, Y=y_) %>% 
  select(PigID, X, Y) %>%
  mutate(type=1) %>% as_tibble()

avail <- avail %>% as_tibble() %>% 
  select(PigID,  x_, y_, case_) %>%
  rename(X=x_, Y=y_, type=case_) %>%
  mutate(type=as.numeric(type))

all <- bind_rows(pigs_res, avail)

write_csv(all, "data/location_data/pigs_used_avail.csv")
write_csv(size, "data/location_data/pigs_est_akde_homeranges.csv")
