library(tidyverse)
library(amt)
library(sf)
library(adehabitatLT)

#### Deer ####
set.seed(1)
# Read in deer location data
deer <- read_csv("data/location_data/deer_filtered.csv") %>%
          mutate(month=as.numeric(format(timestamp, "%m")),
                 year=as.numeric(format(timestamp, "%Y"))) %>%
          dplyr::select(DeerID:study,month,year,timestamp,X,Y) %>%
          unite("DeerID", c(1,5), sep="_") #%>% # separate by year
          # filter(DeerID!="81864_Delta_2022" & DeerID!="27_Central_2020" & DeerID!="340_Central_2019" & DeerID!="27_Central_2019") %>%
          # mutate(month=as.numeric(format(timestamp, "%m")),
          #        day=as.numeric(format(timestamp, "%d"))) %>%
          # filter(!(DeerID=="87924_Delta_2021" & month<6)) %>%
          # filter(!(DeerID=="81864_Delta_2021" & month<8)) %>%
          # mutate(DeerID=case_when(DeerID!="87924_Delta_2022" ~ DeerID,
          #                         DeerID=="87924_Delta_2022" & (month<3 | month>=10) ~ "87924_Delta_2022-1", 
          #                         DeerID=="87924_Delta_2022" & (month>4 & month<9) ~ "87924_Delta_2022-2")) %>%
          # filter(!is.na(DeerID)) %>%
          # dplyr::select(-month, -day)

un.id <- unique(deer$DeerID)


deer %>% group_by(DeerID, month) %>% count() 


# Check how many months of data for each deer
short <- deer %>% mutate(month=format(timestamp, "%m")) %>% group_by(DeerID) %>%
  reframe(n=length(unique(month))) %>%
  filter(n <= 1)

# Remove individuals with < 3 months data
deer <- deer %>% filter(!(DeerID %in% short$DeerID))

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
dist <- numeric()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- deer_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  
  df <- as.data.frame(dat)
  df <- as.ltraj(xy=df[,c("x_","y_")], date=df$t_, id=df$id)

  test <- lavielle(df, Lmin=10, which="R2n", K=3, plotit=TRUE)
  fp <- findpath(test, 3)
  
  ## Regularize trajectories, round values
  refda <- min(df[[1]]$date)
  df_NA <- setNA(df, refda, 240, units="min") 
  df <- sett0(df_NA, refda, 240, units="min")
  
  # Save distance between locations
  dist <- c(dist, df[[1]]$dist)
  
  plot(df[[1]]$R2n, type="l")
  
  
  # Save resampled data
  deer_res <- bind_rows(deer_res, dat)
  
  # Make template raster
  trast <- make_trast(dat, res=30)

  # Create HR
  ## Need to check whether full spatial coverage of available points in HR
  hr_deer <- hr_mcp(dat, level=1)
  hr_deer <- hr_akde(dat, trast=trast, levels=0.95)
  
  # Calculate HR area
  area <- hr_deer %>% hr_area() %>% mutate(areakm2=area/(1000^2)) %>%
              dplyr::select(-level, -area) %>%
              pivot_wider(names_from = "what", values_from="areakm2")
  
  # Generate random points
  r1 <- random_points(hr_deer, n = 10*nrow(dat), type="regular") 
  plot(r1)
  points(dat$x_, dat$y_, col="red", pch=16)

  # Add ID to each location
  r1 <- r1 %>% mutate(DeerID=un.id[i])

  # Save random points
  avail <- bind_rows(avail, r1)
  
  # Save HR area
  area$DeerID <- un.id[i]
  size <- bind_rows(size, area)
  
}

## overall average distance
median(dist, na.rm=TRUE)

# combine all deer points
deer_res <- deer_res %>% rename(DeerID=id, X=x_, Y=y_) %>% 
          dplyr::select(DeerID, X, Y) %>%
          mutate(type=1) %>% as_tibble()

# combine all deer available points
avail <- avail %>% as_tibble() %>% 
            dplyr::select(DeerID,  x_, y_, case_) %>%
            rename(X=x_, Y=y_, type=case_) %>%
            mutate(type=as.numeric(type))

# combine used and available
all <- bind_rows(deer_res, avail)

write_csv(all, "data/location_data/deer_used_avail_all.csv")
write_csv(size, "data/location_data/deer_est_akde_homeranges_all.csv")

# all <- read_csv("data/location_data/deer_used_avail.csv")
# 
# deer_sf <- all %>% filter(type==1) %>% separate("DeerID", into=c("ID", "study", "year"), remove=FALSE) %>%
#               st_as_sf(coords=c("X", "Y"), crs=32616)
# st_write(deer_sf, "data/location_data/filtered_deer_locations.shp")

#### Pigs ####
set.seed(1)
# Read in deer location data
pigs <- read_csv("data/location_data/pigs_filtered.csv") %>%
  mutate(month=as.numeric(format(timestamp, "%m")),
         year=as.numeric(format(timestamp, "%Y"))) %>%
  dplyr::select(PigID:study,month,year,timestamp,X,Y) %>%
  unite("PigID", c(1,5), sep="_") %>% # separate by year
  filter(PigID != "35493_Noxubee_2022")

# Check how many months of data for each pig
short <- pigs %>% mutate(month=format(timestamp, "%m")) %>% group_by(PigID) %>%
  reframe(n=length(unique(month))) %>%
  filter(n < 3)

# Remove individuals with < 3 months data
pigs <- pigs %>% filter(!(PigID %in% short$PigID))

# Create track of all individuals
pigs_tr <- make_track(pigs, X, Y, timestamp, id=PigID, crs=32616)

# Check sampling rate of all individuals
trk_info <- pigs_tr %>% summarize_sampling_rate_many(c("id"), time_unit="min")

range(trk_info$mean)

# Loop over each individual, building home range and extracting points
un.id <- unique(pigs$PigID)
avail <- data.frame()
pigs_res <- data.frame()
size <- data.frame()
dist <- numeric()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- pigs_tr %>% filter(id==un.id[i])

  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  
  df <- as.data.frame(dat)
  df <- as.ltraj(xy=df[,c("x_","y_")], date=df$t_, id=df$id)
  
  ## Regularize trajectories, round values
  refda <- min(df[[1]]$date)
  df_NA <- setNA(df, refda, 240, units="min") # Locations every 3 hours 
  df <- sett0(df_NA, refda, 240, units="min")
  
  # Save distance between locations
  dist <- c(dist, df[[1]]$dist)
  
  # Save resampled data
  pigs_res <- bind_rows(pigs_res, dat)
  
  # Make template raster
  # trast <- make_trast(dat, res=30)
  
  # Create HR
  hr_pigs <- hr_mcp(dat, level=1)
  
  # Calculate HR area
  area <- hr_pigs %>% hr_area() %>% mutate(areakm2=area/(1000^2)) %>%
    dplyr::select(-level, -area) %>%
    pivot_wider(names_from = "what", values_from="areakm2")
  
  # Generate random points
  r1 <- random_points(hr_pigs, n = 10*nrow(dat)) 
  plot(r1, main=un.id[i])
  points(dat$x_, dat$y_, col="red", pch=16)
  
  # Add ID to each location
  r1 <- r1 %>% mutate(PigID=un.id[i])
  
  # Save random points
  avail <- bind_rows(avail, r1)
  
  # Save HR area
  area$PigID <- un.id[i]
  size <- bind_rows(size, area)
}

## overall average distance
median(dist, na.rm=TRUE)
sd(dist, na.rm=TRUE)

pigs_res <- pigs_res %>% rename(PigID=id, X=x_, Y=y_) %>% 
  dplyr::select(PigID, X, Y) %>%
  mutate(type=1) %>% as_tibble()

avail <- avail %>% as_tibble() %>% 
  dplyr::select(PigID,  x_, y_, case_) %>%
  rename(X=x_, Y=y_, type=case_) %>%
  mutate(type=as.numeric(type))

all <- bind_rows(pigs_res, avail)

write_csv(all, "data/location_data/pigs_used_avail.csv")
write_csv(size, "data/location_data/pigs_est_akde_homeranges.csv")


pig_sf <- all %>% filter(type==1) %>% separate("PigID", into=c("ID", "study", "year"), remove=FALSE) %>%
  st_as_sf(coords=c("X", "Y"), crs=32616)
st_write(pig_sf, "data/location_data/filtered_pig_locations.shp")

