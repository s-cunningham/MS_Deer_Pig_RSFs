library(tidyverse)
library(sf)

## Mississippi shapefile
ms <- st_read("data/landscape_data/mississippi_ACEA.shp") 
acea <- st_crs(ms) #crs=9822
ms <- ms %>% st_transform(crs=4326)

# Greater Starkville
strk <- st_read("data/landscape_data/starkville.shp") %>% st_transform(crs=4326)

ms <- rmapshaper::ms_erase(ms, strk)  # Will remove extraneous points in/around Starkville

#### Deer ####
# Read cleaned data
deer1 <- read_csv("data/location_data/Cleaned_Central_DeerGPS_data.csv") %>% 
  rename(id=ID, Temp=CollarTemp, timestamp=date.time) %>% 
  select(id, study, timestamp, Long, Lat, DOP) %>%
  mutate(study="Central")
deer2 <- read_csv("data/location_data/Cleaned_Delta_DeerGPS_data.csv") %>% 
  rename(id=collar_id, timestamp=datetime) %>% 
  select(id, study, timestamp, Long, Lat, DOP) %>%
  mutate(study="Delta")
deer3 <- read_csv("data/location_data/Cleaned_North_DeerGPS_data.csv") %>% 
  rename(timestamp=datetime) %>% 
  select(id, study, timestamp, Long, Lat, DOP) %>%
  mutate(study="North")

deer <- bind_rows(deer1, deer2, deer3)

# Remove locations that are outside of mississippi
deer <- st_as_sf(deer, coords=c("Long", "Lat"), crs=4326) 

deer <- st_intersection(deer, ms)
deer <- deer %>% st_transform(crs=32616)

#### Save shapefiles of each deer study ####

# deersf1 <- deer %>% filter(study=="Central_MS") %>% select(id:DOP) %>% mutate(year=as.numeric(format(timestamp, "%Y")))
# deersf1_1 <- deersf1 %>% filter(year==2017)
# deersf1_2 <- deersf1 %>% filter(year==2018)
# deersf1_3 <- deersf1 %>% filter(year==2019)
# deersf1_4 <- deersf1 %>% filter(year==2020)  
# 
# st_write(deersf1_1, "data/location_data/Central_DeerGPS2017.shp", append=FALSE)
# st_write(deersf1_2, "data/location_data/Central_DeerGPS2018.shp", append=FALSE)
# st_write(deersf1_3, "data/location_data/Central_DeerGPS2019.shp", append=FALSE)
# st_write(deersf1_4, "data/location_data/Central_DeerGPS2020.shp", append=FALSE)
# 
# deersf2 <- deer %>% filter(study=="MS_Delta") %>% select(id:DOP)
# st_write(deersf2, "data/location_data/Delta_DeerGPS.shp", append=FALSE)
# 
# deersf3 <- deer %>% filter(study=="North_MS") %>% select(id:DOP)
# st_write(deersf3, "data/location_data/North_DeerGPS.shp", append=FALSE)

#### Carry on ####

xy <- st_coordinates(deer) %>% as_tibble()
deer <- st_drop_geometry(deer)

deer$X <- xy$X
deer$Y <- xy$Y

deer <- deer %>% unite("DeerID", 1:2, sep="_", remove=FALSE)

# Looks like they have already been filtered to DOP<10
deer <- deer %>% select(DeerID, id, study, timestamp, X, Y, DOP)

npts <- deer %>% group_by(DeerID) %>% count()
hist(npts$n)
sum(npts$n<1000)

# Keep deer with >3 months of data
deer_time <- deer %>% group_by(DeerID) %>%
                reframe(ts=range(timestamp)) %>%
                mutate(time=rep(c("start","end"), 91)) %>%
                pivot_wider(names_from="time", values_from="ts") %>%
                mutate(duration=end-start,
                       months=as.numeric(duration/30))

keep <- deer_time %>% filter(months > 3) %>% select(DeerID) %>% as.vector()
keep <- unlist(unname(keep))

deer <- deer %>% filter(DeerID %in% keep)


write_csv(deer, "data/location_data/deer_filtered.csv")



#### Pigs ####
# Read GPS data
pigs1 <- read_csv("data/location_data/Cleaned_DeltaMS_PigGPS_data.csv") %>% select(-Sex, -Temp) %>% rename(id=collar_id) %>% mutate(study="Delta")
pigs2 <- read_csv("data/location_data/Cleaned_NorthMS_PigGPS_data.csv") %>% select(-Altitude, -Temp, -DOP) %>% mutate(study="Northern")
pigs3 <- read_csv("data/location_data/Cleaned_Noxubee_PigGPS_data.csv") %>% select(-Temp) %>% rename(id=CollarID, datetime=DateTime, Long=x, Lat=y) %>% mutate(study="Noxubee")
pigs4 <- read_csv("data/location_data/Cleaned_EasternMS_PigGPS_data.csv") %>% select(-Altitude, -Temp, -DOP) %>% rename(id=collar_id) %>% mutate(study="Eastern")

# Combine
pigs <- bind_rows(pigs1, pigs2, pigs3, pigs4)

# Remove locations that are outside of mississippi
pigs <- st_as_sf(pigs, coords=c("Long", "Lat"), crs=4326) 

# st_write(pigs, "data/location_data/all_pig_points.shp")

pigs <- st_intersection(pigs, ms)
pigs <- pigs %>% st_transform(crs=acea)






