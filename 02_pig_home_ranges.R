library(tidyverse)
library(amt)
library(sf)
library(adehabitatLT)
library(ctmm)
library(patchwork)

## Set ggplot theme
theme_set(theme_bw())

## Write function to subset and convert to ltraj
# Will need this for segmenting migration
lv_func <- function(x, id_, lmin=30) {
  
  # Subset by ID
  temp <- x %>% 
    # Filter by ID
    filter(id==id_) %>%
    # Select only columns we need to converting back to ltraj
    dplyr::select(id:date) %>%
    # Convert to data.frame
    as.data.frame() %>%
    # Make date a POSIX-formatted column
    mutate(date = as.POSIXct(date, tz="America/Denver"))
  
  # Convert to ltraj
  ltr <- as.ltraj(xy=temp[,c("x","y")], date=temp$date, id=temp$id)
  
  # Run lavielle
  lv <- lavielle(ltr, Lmin=lmin, which="R2n", Kmax=10, plotit=TRUE)
  
  k_try <- chooseseg(lv, S=0.75, output="opt", draw=TRUE)
  
  print("Optimal breaks by chooseseg():")
  print(k_try)
  
  # Plot and segment with select value of breaks
  fp <- findpath(lv, k_try)
  
  # Take input from the console on number of breaks
  print("Enter number of breaks:")
  k <- scan(n = 1)
  
  # Plot and segment with select value of breaks
  fp <- findpath(lv, k)
  
  # Convert to data frame
  df <- ld(fp)
  # remove some columns we don't need (i.e. match pigs_res)
  df <- df %>% as_tibble() %>%
    # reorder columns
    dplyr::select(id, x:rel.angle, burst) %>%
    # Remove extra characters ("Segment.") so that we only have a number in the burst column
    mutate(burst=gsub("Segment.", "", burst)) %>%
    # Combine ID with the segment number
    unite(key, c("id", "burst"), sep="_", remove=FALSE) %>%
    # Reorder columns
    dplyr::select(id, key, x:rel.angle)
  
  # Return segmented data
  return(df)
}

#### Read in pigs location data ####
pigs <- read_csv("data/location_data/pigs_filtered.csv") %>%
  mutate(month=as.numeric(format(timestamp, "%m")),
         year=as.numeric(format(timestamp, "%Y"))) %>%
  dplyr::select(PigID:study,month,year,timestamp,X,Y)# %>%
  # filter(PigID != "35493_Noxubee_2022")

pigs %>% group_by(study) %>% reframe(range=range(timestamp))

# Create track of all individuals
pigs_tr <- make_track(pigs, X, Y, timestamp, id=PigID, crs=32616)

# Check sampling rate of all individuals
trk_info <- pigs_tr %>% summarize_sampling_rate_many(c("id"), time_unit="min")

# List individual IDs
un.id <- unique(pigs$PigID)

## Loop over each ID and: 
# 1. Resample to every 4 hours
# 2. Calculate step lengths (and find median)
# 3. Check patchiness of location data
# 4. Determine whether individuals are moving ranges throughout the year
# 5. Build ctmm home range
# 6. Extract available points from home range


#### 1. and 2. Run loop to resample and calculate step lengths ####
pigs_res <- data.frame()
dist <- numeric()
tr_sum <- data.frame()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- pigs_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  s <- dat %>% summarize_sampling_rate_many(c("id"), time_unit="hour")
  tr_sum <- bind_rows(tr_sum, s)
  
  # Create column for only dates, count how many unique days
  dat <- dat %>% mutate(date=as.Date(format(t_, "%Y-%m-%d"))) 
  
  # How many locations on each day?
  nobs <- dat %>% group_by(date) %>% count() %>% filter(n>=5)
  un.days <- unique(nobs$date)
  
  # Skip to next individual if less than 30 days of data
  if (length(un.days) >= 30) {
    
    # Convert to ltraj
    df <- as.data.frame(dat)
    df <- as.ltraj(xy=df[,c("x_","y_")], date=df$t_, id=df$id)
    
    ## Regularize trajectories, round values
    refda <- min(df[[1]]$date)
    df_NA <- setNA(df, refda, 240, units="min") 
    df <- sett0(df_NA, refda, 240, units="min")
    
    # Save distance between locations
    dist <- c(dist, df[[1]]$dist)
    
    df <- ld(df)
    df <- df %>% as_tibble() %>%
      dplyr::select(id, x:rel.angle)
    
    # Save resampled data
    pigs_res <- bind_rows(pigs_res, df)
  } # if statement
} # i

tr_sum <- tr_sum %>% as_tibble()
tr_sum2 <- tr_sum %>% separate(id, into=c("pigid", "study"), sep="_", remove=FALSE)

#### 4. Determine which animals need to be split ####
pigs_res <- as_tibble(pigs_res)
un.id <- unique(pigs_res$id)

for (i in 1:length(un.id)) {
  
  temp <- pigs_res %>% filter(id==un.id[i])
  
  p1 <- ggplot(temp, aes(x=x, y=y, color=date)) +
    geom_path() + geom_point() +
    ggtitle(un.id[i])
  
  p2 <- ggplot(temp, aes(x=date, y=R2n, color=date)) +
    geom_path()
  
  print(p1 / p2)
  
}
## Look through plots and manually list which IDs probably will need to be segmented, as well as those that might need to be filtered by dates to remove patchy data

to_segment <- c("377123_Delta", "19212_Delta", "19219_Delta", "19223_Delta",
                "152316_Northern", "152317_Northern", "152318_Northern", "26620_Noxubee", "35490_Noxubee", 
                "35493_Noxubee", "30251_Eastern", "26641_Eastern", "26630_Eastern", "26622_Eastern") # "26626_Noxubee", 

# Manually clean up erroneous points
pigs <- pigs_res %>% filter(id!="35489_Noxubee" | (id=="35489_Noxubee" & y<3684000 & x<330000))
pigs <- pigs_res %>% filter(id!="30252_Noxubee" | (id=="30252_Noxubee" & y>3683500 & x>329000))
pigs <- pigs_res %>% filter(id!="152321_Northern" | (id=="152321_Northern" & y<3730000))
pigs <- pigs_res %>% filter(id!="152317_Northern" | (id=="152317_Northern" & x<327500))
pigs <- pigs_res %>% filter(id!="19208_Delta" | (id=="19208_Delta" & x>175000))

# Subset data without pigs that need to be segmented (so that we can add the segments onto it)
pigs <- pigs_res %>% 
  # Filter out NAs
  filter(!is.na(x)) %>%
  # Remove IDs that will be segmented
  filter(!(id %in% to_segment)) %>%
  # Create a 'key' column to match the output from the Lavielle function
  mutate(burst=1) %>%
  unite("key", c("id", "burst"), sep="_", remove=FALSE) %>%
  dplyr::select(id, key, x:rel.angle)

## Need to do this manually, and make sure to add the segmented data to 'pigs'
s <- lv_func(pigs_res, "26622_Eastern", 20)  

# Add segmented data to full dataset
pigs <- bind_rows(pigs, s)

# 377123_Delta : 3
# 19212_Delta : 8 #6
# 19219_Delta : 4
# 19223_Delta : 2
# 152316_Northern : 5
# 152317_Northern : 2
# 152318_Northern : 5
# 26620_Noxubee : 2
# 35490_Noxubee : 2
# 35493_Noxubee : 6
# 30251_Eastern : 3
# 26641_Eastern : 2
# 26630_Eastern : 5
# 26622_Eastern : 2

# Check if there are any IDs in the to_segment vector that are not in the pigs tibble after segmenting and binding
to_segment[!(to_segment %in% pigs$id)]

## Save
write_csv(pigs, "output/segmented_pigs.csv")
# pigs <- read_csv("output/segmented_pigs.csv")

# how many days in each segment?
ndays <- pigs %>% mutate(day=format(date, "%Y-%m-%d")) %>% group_by(key) %>% reframe(ndays=length(unique(day)))
short <- ndays %>% filter(ndays<22)

# Removed anything with less than 3 weeks of data
pigs <- pigs %>% filter(!(key %in% short$key))

#### 5. Loop over individuals and calculate AKDE home range ####

## Set up data
# Convert to lat/long (WGS84)
pigs_sf <- pigs %>% 
  # Remove rows without coordinates
  filter(!is.na(x)) %>% 
  # Convert to sf object
  st_as_sf(coords=c("x", "y"), crs=32616) %>%
  # Reproject
  st_transform(crs=4326)
# save lat/lon
ll <- pigs_sf %>% st_coordinates() %>% as_tibble() %>% 
  rename(location.long=X, location.lat=Y)

pigs <- pigs %>% 
  filter(!is.na(x)) %>% 
  bind_cols(ll)

dat <- pigs %>%      
  # need the lat long, etc. in a specific order 
  dplyr::select(key, location.long, location.lat, date) %>% 
  # also need columns to have a specific name
  rename(individual.local.identifier=key, timestamp=date)

# Create telemetry object
dat <- as.telemetry(dat, timeformat="%Y-%m-%d %H:%M:%S", timezone="America/Chicago")

## Loop over all individuals / individual segments
ids <- unique(pigs$key)

# Create list to save variograms and AKDEs
out <- list()

# Create empty object to store 
size <- data.frame()

# Set progress bar of sanity because this takes long - set it to the number of ids in the dataset (min = 0, max = max number of ids)
pb <- txtProgressBar(min = 0, max = length(ids), style = 3)

# Run loop to fit AKDEs
for (i in 1:length(ids)) {
  
  try({ # Some will likely fail, so prevent the loop from crashing 
    setTxtProgressBar(pb, i)
    
    # Next three lines are to fit variogram and select the best movement model to fit to the data
    var <- variogram(dat[[i]], CI="Gauss") 
    var2 <- variogram.fit(var,  name="GuessBM", interactive=F)
    
    # Automated model selection
    tt <- ctmm.select(dat[[i]], var2, IC="BIC") 
    
    # Fit AKDE
    UD <- akde(dat[[i]], tt, weights=TRUE) 
    
    # Save area of AKDE
    sqkm <- as.data.frame(summary(UD)$CI)
    size <- bind_rows(size, sqkm)
    
    # Plot with points
    plot(dat[[i]], UD)
    
    # Save the variograms and akde in a list  
    out[[i]] <- list(dat[[i]]@info$identity, var, var2, UD) 
  })
}
saveRDS(out, "results/home_ranges/raw_pigs_AKDEs.rds")
# out <- readRDS("results/home_ranges/raw_pigs_AKDEs.rds")

# Take just the AKDE
akde <- lapply(out, function(x) x[[4]]) 

# Save the IDs
ids <- lapply(out, function(x) x[[1]])
ids <- unlist(ids)

# Add IDs to size df
size$id <- ids

# Fix units and reorganzie
size <- size %>%
  # Create a column for the rowname
  rownames_to_column(var="unit") %>%
  # Convert to tibble
  as_tibble() %>%
  # Extract the unit from the rowname
  mutate(unit=str_extract(unit, pattern = "(?<=\\().*(?=\\))")) %>%
  # Convert ha to sq km
  mutate(low=if_else(unit=="hectares",low/100,low),
          est=if_else(unit=="hectares",est/100,est),
          high=if_else(unit=="hectares",high/100,high)) %>%
  # Reorder columns
  dplyr::select(id, low:high) 
# Note: all AKDEs are now in sqare kilometers

write_csv(size, "results/home_ranges/pigs_AKDE_sqkm.csv")

# Convert UD to sf object (multipolygons)
poly <- lapply(akde, as.sf, level.UD=0.99)

# Reproject to UTM 16N
poly <- lapply(poly, st_transform, crs=32616)

# combine into single df instead of list
ud <- do.call(rbind, poly)

# Clean up full set of polygons
ud <- ud %>% 
  # Separate name into important parts (key, 99%, and estimate vs CI)
  separate(1, into=c("key", "level", "which"), sep=" ") %>%
  # Remove rowname
  remove_rownames() %>%
  # Remove the level column, we know they are all 99%
  dplyr::select(-level) %>%
  # Now filter by rows to just the estimate (drop confidence intervals)
  filter(which=="est") %>%
  # Separate key into different pieces
  separate(1, into=c("id", "study", "burst"), sep="_") %>%
  # put id and study back together
  unite("id", 1:2, sep="_")

# Write shapefile
st_write(ud, "output/pigs_AKDE_est.shp", append=FALSE)
  
#### 6. Extract available points from home range ####
ud <- st_read("output/pigs_AKDE_est.shp")
library(terra)

## read a template raster (30 x 30, doesn't really matter which as long as projection is correct & has cells)
cdl <- rast("data/landscape_data/CDL2023_MS50km_30m_nlcdproj.tif")

# project to match raster
ud <- st_transform(ud, st_crs(cdl))

# Create empty tibble for points
avail <- tibble()

# Create empty list to store grid cells
cell_list <- list()

## Loop over each individual
for (i in 1:nrow(ud)) {
  
  # Clip raster (mask & crop extent) to AKDE
  temp <- mask(cdl, vect(ud$geometry[i]))
  temp <- crop(cdl, vect(ud$geometry[i]))
  
  # Convert to points
  pts <- as.points(temp, values=FALSE, na.rm=TRUE, na.all=FALSE)
  
  # Clip points to AKDE
  pts <- crop(pts, vect(ud$geometry[i]))
  
  ## Create empty grid from temp rasters
  # Determine number of cells
  cells <- dim(temp)[1] * dim(temp)[2]
  # create a unique ID for each grid cell
  temp$cell_no <- 1:cells
  # convert to polygons
  pols <- as.polygons(temp[["cell_no"]])
  # Extract the cell numbers that have a point in them
  grid <- extract(pols, pts)[,-1]
  # subset the polygon grid by only the cells that have points in them
  grid <- pols[grid]
  
  # Convert to sf
  grid <- st_as_sf(grid)
  
  # Save as shapefile
  filename <- paste0("output/pigs_avail_cells/", ud$id[i], "_", ud$burst[i], "-avail.shp")
  st_write(grid, filename, append=FALSE)
  
  # convert back to sf
  pts <- st_as_sf(pts)
  
  # convert to tibble and add identifying info, as well as a 0 for available
  pts <- pts %>% 
    st_coordinates() %>% 
    as_tibble() %>%
    mutate(id = ud$id[i],
           burst = ud$burst[i],
           case = 0)
  
  # Add to avail data frame
  avail <- bind_rows(avail, pts)
}
# Save available points
write_csv(avail, "output/pigs_avail_pts.csv")

## Now on to the used locations
pigs_sf <- pigs %>% 
  # Remove rows without coordinates
  filter(!is.na(x)) %>% 
  # Convert to sf object
  st_as_sf(coords=c("x", "y"), crs=32616) %>%
  # remove some columns
  dplyr::select(id, key, geometry) %>%
  # reproject to ud
  st_transform(crs=st_crs(ud))

# Create empty tibble to save new dataset
used <- tibble()
for (i in 1:nrow(ud)) {
  
  # Build key for each individual
  k <- ud$id[i]
  b <- ud$burst[i]
  udkey <- paste0(k, "_", b)
  
  # Filter and convert to SpatVector
  poly <- vect(ud$geometry[i])
  
  # filter to single ID and convert to SpatVector
  temp <- pigs_sf %>% filter(key==udkey) 
  temp <- vect(temp)
  
  # Clip points to AKDE
  temp <- crop(temp, poly)
  
  ## Convert back to tibble
  xy <- temp %>% st_as_sf() %>% st_coordinates() %>% as_tibble()
  
  # remove geometry column
  temp <- temp %>% st_as_sf() %>% st_drop_geometry()
  
  # add xy coordinates back on
  temp <- bind_cols(temp, xy)
  
  # Add a 1 to indicate used points
  temp$case <- 1
  
  # Add to data frame
  used <- bind_rows(used, temp)
}
write_csv(used, "output/pigs_used_locs.csv")

# Organize so used data matches available data
used <- used %>%
  # split key
  separate("key", into=c("rm1", "rm2", "burst"), sep="_") %>%
  # drop extra columns
  dplyr::select(X, Y, id, burst, case)

## Combine used and available points
dat <- bind_rows(used, avail)

# save data
write_csv(dat, "output/pigs_used_avail_locations.csv")

