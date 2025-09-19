library(tidyverse)
library(amt)
library(sf)
library(adehabitatLT)
library(ctmm)
library(patchwork)
library(tidyterra)

## Set ggplot theme
theme_set(theme_bw())

## Write function to subset and convert to ltraj
# Will need this for segmenting migration
lv_func <- function(x, id_, lmin=30) {
  
  # Subset by ID
  temp <- x |> 
    # Filter by ID
    filter(id==id_) |>
    # Select only columns we need to converting back to ltraj
    dplyr::select(id:date) |>
    # Convert to data.frame
    as.data.frame() |>
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
  df <- df |> as_tibble() |>
    # reorder columns
    dplyr::select(id, x:rel.angle, burst) |>
    # Remove extra characters ("Segment.") so that we only have a number in the burst column
    mutate(burst=gsub("Segment.", "", burst)) |>
    # Combine ID with the segment number
    unite(key, c("id", "burst"), sep="_", remove=FALSE) |>
    # Reorder columns
    dplyr::select(id, key, x:rel.angle)
  
  # Return segmented data
  return(df)
}

#### Read in pigs location data ####
pigs <- read_csv("data/location_data/pigs_filtered.csv") |>
  mutate(month=as.numeric(format(timestamp, "%m")),
         year=as.numeric(format(timestamp, "%Y"))) |>
  dplyr::select(PigID:study,month,year,timestamp,X,Y)# |>
# filter(PigID != "35493_Noxubee_2022")

pigs |> group_by(study) |> reframe(range=range(timestamp))

# Create track of all individuals
pigs_tr <- make_track(pigs, X, Y, timestamp, id=PigID, crs=32616)

# Check sampling rate of all individuals
trk_info <- pigs_tr |> summarize_sampling_rate_many(c("id"), time_unit="min")

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
  dat <- pigs_tr |> filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  s <- dat |> summarize_sampling_rate_many(c("id"), time_unit="hour")
  tr_sum <- bind_rows(tr_sum, s)
  
  # Create column for only dates, count how many unique days
  dat <- dat |> mutate(date=as.Date(format(t_, "%Y-%m-%d"))) 
  
  # How many locations on each day?
  nobs <- dat |> group_by(date) |> count() |> filter(n>=5)
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
    df <- df |> as_tibble() |>
      dplyr::select(id, x:rel.angle)
    
    # Save resampled data
    pigs_res <- bind_rows(pigs_res, df)
  } # if statement
} # i

tr_sum <- tr_sum |> as_tibble()
tr_sum2 <- tr_sum |> separate(id, into=c("pigid", "study"), sep="_", remove=FALSE)

pigs <- pigs_res |> 
  # Filter out NAs
  filter(!is.na(x)) 

## Check distribution of how far they move in 4 hours
quantile(pigs_res$dist, probs=c(0.9, 0.95, 0.99), na.rm=TRUE)

# Filter 99th percentiles 
pigs_res <- pigs_res |> filter(dist < (1415.3430*2))

#### 3. Determine which animals need to be split ####
pigs_res <- as_tibble(pigs_res)
un.id <- unique(pigs_res$id)
## Look through plots and manually list which IDs probably will need to be segmented, as well as those that might need to be filtered by dates to remove patchy data
for (i in 1:length(un.id)) {
  
  temp <- pigs_res |> filter(id==un.id[i])
  
  p1 <- ggplot(temp, aes(x=x, y=y, color=date)) +
    geom_path() + geom_point() +
    ggtitle(un.id[i])
  
  p2 <- ggplot(temp, aes(x=date, y=R2n, color=date)) +
    geom_path()
  
  print(p1 / p2)
  
}

pigs_res <- pigs_res |> filter(id!="26622_Eastern" | (id=="26622_Eastern" & x>339000 & y<3739000))
pigs_res <- pigs_res |> filter(id!="30252_Noxubee" | (id=="30252_Noxubee" & x>329000 & y>3683500))
pigs_res <- pigs_res |> filter(id!="26641_Eastern" | (id=="26641_Eastern" & x<340000))
pigs_res <- pigs_res |> filter(id!="26628_Noxubee" | (id=="26628_Noxubee" & x>326000 & x<330000))
pigs_res <- pigs_res |> filter(id!="152320_Northern" | (id=="152320_Northern" & R2n<5E06))
pigs_res <- pigs_res |> filter(id!="19223_Delta" | (id=="19223_Delta" & R2n<5E07))
pigs_res <- pigs_res |> filter(id!="19212_Delta" | (id=="19212_Delta" & R2n<2E08))
pigs_res <- pigs_res |> filter(id!="19210_Delta" | (id=="19210_Delta" & R2n<2E07))
pigs_res <- pigs_res |> filter(id!="19207_Delta" | (id=="19207_Delta" & x<182500))
pigs_res <- pigs_res |> filter(id!="377123_Delta" | (id=="377123_Delta" & x>204500))

pigs_res$burst <- 1
pigs_res <- pigs_res |> dplyr::select(id,burst,x:rel.angle)

## List of IDs that will likely need to be segmented
to_segment <- c("19219_Delta", "19212_Delta")

# Subset data without deer that need to be segmented (so that we can add the segments onto it)
pigs <- pigs_res |> 
  # Filter out NAs
  filter(!is.na(x)) |>
  # Remove IDs that will be segmented
  filter(!(id %in% to_segment)) |>
  # Create a 'key' column to match the output from the Lavielle function
  mutate(burst=1) |>
  unite("key", c("id", "burst"), sep="_", remove=FALSE) |>
  dplyr::select(id, key, x:rel.angle)

## Need to do this manually, and make sure to add the segmented data to 'deer'
s <- lv_func(pigs_res, "19212_Delta", 30)

# Add segmented data to full dataset
pigs <- bind_rows(pigs, s)

# how many days in each segment?
ndays <- pigs |> mutate(day=format(date, "%Y-%m-%d")) |> group_by(id) |> reframe(ndays=length(unique(day)))
short <- ndays |> filter(ndays<30)

# Removed anything with less than 3 weeks of data
pigs <- pigs |> filter(!(id %in% short$id))

# Ending up with some HRs that are likely much too large
temp <- pigs |> filter(key=="19212_Delta_3")
ggplot(temp, aes(x=x, y=y, color=date)) +
  geom_path() + geom_point() 
ggplot(temp, aes(x=date, y=R2n, color=date)) +
  geom_vline(xintercept=as_datetime("2018-08-31 00:00:00")) +
  geom_path()
ggplot(temp) +
  geom_histogram(aes(x=y))

# write_csv(pigs, "output/pigs_filtered_segmented.csv")
pigs <- read_csv("output/pigs_filtered_segmented.csv")
## Set up data
# Convert to lat/long (WGS84)
pigs_sf <- pigs |> 
  # Remove rows without coordinates
  filter(!is.na(x)) |> 
  # Convert to sf object
  st_as_sf(coords=c("x", "y"), crs=32616) |>
  # Reproject
  st_transform(crs=4326)
# save lat/lon
ll <- pigs_sf |> st_coordinates() |> as_tibble() |> 
  rename(location.long=X, location.lat=Y)

pigs <- pigs |> 
  filter(!is.na(x)) |> 
  bind_cols(ll)

# Drop extra forays from 19212_Delta_3
pigs <- pigs |> 
  filter(!(key=="19212_Delta_3" & y<3752500)) 

dat <- pigs |>      
  # need the lat long, etc. in a specific order 
  dplyr::select(key, location.long, location.lat, date) |> 
  # also need columns to have a specific name
  rename(individual.local.identifier=key, timestamp=date)

# Create telemetry object
dat <- as.telemetry(dat, timeformat="%Y-%m-%d %H:%M:%S", timezone="America/Chicago")

# Take ids from telemetry object
ids <- c()
for (i in 1:length(dat)) {
  ids <- c(ids, slot(dat[[i]],"info")$identity)
}

# Create list to save variograms and AKDEs
out <- list()

# Create empty object to store 
size <- data.frame()

# Set progress bar of sanity because this takes long - set it to the number of ids in the dataset (min = 0, max = max number of ids)
pb <- txtProgressBar(min = 0, max = length(ids), style = 3)

out <- readRDS("results/home_ranges/raw_pigs_AKDEs.rds")

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
out <- readRDS("results/home_ranges/raw_pigs_AKDEs.rds")

# Take just the AKDE
akde <- lapply(out, function(x) x[[4]]) 

# Save the IDs
ids <- lapply(out, function(x) x[[1]])
ids <- unlist(ids)

# Add IDs to size df
size$id <- ids

# Fix units and reorganzie
size <- size |>
  # Create a column for the rowname
  rownames_to_column(var="unit") |>
  # Convert to tibble
  as_tibble() |>
  # Extract the unit from the rowname
  mutate(unit=str_extract(unit, pattern = "(?<=\\().*(?=\\))")) |>
  # Convert ha to sq km
  mutate(low=if_else(unit=="hectares",low/100,low),
         est=if_else(unit=="hectares",est/100,est),
         high=if_else(unit=="hectares",high/100,high)) |>
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
ud <- ud |> 
  # Separate name into important parts (key, 99%, and estimate vs CI)
  separate(1, into=c("id", "level", "which"), sep=" ") |>
  # Remove rowname
  remove_rownames() |>
  # Remove the level column, we know they are all 99%
  dplyr::select(-level) |>
  # Now filter by rows to just the estimate (drop confidence intervals)
  filter(which=="est") #|>
# Separate id into different pieces
# separate(1, into=c("id", "study", "burst"), sep="_") |>
# put id and study back together
# unite("id", 1:2, sep="_")

# Write shapefile
st_write(ud, "output/pigs_AKDE_est.shp", append=FALSE)

#### 6. Extract available points from home range ####
ud <- st_read("output/pigs_AKDE_est.shp")
library(terra)

## read a template raster (30 x 30, doesn't really matter which as long as projection is correct & has cells)
cdl <- rast("data/landscape_data/CDL2023_MS50km_30m_nlcdproj.tif")

# aggregate cells
cdl_agg <- aggregate(cdl, fact=7, fun="which.max", cores=2)

# project to match raster
ud <- st_transform(ud, st_crs(cdl))

# Create empty tibble for points
avail <- tibble()

# Create empty list to store grid cells
cell_list <- list()

## Loop over each individual
for (i in 1:nrow(ud)) {
  
  # Clip raster (mask & crop extent) to AKDE
  temp <- mask(cdl, vect(ud$geometry[i])) #_agg 
  temp <- crop(cdl, vect(ud$geometry[i])) #_agg 
  
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
  filename <- paste0("output/pigs_avail_cells/", ud$id[i], "_", ud$burst[i], "-avail_30m.shp")
  # st_write(grid, filename, append=FALSE)
  
  # convert back to sf
  pts <- st_as_sf(pts)
  
  # convert to tibble and add identifying info, as well as a 0 for available
  pts <- pts |> 
    st_coordinates() |> 
    as_tibble() |>
    mutate(id = ud$id[i],
           burst = ud$burst[i],
           case = 0)
  
  # Add to avail data frame
  avail <- bind_rows(avail, pts)
}
# Save available points
write_csv(avail, "output/pigs_avail_pts_30m.csv")
avail <- read_csv("output/pigs_avail_pts_30m.csv")

## Now on to the used locations
pigs_sf <- pigs |> 
  # Remove rows without coordinates
  filter(!is.na(x)) |> 
  # Convert to sf object
  st_as_sf(coords=c("x", "y"), crs=32616) |>
  # remove some columns
  dplyr::select(key, geometry) |>
  # reproject to ud
  st_transform(crs=st_crs(ud))

# Create empty tibble to save new dataset
used <- tibble()
for (i in 1:nrow(ud)) {
  
  # Filter and convert to SpatVector
  poly <- vect(ud$geometry[i])
  
  # filter to single ID and convert to SpatVector
  temp <- pigs_sf |> filter(key==ud$id[i]) 
  temp <- vect(temp)
  
  # Clip points to AKDE
  temp <- crop(temp, poly)
  
  ## Convert back to tibble
  xy <- temp |> st_as_sf() |> st_coordinates() |> as_tibble()
  
  # remove geometry column
  temp <- temp |> st_as_sf() |> st_drop_geometry()
  
  # add xy coordinates back on
  temp <- bind_cols(temp, xy)
  
  # Add a 1 to indicate used points
  temp$case <- 1
  
  # Add to data frame
  used <- bind_rows(used, temp)
}
write_csv(used, "output/pigs_used_locs.csv")
used <- read_csv("output/pigs_used_locs.csv")

# Organize so used data matches available data
used <- used |>
  # split key
  # separate("key", into=c("rm1", "rm2", "burst"), sep="_") |>
  # drop extra columns
  dplyr::select(X, Y, key, case)

## Combine used and available points
avail <- avail |> rename(key=id)
dat <- bind_rows(used, avail)


# Check ratio of used:avail
ratios <- dat |> group_by(key, case) |> count() |>
  pivot_wider(names_from=case, values_from=n) |>
  rename(used=`1`, avail=`0`) |>
  mutate(n_avail_used=avail/used) 

# Drop IDs with HRs that are too big
dat <- dat |> filter(key!="19212_Delta_4")

low <- ratios |> filter(n_avail_used < 1)
un.low <- unique(low$key)

dat <- dat |> filter(key!="19211_Delta_1")

## Per-individual scaling of weights
# total weight per ID
W <- 5000

dat <- dat |>
  left_join(ratios,by="key") |>
  mutate(weight = if_else(case==1, 1, W/avail))

dat <- dat |>
  dplyr::select(-n_avail_used, -pts_needed_for_2x, -additional_pts_needed)

# save data
write_csv(dat, "output/pigs_used_avail_locations_indwt_30m.csv")


