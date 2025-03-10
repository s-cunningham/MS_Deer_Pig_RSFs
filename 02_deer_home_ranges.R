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
lv_func <- function(x, id_) {
  
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
  lv <- lavielle(ltr, Lmin=30, which="R2n", Kmax=10, plotit=TRUE)
  
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
  # remove some columns we don't need (i.e. match deer_res)
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

## Create template raster for generating random points


#### Deer ####
# Read in deer location data
deer <- read_csv("data/location_data/deer_filtered.csv") %>%
          mutate(timestamp=as_datetime(timestamp, tz="America/Chicago")) %>%
          mutate(month=as.numeric(format(timestamp, "%m")),
                 year=as.numeric(format(timestamp, "%Y"))) %>%
          dplyr::select(DeerID:study,month,year,timestamp,X,Y) # %>%
          # unite("DeerID", c(1,5), sep="_") # separate by year

deer %>% group_by(study) %>% reframe(range=range(timestamp))

# Create track of all individuals
deer_tr <- make_track(deer, X, Y, timestamp, id=DeerID, crs=32616)

# Check sampling rate of all individuals
trk_info <- deer_tr %>% summarize_sampling_rate_many(c("id"), time_unit="min")

# List individual IDs
un.id <- unique(deer$DeerID)

## Loop over each ID and: 
# 1. Resample to every 4 hours
# 2. Calculate step lengths (and find median)
# 3. Check patchiness of location data
# 4. Determine whether individuals are moving ranges throughout the year
# 5. Build ctmm home range
# 6. Extract available points from home range


#### 1. and 2. Run loop to resample and calculate step lengths ####
deer_res <- data.frame()
dist <- numeric()
tr_sum <- data.frame()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- deer_tr %>% filter(id==un.id[i])
  
  # Filter all locations to ~4 hours
  dat <- track_resample(dat, rate=hours(4), tolerance=minutes(5))
  s <- dat %>% summarize_sampling_rate_many(c("id"), time_unit="hour")
  tr_sum <- bind_rows(tr_sum, s)
  
  # Create column for only dates, count how many unique days
  dat <- dat %>% mutate(date=as.Date(format(t_, "%Y-%m-%d"))) 

  # How many locations on each day?
  nobs <- dat %>% group_by(date) %>% count() %>% filter(n>=5)
  un.days <- unique(nobs$date)

  # Skip to next individual if less than 90 days of data
  if (length(un.days) >= 90) {
    
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
    deer_res <- bind_rows(deer_res, df)
  } # if statement
} # i

tr_sum <- tr_sum %>% as_tibble()


#### 3. Determine which animals need to be split ####
deer_res <- as_tibble(deer_res)
un.id <- unique(deer_res$id)

for (i in 1:length(un.id)) {
  
  temp <- deer_res %>% filter(id==un.id[i])

  p1 <- ggplot(temp, aes(x=x, y=y, color=date)) +
    geom_path() + geom_point() +
    ggtitle(un.id[i])
  
  p2 <- ggplot(temp, aes(x=date, y=R2n, color=date)) +
         geom_path()
  
  print(p1 / p2)
  
}
## Look through plots and manually list which IDs probably will need to be segmented, as well as those that might need to be filtered by dates to remove patchy data

## IDs that need trimming
temp <- deer_res %>% filter(id=="81864_Delta")
ggplot(temp, aes(x=x, y=y, color=date)) +
  geom_path() + geom_point() 
ggplot(temp, aes(x=date, y=R2n, color=date)) +
  geom_vline(xintercept=as_datetime("2018-08-31 00:00:00")) +
  geom_path()

# 152311_North
deer_res <- deer_res %>% filter((id!="152311_North") | (id=="152311_North" & date <= "2023-12-31 23:59:59"))

# 87924_Delta
deer_res <- deer_res %>% filter((id!="87924_Delta") | (id=="87924_Delta" & (date>="2021-06-30 23:59:59" & date<="2022-09-10 23:59:59")))

# 87899_Delta
deer_res <- deer_res %>% filter((id!="87899_Delta") | (id=="87899_Delta" & date<="2022-12-01 00:00:00"))

# 81864_Delta
deer_res <- deer_res %>% filter((id!="81864_Delta") | (id=="81864_Delta" & (date>="2021-08-25 00:00:00" & date<="2022-02-01 00:00:00")))

# 259_Central
deer_res <- deer_res %>% filter((id!="259_Central") | (id=="259_Central" & date <= "2017-09-15 00:00:00"))

# 252_Central
deer_res <- deer_res %>% filter((id!="252_Central") | (id=="252_Central" & date <= "2017-05-01 00:00:00"))

# 27_Central
deer_res <- deer_res %>% filter((id!="27_Central") | (id=="27_Central" & (date>="2018-08-31 00:00:00" & date<="2019-03-27 23:59:59")))

## List of IDs that will likely need to be segmented
to_segment <- c("152312_North", "152308_North", "87924_Delta", "87919_Delta", "81711_Delta", "97_Central", "90_Central", 
                "47_Central",  "344_Central", "339_Central", "333_Central", "297_Central", "295_Central", "293_Central", 
                "281_Central", "28_Central", "272_Central", "27_Central", "25_Central", "20_Central", "100_Central")

# Subset data without deer that need to be segmented (so that we can add the segments onto it)
deer <- deer_res %>% 
          # Filter out NAs
          filter(!is.na(x)) %>%
          # Remove IDs that will be segmented
          filter(!(id %in% to_segment)) %>%
          # Create a 'key' column to match the output from the Lavielle function
          mutate(burst=1) %>%
          unite("key", c("id", "burst"), sep="_", remove=FALSE) %>%
          dplyr::select(id, key, x:rel.angle)

## try 81711_Delta, 27_Central, and 87924_Delta again
deer <- deer %>% filter(id!="27_Central") #& id!="87924_Delta", id!="81711_Delta" &

## Need to do this manually, and make sure to add the segmented data to 'deer'
s <- lv_func(deer_res, "27_Central")  

# keep only some segments of 87924_Delta
# s <- s %>% filter(key=="87924_Delta_1" | key=="87924_Delta_3" | key=="87924_Delta_5"| key=="87924_Delta_6" | key=="87924_Delta_7")

# Add segmented data to full dataset
deer <- bind_rows(deer, s)

# 152312_North : 4
# 152308_North : 4
# 87924_Delta : 8, keep only 1, 3, 5, 6, 7 
# 87919_Delta : 3
# 81711_Delta : 3 # might need to revise (see note below)
# 97_Central : 3
# 90_Central : 2
# 47_Central : 3
# 344_Central : 7
# 339_Central : 3
# 333_Central : 2
# 297_Central : 8
# 295_Central : 3
# 293_Central : 5
# 281_Central : 3
# 28_Central : 3
# 272_Central : 6
# 27_Central : 4
# 20_Central : 5
# 100_Central : 3
# Run 81711_Delta with 8 breaks...or, stay with the optimal 5 and subset anything that doesn't have 30 or 60 days? What's the average HR crossing time?

# Check if there are any IDs in the to_segment vector that are not in the deer tibble after segmenting and binding
to_segment[!(to_segment %in% deer$id)]

## Save
write_csv(deer, "output/segmented_deer.csv")
deer <- read_csv("output/segmented_deer.csv")

# how many days in each segment?
ndays <- deer %>% mutate(day=format(date, "%Y-%m-%d")) %>% group_by(key) %>% reframe(ndays=length(unique(day)))
short <- ndays %>% filter(ndays <22)

# Removed anything with less than 3 weeks of data
deer <- deer %>% filter(!(key %in% short$key))


#### 4. Loop over individuals and calculate AKDE home range ####

## Set up data
# Convert to lat/long (WGS84)
deer_sf <- deer %>% 
              # Remove rows without coordinates
              filter(!is.na(x)) %>% 
              # Convert to sf object
              st_as_sf(coords=c("x", "y"), crs=32616) %>%
              # Reproject
              st_transform(crs=4326)
# save lat/lon
ll <- deer_sf %>% st_coordinates() %>% as_tibble() %>% 
          rename(location.long=X, location.lat=Y)

deer <- deer %>% 
              filter(!is.na(x)) %>% 
              bind_cols(ll)

dat <- deer %>%      
            # need the lat long, etc. in a specific order 
            dplyr::select(key, location.long, location.lat, date) %>% 
            # also need columns to have a specific name
            rename(individual.local.identifier=key, timestamp=date)

# Create telemetry object
dat <- as.telemetry(dat, timeformat="%Y-%m-%d %H:%M:%S", timezone="America/Chicago")

## Loop over all individuals / individual segments
ids <- unique(deer$key)

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
    var <- variogram(dat[[i]]) 
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
saveRDS(out, "results/raw_deer_AKDEs.rds")

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

write_csv(size, "results/deer_AKDE_sqkm.csv")

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
st_write(ud, "output/deer_AKDE_est.shp")

# ud <- st_read("output/deer_AKDE_est.shp")






# Can plot each UD with:
plot(out[[1]][4])


## Loop over individuals and generate available points
avail <- data.frame()
for (i in 1:length(un.id)) {
  
  
}

