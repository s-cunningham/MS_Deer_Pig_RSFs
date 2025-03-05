library(tidyverse)
library(amt)
library(sf)
library(adehabitatLT)
library(ctmm)
library(patchwork)

theme_set(theme_bw())
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
# 3. Determine whether individuals are moving ranges throughout the year
# 4. Build ctmm home range
# 5. Extract available points from home range


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

  # Skip to next individual if less than 1 month of data
  if (range(dat$t_)[2]-range(dat$t_)[1] >= 30) {
    
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


# test <- lavielle(df, Lmin=10, which="R2n", Kmax=5, plotit=TRUE)
# fp <- findpath(test, 5)
# 


#### 4. Loop over individuals and calculate AKDE home range ####

## Set up data
# Convert to lat/long (WGS84)
deer_sf <- deer_res %>% 
              # Remove rows without coordinates
              filter(!is.na(x)) %>% 
              # Convert to sf object
              st_as_sf(coords=c("x", "y"), crs=32616) %>%
              # Reproject
              st_transform(crs=4326)
# save lat/lon
ll <- deer_sf %>% st_coordinates() %>% as_tibble() %>% 
          rename(location.long=X, location.lat=Y)

deer_res <- deer_res %>% 
              filter(!is.na(x)) %>% 
              bind_cols(ll)

dat <- deer_res %>%      
            # need the lat long, etc. in a specific order 
            dplyr::select(id, location.long, location.lat, date) %>% 
            # also need columns to have a specific name
            rename(individual.local.identifier=id, timestamp=date)

# Create telemetry object
dat <- as.telemetry(dat, timeformat="%Y-%m-%d %H:%M:%S", timezone=attr(deer_res$timestamp, "tzone"), projection=NULL)

## Loop over all individuals

# Create list to save variograms and AKDEs
out <- list()

# Set progress bar of sanity because this takes long - set it to the number of ids in the dataset (min = 0, max = max number of ids)
pb <- txtProgressBar(min = 0, max = length(un.id), style = 3)

# Run loop to fit AKDEs
for (i in 1:length(un.id)) {

  try({ # Some will likely fail, so prevent the loop from crashing 
    setTxtProgressBar(pb, i)
    
    # Next three lines are to fit variogram and select the best movement model to fit to the data
    var <- variogram(dat[[i]]) 
    var2 <- variogram.fit(var,  name="GuessBM", interactive=F)
    # Automated model selection
    tt <- ctmm.select(dat[[i]], var2) 
    
    # Fit AKDE
    UD <- akde(dat[[i]], tt) 
    
    # Plot with points
    plot(dat[[i]], UD)
    
    # Save the variograms and akde in a list  
    out[[i]] <- list(dat[[i]]@info$identity, var, var2, UD) 
  })
}

# Take just the AKDE
akde <- lapply(out, function(x) x[[4]]) #the akde shape
ids <- lapply(out, function(x) x[[1]])
poly <- lapply(akde, function(x) try(SpatialPolygonsDataFrame.UD(x, level.UD=0.999, level=0.99)))



joined<- SpatialPolygons(lapply(poly, function(x){x@polygons[[3]]}), 
                         proj4string=CRS(akde[[3]]@info$projection))

plot(joined@polygons[[1]]@Polygons[[1]]@coords, type="l")
points(dat[[3]]$x, dat[[3]]$y)












# Can plot each UD with:
plot(out[[1]][4])


## Loop over individuals and generate available points
avail <- data.frame()
for (i in 1:length(un.id)) {
  
  
}




#### Pigs ####
# Read in deer location data
pigs <- read_csv("data/location_data/pigs_filtered.csv") %>%
  mutate(month=as.numeric(format(timestamp, "%m")),
         year=as.numeric(format(timestamp, "%Y"))) %>%
  dplyr::select(PigID:study,month,year,timestamp,X,Y) %>%
  unite("PigID", c(1,5), sep="_") %>% # separate by year
  filter(PigID != "35493_Noxubee_2022")
