library(tidyverse)
library(sf)
library(amt)


deer <- read_csv("data/location_data/deer_filtered.csv") %>%
          mutate(month=as.numeric(format(timestamp, "%m")),
                 year=as.numeric(format(timestamp, "%Y"))) %>%
          select(DeerID:study,month,year,timestamp,X,Y)


delta87924 <- deer %>% filter(DeerID=="87924_Delta") %>%
  mutate(minutes=c(NA, diff(timestamp)),
         hours=minutes/60) 
  

ggplot(delta87924) +
  geom_path(aes(x=timestamp, y=X)) +
  geom_point(aes(x=timestamp, y=X, color=hours)) +
  theme_bw()


un.id <- unique(deer$DeerID)
area <- data.frame()
for (i in 1:length(un.id)) {
  
  # Filter to ID
  dat <- deer %>% filter(DeerID==un.id[i])
  
  # Make amt track object
  tr <- make_track(dat, X, Y, timestamp)
  
  # # Make template raster
  trast <- make_trast(tr, res=10)

  # Create HR
  # hr_deer <- hr_kde(tr, trast=trast, level=0.99)
  hr_deer <- hr_mcp(tr, trast=trast)

  tr_mcp <- hr_deer$mcp

  area[i,1] <- hr_deer$mcp$area/1000/1000

  # plot(st_geometry())

  # Add buffer?? See Darlington et al. 2022

  # Generate random points
  # r1 <- random_points(tr, n = 5*nrow(dat)) #, presence = tr
  # plot(r1)
  
}

r1 <- r1 %>% as_tibble() %>% rename(x=x_, y=y_)

ggplot() +
  geom_point(data=r1, aes(x=x, y=y))





