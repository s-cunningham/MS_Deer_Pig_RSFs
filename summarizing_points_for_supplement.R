library(tidyverse)

## Location data wtih dates
deer <- read_csv("output/segmented_deer.csv")

# how many days in each segment?
ndays <- deer |> mutate(day=format(date, "%Y-%m-%d")) |> group_by(key) |> reframe(ndays=length(unique(day)))
short <- ndays |> filter(ndays <30)

# Removed anything with less than 3 weeks of data
deer <- deer |> filter(!(key %in% short$key))

# Summarize number of days, start & end dates
deer <- deer |>
  rename(timestamp=date) |>
  separate("timestamp", into=c("date", "time"), sep=" ") |>
  mutate(date = as.Date(date)) |>
  group_by(id) |> 
  reframe(ndays=length(unique(date)),
          start=min(date),
          end=max(date))

# add column for study
deer <- deer |>
  separate("id", into=c("ID", "study"), sep="_")


# write_csv(pigs, "output/pigs_filtered_segmented.csv")
pigs <- read_csv("output/pigs_filtered_segmented.csv")

# Drop extra forays from 19212_Delta_3
pigs <- pigs |> 
  filter(!(key=="19212_Delta_3" & y<3752500)) 

# Summarize number of days, start & end dates
pigs <- pigs |>
  mutate(date = substring(date, 1, 10)) |>
  mutate(date = as.Date(date)) |>
  group_by(id) |> 
  reframe(ndays=length(unique(date)),
          start=min(date),
          end=max(date))

# add column for study
pigs <- pigs |>
  separate("id", into=c("ID", "study"), sep="_")

## Read 
d_ua <- read_csv("output/deer_used_avail_locations.csv") |>
  group_by(id, case) |>
  count() |>
  pivot_wider(names_from="case", values_from="n") |>
  rename(n_used=`1`, n_avail=`0`)|>
  separate("id", into=c("ID", "study"), sep="_")
p_ua <- read_csv("output/pigs_used_avail_locations.csv") |>
  group_by(id, case) |>
  count() |>
  pivot_wider(names_from="case", values_from="n") |>
  rename(n_used=`1`, n_avail=`0`)|>
  separate("id", into=c("ID", "study"), sep="_")

## Combine

# within species
deer <- deer |>
  left_join(d_ua, by=c("ID", "study")) |>
  # add column for species
  mutate(species = "White-tailed Deer") |>
  # rearrange columns
  select(species, ID, study, ndays:n_used)

pigs <- pigs |> 
  left_join(p_ua, by=c("ID", "study")) |>
  # add column for species
  mutate(species = "Wild Pigs") |>
  # rearrange columns
  select(species, ID, study, ndays:n_used)

# Combine
dat <- bind_rows(deer, pigs)


# write file
write_csv(dat, "output/gps_data_rsf_summary.csv")
