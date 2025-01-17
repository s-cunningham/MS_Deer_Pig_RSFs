library(tidyverse)

# Function for calculating geometric mean
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Function for calculating vector of lambda values
lambda_calc <- function(x) {
  lambda <- numeric()
  for(y in 1:nrow(x)){
    lambda[y] <- x$Average[y+1]/x$Average[y]
  }
  return(lambda)
}

read_sim <- function(file) {
  
  x <- read_tsv(file, skip=14, n_max=11, show_col_types = FALSE) %>%
    # Change obnoxious column name
    rename(allcols=`Minimum    -1 S.D.    Average    +1 S.D.    Maximum`) %>%
    # because there are a different number of spaces, just want all whitespace to become a bar
    mutate(allcols = str_replace_all(allcols, "[[:space:]]+", "|")) %>% 
    # now separate 
    separate_wider_delim(allcols, 
                         names=c("Time", "Minimum", "lSD", "Average", "hSD", "Maximum"),
                         delim="|") %>%
    # convert caracter columns to numeric
    mutate_if(is.character, as.numeric)
  
  return(x)
}

#### PIG MODELS ####
# read files
file1 <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/pigs_core_sim.txt"
p_core <- read_sim(file)

# Calculate lambda
p_core_l <- lambda_calc(p_core)

# Calculate geometric mean lambda
gm_mean(p_core_l)


#### DEER MODELS ####
## Deer core patches, 14.3 deer/km2
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/deer_core_sim.txt"
d_core <- read_sim(file) 

# Calculate lambda
d_core_l <- lambda_calc(d_core)

# Calculate geometric mean lambda
gm_mean(d_core_l)

## Deer core patches, 22 deer/km2
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/deer_core_22deerkm2_sim.txt"
d_core22 <- read_sim(file) 

# Calculate lambda
d_core_l <- lambda_calc(d_core22)

# Calculate geometric mean lambda
gm_mean(d_core_l)


## Deer marginal habitat, 14.3 deer/km2
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/deer_marginal_sim.txt"
d_m <- read_sim(file) 

# Calculate lambda
d_m_l <- lambda_calc(d_m)

# Calculate geometric mean lambda
gm_mean(d_m_l)


## Deer marginal habitat, 22 deer/km2
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/deer_marginal_22deerkm2_sim.txt"
d_m22 <- read_sim(file) 

# Calculate lambda
d_m22_l <- lambda_calc(d_m22)

# Calculate geometric mean lambda
gm_mean(d_m22_l)


## Deer marginal habitat, 14.3 deer/km2
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/deer_highly_marginal_sim.txt"
d_hm <- read_sim(file) 

# Calculate lambda
d_hm_l <- lambda_calc(d_hm)

# Calculate geometric mean lambda
gm_mean(d_hm_l)


## Deer marginal habitat, 22 deer/km2
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/deer_highly_marginal_22deerkm2_sim.txt"
d_hm22 <- read_sim(file) 

# Calculate lambda
d_hm22_l <- lambda_calc(d_hm22)

# Calculate geometric mean lambda
gm_mean(d_hm22_l)


## Read comparison file
file <- "C:/Users/sac793/OneDrive - Mississippi State University/Documents/DeerPigProject/RAMASoutput/trajectory_tables/sens_deer_core_1-5.txt"
