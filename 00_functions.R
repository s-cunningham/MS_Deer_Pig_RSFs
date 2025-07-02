library(tidyverse)
library(R.utils)

## Ramas patch area (RAMAS Spatial Data)
ramas_area <- function(x) {
  
  # Use R.utils package to count how many lines in a file 
  nlines <- as.vector(countLines(x))
  
  # Read in .txt file 
  y <- read_table(x, skip=12, n_max=(nlines-20), col_types="cccccccccccc") %>% suppressWarnings() %>%
    # Fix column names
    rename(Ncells=Area, Area_km2=Area_1, pct_landsc=Area_2, pct_lands=as,
           CoreA=`%`, Edge=`of:`, EdgeA=CoreA., ShapeIdx=Edge, FractDim=`Edge:A`,
           drop1=Shape, drop2=Fract.) %>% 
    # Drop first row
    slice(-1) %>%
    # remove the extra columns
    select(-drop1, -drop2) %>%
    # convert to numeric
    mutate(across(c(Patch, Ncells, Area_km2, CoreA, Edge, EdgeA, ShapeIdx, FractDim), as.numeric)) %>%
    # Replace km2 calculation
    mutate(Area_km2=(Ncells*8100)/1000/1000)
  
  # Split up file name to add column for identifying information (species, bin, density)
  n <- str_split(x, pattern="[:punct:]")
  
  # Species column
  y$species <- n[[1]][3]
  # bin ID
  bintype <- str_detect(n[[1]][4], "\\d") # check if there is a number in string
  y$bin <- ifelse(bintype, as.character(parse_number(n[[1]][4])), n[[1]][4])
  # density
  y$density <- as.numeric(n[[1]][5])
  
  # Rearrange columns
  y <- y %>%
    select(species, density, bin, Patch:FractDim)
  
  return(y)
}

## Ramas carrying capacity (RAMAS Spatial Data)
ramas_K <- function(x) {
  
  # Use R.utils package to count how many lines in a file 
  nlines <- as.vector(countLines(x))
  
  # Read in .txt file 
  y <- read_table(x, skip=12, n_max=(nlines-15), col_types="ddddddddcdd") %>%
    suppressWarnings() %>%
    # Drop first row
    slice(-1) %>% 
    # Fix column names
    rename(TotalHS=Total, AvgHS=Aver., N0=Init., RelFec=Rel.,  
           RelSurv=vital, X=rates, Y=Mean) %>% 
    # remove the extra column
    select(-coord., -Rmax, -RelFec, -RelSurv) %>%
    # Fix X coordinates (drop comma)
    mutate(X=parse_number(X)) 
  
  # Split up file name to add column for identifying information (species, bin, density)
  n <- str_split(x, pattern="[:punct:]")
  
  # Species column
  y$species <- n[[1]][3]
  # bin ID
  # bintype <- str_detect(n[[1]][4], "\\d") & str_detect(n[[1]][4], "CI") # check if there is a number in string
  if (str_detect(n[[1]][4], "CI") & str_detect(n[[1]][4], "\\d")) {
    y$bin <- as.character(parse_number(n[[1]][4]))
    y$value <- str_sub(n[[1]][4], -3, -1)
  } else if (str_detect(n[[1]][4], "\\d") & !str_detect(n[[1]][4], "CI")) {
    # check if there is a digit without CI
    y$bin <- as.character(parse_number(n[[1]][4]))
    y$value <- "mean"
  } else if (!str_detect(n[[1]][4], "C") & !str_detect(n[[1]][4], "\\d")) {
    y$bin <- n[[1]][4]
    y$value <- "mean"
  } else {
    y$bin <- str_sub(n[[1]][4], 1, -4)
    y$value <- str_sub(n[[1]][4], -3, -1)
  }
  
  # y$bin <- ifelse(bintype, as.character(parse_number(n[[1]][4])), n[[1]][4])
  # density
  y$density <- as.numeric(n[[1]][5])
  
  # Rearrange columns
  y <- y %>%
    select(species, density, bin, value, Patch:Y)
  
  return(y)
}

# Get trajectory summary (MP simulation)
ramas_abun <- function(x) {
  
  # Read in .txt file 
  y <- read_table(x, skip=14, n_max=11, col_types= "ddddddd") %>% suppressWarnings()
  # Rename column headers
  names(y) <- c("time", "minimum", "lowerSD", "average", "upperSD", "maximum", "drop")
  y <- y %>% select(-drop)
  
  # Split up file name to add column for identifying information (species, bin, density)
  n <- str_split(x, pattern="[:punct:]")
  
  # Species column
  y$species <- n[[1]][3]
  # bin ID
  # bintype <- str_detect(n[[1]][4], "\\d") & str_detect(n[[1]][4], "CI") # check if there is a number in string
  if (str_detect(n[[1]][4], "CI") & str_detect(n[[1]][4], "\\d")) {
    y$bin <- as.character(parse_number(n[[1]][4]))
    y$value <- str_sub(n[[1]][4], -3, -1)
  } else if (str_detect(n[[1]][4], "\\d") & !str_detect(n[[1]][4], "CI")) {
    # check if there is a digit without CI
    y$bin <- as.character(parse_number(n[[1]][4]))
    y$value <- "mean"
  } else if (!str_detect(n[[1]][4], "C") & !str_detect(n[[1]][4], "\\d")) {
    y$bin <- n[[1]][4]
    y$value <- "mean"
  } else {
    y$bin <- str_sub(n[[1]][4], 1, -4)
    y$value <- str_sub(n[[1]][4], -3, -1)
  }
  # density
  y$density <- as.numeric(n[[1]][5])
  # density dependence?
  # y$densDep <- ifelse(length(n[[1]])==8, "scramble", "exp") 
  
  # Rearrange columns
  y <- y %>%
    select(species, density, bin, value, time:maximum)
  
  return(y)
}


### Calculating lambda from trajectory simulations
# Function for calculating geometric mean
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Function for calculating vector of lambda values
lambda_calc <- function(x) {
  lambda <- numeric()
  for(y in 1:length(x)){
    lambda[y] <- x[y+1]/x[y]
  }
  return(lambda)
}
