
library(tidyverse)
library(R.utils)

### Functions for 07_read_RAMAS_output.R ####
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
  # suitability level
  y$level <- n[[1]][4]
  
  # mean, LCI, UCI
  y$value <- ifelse(length(n[[1]])==7, "mean", n[[1]][6])  

  # density
  y$density <- as.numeric(n[[1]][5])
  
  # Rearrange columns
  y <- y %>%
    select(species, density, level, value, Patch:FractDim)
  
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
  y$suitability <- n[[1]][4]
  # bin ID
  # bintype <- str_detect(n[[1]][4], "\\d") & str_detect(n[[1]][4], "CI") # check if there is a number in string
  if (str_detect(n[[1]][6], "CI")) {
    y$value <- n[[1]][6]
  }  else {
    y$value <- "mean"
  }
  
  # y$bin <- ifelse(bintype, as.character(parse_number(n[[1]][4])), n[[1]][4])
  # density
  y$density <- as.numeric(n[[1]][5])
  
  # Rearrange columns
  y <- y %>%
    select(species, density, suitability, value, Patch:Y)
  
  return(y)
}


### Calculating lambda from trajectory simulations ####
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

### Functions for validation of carrying capacity/MSY exercise ####
simulate_pop <- function(a, A, N0, K_obs, years = 100) {
  N <- matrix(NA, nrow = nrow(A), ncol = years)
  N[, 1] <- N0
  
  for (t in 1:(years - 1)) {
    N_total <- sum(N[, t])
    dd_factor <- 1 / (1 + (N[t - 1] / K)^theta)
    
    A_dd <- A
    A_dd[1, ] <- A[1, ] * dd_factor  # apply DD to fecundity row
    N[, t + 1] <- A_dd %*% N[, t]
  }
  
  final_N <- sum(N[, years])
  return(final_N)
}

## theta-logistic function
# Example theta-logistic projection function
theta_logistic_proj <- function(N0, lambda, K, theta, t_max = 100) {
  N <- numeric(t_max)
  N[1] <- sum(N0)
  for (t in 2:t_max) {
    density_factor <- 1 / (1 + (N[t - 1] / K)^theta)
    N[t] <- N[t - 1] * lambda * density_factor
  }
  return(N)
}

# Estimate with optimize
estimate_theta_opt <- function(N0, lambda, K, t_max = 100) {
  objective_fn <- function(theta) {
    N <- theta_logistic_proj(N0, lambda, K, theta, t_max)
    (tail(N, 1) - K)^2
  }
  result <- optimize(objective_fn, lower = 0.1, upper = 10)
  return(result$minimum)
}


#### Calculate MSY and Objective function for optimizing ####
## Calculate MSY
get_msy <- function(results) {
  
  # Calculate net difference in median Nt - Nt-1
  net <- results$N.median[2:nrow(results)] - results$N.median[1:(nrow(results)-1)]
  
  # Find maximum difference
  peak_idx <- which.max(net)
  
  # Save MSY value
  msy <- net[peak_idx]
  
  return(msy)
}


scale_deer <- function(A_base, surv_scale_fawn, surv_scale_f, surv_scale_m, fec_scale_y, fec_scale_a) {
  
  # Set up deer matrix
  deer.matrix <- matrix(0,12,12)
  
  ## Survival
  # Females
  deer.matrix[2,1] <- A_base[2,1]*surv_scale_fawn
  deer.matrix[3,2] <- A_base[3,2]*surv_scale_f[1]
  deer.matrix[4,3] <- A_base[4,3]*surv_scale_f[2]
  deer.matrix[5,4] <- A_base[5,4]*surv_scale_f[2]
  deer.matrix[6,5] <- A_base[6,5]*surv_scale_f[2]
  deer.matrix[6,6] <- A_base[6,6]*surv_scale_f[2]
  # Males
  deer.matrix[8,7] <- A_base[8,7]*surv_scale_fawn
  deer.matrix[9,8] <- A_base[9,8]*surv_scale_m[1]
  deer.matrix[10,9] <- A_base[10,9]*surv_scale_m[2]
  deer.matrix[11,10] <- A_base[11,10]*surv_scale_m[3]
  deer.matrix[12,11] <- A_base[12,11]*surv_scale_m[4]
  deer.matrix[12,12] <- A_base[12,12]*surv_scale_m[5]
  
  ## Fecundity
  # Maximum fecundity
  R0a <- A_base[1,3]*fec_scale_a
  R0y <- A_base[1,2]*fec_scale_y
  
  # pct_m <- 0.565
  # pct_f <- 0.435
  
  # Calculate for post-breeding census
  # FecundityM <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*pct_m 
  # FecundityF <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*pct_f 
  
  FecundityF <- c(0, R0y, R0a, R0a, R0a, R0a)*0.5
  FecundityM <- c(0, R0y, R0a, R0a, R0a, R0a)*0.5
  
  # Add to matrix
  deer.matrix[1,1:6] <- FecundityF
  deer.matrix[7,1:6] <- FecundityM
  
  return(deer.matrix)
}


deer_pop_proj <- function(A, N0, Kf, theta, Year) {
  
  # Empty array to hold pop count
  deer.array <- matrix(0,nrow=12, ncol=length(Year))
  deer.array[,1] <- N0
  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # save new matrix
    A_dd <- A
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # reduce fecundity
    R0y_dd <- A[1,2] * density_factor 
    R0a_dd <- A[1,3] * density_factor 
    
    # Adjust sex ratio at birth
    A_dd[1,1:6] <- c(0, R0y_dd*A[3,2], R0a_dd*A[4,3], R0a_dd*A[5,4], R0a_dd*A[6,5], R0a_dd*A[6,6]) 
    A_dd[7,1:6] <- c(0, R0y_dd*A[3,2], R0a_dd*A[4,3], R0a_dd*A[5,4], R0a_dd*A[6,5], R0a_dd*A[6,6])
    
    # Calculate new pop size
    deer.array[,y] <- A_dd %*% deer.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
  }
  
  # Add up stage-specific abundance
  N.median <- apply(deer.array,2,sum)
  
  # Put into data frame
  results <- data.frame(Year,N.median)
  
  # Calculate MSY
  msy <- get_msy(results)
  
  return(msy)
}


## Define objective function
objective_fn_deer <- function(params) {
  
  surv_scale_fawn <- params[1]
  surv_scale_f <- params[2:3]
  surv_scale_m <- params[4:8]
  fec_scale_y <- params[9]
  fec_scale_a <- params[10]
  
  # Scale demographic rates
  A_scaled <- scale_deer(A_base, surv_scale_fawn, surv_scale_f, surv_scale_m, fec_scale_y, fec_scale_a)
  
  # Calculate stable stage distribution
  w <- Re(eigen(A_scaled)$vectors[, 1])
  w <- w / sum(w) # normalize to equal 1
  
  # Set up initial popualtion size
  N0 <- w * 0.1 * K  # e.g., start at 50% of K
  
  lambda <- Re(eigen(A_scaled)$values[1])
  
  # Calibrate a for density dependence
  theta <- 3
  
  # Run model and get MSY
  sim_msy <- deer_pop_proj(A_scaled, N0, Kf, theta, Year) 
  
  return((sim_msy - observed_harvest)^2) # Squared error
}




# Scale matrices
scale_pigs <- function(A_base,surv_scale_f,surv_scale_m,fec_scale) {
  
  pig.matrix <- matrix(0,6,6)
  
  ## Survival
  # Females
  pig.matrix[2,1] <- A_base[2,1]*surv_scale_f
  pig.matrix[3,2] <- A_base[3,2]*surv_scale_f
  pig.matrix[4,3] <- A_base[4,3]*surv_scale_f
  # Males
  pig.matrix[5,4] <- A_base[5,4]*surv_scale_m
  pig.matrix[6,5] <- A_base[6,5]*surv_scale_m
  pig.matrix[6,6] <- A_base[6,6]*surv_scale_m
  
  ## Fecundity
  # Maximum fecundity
  R0y <- A_base[1,2]*fec_scale
  R0a <- A_base[1,3]*fec_scale
  
  FecundityF <- c(R0y, R0a)*0.5
  FecundityM <- c(R0y, R0a)*0.5
  
  # Add to matrix
  pig.matrix[1,2:3] <- FecundityF
  pig.matrix[4,2:3] <- FecundityM
  
  return(pig.matrix)
}


## project population
pig_pop_proj <- function(A, N0, Kf, theta, step) {
  
  # Empty array to hold pop count
  pig.array <- matrix(0,nrow=6, ncol=length(step))
  pig.array[,1] <- N0
  
  for (y in 2:length(step)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(pig.array[1:3,y-1])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # save new matrix
    A_dd <- A
    
    # reduce fecundity
    R0j_dd <- A[1,2] * density_factor 
    R0a_dd <- A[1,3] * density_factor 
    
    # Adjust sex ratio at birth
    A_dd[1,2:3] <- c(R0j_dd*A[3,2], R0a_dd*A[3,3]) 
    A_dd[4,2:3] <- c(R0j_dd*A[3,2], R0a_dd*A[3,3])
    
    # Calculate new pop size
    pig.array[,y] <- A_dd %*% pig.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
  }
  
  # Add up stage-specific abundance
  N.median <- apply(pig.array,2,sum)
  
  # Put into data frame
  results <- data.frame(step,N.median)
  
  # Calculate MSY
  msy <- get_msy(results)
  
  return(msy)
}


objective_fn_pigs <- function(params) {
  
  surv_scale_f <- params[1]
  surv_scale_m <- params[2]
  fec_scale <- params[3]
  
  # Scale demographic rates
  A_scaled <- scale_pigs(A_base,surv_scale_f,surv_scale_m,fec_scale)
  
  # Calculate stable stage distribution
  w <- Re(eigen(A_scaled)$vectors[, 1])
  w <- w / sum(w) # normalize to equal 1
  
  # Set up initial popualtion size
  N0 <- w * 0.1 * K  # e.g., start at 50% of K
  
  lambda <- Re(eigen(A_scaled)$values[1])
  
  # Calibrate a for density dependence
  # a <- calibrate_a_opt(A_scaled, N0, K)
  theta <- estimate_theta_opt(N0, lambda, K)
  
  # Run model and get MSY
  sim_msy <- pig_pop_proj(A_scaled, N0, Kf, theta, step) 
  
  return((sim_msy - observed_harvest)^2) # Squared error
}




# Helper function to avoid division by zero
safe_divide <- function(numerator, denominator) {
  if (is.na(denominator) || is.nan(denominator) || denominator <= 0) {
    return(NA)
  } else {
    return(numerator / denominator)
  }
}
