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


### Set up function to run population model with K #### 
run_deer_mod <- function(A_adj, theta, K, Sims=1000, years=50, ev_sd=0.02, harvest=TRUE, theta_mult=2, c_dd=0.2) {
  
  set.seed(1)
  
  ##Stochastistic model
  Year <- 1:years
  
  Kf <- K * 0.6
  
  # Empty array to hold pop count
  deer.array <- array(0,dim=c(12,length(Year),Sims))
  
  # Calculate stable stage distribution
  mnFs <- mean(c(A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6]))
  a2 <- A_adj
  a2[1,1:6] <- a2[1,1:6]*0.5*mnFs
  a2[7,1:6] <- a2[7,1:6]*0.5*mnFs
  
  # Calculate lambda from matrix with split fecundity (like it should be)
  lambda <- Re(eigen(a2)$values[1])
  
  # Set up initial popualtion size
  w <- Re(eigen(a2)$vectors[, 1])
  w <- w / sum(w) # normals to equal 1
  N0 <- w * 1 * K # e.g., start at 50% of K
  
  # Fill starting population
  deer.array[,1,] <- N0
  
  # Create matrix to hold ASR
  asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)
  
  # Set up arrays to save realized demographic rates
  realized_surv <- array(NA, dim=c(12, length(Year), Sims))
  realized_surv[1:12,1,] <- c(A_adj[2,1], A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6],
                              A_adj[8,7], A_adj[9,8], A_adj[10,9], A_adj[11,10], A_adj[12,11], A_adj[12,12])
  
  realized_fecundity <- array(NA, dim=c(4,  length(Year) - 1, Sims)) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,
  
  for (i in 1:Sims) {
    
    
    for (y in 2:length(Year)){
      
      A_s <- A_adj
      
      ## Stochasticity on Survival 
      # Female survival
      for (s in 2:5) {
        # Survive & go
        # Calculate beta distribution parameters
        p <- A_adj[s+1,s]
        
        # prevent extreme values
        p <- min(max(p, 1e-6), 1 - 1e-6)
        
        var <- ev_sd^2
        tmp <- (p * (1 - p) / var) - 1
        
        if (tmp <= 0) {
          # fallback: no stochasticity
          A_s[s+1,s] <- p
        } else {
          alpha <- p * tmp
          beta  <- (1 - p) * tmp
          A_s[s+1,s] <- rbeta(1, alpha, beta)
        }
        
      }
      # Survive & stay
      p <- A_adj[6,6]
      
      # prevent extreme values
      p <- min(max(p, 1e-6), 1 - 1e-6)
      
      var <- ev_sd^2
      tmp <- (p * (1 - p) / var) - 1
      
      if (tmp <= 0) {
        # fallback: no stochasticity
        A_s[6,6] <- p
      } else {
        alpha <- p * tmp
        beta  <- (1 - p) * tmp
        A_s[6,6] <- rbeta(1, alpha, beta)
      }
      
      # Male survival
      for (s in 8:11) {
        # Survive & go
        p <- A_adj[s+1,s]
        
        # prevent extreme values
        p <- min(max(p, 1e-6), 1 - 1e-6)
        
        var <- ev_sd^2
        tmp <- (p * (1 - p) / var) - 1
        
        if (tmp <= 0) {
          # fallback: no stochasticity
          A_s[s+1,s] <- p
        } else {
          alpha <- p * tmp
          beta  <- (1 - p) * tmp
          A_s[s+1,s] <- rbeta(1, alpha, beta)
        }
      }
      # Survive & stay
      p <- A_adj[12,12]
      
      # prevent extreme values
      p <- min(max(p, 1e-6), 1 - 1e-6)
      
      var <- ev_sd^2
      tmp <- (p * (1 - p) / var) - 1
      
      if (tmp <= 0) {
        # fallback: no stochasticity
        A_s[12,12] <- p
      } else {
        alpha <- p * tmp
        beta  <- (1 - p) * tmp
        A_s[12,12] <- rbeta(1, alpha, beta)
      }
      
      # save stochastic matrix
      A_dd <- A_s
      
      # Density dependent adjustment in fecundity
      # What was abundance at time step t-1
      Nf_t <- sum(deer.array[1:6,y-1,i], na.rm=TRUE)
      
      # Check adult age ratio
      female_sum <- sum(deer.array[2:6, y-1, i], na.rm=TRUE)
      male_sum   <- sum(deer.array[8:12, y-1, i], na.rm=TRUE)
      
      if (female_sum > 0) {
        asr <- male_sum / female_sum
      } else {
        asr <- 1
      }
      if (!is.finite(asr)) asr <- 1
      
      # Save sex ratio
      asr_mat[y-1,i] <- asr
      
      # Calculate % males & females
      s_m <- 0.5 + 0.1 * (1 - asr)
      s_m <- pmin(pmax(s_m, 0.4), 0.6)
      
      if (is.na(s_m)) s_m <- 0.5
      s_f <- 1 - s_m
      
      # Calcuate density factor
      density_factor <- 1 / (1 + c_dd * (Nf_t / Kf)^theta)
      density_factor <- pmax(density_factor, 0.2)
      if (!is.finite(density_factor)) density_factor <- 1
      
      # Adjust sex ratio at birth and reduce fecundity according to density
      A_dd[1,1:6] <- c(0, A_s[1,2], A_s[1,3], A_s[1,3], A_s[1,3], A_s[1,3]) * s_f * density_factor
      A_dd[7,1:6] <- c(0, A_s[1,2], A_s[1,3], A_s[1,3], A_s[1,3], A_s[1,3]) * s_m * density_factor
      
      # Stop loop if there are zeros
      if (any(is.na(A_dd))) {
        stop("NA detected in A_dd")
      }
      
      # Save fecundity
      realized_fecundity[1,y-1,i] <- A_dd[1,2] # Females per yearling female
      realized_fecundity[2,y-1,i] <- A_dd[1,3] # Females per adult female
      realized_fecundity[3,y-1,i] <- A_dd[7,2] # Males per yearling female 
      realized_fecundity[4,y-1,i] <- A_dd[7,3] # Males per adult female
      
      # Make sure we have a population to work with
      if (any(is.na(deer.array[, y-1, i]))) {
        stop(paste0("NA detected in deer.array at previous timestep!\n Sim: ",i,"\nYear: ",y ))
      }
      
      ## Harvest deer 
      # Calculate new pop size BEFORE harvest
      N_preharvest <- pmax(A_dd %*% deer.array[,y-1,i], 0) # Make sure to multiply matrix x vector (not vice versa)
      
      # Harvest
      if (harvest) {
        
        # Draw the number of deer to harvest each year
        h <- round(rnorm(1, mean=240000, sd=15000))
        
        doe_h <- h * 0.54
        buck_h <- h - doe_h
        
        doe_r <- c(0.0735, 0.189, 0.1895, 0.183, 0.183, 0.182) * doe_h
        
        # make sure we're not losing males by harvesting more than exist in that age class
        stage_weights <- c(0.04, 0.11, 0.12, 0.23, 0.25, 0.25)
        available <- N_preharvest[7:12]
        
        harvest_alloc <- stage_weights * available
        total_alloc <- sum(harvest_alloc)
        
        if (total_alloc > 0) {
          harvest_alloc <- harvest_alloc / total_alloc
        } else {
          harvest_alloc <- rep(0, length(harvest_alloc))
        }
        
        
        buck_r <- harvest_alloc * buck_h
        
        # Combine for total harvest
        r <- c(doe_r, buck_r)
        
        # make sure we're not taking more animals that exist in a certain age class
        r <- pmin(r, N_preharvest)
        
        # Apply harvest AFTER computing proportion
        N_next <- N_preharvest - r
        # Make sure population is not negative after harvest
        N_next <- pmax(N_next, 0)
        
      } else {
        N_next <- N_preharvest
      }
      
      # Save abundance
      deer.array[,y,i] <- pmax(N_next, 0)
      
      # ---- Realized survival ----
      
      # Female stages
      realized_surv[1, y, i] <- N_next[2] / deer.array[1, y-1, i] # Fawn → yearling (F)
      realized_surv[2, y, i] <- N_next[3] / deer.array[2, y-1, i]  # yearling → 2-year-old  
      realized_surv[3, y, i] <- N_next[4] / deer.array[3, y-1, i]  # 2-year-old → 3-year-old
      realized_surv[4, y, i] <- N_next[5] / deer.array[4, y-1, i]  # 3-year-old → 4-year-old
      
      # Stage 5: comes only from stage 4
      incoming_5 <- deer.array[4, y-1, i]
      realized_surv[5, y, i] <- ifelse(incoming_5 > 0,
                                       N_next[5] / incoming_5,
                                       NA)
      
      # Stage 6: from stage 5 and 6
      incoming_6 <- deer.array[5, y-1, i] + deer.array[6, y-1, i]
      realized_surv[6, y, i] <- ifelse(incoming_6 > 0,
                                       N_next[6] / incoming_6,
                                       NA)
      
      # Male stages
      realized_surv[7, y, i] <- N_next[8] / deer.array[7, y-1, i] # Fawn → yearling (M)
      realized_surv[8, y, i] <- N_next[9] / deer.array[8, y-1, i]   # 8 → 9
      realized_surv[9, y, i] <- N_next[10] / deer.array[9, y-1, i] # 9 → 10
      realized_surv[10, y, i] <- N_next[11] / deer.array[10, y-1, i] # 10 → 11
      
      # Stage 11
      incoming_11 <- deer.array[10, y-1, i]
      realized_surv[11, y, i] <- ifelse(incoming_11 > 1e-6,
                                        N_next[11] / incoming_11,
                                        NA)
      
      # Stage 12
      incoming_12 <- deer.array[11, y-1, i] + deer.array[12, y-1, i]
      realized_surv[12, y, i] <- ifelse(incoming_12 > 1e-6,
                                        N_next[12] / incoming_12,
                                        NA)
      
      # Check buck:doe ratio
      af <- sum(deer.array[2:6,y,i])
      am <- sum(deer.array[8:12,y,i])
      
      if (is.na(af) | sum(af) < 1) {
        warning("Female bottleneck: reproductive females lost!")
      }
      
    }
  }
  
  ## Summarize population trajectory
  N.median <- apply(apply(deer.array,c(2,3),sum),1,median)
  N.10pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.10)
  N.90pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.90)
  
  results <- data.frame(Year,N.median,N.10pct,N.90pct)
  
  ## Summarize recruits
  recruits <- deer.array[1,,] + deer.array[7,,]
  recruits <- rowMeans(recruits)
  
  ## Get realized survival rates
  # Calculate median survival per stage (over all years and sims)
  surv.50pct <- apply(realized_surv, 1, function(x) median(x, na.rm = TRUE))
  surv.10pct <- apply(realized_surv, 1, quantile, probs = 0.10, na.rm=TRUE)
  surv.90pct <- apply(realized_surv, 1, quantile, probs = 0.90, na.rm=TRUE)
  
  r.surv <- data.frame(sex=c(rep("Female",6), rep("Male",6)), 
                       stage=rep(c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"),2),
                       surv.50pct,surv.10pct,surv.90pct,
                       est=c(A_adj[2,1], A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6],
                             A_adj[8,7], A_adj[9,8], A_adj[10,9], A_adj[11,10], A_adj[12,11], A_adj[12,12]),
                       literature=c(A_base[2,1], A_base[3,2], A_base[4,3], A_base[5,4], A_base[6,5], A_base[6,6],
                                    A_base[8,7], A_base[9,8], A_base[10,9], A_base[11,10], A_base[12,11], A_base[12,12]))
  
  r.surv <- r.surv |>
    pivot_longer(cols=c("surv.50pct", "literature"), names_to="source", values_to="survival")
  
  # Set factor levels
  r.surv$stage <- factor(r.surv$stage, levels=c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"),
                         labels=c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"))
  
  r.surv$source <- factor(r.surv$source, levels=c("literature", "surv.50pct"), labels=c("Literature", "Implied"))
  
  ## Get realized fecundities
  # Calculate median survival per stage (over all years and sims)
  median_fec <- apply(realized_fecundity, 1, function(x) median(x, na.rm = TRUE))
  fec.10pct <- apply(realized_fecundity, 1, quantile, probs = 0.10, na.rm=TRUE)
  fec.90pct <- apply(realized_fecundity, 1, quantile, probs = 0.90, na.rm=TRUE)
  
  
  fec <- data.frame(stage=rep(c("Yearling", "Adult"), 2),
                    offspring_sex=rep(c("Female", "Male"), each=2), 
                    realized=median_fec,
                    f10pct=fec.10pct,
                    f90pct=fec.90pct,
                    literature=c(0.57, 0.66, 0.74, 0.85),
                    x=c(0.9, 1.6, 1.1, 1.8))
  fec <- fec |>
    pivot_longer(cols=c("realized", "literature"), names_to="source", values_to="fec")
  fec$source <- factor(fec$source, levels=c("literature", "realized"), labels=c("Literature", "Implied"))
  
  # Summarize sex ratio
  ASR.median <- apply(asr_mat, 1, median)
  ASR.20 <- apply(asr_mat, 1, quantile, probs = 0.20)
  ASR.80 <- apply(asr_mat, 1, quantile, probs = 0.80)
  
  ASR <- data.frame(Year, ASR.median, ASR.20, ASR.80)
  ASR <- ASR[-nrow(ASR),]
  
  # create list for storing results
  results_list <- list()
  results_list$results <- results
  results_list$recruits <- recruits
  results_list$adult_females <- rowMeans(colSums(deer.array[3:6, , ]))
  results_list$three_yr_males <- rowMeans(deer.array[10,,])
  results_list$four_yr_males <- rowMeans(deer.array[11,,])
  results_list$five_yr_males <- rowMeans(deer.array[12,,])
  results_list$last_stage <- sum(deer.array[12, , ] > 0)
  results_list$realized_fecundity <- realized_fecundity
  results_list$realized_surv <- realized_surv
  results_list$r.surv <- r.surv
  results_list$r.fec <- fec
  results_list$ASR <- ASR
  results_list$matrix_lambda <- lambda
  ## Get Lambda over simulation
  results_list$final_lambda <- gm_mean(lambda_calc(results$N.median[(length(Year)/2):length(Year)]))
  
  return(results_list)
  
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


# Calculate MSY and Objective function for optimizing ####
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

## Deer ####
scale_deer <- function(A_base, surv_scale_fawn, surv_scale_f, surv_scale_m, fec_scale) {
  
  # Set up deer matrix
  deer.matrix <- matrix(0,12,12)
  
  ## Survival
  # Females
  deer.matrix[2,1] <- A_base[2,1]*surv_scale_fawn
  deer.matrix[3,2] <- A_base[3,2]*surv_scale_f[1]
  deer.matrix[4,3] <- A_base[4,3]*surv_scale_f[2]
  deer.matrix[5,4] <- A_base[5,4]*surv_scale_f[3]
  deer.matrix[6,5] <- A_base[6,5]*surv_scale_f[4]
  deer.matrix[6,6] <- A_base[6,6]*surv_scale_f[5]
  # Males
  deer.matrix[8,7] <- A_base[8,7]*surv_scale_fawn
  deer.matrix[9,8] <- A_base[9,8]*surv_scale_m[1]
  deer.matrix[10,9] <- A_base[10,9]*surv_scale_m[2]
  deer.matrix[11,10] <- A_base[11,10]*surv_scale_m[3]
  deer.matrix[12,11] <- A_base[12,11]*surv_scale_m[4]
  deer.matrix[12,12] <- A_base[12,12]*surv_scale_m[5]
  
  ## Fecundity
  # Maximum fecundity
  R0a <- A_base[1,3]*fec_scale
  R0y <- A_base[1,2]*fec_scale
  
  # sex ratio at birth
  pct_m <- 0.565
  pct_f <- 0.435
  
  # Calculate for post-breeding census
  FecundityM <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*pct_m
  FecundityF <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*pct_f
  
  # Add to matrix
  deer.matrix[1,1:6] <- FecundityF
  deer.matrix[7,1:6] <- FecundityM
  
  return(deer.matrix)
}

deer_pop_proj <- function(A, N0, Kf, theta, Year, c_dd) {
  
  # Empty array to hold pop count
  deer.array <- matrix(0,nrow=12, ncol=length(Year))
  deer.array[,1] <- N0
  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + c_dd * (Nf_t / Kf)^theta)
    density_factor <- pmax(density_factor, 0.2)
    
    # save new matrix
    A_dd <- A
    
    # reduce fecundity
    R0y_dd <- A[1,2] * density_factor
    R0a_dd <- A[1,3] * density_factor
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1])/sum(deer.array[2:6,y-1], 1e-6)
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,1:6] <- c(0, R0y_dd*A[3,2], R0a_dd*A[4,3], R0a_dd*A[5,4], R0a_dd*A[6,5], R0a_dd*A[6,6])*s_m
    A_dd[7,1:6] <- c(0, R0y_dd*A[3,2], R0a_dd*A[4,3], R0a_dd*A[5,4], R0a_dd*A[6,5], R0a_dd*A[6,6])*s_f
    
    # Calculate new pop size
    deer.array[,y] <- A_dd %*% deer.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
  }
  
  # Add up stage-specific abundance
  N.median <- apply(deer.array,2,sum)
  
  # Put into data frame
  results <- data.frame(Year,N.median)
  
  # Calculate MSY
  msy <- get_msy(results)
  
  # return(msy)
  # Return population trajectory
  return(N.median)
}

## Define objective function
objective_fn_deer <- function(params) {
  
  surv_scale_fawn <- params[1]
  surv_scale_f <- params[2:6]
  surv_scale_m <- params[7:11]
  fec_scale <- params[12]
  
  # Scale demographic rates
  A_scaled <- scale_deer(A_base, surv_scale_fawn, surv_scale_f, surv_scale_m, fec_scale)
  # cap survival (so no chance of >1)
  A_scaled[2:6, 1:6] <- pmin(A_scaled[2:6, 1:6], 0.999)
  A_scaled[8:12, 7:12] <- pmin(A_scaled[8:12, 7:12], 0.999)
  
  # Calculate stable stage distribution
  w <- Re(eigen(A_scaled)$vectors[, 1])
  w <- w / sum(w) # normalize to equal 1
  
  # Set up initial popualtion size
  N0 <- w * 0.01 * K  # e.g., start at 50% of K (or in this case 1% of K)
  
  # Calibrate a for density dependence
  theta <- params[13]
  
  error_theta <- (theta - 1)^2
  
  # --- Population trajectory ---
  sim_N <- deer_pop_proj(A_scaled, N0, Kf, theta, Year, c_dd=c_dd) 
  
  # use last half of time series (avoid transient dynamics)
  tail_N <- tail(sim_N, length(sim_N) / 2)
  
  # simple smooth target: want stable population near K (or at least, the female portion of K)
  error_pop <- mean(((tail_N - Kf) / Kf)^2) + var(tail_N) / Kf^2
  
  # --- MSY --- 
  sim_msy <- get_msy(data.frame(Year, N.median=sim_N))
  
  # Guard against invalid MSY
  if (!is.finite(sim_msy) || sim_msy <= 0) {
    return(1e6)
  }
  
  error_msy <- ((sim_msy - observed_harvest) / observed_harvest)^2
  
  # --- Regularize deviation from baseline (1 = literature values)
  error_reg <- sum((params - 1)^2)
  
  # --- Penalty for slow pop growth rate ----
  lambda <- Re(eigen(A_scaled)$values[1])
  
  error_lambda <- 0
  
  if (lambda < 1.4) {
    error_lambda <- (1.4 - lambda)^2
  } else if (lambda > 1.77) {
    error_lambda <- (lambda - 1.77)^2
  }
  
  # --- Combine ---
  return(error_pop + error_msy + 0.05 * error_reg + 0.1 * error_theta + 50 * error_lambda) 
  
}

## Wild pigs ####
# Scale matrices
scale_pigs <- function(A_base,surv_scale_p,surv_scale_f,surv_scale_m,fec_scale) {
  
  pig.matrix <- matrix(0,6,6)
  
  ## Survival
  # Females
  pig.matrix[2,1] <- A_base[2,1]*surv_scale_p
  pig.matrix[3,2] <- A_base[3,2]*surv_scale_f[1]
  pig.matrix[4,3] <- A_base[4,3]*surv_scale_f[2]
  # Males
  pig.matrix[5,4] <- A_base[5,4]*surv_scale_p
  pig.matrix[6,5] <- A_base[6,5]*surv_scale_m[1]
  pig.matrix[6,6] <- A_base[6,6]*surv_scale_m[2]
  
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
  
  surv_scale_p <- params[1]
  surv_scale_f <- params[2:3]
  surv_scale_m <- params[4:5]
  fec_scale <- params[6]
  
  # Scale demographic rates
  A_scaled <- scale_pigs(A_base,surv_scale_p,surv_scale_f,surv_scale_m,fec_scale)
  
  # Calculate stable stage distribution
  w <- Re(eigen(A_scaled)$vectors[, 1])
  w <- w / sum(w) # normalize to equal 1
  
  # Set up initial population size
  N0 <- w * 0.01 * K 
  
  lambda <- Re(eigen(A_scaled)$values[1])
  
  # Calibrate a for density dependence
  theta <- estimate_theta_opt(N0, lambda, K)
  
  # Run model and get MSY
  sim_msy <- pig_pop_proj(A_scaled, N0, Kf, theta, step) 
  
  return((sim_msy - observed_harvest)^2) # Squared error
}

