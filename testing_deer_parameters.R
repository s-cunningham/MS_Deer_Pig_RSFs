




# Set up grid
# Example ranges (adjust as needed)
female_surv_seq <- seq(0.80, 0.95, by = 0.01)
male_surv_seq   <- seq(0.40, 0.90, by = 0.02)
fecundity_seq   <- seq(1, 1.5, by = 0.1)
harvest_seq     <- seq(100000, 300000, by = 50000)  # total individuals/year

param_grid <- expand.grid(
  female_surv = female_surv_seq,
  male_surv   = male_surv_seq,
  fecundity   = fecundity_seq,
  harvest     = harvest_seq
)


## Set up function to run model


run_pop_proj <- function(female_surv, male_surv, fec_mult, harvest, K) {
  
  # ---------------------
  # PARAMETERS
  # ---------------------
  
  # Time
  Tmax <- 100
  
  # Number of stages
  n_stages <- 6
  
  # Fecundity (base)
  R0a <- 1.8
  R0y <- 1.4
  
  # Stage-specific fecundity
  fec_f <- c(0, R0y, R0a, R0a, R0a, R0a)*f_mult
  
  # Sex ratio at birth: 1.3 males : 1 female
  sr <- 1.3
  pf <- 1 / (1 + sr)
  pm <- sr / (1 + sr)
  

  # ---------------------
  # INITIALIZATION
  # ---------------------
  N_f <- matrix(0, nrow = Tmax + 1, ncol = n_stages)
  N_m <- matrix(0, nrow = Tmax + 1, ncol = n_stages)
  
  # Initial population (example)
  N_f[1, ] <- c(25000, 14000, 10500, 9000, 8000, 7000)
  N_m[1, ] <- c(25000, 14000, 10000, 8000, 7000, 6000)
  
  K_f <- K * 0.6
  
  # ---------------------
  # SIMULATION
  # ---------------------
  for (t in 1:Tmax) {
    
    # Total females at time t (for density dependence)
    N_total_f <- sum(N_f[t, ])
    
    # Density-adjusted fecundity
    DD_scalar <- 1 / (1 + (N_total_f / K_f)^theta)
    births_per_female <- fec_f * DD_scalar
    total_births <- sum(N_f[t, ] * births_per_female)
    
    # Distribute newborns by sex
    N_f[t + 1, 1] <- total_births * pf
    N_m[t + 1, 1] <- total_births * pm
    
    # -------- HARVEST --------
    # Example: proportional harvest by stage (could be custom)
    prop_f <- N_f[t, ] / sum(N_f[t, ])
    prop_m <- N_m[t, ] / sum(N_m[t, ])
    
    harvest_f <- H_f_total * prop_f
    harvest_m <- H_m_total * prop_m
    # -------- SURVIVAL & ADVANCEMENT --------
    
    # Females
    for (s in 2:(n_stages - 1)) {
      N_f[t + 1, s] <- pmax((N_f[t, s - 1] - harvest_f[s - 1]) * S_f[s - 1],0)
    }
    # Final female stage accumulates
    N_f[t + 1, 6] <- pmax((N_f[t, 5] - harvest_f[5]) * S_f[5] +  (N_f[t, 6] - harvest_f[6]) * S_f[6],0)
    
    # Males
    for (s in 2:(n_stages - 1)) {
      N_m[t + 1, s] <- pmax((N_m[t, s - 1] - harvest_m[s - 1]) * S_m[s - 1],0)
    }
    N_m[t + 1, 6] <- pmax((N_m[t, 5] - harvest_m[5]) * S_m[5] + (N_m[t, 6] - harvest_m[6]) * S_m[6], 0)
  } 
}


