library(ggplot2)
library(ggtext)
library(patchwork)
library(truncnorm)
library(patchwork)

source("00_functions.R")

theme_set(theme_bw())

set.seed(1)

# Constants
K <- 1347232
Kf <- K * 0.60
observed_harvest <- 280000
pct_m <- 0.565
pct_f <- 0.435

# Base matrix setup
A_base <- matrix(0, 12, 12)

# Female survival
A_base[2,1] <- 0.70
A_base[3,2] <- 0.952
A_base[4,3] <- 0.914+0.04
A_base[5,4] <- 0.914+0.08
A_base[6,5] <- 0.914+0.04
A_base[6,6] <- 0.914+0.04
fsurv <-c(0.952, 0.914, 0.914, 0.914, 0.914)+0.04

# Male survival
A_base[8,7] <- 0.70
A_base[9,8] <- 0.99
A_base[10,9] <- 0.93
A_base[11,10] <- 0.94
A_base[12,11] <- 0.89
A_base[12,12] <- 0.82

fecundities <- c(1.33+(1-0.952+1), 1.72+(1-0.954+1), 1.83+(1-0.954+1), 1.83+(1-0.954+1), 1.83+(1-0.954+1))
# fecundities <- c(1.5, 1.9, 2, 2, 2)

# Build matrix with Beverton-Holt density dependence
build_matrix <- function(fecundities, A_base) {
  surv <- c(A_base[3,2], A_base[4,3], A_base[5,4], A_base[6,5], A_base[6,6])
  
  Ff <- c(0, fecundities * surv * pct_f)
  Fm <- c(0, fecundities * surv * pct_m)
  
  A <- A_base
  A[1,1:6] <- Ff
  A[7,1:6] <- Fm
  return(A)
}

A <- build_matrix(fecundities, A_base)


# Stable stage distribution
stable_stage <- function(A) {
  w <- Re(eigen(A)$vectors[,1])
  w / sum(w)
}

# Initial abundance by stage
ss <- stable_stage(A)
N0 <- K * 0.01
N0 <- N0*ss

nyears <- 50

# Fill starting population
deer.array <- matrix(0, nrow=12, ncol=nyears)
deer.array[,1] <- N0

theta <- 3

A_dd <- A

# Set up arrays to save realized demographic rates
realized_survival <- matrix(NA, nrow=12, ncol=nyears)
realized_survival[1:12,1] <- c(A[2,1], A[3,2], A[4,3], A[5,4], A[6,5], A[6,6],
                                A[8,7], A[9,8], A[10,9], A[11,10], A[12,11], A[12,12])

realized_fecundity <- matrix(NA, nrow=4, ncol=nyears - 1) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,

n_harvest <- matrix(0, nrow=12, ncol=nyears)

set.seed(1)
for (y in 2:nyears) {
  
  # Female survival
  A_dd[2,1] <- realized_survival[1,y-1]
  A_dd[3,2] <- realized_survival[2,y-1]
  A_dd[4,3] <- realized_survival[3,y-1]
  A_dd[5,4] <- realized_survival[4,y-1]
  A_dd[6,5] <- realized_survival[5,y-1]
  A_dd[6,6] <- realized_survival[6,y-1]
  
  # Male survival
  A_dd[8,7] <- realized_survival[7,y-1]
  A_dd[9,8] <- realized_survival[8,y-1]
  A_dd[10,9] <- realized_survival[9,y-1]
  A_dd[11,10] <- realized_survival[10,y-1]
  A_dd[12,11] <- realized_survival[11,y-1]
  A_dd[12,12] <- realized_survival[12,y-1]
  
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  Nf_t <- sum(deer.array[1:6,y-1], na.rm=TRUE)
  
  # Calcuate density factor
  density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
  
  # Check adult age ratio
  asr <- sum(deer.array[8:12,y-1])/sum(deer.array[2:6,y-1], 1e-6)
  
  # Calculate % males & females
  s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
  s_f <- 1 - s_m
  
  # Adjust sex ratio at birth
  A_dd[1,1:6] <- c(0, fecundities * realized_survival[2:6,y-1] * s_f * density_factor)
  A_dd[7,1:6] <- c(0, fecundities * realized_survival[2:6,y-1] * s_m * density_factor)
  
  # Calculate new pop size
  N_next <- pmax(A_dd %*% deer.array[,y-1],0) # Make sure to multiply matrix x vector (not vice versa)
  
  # Harvest, but only if population is > 1/2K
  # if (sum(N_next) > 0.5) {
  #   # Randomly select a random number to harvest
  #   # h <- round(runif(1, 190000, 260000), digits=0)
  #   # h <- round(rnorm(1, 250000, 15000), digits=0)
  #   h <- 200000
  #   
  #   # Split bucks and does
  #   doe_h <- h * 0.54
  #   buck_h <- h - doe_h
  #   
  #   # Proportional does
  #   doe_r <- c(0.0735, 0.189, 0.1895, 0.183, 0.183, 0.182)*doe_h
  #   
  #   # Proportional bucks
  #   buck_r <- c(0.04, 0.11, 0.12, 0.23, 0.25, 0.25)*buck_h
  #   
  #   # Create vector of individuals to harvest
  #   r <- c(doe_r, buck_r)
  #   
  #   # Remove from population 
  #   N_next <- N_next - r
  #   
  #   # Save number harvested
  #   n_harvest[,y] <- r
  #   
  # }
  
  # Save harvest-adjusted population size
  deer.array[,y] <-  pmax(N_next, 0)
  
  # ---- Realized survival ----
  
  # Female stages
  realized_survival[1, y] <- deer.array[2, y] / deer.array[1, y-1] # Fawn → yearling (F)
  realized_survival[2, y] <- deer.array[3, y] / deer.array[2, y-1]  # yearling → 2-year-old  
  realized_survival[3, y] <- deer.array[4, y] / deer.array[3, y-1]  # 2-year-old → 3-year-old
  realized_survival[4, y] <- deer.array[5, y] / deer.array[4, y-1]  # 3-year-old → 4-year-old
  
  # female ssd
  w_f <- Re(eigen(A_dd)$vectors[, 1])[1:6]
  w_f <- w_f / sum(w_f)
  
  # male ssd
  w_m <- Re(eigen(A_dd)$vectors[, 1])[7:12]
  w_m <- w_m / sum(w_m)
  
  c5_f <- (A_dd[6,5] * w_f[5]) / ((A_dd[6,5] * w_f[5]) + (A_dd[6,6] * w_f[6]))
  c5_m <- (A_dd[12,11] * w_m[5]) / ((A_dd[12,11] * w_m[5]) + (A_dd[12,12] * w_m[6]))
  
  # Stages 5 and 6 are tricky...
  est_from5 <- deer.array[6, y] * c5_f
  est_from6 <- deer.array[6, y] - est_from5
  
  realized_survival[5, y] <- est_from5 / deer.array[5, y-1]
  
  # Stage 6 is tricky: receives from 5 and from itself
  incoming_females <- deer.array[6, y-1] + deer.array[5, y-1]  # crude approx
  if (incoming_females > 0) {
    realized_survival[6, y] <- deer.array[6, y] / incoming_females
  }
  
  # Male stages
  realized_survival[7, y] <- deer.array[8, y] / deer.array[7, y-1] # Fawn → yearling (M)
  realized_survival[8, y] <- deer.array[9, y] / deer.array[8, y-1]   # 8 → 9 (yearling - 2yers)
  realized_survival[9, y] <- deer.array[10, y] / deer.array[9, y-1] # 9 → 10 (2yrs - 3yrs)
  realized_survival[10, y] <- deer.array[11, y] / deer.array[10, y-1] # 10 → 11 (3 yrs - 4 yrs)
  
  # Stages 5 and 6 are tricky...
  est_from11 <- deer.array[12, y] * c5_m
  est_from12 <- deer.array[12, y] - est_from11
  
  realized_survival[11, y] <- est_from11 / deer.array[11, y-1]
  
  # Stage 12 (adult males): 11 → 12 and 12 → 12
  incoming_males <- deer.array[12, y-1] + deer.array[11, y-1]
  if (incoming_males > 0) {
    realized_survival[12, y] <- deer.array[12, y] / incoming_males
  }
  
}

lambda <- Re(eigen(A_dd)$values[1])

## Summarize deer.array so that there is 1 pop size per year
N <- colSums(deer.array)

plot(N, type="l", ylim=c(0, 3000000))
abline(h=K, col="red", lty=2)

# Calculate net difference in median Nt - Nt-1
net <- N[2:length(N)] - N[1:(length(N)-1)]

# Find maximum difference
peak_idx <- which.max(net)

# Save MSY value
msy <- net[peak_idx]
msy

plot(y=net, x=N[2:50], type="l")

