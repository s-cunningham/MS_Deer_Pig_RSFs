

library(truncnorm)
library(ggplot2)
library(ggtext)
library(patchwork)

source("00_functions.R")

theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:100

# Carrying capacity
# K <- 2347197
K <- 1512638

# Set up deer matrix
deer.matrix <- matrix(0,12,12)

## Survival
# Females
deer.matrix[2,1] <- 0.60
deer.matrix[3,2] <- 0.93
deer.matrix[4,3] <- 0.95
deer.matrix[5,4] <- 0.95
deer.matrix[6,5] <- 0.95
deer.matrix[6,6] <- 0.95
# Males
deer.matrix[8,7] <- 0.60
deer.matrix[9,8] <- 0.90
deer.matrix[10,9] <- 0.92
deer.matrix[11,10] <- 0.92
deer.matrix[12,11] <- 0.92
deer.matrix[12,12] <- 0.83

## Fecundity
# Maximum fecundity
R0a <- 2.2#1.8
R0y <- 2#1.4

pct_m <- 0.565
pct_f <- 0.435

# Calculate for post-breeding census
FecundityM <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*pct_m 
FecundityF <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*pct_f 
# Fecundity <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*0.5

# Add to matrix
deer.matrix[1,1:6] <- FecundityF
deer.matrix[7,1:6] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(deer.matrix)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(deer.matrix)$values[1])

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, K)
cat("Estimated theta:", theta, "\n")

Kf <- K * 0.60

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
  A_dd <- deer.matrix
  A_dd[1,2:6] <- A_dd[1,2:6] * density_factor # reduce fecundity
  A_dd[7,2:6] <- A_dd[7,2:6] * density_factor # reduce fecundity
  
  # Calculate new pop size
  deer.array[,y] <- A_dd %*% deer.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
}

N.median <- apply(deer.array,2,sum)

results <- data.frame(Year,N.median)

#Plot population with
ggplot(results) +
  geom_hline(yintercept=K) +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (males + females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


#### Maximum sustained yield ####
net <- results$N.median[2:nrow(results)] - results$N.median[1:(nrow(results)-1)]

peak_idx <- which.max(net)
msy <- net[peak_idx]
pop_at_msy <- results$N.median[peak_idx+1]

df <- data.frame(Nt=results$N.median[2:nrow(results)], net=net)

ggplot(df) +
  coord_cartesian(ylim=c(0,330000)) +
  geom_hline(yintercept=280000, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_vline(xintercept=K) +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))


##### Set up deer matrix ####
orig.deer <- matrix(0,12,12)

## Survival
# Females
orig.deer[2,1] <- 0.52
orig.deer[3,2] <- 0.93
orig.deer[4,3] <- 0.84
orig.deer[5,4] <- 0.84
orig.deer[6,5] <- 0.84
orig.deer[6,6] <- 0.84
# Males
orig.deer[8,7] <- 0.52
orig.deer[9,8] <- 0.82
orig.deer[10,9] <- 0.63
orig.deer[11,10] <- 0.53
orig.deer[12,11] <- 0.44
orig.deer[12,12] <- 0.49

## Fecundity
# Maximum fecundity
R0a <- 1.8
R0y <- 1.4

pct_m <- 0.565
pct_f <- 0.435

# Calculate for post-breeding census
FecundityM <- c(0, R0y*orig.deer[3,2], R0a*orig.deer[4,3], R0a*orig.deer[5,4], R0a*orig.deer[6,5], R0a*orig.deer[6,6])*pct_m 
FecundityF <- c(0, R0y*orig.deer[3,2], R0a*orig.deer[4,3], R0a*orig.deer[5,4], R0a*orig.deer[6,5], R0a*orig.deer[6,6])*pct_f 
# Fecundity <- c(0, R0y*deer.matrix[3,2], R0a*deer.matrix[4,3], R0a*deer.matrix[5,4], R0a*deer.matrix[6,5], R0a*deer.matrix[6,6])*0.5

# Add to matrix
orig.deer[1,1:6] <- FecundityF
orig.deer[7,1:6] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(deer.matrix)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(orig.deer)$values[1])

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, K)
cat("Estimated theta:", theta, "\n")

Kf <- K * 0.60

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
  A_dd <- orig.deer
  A_dd[1,2:6] <- A_dd[1,2:6] * density_factor # reduce fecundity
  A_dd[7,2:6] <- A_dd[7,2:6] * density_factor # reduce fecundity
  
  # Calculate new pop size
  deer.array[,y] <- A_dd %*% deer.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
}

N.median <- apply(deer.array,2,sum)

results <- data.frame(Year,N.median)

#Plot population with
ggplot(results) +
  geom_hline(yintercept=K) +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (males + females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


#### Maximum sustained yield ####
net <- results$N.median[2:nrow(results)] - results$N.median[1:(nrow(results)-1)]

peak_idx <- which.max(net)
msy <- net[peak_idx]
pop_at_msy <- results$N.median[peak_idx+1]

df <- data.frame(Nt=results$N.median[2:nrow(results)], net=net)

ggplot(df) +
  coord_cartesian(ylim=c(0,330000)) +
  geom_hline(yintercept=280000, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_vline(xintercept=K) +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))




















#### Base model (deterministic) ####

A_adj <- deer.matrix

# Empty array to hold pop count
Year <- 1:100
deer.array <- matrix(0,nrow=12, ncol=length(Year))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")

deer.array[,1] <- N0

for (y in 2:length(Year)){
  
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  Nf_t <- sum(deer.array[1:6,y-1])
  
  # Calcuate density factor
  density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
  
  # save new matrix
  A_dd <- A_adj
  A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
  A_dd[7, ] <- A_dd[7, ] * density_factor # reduce fecundity
  
  # Calculate new pop size
  deer.array[,y] <- A_dd %*% deer.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
}


N.median <- apply(deer.array,2,sum)

results <- data.frame(Year,N.median)

#Plot population with
ggplot(results) +
  # coord_cartesian(ylim=c(0, 3200000)) +
  geom_hline(yintercept=K) +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (males + females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

#### Base Model - stochasticity in vital rates ####
set.seed(1)
##Stochastistic model
Year <- 1:20
Sims <- 1000

# Stochasticity on fecundity
f_pct <- 0.01

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Carrying capacity
Ka <- 1345648
Kf <- Ka * 0.6

# Set up initial popualtion size
N0 <- w * 1 * K # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")
# theta <- 5

s_decr <- 0.02

for (i in 1:Sims) {
  
  ## Stochasticity on Survival 
  A_s <- A_adj
  # Female survival
  for (s in 2:5) {
    # Survive & go
    A_s[s+1,s] <- rtruncnorm(1, b=1, mean=A_adj[s+1,s]-s_decr, sd=0.05)
  }
  # Survive & stay
  A_s[6,6] <- rtruncnorm(1, a=0, b=1, mean=A_adj[6,6]-s_decr, sd=0.05)
  # Male survival
  for (s in 8:11) {
    # Survive & go
    A_s[s+1,s] <- rtruncnorm(1, b=1, mean=A_adj[s+1,s]-s_decr, sd=0.05)
  }
  # Survive & stay
  A_s[12,12] <- rtruncnorm(1, b=1, mean=A_adj[12,12]-s_decr, sd=0.05)
  
  ## Stochasticity on fecundity
  R0y_s <- rnorm(1, R0y, R0y*f_pct)
  R0a_s <- rnorm(1, R0a, R0a*f_pct)
  
  for (y in 2:length(Year)){
    
    # save new matrix
    A_dd <- A_s
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)

    # reduce fecundity
    R0y_dd <- R0y_s * density_factor 
    R0a_dd <- R0a_s * density_factor 
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    # cat("Adult Sex Ratio:", asr, "\n")
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,1:6] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_f
    A_dd[7,1:6] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_m
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer (after letting population increase for 10ish years)
    if (y > 0) {
      # Randomly select a random number to harvest
      h <- round(runif(1, 200000, 280000), digits=0)

      # Split bucks and does
      doe_h <- h * 0.54
      buck_h <- h - doe_h

      # Proportional does
      doe_r <- c(0.0735, 0.189, 0.1895, 0.183, 0.183, 0.182)*doe_h

      # Proportional bucks
      buck_r <- c(0.04, 0.11, 0.12, 0.23, 0.25, 0.25)*buck_h

      r <- c(doe_r, buck_r)

      N_next <- N_next - r
    } 
    
    # Save abundance
    deer.array[,y,i] <- pmax(N_next, 0)
    
    # Check buck:doe ratio
    af <- sum(deer.array[2:6,y,i])
    am <- sum(deer.array[8:12,y,i])
    
    if (sum(af) < 1) {
      warning("Female bottleneck: reproductive females lost!")
      # Optionally, set ASR to NA or a capped value
    }
    
  }
}
N.median <- apply(apply(deer.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.80)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  # geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

# Summarize sex ratio
ASR.median <- apply(asr_mat, 1, median)
ASR.20 <- apply(asr_mat, 1, quantile, probs = 0.20)
ASR.80 <- apply(asr_mat, 1, quantile, probs = 0.80)

ASR <- data.frame(Year, ASR.median, ASR.20, ASR.80)
ASR <- ASR[-nrow(ASR),]

# Plot sex ratio
ggplot(ASR) +
  geom_hline(yintercept=1, linetype=2) +
  coord_cartesian(ylim=c(0,1.5)) +
  geom_ribbon(aes(x=Year, ymin=ASR.20, ymax=ASR.80), alpha=0.5) +
  geom_line(aes(x=Year, y=ASR.median))

## Figure for plotting




