library(ggplot2)
library(ggtext)

source("00_functions.R")

theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:100

# Carrying capacity
# K <- 2347197
K <- 1512638
# K <- 1929918

# Set up deer matrix
deer.matrix <- matrix(0,12,12)

## Survival
# Females
deer.matrix[2,1] <- 0.60
deer.matrix[3,2] <- 0.93
deer.matrix[4,3] <- 0.92
deer.matrix[5,4] <- 0.92
deer.matrix[6,5] <- 0.92
deer.matrix[6,6] <- 0.90
# Males
deer.matrix[8,7] <- 0.60
deer.matrix[9,8] <- 0.83
deer.matrix[10,9] <- 0.84
deer.matrix[11,10] <- 0.84
deer.matrix[12,11] <- 0.84
deer.matrix[12,12] <- 0.70

## Fecundity
# Maximum fecundity
R0a <- 2#1.8
R0y <- 1.6#1.4

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
  coord_cartesian(ylim=c(0,300000)) +
  geom_hline(yintercept=250000, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))

#### Start harvesting ####

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


ggmatplot::ggmatplot(t(deer.array), plot_type="line", linewidth=1) 

#### Base Model - stochasticity in vital rates ####
set.seed(1)
##Stochastistic model
Year <- 1:50
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.1
s_pct_m <- 0.1
f_pct <- 0.1

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")

Ka <- 1345648
Ka <- Ka - (Ka*0.15)
Kf <- Ka * 0.6

for (i in 1:Sims) {
  
  ## Stochasticity on Survival (NEED TO FIX FOR DEER)
  A_s <- A_adj
  # Female survival
  for (s in 1:2) {
    # Survive & go
    A_s[s+1,s] <-  + rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_f)
  }
  # Survive & stay
  A_s[3,3] <- rnorm(1, A_adj[3,3], A_adj[3,3]*s_pct_f)
  # Male survival
  for (s in 4:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_m)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, A_adj[6,6], A_adj[6,6]*s_pct_m)
  

  for (y in 2:length(Year)){

    # save new matrix
    A_dd <- A_s
        
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # reduce fecundity
    R0y_dd <- rnorm(1, R0y, R0y*f_pct) * density_factor 
    R0a_dd <- rnorm(1, R0a, R0a*f_pct) * density_factor 
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    # cat("Adult Sex Ratio:", asr, "\n")
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_f
    A_dd[7, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_m
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer (after letting population increase for 10ish years)
    if (y > 1) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 250000, 25000), digits=0)
      
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
N.20pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.10)
N.80pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,3000000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  geom_hline(yintercept=1750000, color="red", linetype=3) +
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
  coord_cartesian(ylim=c(0,1.5)) +
  geom_ribbon(aes(x=Year, ymin=ASR.20, ymax=ASR.80), alpha=0.5) +
  geom_line(aes(x=Year, y=ASR.median))

#### Take more does - stochasticity in vital rates ####
set.seed(1)
##Stochastistic model
Year <- 1:50
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.1
s_pct_m <- 0.1
f_pct <- 0.1

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")

Ka <- 1345648
Ka <- Ka - (Ka*0.15)
Kf <- Ka * 0.6

for (i in 1:Sims) {
  
  ## Stochasticity on Survival (NEED TO FIX FOR DEER)
  A_s <- A_adj
  # Female survival
  for (s in 1:2) {
    # Survive & go
    A_s[s+1,s] <-  + rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_f)
  }
  # Survive & stay
  A_s[3,3] <- rnorm(1, A_adj[3,3], A_adj[3,3]*s_pct_f)
  # Male survival
  for (s in 4:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_m)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, A_adj[6,6], A_adj[6,6]*s_pct_m)
  
  
  for (y in 2:length(Year)){
    
    # save new matrix
    A_dd <- A_s
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # reduce fecundity
    R0y_dd <- rnorm(1, R0y, R0y*f_pct) * density_factor 
    R0a_dd <- rnorm(1, R0a, R0a*f_pct) * density_factor 
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    # cat("Adult Sex Ratio:", asr, "\n")
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_f
    A_dd[7, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_m
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer (after letting population increase for 10ish years)
    if (y > 1) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 250000, 25000), digits=0)
      
      # Split bucks and does
      doe_h <- h * 0.59
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
N.20pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.10)
N.80pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,3000000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  geom_hline(yintercept=1750000, color="red", linetype=3) +
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
  coord_cartesian(ylim=c(0,1.5)) +
  geom_ribbon(aes(x=Year, ymin=ASR.20, ymax=ASR.80), alpha=0.5) +
  geom_line(aes(x=Year, y=ASR.median))

#### Take more bucks - stochasticity in vital rates ####
set.seed(1)
##Stochastistic model
Year <- 1:50
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.1
s_pct_m <- 0.1
f_pct <- 0.1

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")

Ka <- 1345648
Ka <- Ka - (Ka*0.15)
Kf <- Ka * 0.6

for (i in 1:Sims) {
  
  ## Stochasticity on Survival (NEED TO FIX FOR DEER)
  A_s <- A_adj
  # Female survival
  for (s in 1:2) {
    # Survive & go
    A_s[s+1,s] <-  + rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_f)
  }
  # Survive & stay
  A_s[3,3] <- rnorm(1, A_adj[3,3], A_adj[3,3]*s_pct_f)
  # Male survival
  for (s in 4:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_m)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, A_adj[6,6], A_adj[6,6]*s_pct_m)
  
  
  for (y in 2:length(Year)){
    
    # save new matrix
    A_dd <- A_s
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # reduce fecundity
    R0y_dd <- rnorm(1, R0y, R0y*f_pct) * density_factor 
    R0a_dd <- rnorm(1, R0a, R0a*f_pct) * density_factor 
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    # cat("Adult Sex Ratio:", asr, "\n")
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_f
    A_dd[7, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * s_m
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer (after letting population increase for 10ish years)
    if (y > 1) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 250000, 25000), digits=0)
      
      # Split bucks and does
      doe_h <- h * 0.49
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
N.20pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.10)
N.80pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,3000000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  geom_hline(yintercept=1750000, color="red", linetype=3) +
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
  coord_cartesian(ylim=c(0,1.5)) +
  geom_ribbon(aes(x=Year, ymin=ASR.20, ymax=ASR.80), alpha=0.5) +
  geom_line(aes(x=Year, y=ASR.median))



#### Decrease K over time - harvest ~240000 ####
K_t <- seq(K, K - K/4, length.out=length(Year))
# Kf_t <- (Ka/2) * exp(-0.002 * (0:(max(Year) - 1)))
# Kf_t <- (K/2) * exp(-0.03 * (0:(max(Year)/10 - 1)))
# Kf_t <- rep(Kf_t, each=10)

Kf_t <- (K/2) * exp(-0.1 * (0:(max(Year)/50 - 1)))
Kf_t <- rep(Kf_t, each=50)
Kf_t <- Kf * exp(-0.005 * (0:(max(Year) - 1)))



Re(eigen(A_adj)$values[1])

set.seed(1)
##Stochastistic model
Year <- 1:50
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.15
s_pct_m <- 0.15
f_pct <- 0.1

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

theta <- 2.999
# theta <- 2

for (i in 1:Sims) {
  
  ## Stochasticity on Survival (NEED TO FIX FOR DEER)
  A_s <- A_adj
  # Female survival
  for (s in 1:2) {
    # Survive & go
    A_s[s+1,s] <-  + rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_f)
  }
  # Survive & stay
  A_s[3,3] <- rnorm(1, A_adj[3,3], A_adj[3,3]*s_pct_f)
  # Male survival
  for (s in 4:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_m)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, A_adj[6,6], A_adj[6,6]*s_pct_m)
  
  # Build projection
  for (y in 2:length(Year)){
    
    # save new matrix
    A_dd <- A_s
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf_t[y])^theta)
    
    # reduce fecundity
    R0y_dd <- rnorm(1, R0y, R0y*f_pct) * density_factor 
    R0a_dd <- rnorm(1, R0a, R0a*f_pct) * density_factor 
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    cat("Adult Sex Ratio:", asr, "\n")
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    p_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    p_f <- 1 - p_m
    
    # Adjust sex ratio at birth
    A_dd[1, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * p_f
    A_dd[7, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * p_m
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer (after letting population increase for 10ish years)
    if (y > 5) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 240000, 35000), digits=0)
      
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
  geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 



#### Harvest more bucks, K = 1345648 ####

K <- 1345648
Kf <- K * 0.6

set.seed(1)
##Stochastistic model
Year <- 1:100
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.15
s_pct_m <- 0.15
f_pct <- 0.1

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

theta <- 2.999
# theta <- 2

for (i in 1:Sims) {
  
  ## Stochasticity on Survival (NEED TO FIX FOR DEER)
  A_s <- A_adj
  # Female survival
  for (s in 1:2) {
    # Survive & go
    A_s[s+1,s] <-  + rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_f)
  }
  # Survive & stay
  A_s[3,3] <- rnorm(1, A_adj[3,3], A_adj[3,3]*s_pct_f)
  # Male survival
  for (s in 4:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, A_adj[s+1,s], A_adj[s+1,s]*s_pct_m)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, A_adj[6,6], A_adj[6,6]*s_pct_m)
  
  # Build projection
  for (y in 2:length(Year)){
    
    # save new matrix
    A_dd <- A_s
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # reduce fecundity
    R0y_dd <- rnorm(1, R0y, R0y*f_pct) * density_factor 
    R0a_dd <- rnorm(1, R0a, R0a*f_pct) * density_factor 
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    cat("Adult Sex Ratio:", asr, "\n")
    
    # Calculate % males & females
    p_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    p_f <- 1 - p_m
    
    # Adjust sex ratio at birth
    A_dd[1, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * p_f
    A_dd[7, ] <- c(0, R0y_dd, R0a_dd, R0a_dd, R0a_dd, R0a_dd) * p_m
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer (after letting population increase for 10ish years)
    if (y > 5) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 250000, 25000), digits=0)
      
      # Split bucks and does
      doe_h <- h * 0.49
      buck_h <- h - doe_h
      
      # Proportional does
      doe_r <- c(0.0735, 0.189, 0.1895, 0.183, 0.183, 0.182)*doe_h
      
      # Proportional bucks
      buck_r <- c(0.04, 0.18, 0.18, 0.20, 0.20, 0.20)*buck_h
      
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
  geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 









