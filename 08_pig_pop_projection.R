library(ggplot2)

source("00_functions.R")

set.seed(123)
##Stochastistic model
Year <- 1:50
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.05
s_pct_m <- 0.05

# Carrying capacity
K <- 2034396

# Set up pig matrix
pig.array <- array(0,dim=c(6,length(Year),Sims))
pig.matrix <- matrix(0,6,6)

# Fecundity
Fecundity <- c(0.54, 2.24, 2.24)

# Survival
# Females
pig.matrix[1,1:3] <- Fecundity
pig.matrix[4,1:3] <- Fecundity
pig.matrix[2,1] <- 0.3
pig.matrix[3,2] <- 0.35
pig.matrix[3,3] <- 0.35
# Males
pig.matrix[5,4] <- 0.18
pig.matrix[6,5] <- 0.23
pig.matrix[6,6] <- 0.23

# Calculate stable stage distribution
w <- Re(eigen(pig.matrix)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.01 * K  # e.g., start at 50% of K

# Calibrate density dependence parameter
a <- calibrate_a(pig.matrix, N0, K)
cat("Calibrated a =", a, "\n")

# Empty array to hold pop count
pig.array[,1,] <- N0

for(i in 1:Sims){
  
  ## Stochasticity on Survival
  A_s <- pig.matrix
  # Female survival
  for (s in 1:2) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, pig.matrix[s+1,s], s_pct_f)
  }
  # Survive & stay
  A_s[3,3] <- rnorm(1, pig.matrix[3,3], s_pct_f)
  # Male survival
  for (s in 4:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, pig.matrix[s+1,s], s_pct_m)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, pig.matrix[6,6], s_pct_m)

  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(pig.array[,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + a * (N_t / K))

    # save new matrix
    A_dd <- A_s
    A_dd[1, ] <- A_s[1, ] * density_factor # reduce fecundity
    A_dd[4, ] <- A_s[4, ] * density_factor # reduce fecundity
    
    ## Harvest more pigs 
    N <- pig.array[,y-1,i] 
    
    # Calculate new pop size
    pig.array[,y,i]<- A_dd %*% N # Make sure to multiply matrix x vector (not vice versa)
  }
}
N.median <- apply(apply(pig.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.10)
N.80pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  # coord_cartesian(ylim=c(K/2, (K+K/2))) +
  coord_cartesian(ylim=c(0, (K+K/2))) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 



#### Maximum sustained yield ####
net <- results$N.median[2:nrow(results)] - results$N.median[1:(nrow(results)-1)]

peak_idx <- which.max(net)
msy <- net[peak_idx]
pop_at_msy <- results$N.median[peak_idx+1]

plot(results$N.median[2:nrow(results)], net, type="l")
abline(h = msy, lty = 2, col = "blue")  # MSY line
abline(v = pop_at_msy, lty = 2, col = "red")  # N at MSY




