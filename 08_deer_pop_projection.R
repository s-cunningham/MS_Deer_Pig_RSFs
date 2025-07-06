library(ggplot2)


set.seed(123)
##Stochastistic model
Year <- 1:51
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.0001
s_pct_m <- 0.02
a <- 1.35
b <- 1

# Carrying capacity
K <- 600000

# Set up deer matrix
deer.array <- array(0,dim=c(12,length(Year),Sims))
deer.matrix <- matrix(0,12,12)

# Fecundity
Fecundity <- c(0, 0.65, 0.76, 0.76, 0.76, 0.76)

# Survival
# Females
deer.matrix[c(1,7),1:6] <- Fecundity
deer.matrix[2,1] <- 0.52
deer.matrix[3,2] <- 0.93
deer.matrix[4,3] <- 0.84
deer.matrix[5,4] <- 0.84
deer.matrix[6,5] <- 0.84
deer.matrix[6,6] <- 0.84
# Males
deer.matrix[8,7] <- 0.52
deer.matrix[9,8] <- 0.82
deer.matrix[10,9] <- 0.63
deer.matrix[11,10] <- 0.53
deer.matrix[12,11] <- 0.44
deer.matrix[12,12] <- 0.49

# Calculate stable stage distribution
w <- Re(eigen(deer.matrix)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

deer.array[,1,] <- N0

for(i in 1:Sims){
  
  ## Stochasticity on Survival
  A_s <- deer.matrix
  # Female survival
  for (s in 1:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, deer.matrix[s+1,s], s_pct_f)
  }
  # Survive & stay
  A_s[6,6] <- rnorm(1, deer.matrix[6,6], s_pct_f)
  # Male survival
  for (s in 1:5) {
    # Survive & go
    A_s[s+1,s] <- rnorm(1, deer.matrix[s+1,s], s_pct_m)
  }
  # Survive & stay
  A_s[12,12] <- rnorm(1, deer.matrix[12,12], s_pct_m)

  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(deer.array[,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + a * (N_t / K))
    # density_factor <- exp(-b * (N_t / K))

    # save new matrix
    A_dd <- A_s
    A_dd[c(1,7), ] <- A_s[c(1,7), ] * density_factor # reduce fecundity
    
    deer.array[,y,i]<-A_dd%*%deer.array[,y-1,i] #Make sure to multiply matrix x vector (not vice versa)
  }
  
}
N.median<-apply(apply(deer.array,c(2,3),sum),1,median)
N.20pct<-apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct<-apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.80)

results<-data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
