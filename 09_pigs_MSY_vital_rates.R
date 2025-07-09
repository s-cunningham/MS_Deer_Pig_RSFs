library(ggplot2)
library(ggtext)

source("00_functions.R")


theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:50

# Carrying capacity
K <- 2034396

# Set up pig matrix
pig.array <- matrix(0,nrow=6, ncol=length(Year))
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
pig.array[,1] <- N0

for (y in 2:length(Year)){
    
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  N_t <- sum(pig.array[,y-1])
  
  # Calcuate density factor
  # density_factor <- 1 / (1 + a * (N_t / K))
  density_factor <- 1 / (1 + a * (N_t / K))
  
  # save new matrix
  A_dd <- pig.matrix
  A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
  A_dd[4, ] <- A_dd[4, ] * density_factor # reduce fecundity

  # Calculate new pop size
  pig.array[,y]<- A_dd %*% pig.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
}

N.median<-apply(pig.array,2,sum)

results<-data.frame(Year,N.median)

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
  coord_cartesian(ylim=c(0,309000)) +
  geom_hline(yintercept=308887, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))
  
#### Adjust Survival and Fecundity 5-90% ####
adj.f <- seq(0.95, 3.0, by=0.1)
adj.s <- seq(0.75, 1.9, by=0.1)

# Create all combinations of adjustments
adj <- expand.grid(f=adj.f, s=adj.s)

# Empty array to hold pop count
pig.array <- array(0, dim=c(6,length(Year),nrow(adj)))
pig.array[,1,] <- N0

for (i in 1:nrow(adj)) {
  
  A_adj <- pig.matrix
  A_adj[1, ] <- A_adj[1, ] * adj[i,1] 
  A_adj[4, ] <- A_adj[4, ] * adj[i,1] 
  A_adj[2:3,1:3] <- pig.matrix[2:3,1:3] * adj[i,2] 
  A_adj[5:6,4:6] <- pig.matrix[5:6,4:6] * adj[i,2] 
  
  a <- calibrate_a(A_adj, N0, K)
  cat("Calibrated a =", a, "\n")

  print(Re(eigen(A_adj)$values[1]))
  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(pig.array[,y-1,i])

    # Calcuate density factor
    density_factor <- 1 / (1 + a * (N_t / K))
    
    # save new matrix
    A_dd <- A_adj
    A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
    A_dd[4, ] <- A_dd[4, ] * density_factor # reduce fecundity

    # Calculate new pop size
    pig.array[,y,i] <- A_dd %*% pig.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
  }
  
}

both.incr <- apply(pig.array,c(2,3),sum)
both.incr <- as.data.frame(both.incr)
# names(both.incr) <- adj.f

both.incr <- both.incr |>
  as_tibble() |>
  pivot_longer(1:nrow(adj),names_to="increase", values_to="Nt") |>
  arrange(increase)

# Calculate net change
both.net <- both.incr |>
  group_by(increase) |>
  mutate(net=c(NA, diff(Nt))) |>
  # drop NA rows
  filter(!is.na(net))

ggplot(both.net) +
  # coord_cartesian(ylim=c(0,309000)) +
  # geom_vline(xintercept=K) +
  geom_hline(yintercept=308887, linetype=2, color="red") +
  geom_smooth(aes(x=Nt, y=net, group=increase, color=increase), alpha=0.5, linewidth=0.5) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  theme(legend.position="none")


both.net |>
  group_by(increase) |>
  reframe(max_change=max(net)) |>
  slice_max(max_change)

# Which combination had MSY closest to observed harvest
max_pop <- both.net |>
  filter(increase=="V21")

## Plot line with greatest MSY
ggplot(both.net) +
  geom_hline(yintercept=308887, linetype=2, color="red") +
  geom_smooth(aes(x=Nt, y=net, group=increase), color="gray", alpha=0.1, linewidth=0.1) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  geom_smooth(data=max_pop, aes(x=Nt, y=net), linewidth=1, color="#21918c") +
  theme_classic() +
  theme(legend.position="none")


# Population based on matrix with greatest MSY
i <- 21

A_adj <- pig.matrix
A_adj[1, ] <- A_adj[1, ] * adj[i,1] 
A_adj[4, ] <- A_adj[4, ] * adj[i,1] 
A_adj[2:3,1:3] <- pig.matrix[2:3,1:3] * adj[i,2] 
A_adj[5:6,4:6] <- pig.matrix[5:6,4:6] * adj[i,2] 

a <- calibrate_a(A_adj, N0, K)
cat("Calibrated a =", a, "\n")

print(Re(eigen(A_adj)$values[1]))

# Empty array to hold pop count
Year <- 1:30
pig.array <- matrix(0,nrow=6, ncol=length(Year))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

pig.array[,1] <- N0

for (y in 2:length(Year)){
  
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  N_t <- sum(pig.array[,y-1])
  
  # Calcuate density factor
  density_factor <- 1 / (1 + a * (N_t / K))
  
  # save new matrix
  A_dd <- A_adj
  A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
  A_dd[4, ] <- A_dd[4, ] * density_factor # reduce fecundity
  
  # Calculate new pop size
  N_next <- A_dd %*% pig.array[,y-1] # Make sure to multiply matrix x vector (not vice versa)
  
  ## Harvest pigs
  # Randomly select a random number to harvest
  # h <- rnorm(1, 280000, 15000)
  h <- rnorm(1, 300000, 10000)
  print(h)
  r <- (N_next[,1] / sum(N_next[,1]))*h
  
  N_next <- N_next - r
  
  # Save abundance
  pig.array[,y] <- pmax(N_next, 0) 
}

N.median <- apply(pig.array,2,sum)

results <- data.frame(Year,N.median)

#Plot population with
ggplot(results) +
  coord_cartesian(ylim=c(0, 2100000)) +
  geom_hline(yintercept=K) +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (males + females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


#### Add in stochasticity in vital rates ####

##Stochastistic model
Year <- 1:30
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.08
s_pct_m <- 0.08
f_pct <- 0.08

# Empty array to hold pop count
pig.array <- array(0,dim=c(6,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

# Fill starting population
pig.array[,1,] <- N0

set.seed(1)
for (i in 1:Sims) {
  
  ## Stochasticity on Survival
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
  
  ## Some stochasticity in fecundity
  A_s[1,1:3] <- rnorm(3, A_adj[1,1:3], A_adj[1,1:3]*s_pct_m)
  A_s[4,1:3] <- rnorm(3, A_adj[4,1:3], A_adj[4,1:3]*s_pct_m)
  
  # Calibrate a
  a <- calibrate_a(A_s, N0, K)
  cat("Calibrated a =", a, "\n")
  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(pig.array[,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + a * (N_t / K))
    
    # save new matrix
    A_dd <- A_s
    A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
    A_dd[4, ] <- A_dd[4, ] * density_factor # reduce fecundity
    
    # Calculate new pop size
    N_next <- A_dd %*% pig.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest pigs
    # Randomly select a random number to harvest
    # h <- rnorm(1, 280000, 15000)
    h <- rnorm(1, 300000, 10000)
    print(h)
    r <- (N_next[,1] / sum(N_next[,1]))*h
    
    N_next <- N_next - r
    
    # Save abundance
    pig.array[,y,i] <- pmax(N_next, 0) 
  }
}
N.median <- apply(apply(pig.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.80)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0, 2100000)) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

#### Harvest more ####
Year <- 1:30
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.05
s_pct_m <- 0.05
f_pct <- 0.05

# Empty array to hold pop count
pig.array <- array(0,dim=c(6,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

# Pct harvest
h_rate <- matrix(c(0.7,0.5,0.5,0.2,0.5,0.5), nrow=6, ncol=1)

# Save # harvested
n_harvest <- array(0,dim=c(6,length(Year),Sims))

# Fill starting population
pig.array[,1,] <- N0

set.seed(123)
for (i in 1:Sims) {
  
  ## Stochasticity on Survival
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
  
  ## Some stochasticity in fecundity
  A_s[1,1:3] <- rnorm(3, A_adj[1,1:3], A_adj[1,1:3]*s_pct_m)
  A_s[4,1:3] <- rnorm(3, A_adj[4,1:3], A_adj[4,1:3]*s_pct_m)
  
  # Calibrate a
  a <- calibrate_a(A_s, N0, K)

  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(pig.array[,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + a * (N_t / K))
    
    # save new matrix
    A_dd <- A_s
    A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
    A_dd[4, ] <- A_dd[4, ] * density_factor # reduce fecundity
    
    # Calculate new pop size
    N_next <- A_dd %*% pig.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
    
    if (y > 9) {
      ## Harvest pigs
      h <- N_next * h_rate
      N_next <- N_next - h
      
      # Save n harvested
      n_harvest[,y,i] <- h

    } else {
      h <- rnorm(1, 300000, 10000)
      print(h)
      r <- (N_next[,1] / sum(N_next[,1]))*h
      
      N_next <- N_next - r
    }
    # Save abundance
    pig.array[,y,i] <- pmax(N_next, 0) 
  }
}
N.median <- apply(apply(pig.array,c(2,3), sum),1, median)
N.20pct <- apply(apply(pig.array,c(2,3), sum), 1, quantile, probs=0.20)
N.80pct <- apply(apply(pig.array,c(2,3), sum), 1, quantile, probs=0.80)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0, 2100000)) +
  # geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

# median harvest
h.median <- apply(apply(n_harvest,c(2,3),sum),1,median)

