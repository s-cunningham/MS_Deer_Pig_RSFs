library(ggplot2)
library(ggtext)

source("00_functions.R")

theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:100

# Carrying capacity
K <- 3168716

# Set up deer matrix
deer.array <- matrix(0,nrow=12, ncol=length(Year))
deer.matrix <- matrix(0,12,12)

# Fecundity
Fecundity <- c(0, 0.65, 0.76, 0.76, 0.76, 0.76)

# Survival
# Females
deer.matrix[1,1:6] <- Fecundity
deer.matrix[7,1:6] <- Fecundity
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
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.01 * K  # e.g., start at 50% of K

# Calibrate density dependence parameter
a <- calibrate_a(deer.matrix, N0, K)
cat("Calibrated a =", a, "\n")

# a <- 1.63

# Empty array to hold pop count
deer.array[,1] <- N0

lambda <- Re(eigen(deer.matrix)$values[1])  

for (y in 2:length(Year)){
  
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  N_t <- sum(deer.array[,y-1])
  
  # Calcuate density factor
  density_factor <- 1 / (1 + a * (N_t / K))
  # density_factor <- max(1 - (N_t / K)^theta, 0)
  # density_factor <- max(1 - (N_t / K), 0)
  
  # save new matrix
  A_dd <- deer.matrix
  A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
  A_dd[7, ] <- A_dd[7, ] * density_factor # reduce fecundity
  
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
  coord_cartesian(ylim=c(0,221000)) +
  geom_hline(yintercept=220989, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))

#### Adjust Survival and Fecundity ####
adj.f <- seq(1, 2.5, by=0.04)
adj.sm <- seq(0.9, 1.8, by=0.04)
adj.sf <- seq(0.8, 1.2, by=0.04)

# Create all combinations of adjustments
adj <- expand.grid(f=adj.f, sm=adj.sm, sf=adj.sf)

adj <- adj |>
  distinct()

# Empty array to hold pop count
deer.array <- array(0, dim=c(12,length(Year),nrow(adj)))
deer.array[,1,] <- N0

for (i in 1:nrow(adj)) {
  
  A_adj <- deer.matrix
  A_adj[1, ] <- deer.matrix[1, ] * adj[i,1] 
  A_adj[7, ] <- deer.matrix[7, ] * adj[i,1] 
  # Female survival
  for (s in 1:5) {
    A_adj[s+1,s] <- min(deer.matrix[s+1,s]*adj[i,3], 1)
  }
  A_adj[6,6] <- min(deer.matrix[s+1,s]*adj[i,3], 1) 
  # Male survival
  for (s in 1:5+6) {
    # Survive & go
    A_adj[s+1,s] <- min(deer.matrix[s+1,s], 1)
  }
  # Survive & stay
  A_adj[12,12] <- min(deer.matrix[12,12] * adj[i,2], 1)

  # Check lambda
  print(Re(eigen(A_adj)$values[1]))
  
  a <- calibrate_a(A_adj, N0*100, K)
  cat("Calibrated a =", a, "\n")
  
  for (y in 2:length(Year)){
    
    # check buck to doe ratio (see Nagy-Reis et al. 2021 - 1.20)
    # stable stage of the matrix
    w <- Re(eigen(A_adj)$vectors[, 1])
    w <- w / sum(w) # normalize to equal 1
    # sum males
    wf <- sum(w[1:6]*deer.array[1:6,y-1,i])
    wm <- sum(w[7:12]*deer.array[7:12,y-1,i])
    # buck to doe ratio
    bdr <- wf/wm
    
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(deer.array[,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + a * (N_t / K))
    
    # save new matrix
    A_dd <- A_adj
    A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
    A_dd[7, ] <- A_dd[7, ] * density_factor  # reduce fecundity
    
    # Calculate new pop size
    deer.array[,y,i] <- A_dd %*% deer.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
  }

    

}

# Remove rows with NAs


both.incr <- apply(deer.array,c(2,3),sum)
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
  geom_vline(xintercept=K) +
  geom_hline(yintercept=220989, linetype=2, color="red") +
  geom_line(aes(x=Nt, y=net, group=increase, color=increase)) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  theme(legend.position="none")




both.net |>
  group_by(increase) |>
  reframe(max_change=max(net)) |>
  slice_max(max_change)

# Which combination had MSY closest to observed harvest
max_pop <- both.net |>
  filter(increase=="V306")

## Plot line with greatest MSY
ggplot(both.net) +
  geom_hline(yintercept=220989, linetype=2, color="red") +
  geom_smooth(aes(x=Nt, y=net, group=increase), color="gray", alpha=0.1, linewidth=0.1) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  geom_smooth(data=max_pop, aes(x=Nt, y=net), linewidth=1, color="#21918c") +
  theme_classic() +
  theme(legend.position="none")


# Population based on matrix with greatest MSY
i <- 306

A_adj <- deer.matrix
A_adj[1, ] <- A_adj[1, ] * adj[i,1] 
A_adj[7, ] <- A_adj[7, ] * adj[i,1] 
# A_adj[2:6,1:6] <- deer.matrix[2:6,1:6] * adj[i,2] 
A_adj[8:12,7:12] <- deer.matrix[8:12,7:12] * adj[i,2] 
