library(ggplot2)
library(ggtext)

source("00_functions.R")


theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:50

# Carrying capacity
K <- 2050224

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
  

#### Adjust fecundinty only 10-40% ####
adj <- seq(1.05, 3.5, by=0.01)

# Empty array to hold pop count
pig.array <- array(0, dim=c(6,length(Year),length(adj)))
pig.array[,1,] <- N0

for (i in 1:length(adj)) {
  
  A_adj <- pig.matrix
  A_adj[1, ] <- A_adj[1, ] * adj[i] 
  A_adj[4, ] <- A_adj[4, ] * adj[i] 
  
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

fec.incr <- apply(pig.array,c(2,3),sum)
fec.incr <- as.data.frame(fec.incr)
names(fec.incr) <- adj

fec.incr <- fec.incr |>
  as_tibble() |>
  pivot_longer(1:length(adj),names_to="increase", values_to="Nt") |>
  arrange(increase)

# Calculate net change
fec.net <- fec.incr |>
  group_by(increase) |>
  mutate(net=c(NA, diff(Nt))) |>
  # drop NA rows
  filter(!is.na(net))

ggplot(fec.net) +
  # coord_cartesian(ylim=c(0,309000)) +
  geom_vline(xintercept=K) +
  geom_hline(yintercept=308887, linetype=2, color="red") +
  geom_line(aes(x=Nt, y=net, group=increase, color=increase)) +
  theme(legend.position="none")
  


#### Adjust Survival only 5-50% ####
adj <- seq(1.05, 1.5, by=0.01)

# Empty array to hold pop count
pig.array <- array(0, dim=c(6,length(Year),length(adj)))
pig.array[,1,] <- N0

for (i in 1:length(adj)) {
  
  A_adj <- pig.matrix
  A_adj[2:3,1:3] <- A_adj[2:3,1:3] * adj[i] 
  A_adj[5:6,4:6] <- A_adj[5:6,4:6] * adj[i] 
  
  a <- calibrate_a(A_adj, N0, K)
  cat("Calibrated a =", a, "\n")
  
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

surv.incr <- apply(pig.array,c(2,3),sum)
surv.incr <- as.data.frame(surv.incr)
names(surv.incr) <- adj

surv.incr <- surv.incr |>
  as_tibble() |>
  pivot_longer(1:length(adj),names_to="increase", values_to="Nt") |>
  arrange(increase)

# Calculate net change
surv.net <- surv.incr |>
  group_by(increase) |>
  mutate(net=c(NA, diff(Nt))) |>
  # drop NA rows
  filter(!is.na(net))

ggplot(surv.net) +
  # coord_cartesian(ylim=c(0,309000)) +
  geom_hline(yintercept=308887, linetype=2, color="red") +
  geom_line(aes(x=Nt, y=net, group=increase, color=increase)) +
  theme(legend.position="none")


#### Adjust Survival and Fecundity 5-90% ####
adj.f <- seq(0.95, 3, by=0.1)
adj.s <- seq(0.75, 1.5, by=0.1)

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
  geom_vline(xintercept=K) +
  geom_hline(yintercept=308887, linetype=2, color="red") +
  geom_line(aes(x=Nt, y=net, group=increase, color=increase)) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  theme(legend.position="none")


both.net |>
  group_by(increase) |>
  reframe(max_change=max(net)) |>
  slice_max(max_change)


max_pop <- both.net |>
  filter(increase=="V21")


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
  h <- rnorm(1, 280000, 15000)
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

