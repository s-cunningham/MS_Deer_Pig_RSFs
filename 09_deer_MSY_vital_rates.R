library(ggplot2)
library(ggtext)

source("00_functions.R")

theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:100

# Carrying capacity
K <- 2347197
# K <- 3168716
# K <- 2347197 + 500000 + 500000

# Set up deer matrix
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
N0 <- w * 0.1 * K  # e.g., start at 50% of K

# Calibrate density dependence parameter
a <- calibrate_a(deer.matrix, N0, K)
cat("Calibrated a =", a, "\n")

# print(Re(eigen(deer.matrix)$values[1]))

# Empty array to hold pop count
deer.array <- matrix(0,nrow=12, ncol=length(Year))
deer.array[,1] <- N0

lambda <- Re(eigen(deer.matrix)$values[1])  

for (y in 2:length(Year)){
  
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  N_t <- sum(deer.array[,y-1])
  
  # Calcuate density factor
  # density_factor <- 1 / (1 + a * (N_t / K))
  density_factor <- max(1 / (1 + a * (N_t / K)), 0)
  # density_factor <- max(1 - a * (N_t / K), 0)
  
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
  coord_cartesian(ylim=c(0,221000)) +
  geom_hline(yintercept=220989, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))

#### Adjust Survival and Fecundity ####
adj.f <- seq(1, 1.8, by=0.1)  # Fecundity
adj.sF <- seq(0.5, 1.5, by=0.05)  # Fawn survival (male + female)  
adj.sYm <- seq(0.8, 1.2, by=0.05) # Yearling male survival 
adj.sAm3 <- seq(0.8, 1.5, by=0.05)  # 3-yr-old male survival  
adj.sAm <- seq(0.8, 1.6, by=0.05) # adult male survival  
adj.sAf <- seq(0.9, 1.16, by=0.04) # Adult female survival  

# Create all combinations of adjustments
adj <- expand.grid(f=adj.f, sF=adj.sF, sf=adj.sAf, sYm=adj.sYm, sm3=adj.sAm3, sm=adj.sAm)

adj <- adj |>
  distinct()

# Empty array to hold pop count
deer.array <- array(0, dim=c(12,length(Year),nrow(adj)))

# empty vector for matrices that might work
i_vec <- c()
lam_vec <- c()

for (i in 1:nrow(adj)) {
  
  A_adj <- deer.matrix
  # Fecundity
  A_adj[1, ] <- deer.matrix[1, ] * adj[i,1] 
  A_adj[7, ] <- deer.matrix[7, ] * adj[i,1] 
  
  # Fawn survival
  A_adj[2,1] <- deer.matrix[2,1] * adj[i,2]
  A_adj[8,7] <- deer.matrix[2,1] * adj[i,2]
  
  # Adult female survival
  for (s in 3:5) {
    A_adj[s+1,s] <- deer.matrix[s+1,s]*adj[i,3]
  }
  # Survive & stay (oldest females)
  A_adj[6,6] <- deer.matrix[s+1,s]*adj[i,3] 
  # Yearling male survival 
  A_adj[9,8] <- deer.matrix[9,8]*adj[i,4] 
  # 3 yr old males
  A_adj[10,9] <- deer.matrix[10,9]*adj[i,5] 
  # Adult male survival
  for (s in 10:11) {
    A_adj[s+1,s] <- deer.matrix[s+1,s]*adj[i,6]
  }
  # Survive & stay (oldest males)
  A_adj[12,12] <- deer.matrix[12,12] * adj[i,6]

  # Check lambda
  lambda <- Re(eigen(A_adj)$values[1])
  
  if (lambda >1.4 ) {
    print(i)
    
    i_vec <- c(i_vec, i)
    
    lam_vec <- c(lam_vec, lambda)
    
    a <- calibrate_a_opt(A_adj, N0, K)
    # cat("Calibrated a =", a, "\n")
    
    # Calculate stable stage distribution
    w <- Re(eigen(A_adj)$vectors[, 1])
    w <- w / sum(w) # normalize to equal 1
    
    # Set up initial popualtion size
    N0 <- w * 0.1 * K
    
    deer.array[,1,i] <- N0
    
    for (y in 2:length(Year)){
      
      # Density dependent adjustment in fecundity
      # What was abundance at time step t-1
      N_t <- sum(deer.array[,y-1,i])
      
      # Calcuate density factor
      # density_factor <- max(1 - a * (N_t / K), 0)
      density_factor <- max(1 / (1 + a * (N_t / K)), 0)
      
      # save new matrix
      A_dd <- A_adj
      A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
      A_dd[7, ] <- A_dd[7, ] * density_factor  # reduce fecundity
      
      # Calculate new pop size
      deer.array[,y,i] <- A_dd %*% deer.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
    }
    
  }
 
}

# drop slices that don't have i index in i_vec
deer.array <- deer.array[,,i_vec]

both.incr <- apply(deer.array,c(2,3),sum,na.rm=TRUE)
both.incr <- as.data.frame(both.incr)
# names(both.incr) <- adj.f

both.incr <- both.incr |>
  as_tibble() |>
  pivot_longer(1:dim(deer.array)[3],names_to="increase", values_to="Nt") |>
  arrange(increase)

both.incr$matrix <- rep(i_vec, each=100)

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


# Which combination had MSY closest to observed harvest

both.diff <- both.net |>
  mutate(msy_diff=220989-net,
         msy_diff=abs(msy_diff)) |>
  group_by(increase) |>
  reframe(min_change=min(msy_diff)) |>
  slice_min(min_change)


max_pop <- both.net |>
  filter(increase=="V658387")

## Plot line with greatest MSY
ggplot() +
  geom_hline(yintercept=220989, linetype=2, color="red") +
  # geom_smooth(aes(x=Nt, y=net, group=increase), color="gray", alpha=0.1, linewidth=0.1) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  geom_smooth(data=max_pop, aes(x=Nt, y=net), linewidth=1, color="#21918c", se=FALSE) +
  theme_classic() +
  theme(legend.position="none")


# Population based on matrix with greatest MSY
i <- 658387

A_adj <- deer.matrix
# Fecundity
A_adj[1, ] <- deer.matrix[1, ] * adj[i,1] 
A_adj[7, ] <- deer.matrix[7, ] * adj[i,1] 

# Fawn survival
A_adj[2,1] <- deer.matrix[2,1] * adj[i,2]
A_adj[8,7] <- deer.matrix[2,1] * adj[i,2]

# Adult female survival
for (s in 3:5) {
  A_adj[s+1,s] <- deer.matrix[s+1,s]*adj[i,3]
}
# Survive & stay (oldest females)
A_adj[6,6] <- deer.matrix[s+1,s]*adj[i,3] 
# Adult male survival
for (s in 9:11) {
  A_adj[s+1,s] <- deer.matrix[s+1,s]*adj[i,4]
}
# Survive & stay (oldest males)
A_adj[12,12] <- deer.matrix[12,12] * adj[i,4]

print(Re(eigen(A_adj)$values[1]))

