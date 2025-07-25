library(ggplot2)
library(ggtext)

source("00_functions.R")

theme_set(theme_bw())

set.seed(123)
##Stochastistic model
Year <- 1:100

# Carrying capacity
# K <- 2347197*2
K <- 1512638*2
# K <- 2640597

# Set up deer matrix
deer.matrix <- matrix(0,12,12)

# Fecundity
FecundityF <- c(0, 0.566, 0.658, 0.658, 0.658, 0.658)
FecundityM <- c(0, 0.736, 0.854, 0.854, 0.854, 0.854)

# Survival
# Females
deer.matrix[1,1:6] <- FecundityF
deer.matrix[7,1:6] <- FecundityM
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

lambda <- Re(eigen(deer.matrix)$values[1])

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0, lambda, K)
cat("Estimated theta:", theta, "\n")

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
  density_factor <- 1 / (1 + (N_t / K)^theta)
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
  coord_cartesian(ylim=c(0,300000)) +
  geom_hline(yintercept=250000, linetype=2, color="red") +
  geom_hline(yintercept=msy, col="#21918c") +
  geom_line(aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme_classic() +
  theme(axis.title=element_markdown(face="italic"),
        panel.border=element_rect(fill=NA, color="black", linewidth=0.5))

#### Adjust Survival and Fecundity ####
adj.f <- seq(1, 2, by=0.1)  # Fecundity
adj.sF <- seq(0.8, 1.4, by=0.1)  # Fawn survival (male + female)  
adj.sYm <- seq(1, 1.6, by=0.1) # Yearling male survival  sYm=adj.sYm,
adj.sAm3 <- seq(1, 1.6, by=0.1)  # 3-yr-old male survival  
adj.sAm <- seq(1, 2.2, by=0.1) # adult male survival  
adj.sAf <- seq(1, 1.2, by=0.1) # Adult female survival

# Create all combinations of adjustments
adj <- expand.grid(f=adj.f, sF=adj.sF, saf=adj.sAf, sYm=adj.sYm, sm3=adj.sAm3, sm=adj.sAm) #

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
  
  # check buck to doe ratio (see Nagy-Reis et al. 2021 - 1.20)
  # stable stage of the matrix
  # w <- Re(eigen(A_adj)$vectors[, 1])
  # w <- w / sum(w) # normalize to equal 1
  # # sum males
  # wf <- sum(w[3:6]*deer.array[3:6,y-1,i])
  # wm <- sum(w[9:12]*deer.array[9:12,y-1,i])
  # # buck to doe ratio
  # bdr <- wm/wf
  
  # Put survival rates into a vector
  surv <- sum(c(A_adj[2,1], A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6],
                A_adj[8,7], A_adj[9,8], A_adj[10,9], A_adj[11,10], A_adj[12,11], A_adj[12,12]) < 1)
  
  if (lambda > 1.18 & (surv == 12)) {

    i_vec <- c(i_vec, i)
    
    lam_vec <- c(lam_vec, lambda)
    
    # Calibrate theta
    theta <- estimate_theta_opt(N0, lambda, K)

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
      density_factor <- 1 / (1 + (N_t / K)^theta)
      
      # save new matrix
      A_dd <- A_adj
      A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
      A_dd[7, ] <- A_dd[7, ] * density_factor # reduce fecundity
      
      # Calculate new pop size
      deer.array[,y,i] <- A_dd %*% deer.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
    }
    
    # Check buck:doe ratio
    af <- apply(deer.array[3:6,51:100,i], 1, median)
    am <- apply(deer.array[9:12,51:100,i], 1, median)
    bdr <- sum(am)/sum(af)
    print(bdr)
    
    # if (bdr >= 0.7) {
      print(i)
    # }
    
  }
  
}

# drop slices that don't have i index in i_vec
deer.array <- deer.array[,,i_vec]

both.incr <- apply(deer.array,c(2,3),sum)
both.incr <- as.data.frame(both.incr)
# names(both.incr) <- adj.f

both.incr <- both.incr |>
  as_tibble() |>
  pivot_longer(1:length(i_vec),names_to="increase", values_to="Nt") |>
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
  geom_hline(yintercept=238000, linetype=2, color="red") +
  geom_line(aes(x=Nt, y=net, group=increase, color=increase), alpha=0.5, linewidth=0.5) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  labs(x="N<sub>t</sub>", y="N<sub>t-1</sub> - N<sub>t</sub>") +
  theme(legend.position="none",
        axis.title=element_markdown(face="italic"))


# Which combination had MSY closest to observed harvest

both.diff <- both.net |>
  mutate(msy_diff=238000-net,
         msy_diff=abs(msy_diff)) |>
  group_by(increase) |>
  reframe(min_change=min(msy_diff)) |>
  slice_min(min_change)


both.net |>
  group_by(increase) |>
  reframe(max_change=max(net)) |>
  slice_max(max_change)

max_pop <- both.net |>
  filter(increase=="V8425")

## Plot line with greatest MSY
ggplot() +
  geom_vline(xintercept=K) +
  geom_hline(yintercept=238000, linetype=2, color="red") +
  # geom_smooth(aes(x=Nt, y=net, group=increase), color="gray", alpha=0.1, linewidth=0.1) +
  geom_line(data=df, aes(x=Nt, y=net), linewidth=1) +
  geom_line(data=max_pop, aes(x=Nt, y=net), linewidth=1, color="#21918c") + #, se=FALSE
  theme_classic() +
  theme(legend.position="none")


both.incr[both.incr$increase=="V8425",]

# Population based on matrix with greatest MSY
i <- i_vec[8425]
# i <- 310936

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
# A_adj[9,8] <- 0.82
# 3 yr old males
A_adj[10,9] <- deer.matrix[10,9]*adj[i,5]
# A_adj[10,9] <- 0.82
# Adult male survival
for (s in 10:11) {
  A_adj[s+1,s] <- deer.matrix[s+1,s]*adj[i,6]
  # A_adj[s+1,s] <- 0.82
}
# Survive & stay (oldest males)
A_adj[12,12] <- deer.matrix[12,12] * adj[i,6]
# A_adj[12,12] <- 0.82


print(Re(eigen(A_adj)$values[1]))


# Empty array to hold pop count
Year <- 1:100
deer.array <- matrix(0,nrow=12, ncol=length(Year))

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0, lambda, K)
cat("Estimated theta:", theta, "\n")

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K

deer.array[,1] <- N0

for (y in 2:length(Year)){
  
  # Density dependent adjustment in fecundity
  # What was abundance at time step t-1
  N_t <- sum(deer.array[,y-1])
  
  # Calcuate density factor
  density_factor <- 1 / (1 + (N_t / K)^theta)
  
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


ggmatplot::ggmatplot(t(deer.array), plot_type="line") 

#### Add in stochasticity in vital rates ####
set.seed(1)
##Stochastistic model
Year <- 1:100
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.01
s_pct_m <- 0.01
f_pct <- 0.01

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.5 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0


# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0, lambda, K)
cat("Estimated theta:", theta, "\n")
# theta <- 5.75


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

  ## Some stochasticity in fecundity
  A_s[1,1:3] <- rnorm(3, A_adj[1,1:3], A_adj[1,1:3]*s_pct_m)
  A_s[7,1:3] <- rnorm(3, A_adj[7,1:3], A_adj[7,1:3]*s_pct_m)
  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(deer.array[,y-1,i])
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (N_t / K)^theta)
    
    # save new matrix
    A_dd <- A_s
    A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
    A_dd[7, ] <- A_dd[7, ] * density_factor # reduce fecundity
    
    # Calculate new pop size
    N_next <- A_dd %*% deer.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer
    
    # Randomly select a random number to harvest
    h <- round(rnorm(1, 220000, 20000), digits=0)
    
    # Split bucks and does
    doe_h <- h * 0.54
    buck_h <- h - doe_h
    
    # Proportional does
    doe_r <- c(0.0735, 0.189, 0.1895, 0.183, 0.183, 0.182)*doe_h
    
    # Proportional bucks
    buck_r <- c(0.04, 0.11, 0.12, 0.23, 0.25, 0.25)*buck_h
    
    r <- c(doe_r, buck_r)
    
    N_next <- N_next - r
    
    # Save abundance
    deer.array[,y,i] <- pmax(N_next, 0)
    pop.array[,y]  <- A %*% pop.array[,y-1] 
  }
}
N.median <- apply(apply(deer.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.80)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2000000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 




#### Decrease K over time - harvest ~240000 ####

K_t <- K * exp(-0.0001 * (0:(max(Year) - 1)))
# K_t <- K * exp(-0.002 * (0:(max(Year) - 1)))
K_t <- seq(K, K - K/4, length.out=length(Year))

A_adj[10,9] <- 0.86
A_adj[11,10] <- 0.86
A_adj[12,11] <- 0.86
A_adj[12,12] <- 0.86

Re(eigen(A_adj)$values[1])

set.seed(1)
##Stochastistic model
Year <- 1:100
Sims <- 1000

# Variation in Survival
s_pct_f <- 0.01
s_pct_m <- 0.01
f_pct <- 0.01

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
w <- Re(eigen(A_adj)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# Set up initial popualtion size
N0 <- w * 0.75 * K  # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# theta <- 3.5
theta <- 10#4.3

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
  A_s[7,1:3] <- rnorm(3, A_adj[7,1:3], A_adj[7,1:3]*s_pct_m)
  
  for (y in 2:length(Year)){
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    N_t <- sum(deer.array[,y-1,i])
    
    # Calcuate density factor
    # density_factor <- 1 / (1 + (N_t / K)^theta)
    density_factor <- 1 / (1 + (N_t / K_t[y])^theta)
    
    # save new matrix
    A_dd <- A_s
    A_dd[1, ] <- A_dd[1, ] * density_factor # reduce fecundity
    A_dd[7, ] <- A_dd[7, ] * density_factor # reduce fecundity
    
    # print(Re(eigen(A_dd)$values[1]))
    
    # Calculate new pop size
    N_next <- A_dd %*% deer.array[,y-1,i] # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer
    
    # Randomly select a random number to harvest
    h <- round(rnorm(1, 240000, 20000), digits=0)
    
    # Split bucks and does
    doe_h <- h * 0.54
    buck_h <- h - doe_h
    
    # Proportional does
    doe_r <- c(0.0735, 0.189, 0.1895, 0.183, 0.183, 0.182)*doe_h
    
    # Proportional bucks
    buck_r <- c(0.04, 0.11, 0.12, 0.23, 0.25, 0.25)*buck_h
    
    r <- c(doe_r, buck_r)
    
    N_next <- N_next - r
    
    # Save abundance
    deer.array[,y,i] <- pmax(N_next, 0) 
  }
}
N.median <- apply(apply(deer.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.80)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,3500000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 



