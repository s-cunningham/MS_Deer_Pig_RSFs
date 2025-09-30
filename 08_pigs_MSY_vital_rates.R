library(ggplot2)
library(ggtext)
library(patchwork)
library(truncnorm)

source("00_functions.R")

theme_set(theme_bw())

### K at 27 pigs/ km2 - Largest area (Core + Moderate + Marginal) ####

## Stochastistic model
step <- 1:200

# Carrying capacity
K <- 1239278

Kf <- K * 0.5

# Set up pig matrix
A_base <- matrix(0,6,6)

## Survival
# Females
A_base[2,1] <- sqrt(0.3)
A_base[3,2] <- sqrt(0.35)
A_base[3,3] <- sqrt(0.35)
# Males
A_base[5,4] <- sqrt(0.18)
A_base[6,5] <- sqrt(0.23)
A_base[6,6] <- sqrt(0.23)

## Fecundity
# Maximum fecundity (how many piglets per litter?)
R0j <- 6.6*0.75
R0a <- 10.8

FecundityM <- c(R0j, R0a)
FecundityF <- c(R0j, R0a)

# Add to matrix
A_base[1,2:3] <- FecundityF 
A_base[4,2:3] <- FecundityM 

# # Calculate for post-breeding census (assume 1:1 sex ratio)
FecundityM <- c(R0j*A_base[3,2], R0a*A_base[3,3])*0.5
FecundityF <- c(R0j*A_base[3,2], R0a*A_base[3,3])*0.5

A_base[1,2:3] <- FecundityF
A_base[4,2:3] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(A_base)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(A_base)$values[1])

# How many pig are harvested?
observed_harvest <- 150000 

### Optimize survival and fecundities
# set up bounds
pig_lower <- c(1, 1.2, 1.2, 1.2, 1.2, 1.1)
pig_upper <- c(1.5, 1.5, 1.5, 1.6, 1.6, 2.8)

# run optimizer
set.seed(1)
pig_opt <- optim(par=c(runif(1, 1.01, 1.45), runif(4, 1, 1.45), runif(1, 1.01, 2.75)), fn=objective_fn_pigs, method="L-BFGS-B", lower=pig_lower, upper=pig_upper, control = list(trace = 1))

# What are optimized parameter values 
pig_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Female survival
A_adj[2,1] <- A_base[2,1]*pig_opt$par[1]
A_adj[3,2] <- A_base[3,2]*pig_opt$par[2]
A_adj[3,3] <- A_base[3,3]*pig_opt$par[3]

# Male survival
A_adj[5,4] <- A_base[5,4]*pig_opt$par[1]
A_adj[6,5] <- A_base[6,5]*pig_opt$par[4]
A_adj[6,6] <- A_base[6,6]*pig_opt$par[5]

## Fecundity
# Yearlings
A_adj[1,2] <- A_base[1,2]*(pig_opt$par[6])
A_adj[4,2] <- A_base[4,2]*(pig_opt$par[6])
# Adults
A_adj[1,3] <- A_base[1,3]*(pig_opt$par[6])
A_adj[4,3] <- A_base[4,3]*(pig_opt$par[6])

#### Run population model ####

set.seed(1)
##Stochastistic model
Year <- 1:40
Sims <- 1000

# Empty array to hold pop count
pig.array <- array(0,dim=c(6,length(Year),Sims))

# Calculate stable stage distribution
mnFs <- mean(c(A_adj[2,1], A_adj[3,2], A_adj[3,3]))
a2 <- A_adj
a2[1,1:3] <- a2[1,1:3]*0.5*mnFs
a2[4,1:3] <- a2[4,1:3]*0.5*mnFs
w <- Re(eigen(a2)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

lambda <- Re(eigen(a2)$values[1])

# Set up initial popualtion size
N0 <- w * 1 * K#*0.5 # e.g., start at 50% of K

# Fill starting population
pig.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:3], lambda, Kf)
cat("Estimated theta:", theta, "\n")
# theta <- 2.25

# Set up arrays to save realized demographic rates
realized_surv <- array(NA, dim=c(6, length(Year), Sims))
realized_surv[1:6,1,] <- c(A_adj[2,1], A_adj[3,2], A_adj[3,3], A_adj[5,4], A_adj[6,5], A_adj[6,6])

realized_fecundity <- array(NA, dim=c(4,  length(Year) - 1, Sims)) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,

for (i in 1:Sims) {
  
  A_s <- A_adj
  
  ## Stochasticity on Survival 
  # Female survival
  # Survive & go
  A_s[2,1] <- rtruncnorm(1, b=1, mean=A_adj[2,1], sd=0.005)
  A_s[3,2] <- rtruncnorm(1, b=1, mean=A_adj[3,2], sd=0.005)
  # Survive & stay
  A_s[3,3] <- rtruncnorm(1, a=0, b=1, mean=A_adj[3,3], sd=0.005)
  # Male survival
  A_s[5,4] <- rtruncnorm(1, b=1, mean=A_adj[5,4], sd=0.005)
  A_s[6,5] <- rtruncnorm(1, b=1, mean=A_adj[6,5], sd=0.005)
  # Survive & stay
  A_s[6,6] <- rtruncnorm(1, b=1, mean=A_adj[6,6], sd=0.005)
  
  # save new matrix
  A_dd <- A_s
  
  for (y in 2:length(Year)){
    
    ## Adjust fecundity by female survival
    R0j_s <- ifelse(is.na(A_s[1,2]*((realized_surv[2,y-1,i]))), 0, A_s[1,2]*(realized_surv[2,y-1,i]))
    R0a_s <- ifelse(is.na(A_s[1,3]*((realized_surv[3,y-1,i]))), 0, A_s[1,3]*(realized_surv[3,y-1,i]))
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(pig.array[1:3,y-1,i], na.rm=TRUE)
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)*theta)
    
    # Check adult age ratio
    asr <- sum(pig.array[5:6,y-1,i])/sum(pig.array[2:3,y-1,i], 1e-6)
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,2:3] <- c(R0j_s, R0a_s) * s_f * density_factor
    A_dd[4,2:3] <- c(R0j_s, R0a_s) * s_m * density_factor
    
    # Save fecundity
    realized_fecundity[1,y-1,i] <- A_dd[1,2] # Females per yearling female
    realized_fecundity[2,y-1,i] <- A_dd[1,3] # Females per adult female
    
    realized_fecundity[3,y-1,i] <- A_dd[4,2] # Males per yearling female
    realized_fecundity[4,y-1,i] <- A_dd[4,3] # Males per adult female
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% pig.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest pigs 
    if (y > 0 & min(N_next)>0 & sum(!is.na(N_next))==6) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 150000, 15000), digits=0)
      # h <- 100000
      
      sow_h <- h * 0.5
      boar_h <- h - sow_h
      
      # Proportional harvest
      r_f <- (N_next[1:3,1] / sum(N_next[1:3,1]))*sow_h
      r_m <- (N_next[4:6,1] / sum(N_next[4:6,1]))*boar_h
      
      r <- c(r_f, r_m)
      
      # remove
      N_next <- N_next - r
    } 
    
    # Save abundance
    pig.array[,y,i] <- pmax(N_next, 0)
    
    # ---- Realized survival ----
    
    # Female stages
    realized_surv[1, y, i] <- pig.array[2, y, i] / pig.array[1, y-1, i] # Piglet → yearling (F)

    # Stage 2: comes only from stage 1
    incoming_2 <- pig.array[1, y-1, i]
    realized_surv[2, y, i] <- ifelse(incoming_2 > 0,
                                     pig.array[2, y, i] / incoming_2,
                                     NA)
    
    # Stage 3: comes from stage 2 -> 3 and stage 3 -> 3
    incoming_3 <- pig.array[2, y-1, i] + pig.array[3, y-1, i]
    realized_surv[3, y, i] <- ifelse(incoming_3 > 0,
                                     pig.array[3, y, i] / incoming_3,
                                     NA)
    
    # Male stages
    realized_surv[4, y, i] <- pig.array[5, y, i] / pig.array[4, y-1, i] # Piglet → yearling (M)

    # Male stages 5–6 (adults)
    # Stage 5: comes only from stage 4
    incoming_5 <- pig.array[1, y-1, i]
    realized_surv[5, y, i] <- ifelse(incoming_5 > 0,
                                     pig.array[5, y, i] / incoming_5,
                                     NA)
    
    # Stage 6: comes from stage 5 -> 6 and stage 6 -> 6
    incoming_6 <- pig.array[5, y-1, i] + pig.array[6, y-1, i]
    realized_surv[6, y, i] <- ifelse(incoming_6 > 0,
                                     pig.array[6, y, i] / incoming_6,
                                     NA)
    
    # Check buck:doe ratio
    af <- sum(pig.array[2:3,y,i], na.rm=TRUE)
    am <- sum(pig.array[5:6,y,i], na.rm=TRUE)
    
    if (sum(af) < 1) {
      warning("Female bottleneck: reproductive females lost!")
      # Optionally, set ASR to NA or a capped value
    }
    
  }
}
N.median <- apply(apply(pig.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.10)
N.80pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2000000)) +
  geom_hline(yintercept=1000000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000,1000000,1500000, 2000000),
                     labels=c(0, 0.5, 1, 1.5, 2)) +
  scale_x_continuous(breaks=c(0,10,20,30,40), labels=c(0,5,10,15,20), name="Simulation Year") +
  theme_classic() + 
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, color="black")) 

results27 <- results

## Get realized survival rates

# Calculate median survival per stage (over all years and sims)
surv.50pct <- apply(realized_surv , 1, function(x) median(x, na.rm = TRUE))
surv.10pct <- apply(realized_surv, 1, quantile, probs = 0.10, na.rm=TRUE)
surv.90pct <- apply(realized_surv, 1, quantile, probs = 0.90, na.rm=TRUE)

r.surv <- data.frame(sex=c(rep("Female",3), rep("Male",3)), 
                     stage=rep(c("Piglet","Yearling","Adult"),2),
                     surv.50pct,surv.10pct,surv.90pct,
                     literature=c(A_base[2,1], A_base[3,2], A_base[3,3], A_base[5,4], A_base[6,5], A_base[6,6]))
library(dplyr)
r.surv <- r.surv |>
  pivot_longer(cols=c("surv.50pct", "literature"), names_to="source", values_to="survival")

# Square survival to get annual 
r.surv <- r.surv |>
  mutate(survival = survival^2,
         surv.10pct = surv.10pct^2,
         surv.90pct = surv.90pct^2)

# Set factor levels
r.surv$stage <- factor(r.surv$stage, levels=c("Piglet","Yearling","Adult"),
                       labels=c("Piglet","Yearling","Adult"))

r.surv$source <- factor(r.surv$source, levels=c("literature", "surv.50pct"), labels=c("Literature", "Realized"))

ggplot(r.surv) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=stage, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=stage, y=survival, shape=source, color=source), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  scale_shape_manual(values=c(17, 18, 16)) +
  guides(
    shape=guide_legend(position="inside", title="Source"),
    color=guide_legend(position="inside", title="Source")
  ) +
  facet_grid(.~sex) +
  ylab("Survival probability") + xlab("Stage") +
  theme_bw() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        legend.background = element_rect(fill=NA),
        strip.background=element_blank(),
        strip.text=element_text(hjust=0, size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))


r.surv27 <- r.surv


# Calculate median survival per stage (over all years and sims)
median_fec <- apply(realized_fecundity, 1, function(x) median(x, na.rm = TRUE))
fec.10pct <- apply(realized_fecundity, 1, quantile, probs = 0.10, na.rm=TRUE)
fec.90pct <- apply(realized_fecundity, 1, quantile, probs = 0.90, na.rm=TRUE)

fec27 <- data.frame(stage=rep(c("Piglet", "Yearling + Adult"), 2),
                    offspring_sex=rep(c("Female", "Male"), each=2), 
                    realized=median_fec,
                    f10pct=fec.10pct,
                    f90pct=fec.90pct,
                    literature=c(0.54, 2.24, 0.54, 2.24),  # annual 
                    x=c(0.9, 1.6, 1.1, 1.8)) |>
  mutate(realized=realized*2, f10pct=f10pct*2, f90pct=f90pct*2)

fec27 <- fec27 |>
  pivot_longer(cols=c("realized", "literature"), names_to="source", values_to="fec")
fec27$source <- factor(fec27$source, levels=c("literature", "realized"), labels=c("Literature", "Realized"))

# drop males
# fec27 <- fec27 |> 
#   select(-offspring_sex, -x) |>
#   distinct() 


### K at 8 pigs/ km2 - full state ####

## Stochastistic model
step <- 1:200

# Carrying capacity
K <- 367822

Kf <- K * 0.5

# Set up pig matrix
A_base <- matrix(0,6,6)

# Fecundity
# Fecundity <- c(0.54, 2.24, 2.24)

## Survival
# Females
A_base[2,1] <- sqrt(0.3)
A_base[3,2] <- sqrt(0.35)
A_base[3,3] <- sqrt(0.35)
# Males
A_base[5,4] <- sqrt(0.18)
A_base[6,5] <- sqrt(0.23)
A_base[6,6] <- sqrt(0.23)

## Fecundity
# Maximum fecundity
R0j <- 6.6*0.75 
R0a <- 10.8

FecundityM <- c(R0j, R0a)
FecundityF <- c(R0j, R0a)

# Add to matrix
A_base[1,2:3] <- FecundityF
A_base[4,2:3] <- FecundityM

pig.matrix <- A_base
# 1.3:1 males:females
pct_m <- 0.5
pct_f <- 0.5

# # Calculate for post-breeding census
FecundityM <- c(R0j*A_base[3,2], R0a*A_base[3,3])*pct_m
FecundityF <- c(R0j*A_base[3,2], R0a*A_base[3,3])*pct_f

pig.matrix[1,2:3] <- FecundityF
pig.matrix[4,2:3] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(pig.matrix)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(pig.matrix)$values[1])

# How many pig are harvested?
observed_harvest <- 150000 

### Optimize survival and fecundities
# set up bounds
pig_lower <- c(1, 1.2, 1.2, 1.2, 1.2, 1.1)
pig_upper <- c(1.5, 1.9, 1.9, 2, 2, 2.8)

# run optimizer
set.seed(1)
pig_opt <- optim(par=c(runif(1, 1.01, 1.45), runif(4, 1, 1.95), runif(1, 1.01, 2.75)), fn=objective_fn_pigs, method="L-BFGS-B", lower=pig_lower, upper=pig_upper, control = list(trace = 1))

# What are optimized parameter values 
pig_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Female survival
A_adj[2,1] <- A_base[2,1]*pig_opt$par[1]
A_adj[3,2] <- A_base[3,2]*pig_opt$par[1]#+0.05
A_adj[3,3] <- A_base[3,3]*pig_opt$par[1]#+0.05

# Male survival
A_adj[5,4] <- A_base[5,4]*pig_opt$par[2]#-0.08
A_adj[6,5] <- A_base[6,5]*pig_opt$par[2]#-0.08
A_adj[6,6] <- A_base[6,6]*pig_opt$par[2]#-0.08

## Fecundity
# Yearlings
A_adj[1,2] <- A_base[1,2]*pig_opt$par[3]
A_adj[4,2] <- A_base[4,2]*pig_opt$par[3]
# Adults
A_adj[1,3] <- A_base[1,3]*pig_opt$par[3]
A_adj[4,3] <- A_base[4,3]*pig_opt$par[3]

#### Run population model ####

set.seed(1)
##Stochastistic model
Year <- 1:40
Sims <- 1000

# Empty array to hold pop count
pig.array <- array(0,dim=c(6,length(Year),Sims))

# Calculate stable stage distribution
mnFs <- mean(c(A_adj[2,1], A_adj[3,2], A_adj[3,3]))
a2 <- A_adj
a2[1,1:3] <- a2[1,1:3]*0.5*mnFs
a2[4,1:3] <- a2[4,1:3]*0.5*mnFs
w <- Re(eigen(a2)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# female ssd
w_f <- Re(eigen(a2)$vectors[, 1])[1:3]
w_f <- w_f / sum(w_f)

# male ssd
w_m <- Re(eigen(a2)$vectors[, 1])[4:6]
w_m <- w_m / sum(w_m)

lambda <- Re(eigen(a2)$values[1])

# Set up initial popualtion size
N0 <- w * 1 * K # e.g., start at 50% of K

# Fill starting population
pig.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:3], lambda, Kf)
cat("Estimated theta:", theta, "\n")
theta <- 0.8

# Set up arrays to save realized demographic rates
realized_surv <- array(NA, dim=c(6, length(Year), Sims))
realized_surv[1:6,1,] <- c(A_adj[2,1], A_adj[3,2], A_adj[3,3], A_adj[5,4], A_adj[6,5], A_adj[6,6])

realized_fecundity <- array(NA, dim=c(4,  length(Year) - 1, Sims)) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,

for (i in 1:Sims) {
  
  A_s <- A_adj
  
  ## Stochasticity on Survival 
  # Female survival
  # Survive & go
  A_s[2,1] <- rtruncnorm(1, b=1, mean=A_adj[2,1], sd=0.01)
  A_s[3,2] <- rtruncnorm(1, b=1, mean=A_adj[3,2], sd=0.01)
  # Survive & stay
  A_s[3,3] <- rtruncnorm(1, a=0, b=1, mean=A_adj[3,3], sd=0.01)
  # Male survival
  A_s[5,4] <- rtruncnorm(1, b=1, mean=A_adj[5,4], sd=0.01)
  A_s[6,5] <- rtruncnorm(1, b=1, mean=A_adj[6,5], sd=0.01)
  # Survive & stay
  A_s[6,6] <- rtruncnorm(1, b=1, mean=A_adj[6,6], sd=0.01)
  
  # save new matrix
  A_dd <- A_s
  
  for (y in 2:length(Year)){
    
    ## Adjust fecundity by female survival
    R0j_s <- ifelse(is.na(A_s[1,2]*(realized_surv[3,y-1,i])), 0, A_s[1,2]*(realized_surv[3,y-1,i]))
    R0a_s <- ifelse(is.na(A_s[1,3]*(realized_surv[4,y-1,i])), 0, A_s[1,3]*(realized_surv[4,y-1,i]))
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(pig.array[1:3,y-1,i], na.rm=TRUE)
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)*theta)
    
    # Check adult age ratio
    asr <- sum(pig.array[5:6,y-1,i])/sum(pig.array[2:3,y-1,i], 1e-6)
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,2:3] <- c(R0j_s, R0a_s) * s_f * density_factor
    A_dd[4,2:3] <- c(R0j_s, R0a_s) * s_m * density_factor
    
    # Save fecundity
    realized_fecundity[1,y-1,i] <- A_dd[1,2] # Females per yearling female
    realized_fecundity[2,y-1,i] <- A_dd[1,3] # Females per adult female
    
    realized_fecundity[3,y-1,i] <- A_dd[4,2] # Males per yearling female
    realized_fecundity[4,y-1,i] <- A_dd[4,3] # Males per adult female
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% pig.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest pigs 
    if (y > 0 & min(N_next)>0 & sum(!is.na(N_next))==6) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 150000, 15000), digits=0)
      
      sow_h <- h * 0.45
      boar_h <- h - sow_h
      
      # Proportional harvest
      r_f <- (N_next[1:3,1] / sum(N_next[1:3,1]))*sow_h
      r_m <- (N_next[4:6,1] / sum(N_next[4:6,1]))*boar_h
      
      r <- c(r_f, r_m)
      
      # remove
      N_next <- N_next - r
    } 
    
    # Save abundance
    pig.array[,y,i] <- pmax(N_next, 0)
    
    # ---- Realized survival ----
    
    # Female stages
    realized_surv[1, y, i] <- pig.array[2, y, i] / pig.array[1, y-1, i] # Piglet → yearling (F)
    
    # Stages 2 and 3 are trickier
    incoming_2 <- pig.array[1, y-1, i] + pig.array[2, y-1, i]
    if (incoming_2 > 0) {
      realized_surv[2, y, i] <- pig.array[2, y, i] / incoming_2
    } else {
      realized_surv[2, y, i] <- NA
    }
    
    # Stage 3 receives from stage 2 → 3 and stage 3 → 3
    incoming_3 <- pig.array[2, y-1, i] + pig.array[3, y-1, i]
    if (incoming_3 > 0) {
      realized_surv[3, y, i] <- pig.array[3, y, i] / incoming_3
    } else {
      realized_surv[3, y, i] <- NA
    }
    
    # Male stages
    realized_surv[4, y, i] <- pig.array[5, y, i] / pig.array[4, y-1, i] # Piglet → yearling (M)
    
    # Male stages 5–6 (adults)
    # Stage 5 receives from stage 4 → 11 
    incoming_5 <- pig.array[4, y-1, i] + pig.array[5, y-1, i]
    if (incoming_5 > 0) {
      realized_surv[5, y, i] <- pig.array[5, y, i] / incoming_5
    } else {
      realized_surv[5, y, i] <- NA
    }
    
    # Stage 12 receives from stage 5 → 6 and stage 6 → 6
    incoming_6 <- pig.array[5, y-1, i] + pig.array[6, y-1, i]
    if (incoming_6 > 0) {
      realized_surv[6, y, i] <- pig.array[6, y, i] / incoming_6
    } else {
      realized_surv[6, y, i] <- NA
    }
    
    # Check buck:doe ratio
    af <- sum(pig.array[2:3,y,i], na.rm=TRUE)
    am <- sum(pig.array[5:6,y,i], na.rm=TRUE)
    
    if (sum(af) < 1) {
      warning("Female bottleneck: reproductive females lost!")
      # Optionally, set ASR to NA or a capped value
    }
    
  }
}
N.median <- apply(apply(pig.array,c(2,3),sum),1,median)
N.20pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.10)
N.80pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2000000)) +
  geom_hline(yintercept=1000000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000,1000000,1500000,2000000),
                     labels=c(0,0.5, 1,1.5, 2)) +
  scale_x_continuous(breaks=c(0,10,20,30,40), labels=c(0,5,10,15,20), name="Simulation Year") +
  theme_classic() + 
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, color="black")) 

results8 <- results

## Get realized survival rates
# subset realized survival to last 10 years
last_10_years <- (max(Year)-10):(max(Year)-1)
realized_surv_subset <- realized_surv[, last_10_years, ]

# Calculate median survival per stage (over all years and sims)
surv.50pct <- apply(realized_surv , 1, function(x) median(x, na.rm = TRUE))
surv.10pct <- apply(realized_surv, 1, quantile, probs = 0.10, na.rm=TRUE)
surv.90pct <- apply(realized_surv, 1, quantile, probs = 0.90, na.rm=TRUE)

r.surv <- data.frame(sex=c(rep("Female",3), rep("Male",3)), 
                     stage=rep(c("Piglet","Yearling","Adult"),2),
                     surv.50pct,surv.10pct,surv.90pct,
                     est=c(A_adj[2,1], A_adj[3,2], A_adj[3,3], A_adj[5,4], A_adj[6,5], A_adj[6,6]),
                     literature=c(A_base[2,1], A_base[3,2], A_base[3,3], A_base[5,4], A_base[6,5], A_base[6,6]))

r.surv <- r.surv |>
  pivot_longer(cols=c("surv.50pct", "est", "literature"), names_to="source", values_to="survival") |>
  mutate(survival = survival^2,
         surv.10pct = surv.10pct^2,
         surv.90pct = surv.90pct^2)


# Set factor levels
r.surv$stage <- factor(r.surv$stage, levels=c("Piglet","Yearling","Adult"),
                       labels=c("Piglet","Yearling","Adult"))

r.surv$source <- factor(r.surv$source, levels=c("literature", "est", "surv.50pct"), labels=c("Literature", "Optimized", "Implied"))

ggplot(r.surv) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=stage, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=stage, y=survival, shape=source, color=source), size=3) +
  scale_color_manual(values=c("#5ec962", "#440154", "#21918c")) +
  scale_shape_manual(values=c(17, 18, 16)) +
  guides(
    shape=guide_legend(position="inside", title="Source"),
    color=guide_legend(position="inside", title="Source")
  ) +
  facet_grid(.~sex) +
  ylab("Survival probability") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        legend.background = element_rect(fill=NA),
        strip.background=element_blank(),
        strip.text=element_text(hjust=0))

r.surv8 <- r.surv


# Calculate median survival per stage (over all years and sims)
median_fec <- apply(realized_fecundity, 1, function(x) median(x, na.rm = TRUE))
fec.10pct <- apply(realized_fecundity, 1, quantile, probs = 0.10, na.rm=TRUE)
fec.90pct <- apply(realized_fecundity, 1, quantile, probs = 0.90, na.rm=TRUE)

fec8 <- data.frame(stage=rep(c("Piglet", "Yearling + Adult"), 2),
                    offspring_sex=rep(c("Female", "Male"), each=2), 
                    realized=median_fec,
                    f10pct=fec.10pct,
                    f90pct=fec.90pct,
                    literature=c(0.54, 2.24, 0.54, 2.24),  # annual 
                    x=c(0.9, 1.6, 1.1, 1.8)) |>
  mutate(realized=realized*2, f10pct=f10pct*2, f90pct=f90pct*2)

fec8 <- fec8 |>
  pivot_longer(cols=c("realized", "literature"), names_to="source", values_to="fec")
fec8$source <- factor(fec8$source, levels=c("literature", "realized"), labels=c("Literature", "Realized"))

# drop males
fec8 <- fec8 |> 
  select(-offspring_sex, -x) |>
  distinct() 


#### Combine and plot for manuscript ####


results27$density <- "27 pigs km<sup>-2</sup>"
results8$density <- "8 pigs km<sup>-2</sup>"

results <- bind_rows(results27,results8)

results <- results |>
  as_tibble() |>
  rename(N.10pct=N.20pct, N.90pct=N.80pct)

write_csv(results, "results/pigs_population_projections.csv")


r.surv27$density <- "27 pigs km<sup>-2</sup>"
r.surv8$density <- "8 pigs km<sup>-2</sup>"


surv <- bind_rows(r.surv27, r.surv8)

surv <- surv |>
  filter(source != "Optimized")

write_csv(surv, "results/pigs_survival.csv")


write_csv(fec27, "results/pig_fecundity.csv")






res_plot <- ggplot(results) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1000000, color="black", linetype=3) +
  geom_hline(yintercept=2034396, linetype=2, color="#ed7953") +
  geom_hline(yintercept=603816, linetype=2, color="#0d0887") +
  geom_ribbon(aes(x=Year, ymin=N.10pct, ymax=N.90pct, fill=density, group=density), alpha=.2) +
  geom_line(aes(x=Year, y=N.median, color=density, group=density), alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000, 1000000, 1500000, 2000000, 2500000),
                     labels=c(0, 0.5, 1, 1.5, 2, 2.5)) +
  scale_x_continuous(breaks=c(0,10,20,30,40), labels=c(0,5,10,15,20)) +
  scale_color_manual(values=c("#ed7953","#0d0887")) +
  scale_fill_manual(values=c("#ed7953","#0d0887")) +
  guides(
    fill=guide_legend(position="inside", title="Density"),
    color=guide_legend(position="inside", title="Density")
  ) +
  xlab("Simulation Year") +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
        legend.position.inside=c(0.96,0.05),
        legend.justification=c(1,0),
        # legend.background=element_rect(fill=NA),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12)) 

## Survival plot
r.surv27 <- r.surv27 |>
  mutate(x=case_when(stage=="Piglet" ~ 0.5,
                     stage=="Yearling" ~ 0.6,
                     stage=="Adult" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01))


s_plot <- ggplot(r.surv27_all) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=x, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=x, y=survival, color=source, shape=sex), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  scale_shape_manual(values=c(17, 16)) +
  guides(
    color="none", #guide_legend(position="inside", title="Source"),
    shape=guide_legend(position="inside", title="Sex")
  ) +
  scale_x_continuous(breaks=c(0.5,0.6,0.7), labels=c("Piglet", "Yearling", "Adult"), expand = expansion(mult = 0.25, add = 0)) +
  ylab("Survival probability") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        legend.background = element_rect(fill=NA),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))

## Fecundity
f_plot <- ggplot(fec27) +
  coord_cartesian(ylim=c(0,2.5)) +
  geom_segment(aes(x=stage, y=f10pct, yend=f90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=stage, y=fec, group=source, color=source), size=3) +
  scale_color_manual(values=c("#5ec962", "#440154", "#21918c")) +
  guides(
    color=guide_legend(position="inside", title="Source")
  ) +
  ylab("Fecundity") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        legend.background = element_rect(fill=NA),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))

(res_plot / (f_plot | s_plot)) + plot_annotation(tag_levels = 'A') #+ plot_layout(heights=c(3,2))
