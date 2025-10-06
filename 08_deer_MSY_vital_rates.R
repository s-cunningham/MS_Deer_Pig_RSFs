library(ggplot2)
library(ggtext)
library(patchwork)
library(truncnorm)

source("00_functions.R")

theme_set(theme_bw())

set.seed(1)

#### Deer density 14.3 / km2 ####
# Years in simulation
Year <- 1:100

# Carrying capacity 14.3 deer/km2
K <- 1347232

# Set up deer matrix
A_base <- matrix(0,12,12)

## Survival
# Females
A_base[2,1] <- 0.52
A_base[3,2] <- 0.93
A_base[4,3] <- 0.84
A_base[5,4] <- 0.84
A_base[6,5] <- 0.84
A_base[6,6] <- 0.84
# Males
A_base[8,7] <- 0.52
A_base[9,8] <- 0.82
A_base[10,9] <- 0.63
A_base[11,10] <- 0.53
A_base[12,11] <- 0.44
A_base[12,12] <- 0.49

## Fecundity
# Maximum fecundity
R0a <- 1.8
R0y <- 1.4

FecundityM <- c(0, R0y, R0a, R0a, R0a, R0a)
FecundityF <- c(0, R0y, R0a, R0a, R0a, R0a)

# Add to matrix
A_base[1,1:6] <- FecundityF
A_base[7,1:6] <- FecundityM


deer.matrix <- A_base
# 1.3:1 males:females
pct_m <- 0.565
pct_f <- 0.435

# # Calculate for post-breeding census
FecundityM <- c(0, R0y*A_base[3,2], R0a*A_base[4,3], R0a*A_base[5,4], R0a*A_base[6,5], R0a*A_base[6,6])*pct_m
FecundityF <- c(0, R0y*A_base[3,2], R0a*A_base[4,3], R0a*A_base[5,4], R0a*A_base[6,5], R0a*A_base[6,6])*pct_f

deer.matrix[1,1:6] <- FecundityF
deer.matrix[7,1:6] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(deer.matrix)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(deer.matrix)$values[1])

Kf <- K * 0.60

# How many deer are harvested?
observed_harvest <- 280000 

### Optimize survival and fecundities
# set up bounds
deer_lower <- c(1, 1, 1, 1, 1, 1, 1.2, 1.4, 1.7, 1.7, 1, 1)
deer_upper <- c(1.25, 1.05 ,1.05, 1.05, 1.2, 1.15, 1.5, 1.8, 2.1, 1.95, 1.3, 1.3)

# run optimizer
set.seed(1)
deer_opt <- optim(par=runif(12, 1, 1.05), fn=objective_fn_deer, method="L-BFGS-B", lower=deer_lower, upper=deer_upper, control = list(trace = 1))

# What are optimized parameter values 
deer_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Fawn survival
A_adj[2,1] <- A_base[2,1]*deer_opt$par[1]
A_adj[8,7] <- A_base[8,7]*deer_opt$par[1]

# Female survival
A_adj[3,2] <- A_base[3,2]*deer_opt$par[2]
A_adj[4,3] <- A_base[4,3]*deer_opt$par[3]
A_adj[5,4] <- A_base[5,4]*deer_opt$par[4]
A_adj[6,5] <- A_base[6,5]*deer_opt$par[5]
A_adj[6,6] <- A_base[6,6]*deer_opt$par[6]

# Male survival
A_adj[9,8] <- A_base[9,8]*deer_opt$par[7]#0.95#
A_adj[10,9] <- A_base[10,9]*deer_opt$par[8]#0.93#
A_adj[11,10] <- A_base[11,10]*deer_opt$par[9] #0.94#
A_adj[12,11] <- A_base[12,11]*deer_opt$par[10] #0.89#
A_adj[12,12] <- A_base[12,12]*deer_opt$par[11] #0.82#

## Fecundity
# Yearlings
A_adj[1,2] <- A_base[1,2]*deer_opt$par[12]
A_adj[7,2] <- A_base[7,2]*deer_opt$par[12]
# Adults
A_adj[1,3:6] <- A_base[1,3:6]*deer_opt$par[12]
A_adj[7,3:6] <- A_base[7,3:6]*deer_opt$par[12]

# if >1, set to 1
A_adj[6, 5] <- 1
A_adj[10, 9] <- 1
A_adj[11, 10] <- 1



#### Run population model ####
set.seed(1)
##Stochastistic model
Year <- 1:20
Sims <- 1000

# Stochasticity on fecundity
f_pct <- 0.01

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
mnFs <- mean(c(A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6]))
a2 <- A_adj
a2[1,1:6] <- a2[1,1:6]*0.5*mnFs
a2[7,1:6] <- a2[7,1:6]*0.5*mnFs
w <- Re(eigen(a2)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

lambda <- Re(eigen(a2)$values[1])

# Set up initial popualtion size
N0 <- w * 1 * K # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")

# Set up arrays to save realized demographic rates
realized_surv <- array(NA, dim=c(12, length(Year), Sims))
realized_surv[1:12,1,] <- c(A_adj[2,1], A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6],
                                A_adj[8,7], A_adj[9,8], A_adj[10,9], A_adj[11,10], A_adj[12,11], A_adj[12,12])

realized_fecundity <- array(NA, dim=c(4,  length(Year) - 1, Sims)) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,
set.seed(1)
for (i in 1:Sims) {
  
  A_s <- A_adj
  
  ## Stochasticity on Survival 
  # Female survival
  for (s in 2:5) {
    # Survive & go
    A_s[s+1,s] <- rtruncnorm(1, b=1, mean=A_adj[s+1,s], sd=0.02)
  }
  # Survive & stay
  A_s[6,6] <- rtruncnorm(1, a=0, b=1, mean=A_adj[6,6], sd=0.02)
  # Male survival
  for (s in 8:11) {
    # Survive & go
    A_s[s+1,s] <- rtruncnorm(1, b=1, mean=A_adj[s+1,s], sd=0.02)
  }
  # Survive & stay
  A_s[12,12] <- rtruncnorm(1, b=1, mean=A_adj[12,12], sd=0.02)
  
  # save new matrix
  A_dd <- A_s
  
  for (y in 2:length(Year)){
    
    ## Adjust fecundity by female survival
    R0y_s <- A_s[1,2]*(realized_surv[2,y-1,i])
    R0a_s3 <- A_s[1,3]*(realized_surv[3,y-1,i])
    R0a_s4 <- A_s[1,3]*(realized_surv[4,y-1,i])
    R0a_s5 <- A_s[1,3]*(realized_surv[5,y-1,i])
    R0a_s6 <- A_s[1,3]*(realized_surv[6,y-1,i])
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i], na.rm=TRUE)
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)

    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)

    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,1:6] <- c(0, R0y_s, R0a_s3, R0a_s4, R0a_s5, R0a_s6) * s_f * density_factor
    A_dd[7,1:6] <- c(0, R0y_s, R0a_s3, R0a_s4, R0a_s5, R0a_s6) * s_m * density_factor
    
    # Save fecundity
    realized_fecundity[1,y-1,i] <- A_dd[1,2] # Females per yearling female
    realized_fecundity[2,y-1,i] <- A_dd[1,3] # Females per adult female
    realized_fecundity[3,y-1,i] <- A_dd[7,2] # Males per yearling female 
    realized_fecundity[4,y-1,i] <- A_dd[7,3] # Males per adult female
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer 
    if (y > 0) {
      # Randomly select a random number to harvest
      # h <- round(runif(1, 190000, 280000), digits=0)
      h <- round(rnorm(1, 240000, 15000),digits=0)
      # h=200000
      
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
    
    # ---- Realized survival ----
    
    # Female stages
    realized_surv[1, y, i] <- deer.array[2, y, i] / deer.array[1, y-1, i] # Fawn → yearling (F)
    realized_surv[2, y, i] <- deer.array[3, y, i] / deer.array[2, y-1, i]  # yearling → 2-year-old  
    realized_surv[3, y, i] <- deer.array[4, y, i] / deer.array[3, y-1, i]  # 2-year-old → 3-year-old
    realized_surv[4, y, i] <- deer.array[5, y, i] / deer.array[4, y-1, i]  # 3-year-old → 4-year-old
    
    # Stage 5: comes only from stage 4
    incoming_5 <- deer.array[4, y-1, i]
    realized_surv[5, y, i] <- ifelse(incoming_5 > 0,
                                      deer.array[5, y, i] / incoming_5,
                                      NA)
    
    # Stage 6: comes from stage 5 -> 6 and stage 6 -> 6
    incoming_6 <- deer.array[5, y-1, i] + deer.array[6, y-1, i]
    realized_surv[6, y, i] <- ifelse(incoming_6 > 0,
                                      deer.array[6, y, i] / incoming_6,
                                      NA)
    
    # Male stages
    realized_surv[7, y, i] <- deer.array[8, y, i] / deer.array[7, y-1, i] # Fawn → yearling (M)
    realized_surv[8, y, i] <- deer.array[9, y, i] / deer.array[8, y-1, i]   # 8 → 9
    realized_surv[9, y, i] <- deer.array[10, y, i] / deer.array[9, y-1, i] # 9 → 10
    realized_surv[10, y, i] <- deer.array[11, y, i] / deer.array[10, y-1, i] # 10 → 11
    
    # Stage 11: comes only from stage 10
    incoming_11 <- deer.array[10, y-1, i]
    realized_surv[11, y, i] <- ifelse(incoming_11 > 0,
                                      deer.array[11, y, i] / incoming_11,
                                      NA)
    
    # Stage 12: comes from stage 11 -> 12 and stage 12 -> 12
    incoming_12 <- deer.array[11, y-1, i] + deer.array[12, y-1, i]
    realized_surv[12, y, i] <- ifelse(incoming_12 > 0,
                                      deer.array[12, y, i] / incoming_12,
                                      NA)
    
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
N.10pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.10)
N.90pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.10pct,N.90pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  # geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.10pct, ymax=N.90pct), alpha=.2,fill="#21918c") +
  geom_line(aes(x=Year, y=N.median),colour="#21918c",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000, 1000000, 1500000, 2000000, 2500000),
                     labels=c(0, 0.5, 1, 1.5, 2, 2.5)) +
  theme_bw() +
  xlab("Simulation Year") +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
        legend.position.inside=c(0,0),
        legend.justification=c(0,0),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) 

results14 <- results

## Get realized survival rates

# Calculate median survival per stage (over all years and sims)
surv.50pct <- apply(realized_surv, 1, function(x) median(x, na.rm = TRUE))
surv.10pct <- apply(realized_surv, 1, quantile, probs = 0.10, na.rm=TRUE)
surv.90pct <- apply(realized_surv, 1, quantile, probs = 0.90, na.rm=TRUE)

r.surv <- data.frame(sex=c(rep("Female",6), rep("Male",6)), 
                     stage=rep(c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"),2),
                     surv.50pct,surv.10pct,surv.90pct,
                     est=c(A_adj[2,1], A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6],
                           A_adj[8,7], A_adj[9,8], A_adj[10,9], A_adj[11,10], A_adj[12,11], A_adj[12,12]),
                     literature=c(A_base[2,1], A_base[3,2], A_base[4,3], A_base[5,4], A_base[6,5], A_base[6,6],
                                  A_base[8,7], A_base[9,8], A_base[10,9], A_base[11,10], A_base[12,11], A_base[12,12]))
library(dplyr)
r.surv <- r.surv |>
  pivot_longer(cols=c("surv.50pct", "literature"), names_to="source", values_to="survival")


# Set factor levels
r.surv$stage <- factor(r.surv$stage, levels=c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"),
                       labels=c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"))

r.surv$source <- factor(r.surv$source, levels=c("literature", "surv.50pct"), labels=c("Literature", "Implied"))

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
  theme(#panel.border=element_rect(fill=NA, color="black"),
        legend.position.inside = c(0,0),
        legend.justification=c(0,0),
        legend.background = element_rect(fill=NA),
        strip.background=element_blank(),
        strip.text=element_text(hjust=0, size=14),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(angle=25, vjust=0.8, hjust=0.8, size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))

surv14 <- r.surv

## Get realized fecundities

# Calculate median survival per stage (over all years and sims)
median_fec <- apply(realized_fecundity, 1, function(x) median(x, na.rm = TRUE))
fec.10pct <- apply(realized_fecundity, 1, quantile, probs = 0.10, na.rm=TRUE)
fec.90pct <- apply(realized_fecundity, 1, quantile, probs = 0.90, na.rm=TRUE)


fec14 <- data.frame(stage=rep(c("Yearling", "Adult"), 2),
                    offspring_sex=rep(c("Female", "Male"), each=2), 
                    realized=median_fec,
                    f10pct=fec.10pct,
                    f90pct=fec.90pct,
                    literature=c(0.57, 0.66, 0.74, 0.85),
                    x=c(0.9, 1.6, 1.1, 1.8))

fec14 <- fec14 |>
  pivot_longer(cols=c("realized", "literature"), names_to="source", values_to="fec")
fec14$source <- factor(fec14$source, levels=c("literature", "realized"), labels=c("Literature", "Implied"))

# Summarize sex ratio
ASR.median <- apply(asr_mat, 1, median)
ASR.20 <- apply(asr_mat, 1, quantile, probs = 0.20)
ASR.80 <- apply(asr_mat, 1, quantile, probs = 0.80)

ASR <- data.frame(Year, ASR.median, ASR.20, ASR.80)
ASR <- ASR[-nrow(ASR),]

# Plot sex ratio
ggplot(ASR) +
  geom_hline(yintercept=1, linetype=2) +
  geom_hline(yintercept=0.4, linetype=2, color="red") +
  coord_cartesian(ylim=c(0,1.5)) +
  geom_ribbon(aes(x=Year, ymin=ASR.20, ymax=ASR.80), alpha=0.5) +
  geom_line(aes(x=Year, y=ASR.median)) +
  ylab("ASR (males/females)") +
  theme_classic()

## Get Lambda over simulation
gm_mean(lambda_calc(results14$N.median))

#### Deer density 22 / km2 ####

set.seed(1)

# Years in simulation
Year <- 1:100

# Carrying capacity
K <- 2090532

# Set up deer matrix
A_base <- matrix(0,12,12)

## Survival
# Females
A_base[2,1] <- 0.52
A_base[3,2] <- 0.93
A_base[4,3] <- 0.84
A_base[5,4] <- 0.84
A_base[6,5] <- 0.84
A_base[6,6] <- 0.84
# Males
A_base[8,7] <- 0.52
A_base[9,8] <- 0.82
A_base[10,9] <- 0.63
A_base[11,10] <- 0.53
A_base[12,11] <- 0.44
A_base[12,12] <- 0.49

## Fecundity
# Maximum fecundity
R0a <- 1.8
R0y <- 1.4

FecundityM <- c(0, R0y, R0a, R0a, R0a, R0a)
FecundityF <- c(0, R0y, R0a, R0a, R0a, R0a)

# Add to matrix
A_base[1,1:6] <- FecundityF
A_base[7,1:6] <- FecundityM


deer.matrix <- A_base
# 1.3:1 males:females
pct_m <- 0.565
pct_f <- 0.435

# # Calculate for post-breeding census
FecundityM <- c(0, R0y*A_base[3,2], R0a*A_base[4,3], R0a*A_base[5,4], R0a*A_base[6,5], R0a*A_base[6,6])*pct_m
FecundityF <- c(0, R0y*A_base[3,2], R0a*A_base[4,3], R0a*A_base[5,4], R0a*A_base[6,5], R0a*A_base[6,6])*pct_f

deer.matrix[1,1:6] <- FecundityF
deer.matrix[7,1:6] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(deer.matrix)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(deer.matrix)$values[1])

Kf <- K * 0.60

# How many deer are harvested?
observed_harvest <- 290000 

### Optimize survival and fecundities
# set up bounds
deer_lower <- c(1, 1, 1, 1, 1, 1, 1.2, 1.4, 1.7, 1.7, 1, 1)
deer_upper <- c(1.25, 1.05 ,1.05, 1.05, 1.2, 1.15, 1.5, 1.8, 2.1, 1.95, 1.3, 1.3)

# run optimizer
set.seed(1)
deer_opt <- optim(par=runif(12, 1, 1.05), fn=objective_fn_deer, method="L-BFGS-B", lower=deer_lower, upper=deer_upper, control = list(trace = 1))

# What are optimized parameter values 
deer_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Fawn survival
A_adj[2,1] <- A_base[2,1]*deer_opt$par[1]
A_adj[8,7] <- A_base[8,7]*deer_opt$par[1]

# Female survival
A_adj[3,2] <- A_base[3,2]*deer_opt$par[2]
A_adj[4,3] <- A_base[4,3]*deer_opt$par[3]
A_adj[5,4] <- A_base[5,4]*deer_opt$par[4]
A_adj[6,5] <- A_base[6,5]*deer_opt$par[5]
A_adj[6,6] <- A_base[6,6]*deer_opt$par[6]

# Male survival
A_adj[9,8] <- A_base[9,8]*deer_opt$par[7]
A_adj[10,9] <- A_base[10,9]*deer_opt$par[8]
A_adj[11,10] <- A_base[11,10]*deer_opt$par[9]
A_adj[12,11] <- A_base[12,11]*deer_opt$par[10]
A_adj[12,12] <- A_base[12,12]*deer_opt$par[11]

## Fecundity
# Yearlings
A_adj[1,2] <- A_base[1,2]*deer_opt$par[12]
A_adj[7,2] <- A_base[7,2]*deer_opt$par[12]
# Adults
A_adj[1,3:6] <- A_base[1,3:6]*deer_opt$par[12]
A_adj[7,3:6] <- A_base[7,3:6]*deer_opt$par[12]

# if >1, set to 1
A_adj[6, 5] <- 1

#### Run population model ####

set.seed(1)
##Stochastistic model
Year <- 1:20
Sims <- 1000

# Stochasticity on fecundity
f_pct <- 0.01

# Empty array to hold pop count
deer.array <- array(0,dim=c(12,length(Year),Sims))

# Calculate stable stage distribution
mnFs <- mean(c(A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6]))
a2 <- A_adj
a2[1,1:6] <- a2[1,1:6]*0.5*mnFs
a2[7,1:6] <- a2[7,1:6]*0.5*mnFs
w <- Re(eigen(a2)$vectors[, 1])
w <- w / sum(w) # normals to equal 1

# female ssd
w_f <- Re(eigen(a2)$vectors[, 1])[1:6]
w_f <- w_f / sum(w_f)

# male ssd
w_m <- Re(eigen(a2)$vectors[, 1])[7:12]
w_m <- w_m / sum(w_m)

lambda <- Re(eigen(a2)$values[1])

# Set up initial popualtion size
N0 <- w * 1 * K # e.g., start at 50% of K

# Fill starting population
deer.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
theta <- estimate_theta_opt(N0[1:6], lambda, Kf)
cat("Estimated theta:", theta, "\n")

# Set up arrays to save realized demographic rates
realized_surv <- array(NA, dim=c(12, length(Year), Sims))
realized_surv[1:12,1,] <- c(A_adj[2,1], A_adj[3,2], A_adj[4,3], A_adj[5,4], A_adj[6,5], A_adj[6,6],
                                A_adj[8,7], A_adj[9,8], A_adj[10,9], A_adj[11,10], A_adj[12,11], A_adj[12,12])

realized_fecundity <- array(NA, dim=c(4,  length(Year) - 1, Sims)) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,

for (i in 1:Sims) {
  
  A_s <- A_adj
  
  ## Stochasticity on Survival 
  # Female survival
  for (s in 2:5) {
    # Survive & go
    A_s[s+1,s] <- rtruncnorm(1, b=1, mean=A_adj[s+1,s], sd=0.02)
  }
  # Survive & stay
  A_s[6,6] <- rtruncnorm(1, a=0, b=1, mean=A_adj[6,6], sd=0.02)
  # Male survival
  for (s in 8:11) {
    # Survive & go
    A_s[s+1,s] <- rtruncnorm(1, b=1, mean=A_adj[s+1,s], sd=0.02)
  }
  # Survive & stay
  A_s[12,12] <- rtruncnorm(1, b=1, mean=A_adj[12,12], sd=0.02)
  
  # save new matrix
  A_dd <- A_s
  
  for (y in 2:length(Year)){
    
    ## Adjust fecundity by female survival
    R0y_s <- A_s[1,2]*(realized_surv[2,y-1,i])
    R0a_s3 <- A_s[1,3]*(realized_surv[3,y-1,i])
    R0a_s4 <- A_s[1,3]*(realized_surv[4,y-1,i])
    R0a_s5 <- A_s[1,3]*(realized_surv[5,y-1,i])
    R0a_s6 <- A_s[1,3]*(realized_surv[6,y-1,i])
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(deer.array[1:6,y-1,i], na.rm=TRUE)
    
    # Calcuate density factor
    density_factor <- 1 / (1 + (Nf_t / Kf)^theta)
    
    # Check adult age ratio
    asr <- sum(deer.array[8:12,y-1,i])/sum(deer.array[2:6,y-1,i], 1e-6)
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,1:6] <- c(0, R0y_s, R0a_s3, R0a_s4, R0a_s5, R0a_s6) * s_f * density_factor
    A_dd[7,1:6] <- c(0, R0y_s, R0a_s3, R0a_s4, R0a_s5, R0a_s6) * s_m * density_factor
    
    # Save fecundity
    realized_fecundity[1,y-1,i] <- A_dd[1,2] # Females per yearling female
    realized_fecundity[2,y-1,i] <- A_dd[1,3] # Females per adult female
    realized_fecundity[3,y-1,i] <- A_dd[7,2] # Males per yearling female 
    realized_fecundity[4,y-1,i] <- A_dd[7,3] # Males per adult female
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% deer.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest deer 
    if (y > 0) {
      # Randomly select a random number to harvest
      # h <- round(runif(1, 190000, 280000), digits=0)
      h <- round(rnorm(1, 240000, 15000),digits=0)
      # h=200000
      
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
    
    # ---- Realized survival ----
    
    # Female stages
    realized_surv[1, y, i] <- deer.array[2, y, i] / deer.array[1, y-1, i] # Fawn → yearling (F)
    realized_surv[2, y, i] <- deer.array[3, y, i] / deer.array[2, y-1, i]  # yearling → 2-year-old  
    realized_surv[3, y, i] <- deer.array[4, y, i] / deer.array[3, y-1, i]  # 2-year-old → 3-year-old
    realized_surv[4, y, i] <- deer.array[5, y, i] / deer.array[4, y-1, i]  # 3-year-old → 4-year-old
    
    # Stage 5: comes only from stage 4
    incoming_5 <- deer.array[4, y-1, i]
    realized_surv[5, y, i] <- ifelse(incoming_5 > 0,
                                     deer.array[5, y, i] / incoming_5,
                                     NA)
    
    # Stage 6: comes from stage 5 -> 6 and stage 6 -> 6
    incoming_6 <- deer.array[5, y-1, i] + deer.array[6, y-1, i]
    realized_surv[6, y, i] <- ifelse(incoming_6 > 0,
                                     deer.array[6, y, i] / incoming_6,
                                     NA)
    
    # Male stages
    realized_surv[7, y, i] <- deer.array[8, y, i] / deer.array[7, y-1, i] # Fawn → yearling (M)
    realized_surv[8, y, i] <- deer.array[9, y, i] / deer.array[8, y-1, i]   # 8 → 9
    realized_surv[9, y, i] <- deer.array[10, y, i] / deer.array[9, y-1, i] # 9 → 10
    realized_surv[10, y, i] <- deer.array[11, y, i] / deer.array[10, y-1, i] # 10 → 11
    
    # Stage 11: comes only from stage 10
    incoming_11 <- deer.array[10, y-1, i]
    realized_surv[11, y, i] <- ifelse(incoming_11 > 0,
                                      deer.array[11, y, i] / incoming_11,
                                      NA)
    
    # Stage 12: comes from stage 11 -> 12 and stage 12 -> 12
    incoming_12 <- deer.array[11, y-1, i] + deer.array[12, y-1, i]
    realized_surv[12, y, i] <- ifelse(incoming_12 > 0,
                                      deer.array[12, y, i] / incoming_12,
                                      NA)
    
    # Check buck:doe ratio
    af <- sum(deer.array[2:6,y,i])
    am <- sum(deer.array[8:12,y,i])
    
    if (sum(af) < 1) {
      warning("Female bottleneck: reproductive females lost!")
    }
    
  }
}
N.median <- apply(apply(deer.array,c(2,3),sum),1,median)
N.10pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.10)
N.90pct <- apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.90)

results <- data.frame(Year,N.median,N.10pct,N.90pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  # geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.10pct, ymax=N.90pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

results22 <- results


## Get Lambda over simulation
gm_mean(lambda_calc(results22$N.median))

## Get realized survival rates
# subset realized survival to last 10 years
last_10_years <- (max(Year)-10):(max(Year)-1)
realized_surv_subset <- realized_surv[, last_10_years, ]

# Calculate median survival per stage (over all years and sims)
surv.50pct <- apply(realized_surv, 1, function(x) median(x, na.rm = TRUE))
surv.10pct <- apply(realized_surv, 1, quantile, probs = 0.10, na.rm=TRUE)
surv.90pct <- apply(realized_surv, 1, quantile, probs = 0.90, na.rm=TRUE)

r.surv <- data.frame(sex=c(rep("Female",6), rep("Male",6)), 
                     stage=rep(c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"),2),
                     surv.50pct,surv.10pct,surv.90pct,
                     literature=c(A_base[2,1], A_base[3,2], A_base[4,3], A_base[5,4], A_base[6,5], A_base[6,6],
                                  A_base[8,7], A_base[9,8], A_base[10,9], A_base[11,10], A_base[12,11], A_base[12,12]))
library(dplyr)
r.surv <- r.surv |>
  pivot_longer(cols=c("surv.50pct", "literature"), names_to="source", values_to="survival")


# Set factor levels
r.surv$stage <- factor(r.surv$stage, levels=c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"),
                       labels=c("Fawn","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"))

r.surv$source <- factor(r.surv$source, levels=c("literature", "surv.50pct"), labels=c("Literature",  "Implied"))

ggplot(r.surv) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=stage, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=stage, y=survival, shape=source, color=source), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  scale_shape_manual(values=c(18, 16)) +
  guides(
    shape=guide_legend(position="inside", title="Source"),
    color=guide_legend(position="inside", title="Source")
  ) +
  facet_grid(.~sex) +
  ylab("Survival probability") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,0),
        legend.justification=c(0,0),
        legend.background = element_rect(fill=NA),
        strip.background=element_blank(),
        strip.text=element_text(hjust=0))

surv22 <- r.surv

surv22 <- surv22 |>
  mutate(x=case_when(stage=="Fawn" ~ 0.2,
                     stage=="Yearling" ~ 0.3,
                     stage=="2-year-old" ~ 0.4,
                     stage=="3-year-old" ~ 0.5,
                     stage=="4-year-old" ~ 0.6,
                     stage=="5+ years-old" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01))

s22_plot <- ggplot(surv22) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=x, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=x, y=survival, color=source, shape=sex), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  scale_shape_manual(values=c(17, 16)) +
  guides(
    color=guide_legend(position="inside", title="Source"),
    shape=guide_legend(position="inside", title="Sex")
  ) +
  scale_x_continuous(breaks=seq(0.2,0.7, by=0.1), 
                     labels=c("Fawn", "Yearling", "2 years", "3 years", "4 years", "5+ years"), 
                     expand = expansion(mult = 0.15, add = 0)) +
  ylab("Survival probability") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,0),
        legend.justification=c(0,0),
        legend.background = element_rect(fill=NA),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))

## Get realized fecundities

# Calculate median survival per stage (over all years and sims)
median_fec <- apply(realized_fecundity, 1, function(x) median(x, na.rm = TRUE))
fec.10pct <- apply(realized_fecundity, 1, quantile, probs = 0.10, na.rm=TRUE)
fec.90pct <- apply(realized_fecundity, 1, quantile, probs = 0.90, na.rm=TRUE)


fec22 <- data.frame(stage=rep(c("Yearling", "Adult"), 2),
                    offspring_sex=rep(c("Female", "Male"), each=2), 
                    realized=median_fec,
                    f10pct=fec.10pct,
                    f90pct=fec.90pct,
                    literature=c(0.57, 0.66, 0.74, 0.85),
                    x=c(0.9, 1.6, 1.1, 1.8))

fec22 <- fec22 |>
  pivot_longer(cols=c("realized", "literature"), names_to="source", values_to="fec")
fec22$source <- factor(fec22$source, levels=c("literature", "realized"), labels=c("Literature", "Implied"))


f22_plot <- ggplot(fec22) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,1.4)) +
  geom_segment(aes(x=x, y=f10pct, yend=f90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=x, y=fec, group=source, color=source, shape=offspring_sex), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  scale_shape_manual(values=c(17, 16)) +
  scale_x_continuous(breaks=c(1,1.7), 
                     labels=c("Yearling", "Adult"), 
                     expand = expansion(mult = 0.15, add = 0)) +
  guides(
    color=guide_legend(position="inside", title="Source"),
    shape=guide_legend(position="inside", title="Offspring Sex")
  ) +
  ylab("Fecundity") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        # legend.direction="horizontal",
        legend.background = element_rect(fill=NA),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11))

library(patchwork)

s22_plot + f22_plot + plot_layout(widths = c(1.8, 1))


#### Save results from both densities ####

# Combine population projections
results14$density <- "14.3 deer km<sup>-2</sup>"
results22$density <- "22 deer km<sup>-2</sup>"
results <- bind_rows(results14, results22)

write_csv(results, "results/deer_population_projections.csv")

# Combine survival
surv14$density <- "14.3 deer km<sup>-2</sup>"
surv22$density <- "22 deer km<sup>-2</sup>"

survival <- bind_rows(surv14, surv22)

write_csv(survival, "results/deer_survival.csv")

# Combine fecundity
fec14$density <- "14.3 deer km<sup>-2</sup>"
fec22$density <- "22 deer km<sup>-2</sup>"

fecundity <- bind_rows(fec14, fec22)

write_csv(fecundity, "results/deer_fecundity.csv")













res_plot <- ggplot(results14) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1610000, color="black", linetype=3) +
  geom_hline(yintercept=1512638, linetype=2, color="#ed7953") +
  geom_hline(yintercept=2347197, linetype=2, color="#0d0887") +
  geom_ribbon(aes(x=Year, ymin=N.10pct, ymax=N.90pct, fill=density, group=density), alpha=.2) +
  geom_line(aes(x=Year, y=N.median, color=density, group=density), alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000, 1000000, 1500000, 2000000, 2500000),
                     labels=c(0, 0.5, 1, 1.5, 2, 2.5))+
  scale_color_manual(values=c("#ed7953","#0d0887")) +
  scale_fill_manual(values=c("#ed7953","#0d0887")) +
  guides(
    fill=guide_legend(position="inside", title="Density"),
    color=guide_legend(position="inside", title="Density")
  ) +
  xlab("Simulation Year") +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    legend.position.inside=c(0,0),
    legend.justification=c(0,0),
    legend.background=element_rect(fill=NA),
    legend.text=element_markdown(size=11),
    legend.title=element_text(size=12),
    axis.text=element_text(size=11),
    axis.title=element_text(size=12)) 

s_text <- data.frame(x=rep("Fawn",2), y=1, sex=c("Female", "Male"), label=c("Female", "Male"))

s_plot <- ggplot(r.surv14) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=stage, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=stage, y=survival, color=source), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  # scale_shape_manual(values=c(17, 18, 16)) +
  guides(
    # shape=guide_legend(position="inside", title="Source"),
    color=guide_legend(position="inside", title="Source")
  ) +
  geom_text(data=s_text, aes(x=x, y=1, label=label)) +
  facet_grid(sex~.) +
  ylab("Survival probability") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,0),
        legend.justification=c(0,0),
        legend.background = element_rect(fill=NA),
        strip.background=element_blank(),
        axis.text.x=element_text(angle=25, vjust=0.8, hjust=0.8, size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=12),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12),
        strip.text=element_blank())

## Survival plot
r.surv27_all <- r.surv27_all |>
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





 f_plot <- ggplot(fec14) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,1.4)) +
  geom_segment(aes(x=x, y=f10pct, yend=f90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=x, y=fec, group=source, color=source, shape=offspring_sex), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  scale_shape_manual(values=c(17, 16)) +
  scale_x_continuous(breaks=c(1,1.7), labels=c("Yearling", "Adult")) +
  guides(
    color=guide_legend(position="inside", title="Source"),
    shape=guide_legend(position="inside", title="Offspring Sex")
  ) +
  ylab("Fecundity") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        # legend.direction="horizontal",
        legend.background = element_rect(fill=NA),
        legend.text=element_markdown(size=12),
        legend.title=element_text(size=14),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12))

(free(res_plot / f_plot) | s_plot) + plot_annotation(tag_levels = 'A') + plot_layout(widths=c(1.2,1))
