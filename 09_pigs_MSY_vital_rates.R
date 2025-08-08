library(ggplot2)
library(ggtext)
library(patchwork)
library(truncnorm)

source("00_functions.R")

theme_set(theme_bw())

### K at 27 pigs/ km2 ####

## Stochastistic model
Year <- 1:50

# Carrying capacity
K <- 2034396 

Kf <- K * 0.5

# Set up pig matrix
A_base <- matrix(0,6,6)

# Fecundity
# Fecundity <- c(0.54, 2.24, 2.24)

## Survival
# Females
A_base[2,1] <- 0.3
A_base[3,2] <- 0.35
A_base[3,3] <- 0.35
# Males
A_base[5,4] <- 0.18
A_base[6,5] <- 0.23
A_base[6,6] <- 0.23

## Fecundity
# Maximum fecundity
R0p <- 4.8*0.75 
R0y <- 6.4*2
R0a <- 6.4*2

FecundityM <- c(R0p, R0y, R0a)
FecundityF <- c(R0p, R0y, R0a)

# Add to matrix
A_base[1,1:3] <- FecundityF
A_base[4,1:3] <- FecundityM

pig.matrix <- A_base
# 1.3:1 males:females
pct_m <- 0.5
pct_f <- 0.5

# # Calculate for post-breeding census
FecundityM <- c(R0p*A_base[2,1], R0y*A_base[3,2], R0a*A_base[3,3])*pct_m
FecundityF <- c(R0p*A_base[5,4], R0y*A_base[3,2], R0a*A_base[3,3])*pct_f

pig.matrix[1,1:3] <- FecundityF
pig.matrix[4,1:3] <- FecundityM

# Calculate stable stage distribution
w <- Re(eigen(pig.matrix)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(pig.matrix)$values[1])

# How many pig are harvested?
observed_harvest <- 350000 

### Optimize survival and fecundities
# set up bounds
pig_lower <- c(0.75,0.75,0.75,0.75,.95,0.95)
pig_upper <- c(3,2.8,3,3,3.25,3.25)

# run optimizer
set.seed(1)
pig_opt <- optim(par=runif(8, 1.01, 2.4), fn=objective_fn_pigs, method="L-BFGS-B", lower=pig_lower, upper=pig_upper, control = list(trace = 1))

# What are optimized parameter values 
pig_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Female survival
A_adj[2,1] <- A_base[2,1]*pig_opt$par[1]
A_adj[3,2] <- A_base[3,2]*pig_opt$par[2]
A_adj[3,3] <- A_base[3,3]*pig_opt$par[3]

# Male survival
A_adj[5,4] <- A_base[5,4]*pig_opt$par[4]
A_adj[6,5] <- A_base[6,5]*pig_opt$par[5]
A_adj[6,6] <- A_base[6,6]*pig_opt$par[6]

## Fecundity
# Yearlings
A_adj[1,1] <- A_base[1,1]*pig_opt$par[7]
A_adj[4,1] <- A_base[4,1]*pig_opt$par[7]
# Adults
A_adj[1,2:3] <- A_base[1,2:3]*pig_opt$par[8]
A_adj[4,2:3] <- A_base[4,2:3]*pig_opt$par[8]



#### Run population model ####


set.seed(1)
##Stochastistic model
Year <- 1:20
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

c5_f <- (a2[3,2] * w_f[2]) / ((a2[3,2] * w_f[2]) + (a2[3,3] * w_f[3]))
c5_m <- (a2[6,5] * w_m[2]) / ((a2[6,5] * w_m[2]) + (a2[6,6] * w_m[3]))

lambda <- Re(eigen(a2)$values[1])

# Set up initial popualtion size
N0 <- w * 0.75 * K # e.g., start at 50% of K

# Fill starting population
pig.array[,1,] <- N0

# Create matrix to hold ASR
asr_mat <- matrix(0, nrow=length(Year), ncol=Sims)

# Calibrate density dependence parameter
a <- calibrate_a_opt(a2, N0, K)
cat("Calibrated a =", a, "\n")
a <- 3

# Set up arrays to save realized demographic rates
realized_survival <- array(NA, dim=c(6, length(Year), Sims))
realized_survival[1:6,1,] <- c(A_adj[2,1], A_adj[3,2], A_adj[3,3], A_adj[5,4], A_adj[6,5], A_adj[6,6])

realized_fecundity <- array(NA, dim=c(6,  length(Year) - 1, Sims)) # row 1: yearling F, row 2: adult F, row 3: yearling M, row 4: adult M,

for (i in 1:Sims) {
  
  A_s <- A_adj
  
  ## Stochasticity on Survival 
  # Female survival
  # Survive & go
  A_s[2,1] <- rtruncnorm(1, b=1, mean=A_adj[2,1], sd=0.02)
  A_s[3,2] <- rtruncnorm(1, b=1, mean=A_adj[3,2], sd=0.02)
  # Survive & stay
  A_s[3,3] <- rtruncnorm(1, a=0, b=1, mean=A_adj[6,6], sd=0.02)
  # Male survival
  A_s[5,4] <- rtruncnorm(1, b=1, mean=A_adj[5,4], sd=0.02)
  A_s[6,5] <- rtruncnorm(1, b=1, mean=A_adj[6,5], sd=0.02)
  # Survive & stay
  A_s[6,6] <- rtruncnorm(1, b=1, mean=A_adj[6,6], sd=0.02)
  
  # save new matrix
  A_dd <- A_s
  
  for (y in 2:length(Year)){
    
    ## Adjust fecundity by female survival
    R0p_s <- A_s[1,1]*(realized_survival[2,y-1,i])
    R0y_s <- A_s[1,2]*(realized_survival[3,y-1,i])
    R0a_s <- A_s[1,3]*(realized_survival[4,y-1,i])
    
    # Density dependent adjustment in fecundity
    # What was abundance at time step t-1
    Nf_t <- sum(pig.array[1:3,y-1,i], na.rm=TRUE)
    
    # Calcuate density factor
    # density_factor <- 1 / (1 + a * (Nf_t / Kf))
    density_factor <- 1 / (1 + (Nf_t / Kf)*theta)
    
    # Check adult age ratio
    asr <- sum(pig.array[5:6,y-1,i])/sum(pig.array[2:3,y-1,i], 1e-6)
    
    # Save sex ratio
    asr_mat[y-1,i] <- asr
    
    # Calculate % males & females
    s_m <- pmin(pmax(0.5 + 0.1 * (1 - asr), 0.4), 0.6)
    s_f <- 1 - s_m
    
    # Adjust sex ratio at birth
    A_dd[1,1:3] <- c(R0p_s, R0y_s, R0a_s) * s_f * density_factor
    A_dd[4,1:3] <- c(R0p_s, R0y_s, R0a_s) * s_m * density_factor
    
    # Save fecundity
    realized_fecundity[1,y-1,i] <- A_dd[1,1] # Females per yearling female
    realized_fecundity[2,y-1,i] <- A_dd[1,2] # Females per yearlingfemale
    realized_fecundity[3,y-1,i] <- A_dd[1,3] # Females per adult female
    
    realized_fecundity[4,y-1,i] <- A_dd[4,1] # Males per yearling female 
    realized_fecundity[5,y-1,i] <- A_dd[4,2] # Males per yearling female
    realized_fecundity[6,y-1,i] <- A_dd[4,3] # Males per adult female
    
    # Calculate new pop size
    N_next <- pmax(A_dd %*% pig.array[,y-1,i],0) # Make sure to multiply matrix x vector (not vice versa)
    
    ## Harvest pigs 
    if (y > 0 & min(N_next)>0 & sum(!is.na(N_next))==6) {
      # Randomly select a random number to harvest
      h <- round(rnorm(1, 300000, 20000), digits=0)
      
      # Proportional harvest
      r <- (N_next[,1] / sum(N_next[,1]))*h
      
      # remove
      N_next <- N_next - r
    } 
    
    # Save abundance
    pig.array[,y,i] <- pmax(N_next, 0)
    
    # ---- Realized survival ----
    
    # Female stages
    realized_survival[1, y, i] <- safe_divide(pig.array[2, y, i], pig.array[1, y-1, i]) # Piglet → yearling (F)
    
    # female ssd
    w_f <- Re(eigen(A_dd)$vectors[, 1])[1:3]
    w_f <- w_f / sum(w_f)
    
    # male ssd
    w_m <- Re(eigen(A_dd)$vectors[, 1])[4:6]
    w_m <- w_m / sum(w_m)
    
    c5_f <- (A_dd[3,2] * w_f[2]) / ((A_dd[3,2] * w_f[2]) + (A_dd[3,3] * w_f[3]))
    c5_m <- (A_dd[6,5] * w_m[2]) / ((A_dd[6,5] * w_m[2]) + (A_dd[6,6] * w_m[3]))
    
    # Stages 5 and 6 are tricky...
    est_from5 <- pig.array[3, y, i] * c5_f
    est_from6 <- pig.array[3, y, i] - est_from5
    
    realized_survival[2, y, i] <- safe_divide(est_from5, pig.array[2, y-1, i])
    
    # Stage 6 is tricky: receives from 5 and from itself
    incoming_females <- pig.array[3, y-1, i] + pig.array[2, y-1, i]  # crude approx
    
    if (incoming_females > 0) {
      realized_survival[3, y, i] <- pig.array[3, y, i] / incoming_females
    } 
    
    # Male stages
    realized_survival[4, y, i] <- safe_divide(pig.array[5, y, i], pig.array[4, y-1, i]) # Piglet → yearling (M)
    
    # Stages 5 and 6 are tricky...
    est_from11 <- pig.array[6, y, i] * c5_m
    est_from12 <- pig.array[6, y, i] - est_from11
    
    realized_survival[5, y, i] <- safe_divide(est_from11, pig.array[5, y-1, i])
    
    # Stage 12 (adult males): 11 → 12 and 12 → 12
    incoming_males <- pig.array[6, y-1, i] + pig.array[5, y-1, i]
    if (incoming_males > 0) {
      realized_survival[6, y, i] <- pig.array[6, y, i] / incoming_males
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
N.20pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct <- apply(apply(pig.array,c(2,3),sum),1,quantile, probs = 0.80)

results <- data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1000000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 



## Get realized survival rates
# subset realized survival to last 10 years
last_10_years <- (max(Year)-10):(max(Year)-1)
realized_survival_subset <- realized_survival[, last_10_years, ]

# Calculate median survival per stage (over all years and sims)
surv.50pct <- apply(realized_survival , 1, function(x) median(x, na.rm = TRUE))
surv.10pct <- apply(realized_survival, 1, quantile, probs = 0.10, na.rm=TRUE)
surv.90pct <- apply(realized_survival, 1, quantile, probs = 0.90, na.rm=TRUE)

r.surv <- data.frame(sex=c(rep("Female",6), rep("Male",6)), 
                     stage=rep(c("Piglet","Yearling","Adult"),2),
                     surv.50pct,surv.10pct,surv.90pct,
                     est=c(A_adj[2,1], A_adj[3,2], A_adj[3,3], A_adj[5,4], A_adj[6,5], A_adj[6,6]),
                     literature=c(A_base[2,1], A_base[3,2], A_base[3,3], A_base[5,4], A_base[6,5], A_base[6,6]))
library(dplyr)
r.surv <- r.surv |>
  pivot_longer(cols=c("surv.50pct", "est", "literature"), names_to="source", values_to="survival")


# Set factor levels
r.surv$stage <- factor(r.surv$stage, levels=c("Piglet","Yearling","Adult"),
                       labels=c("Piglet","Yearling","Adult"))

r.surv$source <- factor(r.surv$source, levels=c("literature", "est", "surv.50pct"), labels=c("Literature", "Optimized", "Realized"))

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
