library(ggplot2)
library(ggtext)
library(patchwork)

source("00_functions.R")

theme_set(theme_bw())

### K at 27 pigs/ km2 - Largest area (Core + Moderate + Marginal) ####

## Stochastistic model
step <- 1:200

# Carrying capacity
K <- 1239278

Kf <- K * 0.5 # Carrying capacity of only females

# Set up pig matrix
A_base <- matrix(0,6,6)

## Survival
# Females
A_base[2,1] <- sqrt(0.35)  # Piglet survival
A_base[3,2] <- sqrt(0.35) # Subadult survival
A_base[3,3] <- sqrt(0.35) # Adult survival
# Males
A_base[5,4] <- sqrt(0.35) # Piglet survival
A_base[6,5] <- sqrt(0.23) # Subadult survival
A_base[6,6] <- sqrt(0.23) # Adult survival

## Fecundity
# Maximum fecundity (how many piglets per litter?) x how many females produce at a given time step
R0j <- 4.8*0.75
R0a <- 6.4

FecundityM <- c(R0j, R0a)
FecundityF <- c(R0j, R0a)

# Add to matrix
A_base[1,2:3] <- FecundityF 
A_base[4,2:3] <- FecundityM 

# Calculate stable stage distribution
w <- Re(eigen(A_base)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(A_base)$values[1])

# How many pig are removed?
observed_harvest <- 150000 

### Optimize survival and fecundities
# Caps on survival (maximum values)
female_cap <- sqrt(0.65)
male_cap <- sqrt(0.60)

female_base <- c(sqrt(0.35), sqrt(0.35))
female_upper <- female_cap / female_base

male_base <- c(sqrt(0.23), sqrt(0.23))
male_upper <- male_cap / male_base

# Caps on survival (minimum values)
female_min <- sqrt(0.25)
male_min <- sqrt(0.25)

female_lower <- female_min / female_base
male_lower <- male_min / male_base

# set up bounds (Fawn survival (1), female survival (2), male survival (2), fecundity (1), theta(1))
pigs_lower <- c(0.3,          # Piglet survival (0-6 months)
                0.3,          # Piglet survival (6-12 months)
                female_lower, # Female survival
                male_lower,   # Male survival
                1.2,          # Fecundity
                0.001)        # Theta

pigs_upper <- c(1.2,          # Piglet survival (0-6 months)
                1.4,          # Piglet survival (6-12 months)
                female_upper, # female survival
                male_upper,   # male survival
                2,          # Fecundity
                1)            # Theta

c_dd <- 0.5

# run optimizer
set.seed(1)
pig_opt <- optim(par=runif(8, 0.9, 1.2), fn=objective_fn_pigs, method="L-BFGS-B", lower=pigs_lower, upper=pigs_upper, control=list(trace=1, maxit=1000))

# What are optimized parameter values 
pig_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Female survival
A_adj[2,1] <- A_base[2,1]*pig_opt$par[2]
A_adj[3,2] <- A_base[3,2]*pig_opt$par[3]
A_adj[3,3] <- A_base[3,3]*pig_opt$par[4]

# Male survival
A_adj[5,4] <- A_base[5,4]*pig_opt$par[2]
A_adj[6,5] <- A_base[6,5]*pig_opt$par[5]
A_adj[6,6] <- A_base[6,6]*pig_opt$par[6]

## Fecundity
# Yearlings
A_adj[1,2] <- A_base[1,2]*(pig_opt$par[7])
A_adj[4,2] <- A_base[4,2]*(pig_opt$par[7])
# Adults
A_adj[1,3] <- A_base[1,3]*(pig_opt$par[7])
A_adj[4,3] <- A_base[4,3]*(pig_opt$par[7])


a2 <- A_adj
a2[1, ] <- a2[1,]/2
a2[4, ] <- a2[4,]/2
Re(eigen(A_base)$values[1])^2

#### Run population model ####
set.seed(1)
res27 <- run_pig_mod(A_adj, neonate_surv=pig_opt$par[1], theta=pig_opt$par[8], K=1239278, Sims=1000, steps=200, ev_sd=0.05, harvest=TRUE, c_dd=c_dd) 

# Lambda from adjusted matrix
res27$matrix_lambda

# Lambda from final matrix
res27$final_lambda

#Plot population projection
ggplot(res27$results) +
  coord_cartesian(ylim=c(0,4000000)) +
  geom_hline(yintercept=1000000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.10pct, ymax=N.90pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)",
                     breaks=c(0,500000,1000000,1500000, 2000000, 2500000, 3000000, 3500000, 4000000),
                     labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
  scale_x_continuous(breaks=c(0,10,20,30,40), labels=c(0,5,10,15,20), name="Simulation Year") +
  theme_classic() + 
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, color="black")) 

ggplot(res27$r.surv) +
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

ggplot(res27$r.fec) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,4)) +
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
        legend.direction="horizontal",
        legend.background = element_rect(fill=NA),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11))



### K at 8 pigs/ km2 - full state ####
set.seed(1)
## Stochastistic model
step <- 1:200

# Carrying capacity
K <- 367822

Kf <- K * 0.5

# Set up pig matrix
A_base <- matrix(0,6,6)

## Survival
# Females
A_base[2,1] <- sqrt(0.35)  # Piglet survival
A_base[3,2] <- sqrt(0.35) # Subadult survival
A_base[3,3] <- sqrt(0.35) # Adult survival
# Males
A_base[5,4] <- sqrt(0.23) # Piglet survival
A_base[6,5] <- sqrt(0.23) # Subadult survival
A_base[6,6] <- sqrt(0.23) # Adult survival

## Fecundity
# Maximum fecundity (how many piglets per litter?) x how many females produce at a given time step
R0j <- 4.8*0.75
R0a <- 6.4

FecundityM <- c(R0j, R0a)
FecundityF <- c(R0j, R0a)

# Add to matrix
A_base[1,2:3] <- FecundityF 
A_base[4,2:3] <- FecundityM 

# Calculate stable stage distribution
w <- Re(eigen(A_base)$vectors[, 1])
w <- w / sum(w) # normalize to equal 1

# Set up initial popualtion size
N0 <- w * 0.1 * K  # e.g., start at 50% of K

lambda <- Re(eigen(A_base)$values[1])

# How many pig are removed?
observed_harvest <- 150000 

### Optimize survival and fecundities
# Caps on survival (maximum values)
female_cap <- sqrt(0.60)
male_cap <- sqrt(0.55)

female_base <- c(sqrt(0.35), sqrt(0.35))
female_upper <- female_cap / female_base

male_base <- c(sqrt(0.23), sqrt(0.23))
male_upper <- male_cap / male_base

# Caps on survival (minimum values)
female_min <- sqrt(0.30)
male_min <- sqrt(0.25)

female_lower <- female_min / female_base
male_lower <- male_min / male_base

# set up bounds (Fawn survival (1), female survival (2), male survival (2), fecundity (1), theta(1))
pigs_lower <- c(0.2,          # Piglet survival (0-6 months)
                0.3,          # Piglet survival (6-12 months)
                female_lower, # Female survival
                male_lower,   # Male survival
                1.2,          # Fecundity
                0.001)        # Theta

pigs_upper <- c(1.6,          # Piglet survival (0-6 months)
                1.2,          # Piglet survival (6-12 months)
                female_upper, # female survival
                male_upper,   # male survival
                2,          # Fecundity
                1)            # Theta

c_dd <- 0.12

# run optimizer
set.seed(1)
pig_opt <- optim(par=runif(8, 0.9, 1.2), fn=objective_fn_pigs, method="L-BFGS-B", lower=pigs_lower, upper=pigs_upper, control=list(trace=1, maxit=1000))

# What are optimized parameter values 
pig_opt$par

#### Adjusted matrix ####
A_adj <- A_base

# Female survival
A_adj[2,1] <- A_base[2,1]*pig_opt$par[2]
A_adj[3,2] <- A_base[3,2]*pig_opt$par[3]
A_adj[3,3] <- A_base[3,3]*pig_opt$par[4]

# Male survival
A_adj[5,4] <- A_base[5,4]*pig_opt$par[2]
A_adj[6,5] <- A_base[6,5]*pig_opt$par[5]
A_adj[6,6] <- A_base[6,6]*pig_opt$par[6]

## Fecundity
# Yearlings
A_adj[1,2] <- A_base[1,2]*(pig_opt$par[7])
A_adj[4,2] <- A_base[4,2]*(pig_opt$par[7])
# Adults
A_adj[1,3] <- A_base[1,3]*(pig_opt$par[7])
A_adj[4,3] <- A_base[4,3]*(pig_opt$par[7])

a2 <- A_adj
a2[1, ] <- a2[1,]/2
a2[4, ] <- a2[4,]/2
Re(eigen(A_base)$values[1])^2


#### Run population model ####
set.seed(1)
res8 <- run_pig_mod(A_adj, neonate_surv=pig_opt$par[1], theta=pig_opt$par[7], K=367822, Sims=1000, steps=200, ev_sd=0.05, harvest=TRUE, c_dd=c_dd) 

# Lambda from adjusted matrix
res8$matrix_lambda

# Lambda from final matrix
res8$final_lambda

# Plot population projection
ggplot(res8$results) +
  coord_cartesian(ylim=c(0,2000000)) +
  geom_hline(yintercept=1000000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.10pct, ymax=N.90pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  theme_classic() + 
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, color="black")) 

ggplot(res8$r.surv) +
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


ggplot(res8$r.fec) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,4)) +
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
        legend.direction="horizontal",
        legend.background = element_rect(fill=NA),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11))


#### Combine and plot for manuscript ####
res27$results$density <- "27 pigs km<sup>-2</sup>"
res8$results$density <- "8 pigs km<sup>-2</sup>"

# Population projections
results <- bind_rows(res27$results,res8$results)

write_csv(results, "results/pigs_population_projections.csv")

# Survival
# Combine survival
surv8 <- res8$r.surv
surv27 <- res27$r.surv

surv8$density <- "8 pigs km<sup>-2</sup>"
surv27$density <- "27 pigs km<sup>-2</sup>"

surv <- bind_rows(surv8,surv27)

write_csv(surv, "results/pigs_survival.csv")

# Fecundity
fec8 <- res8$r.fec
fec27 <- res27$r.fec

fec8$density <- "8 pigs km<sup>-2</sup>"
fec27$density <- "27 pigs km<sup>-2</sup>"

fec <- bind_rows(fec8,fec27)

write_csv(fec, "results/pig_fecundity.csv")



nsurv_8f <- res8$neonate_surv |>
  mutate(density = "8 pigs km<sup>-2</sup>",
         stage="Neonate",
         sex="Female")
nsurv_8m <- res8$neonate_surv |>
  mutate(density = "8 pigs km<sup>-2</sup>",
         stage="Neonate",
         sex="Male")


nsurv_27f <- res27$neonate_surv|>
  mutate(density = "27 pigs km<sup>-2</sup>",
         stage="Neonate",
         sex="Female")
nsurv_27m <- res27$neonate_surv|>
  mutate(density = "27 pigs km<sup>-2</sup>",
         stage="Neonate",
         sex="Male")

nsurv <- bind_rows(nsurv_8f, nsurv_8m, nsurv_27f, nsurv_27m)
write_csv(nsurv, "results/pigs_neonate_survival.csv")

