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

Kf <- K * 0.5 # Carrying capacity of only females

# Set up pig matrix
A_base <- matrix(0,6,6)

## Survival
# Females
A_base[2,1] <- sqrt(0.44)  # Piglet survival
A_base[3,2] <- sqrt(0.35) # Subadult survival
A_base[3,3] <- sqrt(0.35) # Adult survival
# Males
A_base[5,4] <- sqrt(0.44) # Piglet survival
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
female_cap <- sqrt(0.55)
male_cap <- sqrt(0.55)

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
pigs_lower <- c(0.3,          # Piglet survival
                female_lower, # Female survival
                male_lower,   # Male survival
                1,          # Fecundity
                0.001)        # Theta

pigs_upper <- c(1.6,          # Piglet survival
                female_upper, # female survival
                male_upper,   # male survival
                2.8,          # Fecundity
                1.5)            # Theta

c_dd <- 1.1

# run optimizer
set.seed(1)
pig_opt <- optim(par=runif(7, 0.9, 1.2), fn=objective_fn_pigs, method="L-BFGS-B", lower=pigs_lower, upper=pigs_upper, control=list(trace=1, maxit=1000))

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


a2 <- A_adj
a2[1, ] <- a2[1,]/2
a2[4, ] <- a2[4,]/2
Re(eigen(A_base)$values[1])^2

#### Run population model ####
set.seed(1)
res27 <- run_pig_mod(A_adj, theta=pig_opt$par[7], K=1239278, Sims=1000, steps=100, ev_sd=0.05, harvest=TRUE, c_dd=c_dd) 

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
A_base[2,1] <- sqrt(0.44)  # Piglet survival
A_base[3,2] <- sqrt(0.35) # Subadult survival
A_base[3,3] <- sqrt(0.35) # Adult survival
# Males
A_base[5,4] <- sqrt(0.44) # Piglet survival
A_base[6,5] <- sqrt(0.23) # Subadult survival
A_base[6,6] <- sqrt(0.23) # Adult survival

## Fecundity
# Maximum fecundity (how many piglets per litter?) x how many females produce at a given time step
R0j <- 4.8*0.75
R0a <- 6.4*0.95

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
female_cap <- sqrt(0.55)
male_cap <- sqrt(0.55)

female_base <- c(sqrt(0.35), sqrt(0.35))
female_upper <- female_cap / female_base

male_base <- c(sqrt(0.23), sqrt(0.23))
male_upper <- male_cap / male_base

# Caps on survival (minimum values)
female_min <- sqrt(0.4)
male_min <- sqrt(0.4)

female_lower <- female_min / female_base
male_lower <- male_min / male_base

# set up bounds (Fawn survival (1), female survival (2), male survival (2), fecundity (1), theta(1))
pigs_lower <- c(0.3,          # Piglet survival
                female_lower, # Female survival
                male_lower,   # Male survival
                0.8,          # Fecundity
                0.001)        # Theta

pigs_upper <- c(1.6,          # Piglet survival
                female_upper, # female survival
                male_upper,   # male survival
                2.6,          # Fecundity
                1.5)            # Theta

c_dd <- 2 #1.9

# run optimizer
set.seed(1)
pig_opt <- optim(par=runif(7, 0.9, 1.2), fn=objective_fn_pigs, method="L-BFGS-B", lower=pigs_lower, upper=pigs_upper, control=list(trace=1, maxit=1000))

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

a2 <- A_adj
a2[1, ] <- a2[1,]/2
a2[4, ] <- a2[4,]/2
Re(eigen(A_base)$values[1])^2

#### Run population model ####
set.seed(1)
res8 <- run_pig_mod(A_adj, theta=pig_opt$par[7], K=367822, Sims=1000, steps=100, ev_sd=0.05, harvest=TRUE, c_dd=c_dd) 

# Lambda from adjusted matrix
res8$matrix_lambda

# Lambda from final matrix
res8$final_lambda

# Plot population projection
ggplot(res8$results) +
  # coord_cartesian(ylim=c(0,2000000)) +
  geom_hline(yintercept=1000000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.10pct, ymax=N.90pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  # scale_y_continuous(name="Abundance (in millions)", 
  #                    breaks=c(0,500000,1000000,1500000,2000000),
  #                    labels=c(0,0.5, 1,1.5, 2)) +
  # scale_x_continuous(breaks=c(0,10,20,30,40), labels=c(0,5,10,15,20), name="Simulation Year") +
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
