library(ggplot2)
library(dplyr)
library(ggtext)
library(patchwork)

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
A_base[2,1] <- 0.27
A_base[3,2] <- 0.93
A_base[4,3] <- 0.84
A_base[5,4] <- 0.84
A_base[6,5] <- 0.84
A_base[6,6] <- 0.84
# Males
A_base[8,7] <- 0.27
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
# 1.3:1 males:females - splits offspring into males & females
pct_m <- 0.565
pct_f <- 0.435

# # Calculate for post-breeding census (need to multiply fecundity by doe survival)
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

# Caps on survival (maximum values)
female_cap <- 0.97
male_cap <- 0.93

female_base <- c(0.93, 0.84, 0.84, 0.84, 0.84)
female_upper <- female_cap / female_base

male_base <- c(0.82, 0.63, 0.53, 0.44, 0.49)
male_upper <- male_cap / male_base

# Caps on survival (minimum values)
female_min <- 0.75
male_min <- 0.4

female_lower <- female_min / female_base
male_lower <- male_min / male_base

# set up bounds (Fawn survival (1), female survival (5), male survival (5), fecundity (1), theta(1))
deer_lower <- c(0.6, # Fawn survival
                female_lower, # Female survival
                male_lower, # Male survival
                0.75,  # Fecundity
                3) # Theta

deer_upper <- c(2.2,          # Fawn survival
                female_upper, # female survival
                male_upper,   # male survival
                1.5,          # Fecundity
                5)            # Theta

theta_mult <- 2
c_dd <- 0.20

# run optimizer
set.seed(1)
deer_opt <- optim(par=rep(1.4, 13), fn=objective_fn_deer, method="L-BFGS-B", lower=deer_lower, upper=deer_upper, control = list(trace = 1, maxit=1000))

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


#### Run population model ####
set.seed(1)
res14 <- run_deer_mod(A_adj, theta=deer_opt$par[13], K=1347232, Sims=1000, years=30, ev_sd=0.02, harvest=TRUE, theta_mult=theta_mult, c_dd=c_dd) 

plot(res14$three_yr_males,type="l", col="purple", ylim=c(0, 100000))
lines(res14$four_yr_males, col="blue")
lines(res14$five_yr_males, col="darkgreen")

# Lambda from adjusted matrix
res14$matrix_lambda

# Lambda from final matrix
res14$final_lambda

# Plot population projection
ggplot(res14$results) +
  coord_cartesian(ylim=c(0,2500000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +  # Estimated current population size
  # geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) + # Carrying capacity
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

ggplot(res14$r.surv) +
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

ggplot(res14$r.fec) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,2)) +
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

# Plot sex ratio
ggplot(res14$ASR) +
  geom_hline(yintercept=1, linetype=2) +
  geom_hline(yintercept=0.4, linetype=2, color="red") +
  coord_cartesian(ylim=c(0,1.5)) +
  geom_ribbon(aes(x=Year, ymin=ASR.20, ymax=ASR.80), alpha=0.5) +
  geom_line(aes(x=Year, y=ASR.median)) +
  ylab("ASR (males/females)") +
  theme_classic()

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
A_base[2,1] <- 0.23
A_base[3,2] <- 0.93
A_base[4,3] <- 0.84
A_base[5,4] <- 0.84
A_base[6,5] <- 0.84
A_base[6,6] <- 0.84
# Males
A_base[8,7] <- 0.23
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
# 1.3:1 males:females - splits offspring into males & females
pct_m <- 0.565
pct_f <- 0.435

# # Calculate for post-breeding census (need to multiply fecundity by doe survival)
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
# set up bounds (Fawn survival (1), female survival (5), male survival (5), fecundity (1), theta(1))
deer_lower <- c(0.95, # Fawn survival
                0.85, # Female survival
                0.85,
                0.85,
                0.85, 
                0.85, 
                0.99, # Male survival
                0.95,
                0.95,
                0.95,
                0.95, 
                0.75,  # Fecundity
                0.9) # Theta

deer_upper <- c(3.5, # Fawn survival
                1.05, # Female survival
                1.05, 
                1.05, 
                1.05,
                1.05, 
                1.05, # Male survival
                1.5, 
                1.8,
                1.9, 
                1.9,
                2.5, # Fecundity
                6) # Theta

# run optimizer
set.seed(1)
deer_opt <- optim(par=rep(1.2, 13), fn=objective_fn_deer, method="L-BFGS-B", lower=deer_lower, upper=deer_upper, control = list(trace = 1))

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
res22 <- run_deer_mod(A_adj, theta=deer_opt$par[13], K=2090532, Sims=1000, years=50) 

#Plot population projection
ggplot(res22$results) +
  coord_cartesian(ylim=c(0,2800000)) +
  geom_hline(yintercept=1610000, color="red", linetype=3) +
  # geom_hline(yintercept=1750000, color="red", linetype=3) +
  geom_hline(yintercept=K) +
  geom_ribbon(aes(x=Year,ymin=N.10pct, ymax=N.90pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (total population)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

ggplot(res22$r.surv) +
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

ggplot(res22$r.fec) +
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



#### Save results from both densities ####
results14 <- res14$results
results22 <- res22$results

# Combine population projections
results14$density <- "14.3 deer km<sup>-2</sup>"
results22$density <- "22 deer km<sup>-2</sup>"
results <- bind_rows(results14, results22)

write_csv(results, "results/deer_population_projections.csv")

# Combine survival
surv14 <- res14$r.surv
surv22 <- res22$r.surv

surv14$density <- "14.3 deer km<sup>-2</sup>"
surv22$density <- "22 deer km<sup>-2</sup>"

survival <- bind_rows(surv14, surv22)

write_csv(survival, "results/deer_survival.csv")

# Combine fecundity
fec14 <- res14$r.fec
fec22 <- res22$r.fec

fec14$density <- "14.3 deer km<sup>-2</sup>"
fec22$density <- "22 deer km<sup>-2</sup>"

fecundity <- bind_rows(fec14, fec22)

write_csv(fecundity, "results/deer_fecundity.csv")







### --- Manuscript figure ---------


res_plot <- ggplot(results) +
  coord_cartesian(ylim=c(0,3000000)) +
  geom_hline(yintercept=1610000, color="black", linetype=3) +
  geom_hline(yintercept=1347232, linetype=2, color="#ed7953") +
  geom_hline(yintercept=2090532, linetype=2, color="#0d0887") +
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

s_plot <- ggplot(survival) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=stage, y=surv.10pct, yend=surv.90pct), color="#21918c", linewidth=1) +
  geom_point(aes(x=stage, y=survival, color=source), size=3) +
  scale_color_manual(values=c("#440154", "#21918c")) +
  # scale_shape_manual(values=c(17, 18, 16)) +
  guides(
    # shape=guide_legend(position="inside", title="Source"),
    color=guide_legend(position="inside", title="Source")
  ) +
  # geom_text(data=s_text, aes(x=x, y=1, label=label)) +
  facet_grid(sex~density) +
  ylab("Survival probability") + xlab("Stage") +
  theme_classic() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        legend.position.inside = c(1,0),
        legend.justification=c(1,0),
        legend.background = element_rect(fill=NA),
        strip.background=element_blank(),
        axis.text.x=element_text(angle=25, vjust=0.8, hjust=0.8, size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=12),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12))#,
        # strip.text=element_blank())

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
