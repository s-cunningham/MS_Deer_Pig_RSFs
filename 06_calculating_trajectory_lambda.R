library(tidyverse)
library(patchwork)
library(tagger)

# Function for calculating geometric mean
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Function for calculating vector of lambda values
lambda_calc <- function(x) {
  lambda <- numeric()
  for(y in 1:nrow(x)){
    lambda[y] <- x$Average[y+1]/x$Average[y]
  }
  return(lambda)
}

read_sim <- function(file, skip_rows=14) {

  x <- read_tsv(file, skip=skip_rows, n_max=11, show_col_types = FALSE) %>%
    # Change obnoxious column name
    rename(any_of(c(allcols="Minimum      Lower       Mean      Upper    Maximum", 
                    allcols="Minimum    -1 S.D.    Average    +1 S.D.    Maximum"))) %>%
    # because there are a different number of spaces, just want all whitespace to become a bar
    mutate(allcols = str_replace_all(allcols, "[[:space:]]+", "|")) %>%
    # now separate 
    separate_wider_delim(allcols,
                         names=c("Time", "Minimum", "lSD", "Average", "hSD", "Maximum"),
                         delim="|") %>%
    # convert caracter columns to numeric
    mutate_if(is.character, as.numeric)
  
  return(x)
}

### Read simulation results
sim_file <- list.files(path="results/simulations/", pattern=".txt", full.names=TRUE)

sims <- data.frame()
for (i in 1:length(sim_file)) {
  
  # Read File
  temp <- read_sim(sim_file[i])
  temp <- temp %>% mutate(sd=((hSD-Average) + (Average=lSD))/2)
  
  # Standard deviation lambda
  lSD <- numeric()
  uSD <- numeric()
  for(y in 1:nrow(temp)){
    lSD[y] <- temp$lSD[y+1]/temp$lSD[y]
    uSD[y] <- temp$hSD[y+1]/temp$hSD[y]
  }
  
  # Calculate geomlSD# Calculate geometric mean lambda
  lmda <- gm_mean(lambda_calc(temp))
  
  # Combine with filename
  res <- data.frame(file=str_split(sim_file[i], pattern="/")[[1]][3], lambda=lmda,
                    lower=gm_mean(lSD), upper=gm_mean(uSD))
  
  # Add to results df
  sims <- bind_rows(sims, res)
  
}
sims <- sims %>% as_tibble %>% 
  select(file, lower, lambda, upper) %>%
  mutate(file=gsub(".txt", "", x=file)) %>%
  mutate(species=str_split_i(file, "_", i=1),
         patch=str_split_i(file, "_", i=3),
         density=str_split_i(file, "_", i=2)) %>%
  mutate(density=if_else(density=="14", "14.3",density)) %>%
  rename(density_int=density) %>%
  mutate(density=if_else(density_int=="14.3" | density_int=="8", "Low", "High")) %>%
  mutate(level="0") %>%
  select(species, patch, density, level, lower:upper) %>%
  mutate(patch=case_when(patch=="c" ~ "core",
                         patch=="hm" ~ "highlymarginal",
                         patch=="m" ~ "marginal"))

#### List trajectory files for Rmax sensitivity analysis ####
files <- list.files(path="results/sensitivity_analysis/Rmax_sensitivity/", pattern=".txt", full.names=TRUE)

results <- data.frame()
for (i in 1:length(files)) {
  
  # Read File
  temp <- read_sim(files[i], skip_rows=15)
  
  # Standard deviation lambda
  lSD <- numeric()
  uSD <- numeric()
  for(y in 1:nrow(temp)){
    lSD[y] <- temp$lSD[y+1]/temp$lSD[y]
    uSD[y] <- temp$hSD[y+1]/temp$hSD[y]
  }
  
  # Calculate geometric mean lambda
  lmda <- gm_mean(lambda_calc(temp))
  sqrt(prod(lmda,na.rm=TRUE))
  
  # Combine with filename
  res <- data.frame(file=str_split(files[i], pattern="/")[[1]][4], lambda=lmda,
                    lower=gm_mean(lSD), upper=gm_mean(uSD))
  
  # Add to results df
  results <- bind_rows(results, res)
  
}
results

results <- results %>% as_tibble %>% 
               select(file, lower, lambda, upper) %>%
               mutate(file=gsub(".txt", "", x=file),
                      file=gsub("rmax", "", x=file)) %>%
               mutate(species=str_split_i(file, "_", i=1),
                      patch=str_split_i(file, "_", i=3),
                      density=str_split_i(file, "_", i=2),
                      level=str_split_i(file, "_", i=5)) %>%
                mutate(density=if_else(density=="14", "14.3",density)) %>%
                rename(density_int=density) %>%
                mutate(density=if_else(density_int=="14.3" | density_int=="8", "Low", "High")) %>%
                select(species, patch, density, level, lower:upper) %>%
                mutate(patch=case_when(patch=="c" ~ "core",
                         patch=="hm" ~ "highlymarginal",
                         patch=="m" ~ "marginal"))

results <- bind_rows(sims, results)


results$level <- factor(results$level, levels=c("-40","-30","-20","-10","0", "10","20","30","40"))
results$species <- factor(results$species, levels=c("deer", "pigs"), labels=c("White-tailed Deer", "Wild Pigs"))
# results$density <- factor(results$density, levels=c("Low", "High"), labels=c("Low Density", "High Density"))
results$patch <- factor(results$patch, levels=c("core", "marginal"), labels=c("Core", "Core+Marginal"))

results <- results %>% filter(patch!="highlymarginal" & density=="Low")

ggplot(results) +
  # geom_vline(xintercept=5, color="gray90") +
  geom_hline(yintercept=1, color="gray5") +
  geom_segment(aes(x=level, y=lower, yend=upper, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  geom_point(aes(x=level, y=lambda, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  scale_color_manual(values=c("#0d0887", "#cc4778"), name="Patch Type") +
  facet_grid(.~species) +
  ylab(expression("Finite rate of increase ("*lambda*")")) + xlab("Change in maximum growth rate (%)") +
  guides(color=guide_legend(position = "inside")) +
  theme_classic() +
  # tag_facets() +
  theme(strip.background=element_rect(color=NA,fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.line=element_line(linewidth=0),
        # legend.position.inside=c(0.01,0.99),
        legend.position.inside=c(0.3,0.15),
        legend.justification=c(0,1),
        strip.text=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11))
ggsave(file="figs/lambda_rmax_sensitivity.svg")
# Saving 9.42 x 5.22 in image
# Saving 8.5 x 4.97 in image

#### List simulation files for K sensitivity analysis ####
files <- list.files(path="results/sensitivity_analysis/carrying_capacity_sensitivity/", pattern=".txt", full.names=TRUE)

resultsK <- data.frame()
for (i in 1:length(files)) {
  
  # Read File
  temp <- read_sim(files[i], skip_rows=15)
  
  # Standard deviation lambda
  lSD <- numeric()
  uSD <- numeric()
  for(y in 1:nrow(temp)){
    lSD[y] <- temp$lSD[y+1]/temp$lSD[y]
    uSD[y] <- temp$hSD[y+1]/temp$hSD[y]
  }
  
  # Calculate geometric mean lambda
  lmda <- gm_mean(lambda_calc(temp))
  
  # Combine with filename
  res <- data.frame(file=str_split(files[i], pattern="/")[[1]][4], lambda=lmda,
                    lower=gm_mean(lSD), upper=gm_mean(uSD))
  
  # Add to results df
  resultsK <- bind_rows(resultsK, res)
  
}
resultsK

resultsK <- resultsK %>% as_tibble %>% 
  select(file, lower, lambda, upper) %>%
  mutate(file=gsub(".txt", "", x=file),
         file=gsub("K", "", x=file)) %>%
  mutate(species=str_split_i(file, "_", i=1),
         patch=str_split_i(file, "_", i=3),
         density=str_split_i(file, "_", i=2),
         level=str_split_i(file, "_", i=5)) %>%
  mutate(density=if_else(density=="14", "14.3",density)) %>%
  rename(density_int=density) %>%
  mutate(density=if_else(density_int=="14.3" | density_int=="8", "Low", "High")) %>%
  select(species, patch, density, level, lower:upper) %>%
  mutate(patch=case_when(patch=="c" ~ "core",
                         patch=="hm" ~ "highlymarginal",
                         patch=="m" ~ "marginal"))

resultsK <- bind_rows(sims, resultsK)

resultsK$level <- factor(resultsK$level, levels=c("-40","-30","-20","-10","0", "10","20","30","40"))
resultsK$species <- factor(resultsK$species, levels=c("deer", "pigs"), labels=c("White-tailed Deer", "Wild Pigs"))
# results$density <- factor(results$density, levels=c("Low", "High"), labels=c("Low Density", "High Density"))
resultsK$patch <- factor(resultsK$patch, levels=c("core", "marginal"), labels=c("Core", "Core+Marginal"))

resultsK <- resultsK %>% filter(patch!="highlymarginal" & density=="Low")

ggplot(resultsK) +
  # geom_vline(xintercept=5, color="gray90") +
  geom_hline(yintercept=1, color="gray5") +
  geom_segment(aes(x=level, y=lower, yend=upper, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  geom_point(aes(x=level, y=lambda, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  scale_color_manual(values=c("#0d0887", "#cc4778"), name="Patch Type") +
  facet_grid(.~species) +
  ylab(expression("Finite rate of increase ("*lambda*")")) + xlab("Change in carrying capacity (%)") +
  guides(color=guide_legend(position = "inside")) +
  theme_classic() +
  theme(strip.background=element_rect(color=NA,fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.line=element_line(linewidth=0),
        legend.position.inside=c(0.01,0.99),
        legend.justification=c(0,1),
        strip.text=element_text(size=11),
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11))

ggsave(file="figs/lambda_K_sensitivity.svg")
# Saving 9.42 x 5.22 in image














### Separate species, plot by density
deer <- results %>% filter(species=="White-tailed Deer") 
pigs <- results %>% filter(species=="Wild Pigs")

deer_labeller <- as_labeller(c(Low="14.3~deer/km^2", High="22~deer/km^2"), default=label_parsed)
deer$density <- factor(deer$density, levels=c("Low", "High"))

pigs_labeller <- as_labeller(c(Low="8~pigs/km^2", High="27~pigs/km^2"), default=label_parsed)
pigs$density <- factor(pigs$density, levels=c("Low", "High"))

deer_plot <- ggplot(deer) +
  coord_cartesian(ylim=c(0.95, 1.08)) +
  geom_hline(yintercept=1, color="gray") +
  geom_segment(aes(x=level, y=lower, yend=upper, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  geom_point(aes(x=level, y=lambda, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  scale_color_manual(values=c("#0d0887", "#cc4778"), name="Patch Type") +
  facet_wrap(~density, labeller=deer_labeller) +
  ylab(expression("Finite rate of increase ("*lambda*")")) + xlab("Change in maximum growth rate (%)") +
  guides(color=guide_legend(position = "inside")) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=1),
        strip.background=element_rect(color=NA,fill=NA),
        strip.text=element_text(size=11),
        axis.line=element_line(linewidth=0),
        legend.position.inside=c(0.01,0.99),
        legend.justification=c(0,1),
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.text=element_text(size=11))

pig_plot <- ggplot(pigs) +
  coord_cartesian(ylim=c(0.95, 1.08)) +
  geom_hline(yintercept=1, color="gray") +
  geom_segment(aes(x=level, y=lower, yend=upper, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  geom_point(aes(x=level, y=lambda, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  scale_color_manual(values=c("#0d0887", "#cc4778")) +
  facet_wrap(~density, labeller=pigs_labeller) +
  ylab(expression("Finite rate of increase ("*lambda*")")) + xlab("Change in maximum growth rate (%)") +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=1),
        strip.background=element_rect(color=NA,fill=NA),
        strip.text=element_text(size=11),
        axis.line=element_line(linewidth=0),
        legend.position="none",
        legend.justification=c(0,1),
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11))

deer_plot / pig_plot + plot_annotation(tag_levels="A", tag_prefix="(", tag_suffix=")") 


###  Geometric mean the other way
# lambda <- numeric()
# for(y in 1:nrow(temp)){
#   lambda[y] <- temp$Average[y+1]/temp$Average[y]
# }
# 
# nthroot <- function(x,n) {
#   (abs(x)^(1/n))*sign(x)
# }
# nthroot(prod(lambda, na.rm=TRUE),(length(lambda)-1))

