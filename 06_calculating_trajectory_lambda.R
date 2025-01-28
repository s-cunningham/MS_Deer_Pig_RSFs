library(tidyverse)
library(patchwork)

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
  mutate(file=gsub(".txt", "", x=file),
         file=gsub("rmax", "", x=file)) %>%
  mutate(species=str_split_i(file, "_", i=1),
         patch=str_split_i(file, "_", i=2),
         density=str_split_i(file, "_", i=3)) %>%
  mutate(density=if_else(density=="14", "14.3",density)) %>%
  rename(density_int=density) %>%
  mutate(density=if_else(density_int=="14.3" | density_int=="8", "Low", "High")) %>%
  mutate(level="0") %>%
  select(species, patch, density, level, lower:upper) #%>%
  # mutate(sd=((upper-lambda)+(lambda-lower))/2,
  #        se=sd/sqrt(5000),
  #        lCI=lambda-(1.96*se),
  #        uCI=lambda+(1.96*se))

#### List trajectory files ####
files <- list.files(path="results/sensitivity_analysis/", pattern=".txt", full.names=TRUE)

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
  
  # Combine with filename
  res <- data.frame(file=str_split(files[i], pattern="/")[[1]][3], lambda=lmda,
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
                      patch=str_split_i(file, "_", i=2),
                      density=str_split_i(file, "_", i=3),
                      level=str_split_i(file, "_", i=4)) %>%
                mutate(density=if_else(density=="14", "14.3",density)) %>%
                rename(density_int=density) %>%
                mutate(density=if_else(density_int=="14.3" | density_int=="8", "Low", "High")) %>%
                select(species, patch, density, level, lower:upper)

results <- bind_rows(sims, results)


results$level <- factor(results$level, levels=c("-40","-30","-20","-10","0", "10","20","30","40"))
results$species <- factor(results$species, levels=c("deer", "pigs"), labels=c("White-tailed Deer", "Wild Pigs"))
# results$density <- factor(results$density, levels=c("Low", "High"), labels=c("Low Density", "High Density"))
results$patch <- factor(results$patch, levels=c("core", "marginal"), labels=c("Core", "Core+Marginal"))

results <- results %>% filter(patch!="highlymarginal" & density=="Low")

ggplot(results) +
  geom_vline(xintercept=5, color="gray90") +
  geom_hline(yintercept=1, color="gray5") +
  geom_segment(aes(x=level, y=lower, yend=upper, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  geom_point(aes(x=level, y=lambda, group=patch, color=patch), position=position_dodge(width = 0.5)) +
  scale_color_manual(values=c("#0d0887", "#cc4778"), name="Patch Type") +
  facet_grid(.~species) +
  ylab(expression("Finite rate of increase ("*lambda*")")) + xlab("Change in maximum growth rate (%)") +
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


## N0 models (pigs marginal)

N0_files <- list.files(path="results/initial_abundance_sensitivity/", pattern=".txt", full.names=TRUE)
results <- data.frame()
for (i in 1:length(N0_files)) {
  
  # Read File
  temp <- read_sim(N0_files[i], skip_rows=15)
  
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
  res <- data.frame(file=str_split(N0_files[i], pattern="/")[[1]][3], lambda=lmda,
                    lower=gm_mean(lSD), upper=gm_mean(uSD))
  
  # Add to results df
  results <- bind_rows(results, res)
  
}
results



