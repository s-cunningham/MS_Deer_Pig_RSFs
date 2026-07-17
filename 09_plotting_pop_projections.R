
library(tidyverse)
library(ggtext)
library(patchwork)

## Read demographic rates from literature ####
rates <- read_csv("data/literature_demographic_rates.csv") |>
  mutate(sex = ifelse(sex=="female", "Female", "Male"))

## Deer ---------------------------------------------------------------------------------------------------------
d_rates_s <- rates |>
  filter(param=="survival" & species=="deer")

d_rates_s$age <- factor(d_rates_s$age, 
                    levels=c("fawn","yearling","2-year-old","3-year-old","4-year-old","5+ years-old"), 
                    labels=c("Neonate","Yearling","2-year-old","3-year-old","4-year-old","5+ years-old"))

d_rates_s <- d_rates_s |>
  rename(stage = age) |>
  mutate(stage = as.character(stage)) |>
  mutate(x=case_when(stage=="Neonate" ~ 0.1,
                     stage=="Fawn" ~ 0.2,
                     stage=="Yearling" ~ 0.3,
                     stage=="2-year-old" ~ 0.4,
                     stage=="3-year-old" ~ 0.5,
                     stage=="4-year-old" ~ 0.6,
                     stage=="5+ years-old" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01))

d_rates_r <- rates |>
  filter(param=="fecundity" & species=="deer")|>
  rename(rate = survival)

## Pigs ------------------------------------------------------------------------------------------------------------
p_rates_s <- rates |>
  filter(param=="survival" & species=="pigs") |>
  mutate(age = ifelse(age=="Yearling", "Subadult", age))

p_rates_s$age <- factor(p_rates_s$age, 
                        levels=c("Neonate", "Piglet", "Subadult", "Adult"),
                        labels=c("Neonate", "Piglet", "Subadult", "Adult"))

p_rates_s <- p_rates_s |>
  rename(stage = age) |>
  mutate(stage = as.character(stage)) |>
  mutate(x=case_when(stage=="Neonate" ~ 0.4,
                     stage=="Piglet" ~ 0.5,
                     stage=="Subadult" ~ 0.6,
                     stage=="Adult" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01)) 

p_rates_r <- rates |>
  filter(param=="fecundity" & species=="pigs") |>
  rename(rate = survival) |>
  rename(stage=age, fec=rate) |>
  mutate(stage = ifelse(stage=="Yearling", "Subadult", stage))



## Read results ####

# Deer
deer <- read_csv("results/deer_population_projections.csv")
deer_s <- read_csv("results/deer_survival.csv") |>
  filter(source=="Implied") |>
  select(-source, -est)
deer_f <- read_csv("results/deer_fecundity.csv")
deer_neo_s <- read_csv("results/deer_neonate_survival.csv") |>
  rename(survival=surv.50pct)

deer_s <- bind_rows(deer_s, deer_neo_s)

# Pigs
pigs <- read_csv("results/pigs_population_projections.csv")
pigs_s <- read_csv("results/pigs_survival.csv") |>
  filter(source=="Implied") |>
  select(-source)
pigs_f <- read_csv("results/pig_fecundity.csv") |>
  filter(source=="Implied")
pigs_neo_s <- read_csv("results/pigs_neonate_survival.csv") |>
  rename(survival=surv.50pct)
pigs_s <- bind_rows(pigs_s, pigs_neo_s)

## Plot ####

### Deer ####
## Population projection
deer_pop <- ggplot(deer) +
  coord_cartesian(ylim=c(0,3000000)) +
  geom_hline(yintercept=1610000, color="gray", linewidth=1) +
  geom_hline(yintercept=1347232, linetype=2, color="#ed7953", linewidth=1, alpha=0.6) +
  geom_hline(yintercept=2090532, linetype=2, color="#9c179e", linewidth=1, alpha=0.6) +
  geom_ribbon(aes(x=Year, ymin=N.10pct, ymax=N.90pct, fill=density, group=density), alpha=0.2) +
  geom_line(aes(x=Year, y=N.median, color=density, group=density), alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000, 1000000, 1500000, 2000000, 2500000, 3000000),
                     labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3))+
  scale_color_manual(values=c("#ed7953","#9c179e")) +
  scale_fill_manual(values=c("#ed7953","#9c179e")) +
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

## Survival
deer_s <- deer_s |>
  mutate(x=case_when(stage=="Neonate" ~ 0.1,
                     stage=="Fawn" ~ 0.2,
                     stage=="Yearling" ~ 0.3,
                     stage=="2-year-old" ~ 0.4,
                     stage=="3-year-old" ~ 0.5,
                     stage=="4-year-old" ~ 0.6,
                     stage=="5+ years-old" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01)) #|>
  # filter(density=="14.3 deer km<sup>-2</sup>")
  # filter(stage!="Fawn")

## Plot 14 deer/km2
ds_plot <- ggplot() +
  coord_cartesian(ylim=c(0,1)) +
  geom_point(data=d_rates_s, aes(x=x, y=survival), shape=21, color="gray10", fill="gray",  size=3, alpha=0.45) +
  geom_segment(data=deer_s, aes(x=x, y=surv.10pct, yend=surv.90pct, color=sex), linewidth=1) +
  geom_point(data=deer_s, aes(x=x, y=survival, color=sex), size=3) +
  guides(
    color=guide_legend(position="inside", title="Sex"),
  ) +
  scale_color_manual(values=c("#21918c", "#440154")) +
  scale_x_continuous(breaks=seq(0.1,0.7, by=0.1), 
                     labels=c("Neonate*", "Fawn", "Yearling", "2 years", "3 years", "4 years", "5+ years"), 
                     expand = expansion(mult = 0.15, add = 0)) +
  ylab("Survival probability") + xlab("Stage") +
  facet_grid(.~density) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position.inside=c(0.15, 0.3),
        legend.justification = c(0,1),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_markdown(size=12),
        legend.text=element_text(size=11),
        legend.title=element_blank())

## Fecundity
deer_f <- deer_f |>
  # replace 'Realized' with 'Implied'
  filter(source=="Implied") |>
  arrange(stage) |>
  group_by(density, stage) |>
  reframe(fec=sum(fec),
          f10pct=sum(f10pct),
          f90pct=sum(f90pct)) |>
  mutate(x=case_when(stage=="Yearling" ~ 0.1,
                     stage=="Adult" ~ 0.2))

deer_f$stage <- factor(deer_f$stage, levels=c("Yearling", "Adult"), labels=c("Yearling", "Adult\n(2-5+ years)"))

d_rates_r <- d_rates_r |>
  rename(stage=age)

d_rates_r$stage <- factor(d_rates_r$stage, levels=c("Yearling", "Adult"), labels=c("Yearling", "Adult\n(2-5+ years)"))


d_rates_r <- d_rates_r |>
  group_by(title, stage) |>
  reframe(fawns = sum(rate))


df_plot <- ggplot(deer_f) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,2)) +
  geom_point(data=d_rates_r, aes(x=stage, y=fawns), shape=21, color="gray10", fill="gray",  size=3, alpha=0.45) +
  geom_segment(aes(x=stage, y=f10pct, yend=f90pct), linewidth=1) +
  geom_point(aes(x=stage, y=fec), size=3) +
  facet_grid(.~density) +
  ylab("Fecundity (fawns/doe)") + xlab("Stage") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_markdown(size=12),
        legend.text=element_text(size=11),
        legend.title=element_blank())



(deer_pop | df_plot) / ds_plot + plot_annotation(tag_levels='a', tag_prefix = "(", tag_suffix=")")


ggsave("figs/deer_implied_rates_6.svg", height=8, width=10, units="in") #Saving 10 x 8 in image

### Pigs ####

## Population projection
pigs_pop <- ggplot(pigs) +
  coord_cartesian(ylim=c(0,3000000)) +
  geom_hline(yintercept=1000000, color="gray", linewidth=1) +
  geom_hline(yintercept=367822, linetype=2, color="#9c179e", linewidth=1) +
  geom_hline(yintercept=1239278, linetype=2, color="#ed7953", linewidth=1) +
  geom_ribbon(aes(x=Year, ymin=N.10pct, ymax=N.90pct, fill=density, group=density), alpha=.2) +
  geom_line(aes(x=Year, y=N.median, color=density, group=density), alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000, 1000000, 1500000, 2000000, 2500000, 3000000),
                     labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  scale_x_continuous(breaks=seq(0, 200, by=20), labels=seq(0,100, by=10)) +
  scale_color_manual(values=c("#ed7953","#9c179e")) +
  scale_fill_manual(values=c("#ed7953","#9c179e")) +
  guides(
    fill=guide_legend(position="inside", title="Density"),
    color=guide_legend(position="inside", title="Density")
  ) +
  xlab("Simulation Year") +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
        legend.position.inside=c(0,1),
        legend.justification=c(0,1),
        legend.background=element_rect(fill=NA),
        legend.text=element_markdown(size=11),
        legend.title=element_text(size=12),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12)) 

## Survival
pigs_s <- pigs_s |>
  mutate(x=case_when(stage=="Neonate" ~ 0.4,
                     stage=="Piglet" ~ 0.5,
                     stage=="Yearling" ~ 0.6,
                     stage=="Adult" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01))

pigs_s$density <- factor(pigs_s$density, levels=c("8 pigs km<sup>-2</sup>","27 pigs km<sup>-2</sup>"))


ps_plot <- ggplot() +
  coord_cartesian(ylim=c(0,1)) +
  geom_point(data=p_rates_s, aes(x=x, y=survival), shape=21, color="gray10", fill="gray",  size=3, alpha=0.45) +
  geom_segment(data=pigs_s, aes(x=x, y=surv.10pct, yend=surv.90pct, color=sex), linewidth=1) +
  geom_point(data=pigs_s, aes(x=x, y=survival, color=sex), size=3) +
  guides(
    color=guide_legend(position="inside", title="Sex"),
  ) +
  scale_color_manual(values=c("#21918c", "#440154")) +
  scale_x_continuous(breaks=seq(0.4,0.7, by=0.1),
                     labels=c("Neonate*", "Piglet", "Subadult", "Adult"),
                     expand = expansion(mult = 0.15, add = 0)) +
  ylab("Survival probability") + xlab("Stage") +
  facet_grid(.~density) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position.inside=c(0.01, 0.98),
        legend.justification = c(0,1),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_markdown(size=12),
        legend.text=element_text(size=11),
        legend.title=element_blank())


## Fecundity
pigs_f <- pigs_f |>
  # replace 'Realized' with 'Implied'
  filter(source=="Implied") |>
  arrange(stage) |>
  group_by(density, stage) |>
  reframe(fec=sum(fec),
          f10pct=sum(f10pct),
          f90pct=sum(f90pct)) |>
  mutate(x=case_when(stage=="Yearling" ~ 0.1,
                     stage=="Adult" ~ 0.2))

pigs_f$stage <- factor(pigs_f$stage, levels=c("Yearling", "Adult"), labels=c("Subadult", "Adult\n(1+ years)"))
pigs_f$density <- factor(pigs_f$density, levels=c("8 pigs km<sup>-2</sup>","27 pigs km<sup>-2</sup>"))


p_rates_r$stage <- factor(p_rates_r$stage, levels=c("Subadult", "Adult"), labels=c("Subadult", "Adult\n(1+ years)"))

pf_plot <- ggplot(pigs_f) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,10)) +
  geom_point(data=p_rates_r, aes(x=stage, y=fec), shape=21, color="gray10", fill="gray",  size=3, alpha=0.45) +
  geom_segment(aes(x=stage, y=f10pct, yend=f90pct), linewidth=1) +
  geom_point(aes(x=stage, y=fec), size=3) +
  facet_grid(.~density) +
  ylab("Fecundity (piglets/sow)") + xlab("Stage") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_markdown(size=12),
        legend.text=element_text(size=11),
        legend.title=element_blank())



(pigs_pop | pf_plot) / ps_plot + plot_annotation(tag_levels='a', tag_prefix = "(", tag_suffix=")")


ggsave("figs/pigs_implied_rates_6.svg", height=8, width=10, units="in") #Saving 10 x 8 in image

