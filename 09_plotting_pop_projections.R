
library(tidyverse)
library(ggtext)
library(patchwork)

## Read results ####

# Deer
deer <- read_csv("results/deer_population_projections.csv")
deer_s <- read_csv("results/deer_survival.csv") 
deer_f <- read_csv("results/deer_fecundity.csv")

# Pigs
pigs <- read_csv("results/pigs_population_projections.csv")
pigs_s <- read_csv("results/pigs_survival.csv")
pigs_f <- read_csv("results/pig_fecundity.csv")

## Plot ####

### Deer ####
## Population projection
deer_pop <- ggplot(deer) +
  coord_cartesian(ylim=c(0,3000000)) +
  # geom_hline(yintercept=1610000, color="gray", linewidth=1) +
  geom_hline(yintercept=1347232, linetype=2, color="#ed7953", linewidth=1) +
  geom_hline(yintercept=2090532, linetype=2, color="#9c179e", linewidth=1) +
  geom_ribbon(aes(x=Year, ymin=N.10pct, ymax=N.90pct, fill=density, group=density), alpha=.2) +
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
  mutate(x=case_when(stage=="Fawn" ~ 0.2,
                     stage=="Yearling" ~ 0.3,
                     stage=="2-year-old" ~ 0.4,
                     stage=="3-year-old" ~ 0.5,
                     stage=="4-year-old" ~ 0.6,
                     stage=="5+ years-old" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01))

# Show only 1 density
deer_s <- deer_s |>
  filter(density=="14.3 deer km<sup>-2</sup>")


ds_plot <- ggplot(deer_s) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=x, y=surv.10pct, yend=surv.90pct, color=sex), linewidth=1) +
  geom_point(aes(x=x, y=survival, shape=source, color=sex), size=3) +
  scale_color_manual(values=c("#21918c", "#440154")) +
  scale_shape_manual(values=c(17, 16)) +
  guides(
    color=guide_legend(position="inside", title="Sex"),
    shape=guide_legend(position="inside", title="Source")
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

## Fecundity

# Filter to 14.3 deer
deer_f <- deer_f |>
  # replace 'Realized' with 'Implied'
  mutate(source = if_else(source=="Realized", "Implied", source)) |>
  # drop 22 deer/km2
  filter(density=="14.3 deer km<sup>-2</sup>")


df_plot <- ggplot(deer_f) +
  coord_cartesian(xlim=c(0.7, 2), ylim=c(0,1.4)) +
  geom_segment(aes(x=x, y=f10pct, yend=f90pct, color=offspring_sex), linewidth=1) +
  geom_point(aes(x=x, y=fec, group=source, color=offspring_sex, shape=source), size=3) +
  scale_color_manual(values=c("#21918c", "#440154")) +
  scale_shape_manual(values=c(17, 16)) +
  scale_x_continuous(breaks=c(1,1.7), 
                     labels=c("Yearling", "Adult"), 
                     expand = expansion(mult = 0.15, add = 0)) +
  guides(
    color=guide_legend(position="inside", title="Offspring Sex"),
    shape=guide_legend(position="inside", title="Source")
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


deer_plots <- deer_pop + ds_plot + df_plot + plot_layout(widths = c(1.2, 1.6, 1))


### Pigs ####

## Population projection
pigs_pop <- ggplot(pigs) +
  coord_cartesian(ylim=c(0,2000000)) +
  # geom_hline(yintercept=1000000, color="gray", linewidth=1) +
  geom_hline(yintercept=367822, linetype=2, color="#9c179e", linewidth=1) +
  geom_hline(yintercept=1239278, linetype=2, color="#ed7953", linewidth=1) +
  geom_ribbon(aes(x=Year, ymin=N.10pct, ymax=N.90pct, fill=density, group=density), alpha=.2) +
  geom_line(aes(x=Year, y=N.median, color=density, group=density), alpha=1,linewidth=1) +
  scale_y_continuous(name="Abundance (in millions)", 
                     breaks=c(0,500000, 1000000, 1500000, 2000000),
                     labels=c(0, 0.5, 1, 1.5, 2)) +
  scale_x_continuous(breaks=c(10,20,30,40), labels=c(5,10,15,20)) +
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
  mutate(x=case_when(stage=="Piglet" ~ 0.5,
                     stage=="Yearling" ~ 0.6,
                     stage=="Adult" ~ 0.7)) |>
  mutate(x = if_else(sex=="Female", x - 0.01, x + 0.01))

# Show only 1 density
pigs_s <- pigs_s |>
  filter(density=="27 pigs km<sup>-2</sup>") |>
  # replace 'Realized' with 'Implied'
  mutate(source = if_else(source=="Realized", "Implied", source))


ps_plot <- ggplot(pigs_s) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=x, y=surv.10pct, yend=surv.90pct, color=sex), linewidth=1) +
  geom_point(aes(x=x, y=survival, color=sex, shape=source), size=3) +
  scale_color_manual(values=c("#21918c", "#440154")) +
  scale_shape_manual(values=c(17, 16)) +
  guides(
    color=guide_legend(position="inside", title="Source"),
    shape=guide_legend(position="inside", title="Sex")
  ) +
  scale_x_continuous(breaks=seq(0.5,0.7, by=0.1), 
                     labels=c("Piglet", "Yearling", "Adult"), 
                     expand = expansion(mult = 0.15, add = 0)) +
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

# Filter to 14.3 deer
pigs_f <- pigs_f |>
  # replace 'Realized' with 'Implied'
  mutate(source = if_else(source=="Realized", "Implied", source)) |>
  mutate(x=case_when(stage=="Piglet" ~ 0.5,
                     stage=="Yearling + Adult" ~ 0.6)) |>
  mutate(x = if_else(offspring_sex=="Female", x - 0.01, x + 0.01))



pf_plot <- ggplot(pigs_f) +
  coord_cartesian(ylim=c(0,2.5)) +
  geom_segment(aes(x=x, y=f10pct, yend=f90pct, color=offspring_sex), linewidth=1) +
  geom_point(aes(x=x, y=fec, group=source, color=offspring_sex, shape=source), size=3) +
  scale_color_manual(values=c("#21918c", "#440154")) +
  scale_x_continuous(breaks=seq(0.5,0.6, by=0.1), 
                     labels=c("Piglet", "Yearling + Adult"), 
                     expand = expansion(mult = 0.30, add = 0)) +
  guides(
    color=guide_legend(position="inside", title="Offspring Sex"),
    shape=guide_legend(position="inside", title="Source")
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


pig_plots <- pigs_pop + ps_plot + pf_plot 


#### Plot them all 


deer_plots / pig_plots + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("figs/figure5_pop_sim.pdf")
# Saving 12.4 x 7.62 in image


deer_s |> filter(source=="Implied")

pigs_s |> filter(source=="Implied")
