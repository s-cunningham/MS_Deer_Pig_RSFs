library(tidyverse)
library(broom.mixed)
library(patchwork)

## Load RDS of model
deer_rsf <- readRDS("data/deer_rsf_model.RDS")
m1_deer <- broom.mixed::tidy(deer_rsf, effects="fixed", exponentiate=FALSE) %>% filter(term!="(Intercept)")
m1_deer$term <- factor(m1_deer$term, levels=rev(c("pctBottomland", "pctDecid", "pctCrops", "pctHerbaceous", "pctEvergreen", "pctWater")), 
                       labels=rev(c("Bottomland\nHardwoods", "Deciduous\nForest", "Crops", "Herbaceous", "Evergreen\nForest", "Water")))

pigs_rsf <- readRDS("data/pigs_rsf_model.RDS")
m1_pigs <- broom.mixed::tidy(pigs_rsf, effects="fixed", exponentiate=FALSE) %>% filter(term!="(Intercept)")
m1_pigs$term <- factor(m1_pigs$term, levels=rev(c("bottomlandhw","foodcrops",  "uplandforest", "streams", "roads")), 
                       labels=rev(c("Bottomland\nHardwoods", "Palatable\nCrops", "Upland\nForest", "Streams", "Roads")))

## Create sequence for gradient
deer <- m1_deer %>% select(term, estimate) %>%
  group_by(term) %>%
  mutate(estimate = list(seq(0, estimate, length.out = 100))) %>% # create sequence
  unnest(estimate)

pigs <- m1_pigs %>% select(term, estimate) %>%
  group_by(term) %>%
  mutate(estimate = list(seq(0, estimate, length.out = 100))) %>% # create sequence
  unnest(estimate)


deer_coef <- ggplot() +
  # coord_flip() +
  geom_hline(yintercept=0, linewidth=0.2) +
  geom_line(data=deer, aes(x=term, y=estimate, colour=estimate), linewidth = 15) +
  scale_color_viridis_c(option = "plasma", direction=-1, guide="none") +
  ylab("Selection coefficient") + xlab("Percentage of") +
  geom_errorbar(data=m1_deer, aes(x=term, ymin=estimate-std.error, ymax=estimate+std.error), width=.1, color="gray60") +
  theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title=element_text(size=12),
        axis.text=element_text(size=11,angle = 90, vjust = 0.5, hjust=1))

pigs_coef <- ggplot() +
  # coord_flip() +
  geom_hline(yintercept=0, linewidth=0.2) +
  geom_line(data=pigs, aes(x=term, y=estimate, colour=estimate), linewidth = 15) +
  scale_color_viridis_c(option = "plasma", direction=1, end=0.5, guide="none") + #
  ylab("Selection coefficient") + xlab("Distance to") +
  geom_errorbar(data=m1_pigs, aes(x=term, ymin=estimate-std.error, ymax=estimate+std.error), width=.1, color="gray60") +
  theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text=element_text(size=11,angle = 90, vjust = 0.5, hjust=1))

layout <- "
AAAAAA
AAAAAA
AAAAAA
AAAAAA
BBBBBB
"

(deer_rsf + pigs_rsf) / (deer_coef + pigs_coef) +
  plot_annotation(tag_levels = 'a', tag_prefix="(", tag_suffix=")") + plot_layout(design = layout)
  # plot_layout(heights = c(8, 1))

ggsave("figs/figure_3_combined_730x1000.svg", width=7.7, height=10)
