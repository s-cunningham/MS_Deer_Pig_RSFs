library(tidyverse)
library(patchwork)

# Read data
deer <- read_csv("output/deer_functional.csv") %>%
  # add column for species
  mutate(species = "White-tailed Deer")
pigs <- read_csv("output/pigs_functional.csv") %>%
  # add column for species
  mutate(species = "Wild Pigs")


# Create factor levels for covariates
deer$value <- factor(deer$value, levels=c("Hardwoods", "Graminoids", "Shrubs", "Developed", "Crops","Water", "Water2"),
                     labels=c("Hardwoods", "Graminoids", "Shrubs", "Developed","Crops", "Water", expression(Water^2)))
pigs$value <- factor(pigs$value, levels=c("Hardwoods", "Graminoids", "Shrubs", "Developed","Water", "Water2"),
                     labels=c("Hardwoods", "Graminoids", "Shrubs", "Developed","Water", expression(Water^2)))

# Split into distance and pct cover data frames
# Filter by covariates that have water
deer_dist <- deer %>%
  filter(str_detect(value, "Water"))
pigs_dist <- pigs %>%
  filter(str_detect(value, "Water"))

# Filter by covariates that do not have water
deer_pct <- deer %>%
  filter(!str_detect(value, "Water"))
pigs_pct <- pigs%>%
  filter(!str_detect(value, "Water"))

#### Make plots ####
# White-tailed deer
dpct <- ggplot(deer_pct) +
  geom_hline(yintercept=0, linewidth=0.5, color="gray60") +
  geom_point(aes(x=freq, y=beta), alpha=0.5) +
  facet_wrap(vars(value), scales="free", ncol=5, strip.position="bottom", labeller=label_parsed) +
  labs(y = expression(beta)) +
  xlab("Percent landcover in home range") +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.y=element_text(angle = 0, vjust = 0.5, size=16),
        axis.text=element_text(size=10),
        strip.placement="outside",
        strip.background=element_rect(color=NA, fill=NA))

ddist <- ggplot(deer_dist) +
  geom_hline(yintercept=0, linewidth=0.5, color="gray60") +
  geom_point(aes(x=freq, y=beta), alpha=0.5) +
  facet_wrap(vars(value), scales="free", ncol=2, strip.position="bottom", labeller=label_parsed) +
  xlab("Average distance") +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.y=element_blank(),
        axis.text=element_text(size=10),
        strip.placement="outside",
        strip.background=element_rect(color=NA, fill=NA))


# Set up layout
layout <- "
AAAAABB
"
# Combine figures
deer_plt <- dpct + ddist + plot_layout(design = layout)


# Wild pigs
ppct <- ggplot(pigs_pct) +
  geom_hline(yintercept=0, linewidth=0.5, color="gray60") +
  geom_point(aes(x=freq, y=beta), alpha=0.5) +
  facet_wrap(vars(value), scales="free", ncol=4, strip.position="bottom", labeller=label_parsed) +
  labs(y = expression(beta)) +
  xlab("Percent landcover in home range") +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.y=element_text(angle = 0, vjust = 0.5, size=16),
        axis.text=element_text(size=10),
        strip.placement="outside",
        strip.background=element_rect(color=NA, fill=NA))


pdist <- ggplot(pigs_dist) +
  geom_hline(yintercept=0, linewidth=0.5, color="gray60") +
  geom_point(aes(x=freq, y=beta), alpha=0.5) +
  facet_wrap(vars(value), scales="free", ncol=2, strip.position="bottom", labeller=label_parsed) +
  xlab("Average distance") +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5),
        axis.title.y=element_blank(),
        axis.text=element_text(size=10),
        strip.placement="outside",
        strip.background=element_rect(color=NA, fill=NA))

# Set up layout
layout <- "
AAAABCC
"
# Combine figures
pigs_plt <- ppct + plot_spacer() + pdist + plot_layout(design = layout)


deer_plt / pigs_plt

# save figure
ggsave("figs/functional_response_plot.svg")
# Saving 13.1 x 5.89 in image
