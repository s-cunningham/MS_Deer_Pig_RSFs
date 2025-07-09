library(tidyverse)

# Create data frame for suitability and cells
hab <- data.frame(species=rep(c("deer", "pigs"), each=3),
                  type=rep(c("core","marginal","verymarginal"), 2),
                  ncells=c(5436507,5409943,3665703,6346538,6340520,1906219),
                  totalHS=c(5301442,10316653,13039984,4484418,8434372,9289481),
                  avgHS=c(0.975,0.927,0.743,0.702,0.619,0.515))

# Densities

hab$totalHS[hab$species=="deer" & hab$type=="verymarginal"] * 30 * 0.0081
