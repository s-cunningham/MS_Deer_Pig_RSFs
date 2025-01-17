library(tidyverse)


Year<-1:10

###################################################################
## Deterministic model
#Matrix Calculations

pig.array <- matrix(0,nrow=6, ncol=length(Year))
pig.matrix<-matrix(0,6,6)

Fecundity<-c(0.54, 2.24, 2.24)
N.0 <- c(10000,10000,10000,10000,10000,10000)

pig.matrix[1,1:3] <- Fecundity
pig.matrix[4,1:3] <- Fecundity
pig.matrix[2,1] <- 0.3
pig.matrix[3,2] <- 0.35
pig.matrix[3,3] <- 0.35
pig.matrix[5,4] <- 0.18
pig.matrix[6,5] <- 0.23
pig.matrix[6,6] <- 0.23

pig.array[,1] <- N.0

for (y in 2:length(Year)){
  pig.array[,y]<-pig.matrix%*%pig.array[,y-1] #Make sure to multiply matrix x vector (not vice versa)
}
N.median<-apply(pig.array,2,sum)

results<-data.frame(Year,N.median)

#Plot population with
ggplot(results) +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (males + females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

pig.array <- as.data.frame(pig.array)
names(pig.array) <- Year

pig.array$stage <- c("PigletF", "YearlingF", "AdultF", "PigletM", "YearlingM", "AdulM")

pig.array <- pig.array %>% select(stage, 1:max(Year)) %>% pivot_longer(2:(1+max(Year)), names_to="year", values_to="N") %>% as_tibble() %>%
                mutate(year=as.numeric(year))


ggplot(pig.array) +
  geom_line(aes(x=year, y=N, group=stage, color=stage))

# Calculate lambda
lambda <- numeric()
for(y in 1:nrow(results)){
  lambda[y] <- results$N.median[y+1]/results$N.median[y]
}
mean(lambda, na.rm=TRUE) ## arithmatic mean

# geometric mean
prod(lambda[1:9])


