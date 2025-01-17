library(ggplot2)


Year<-1:10

###################################################################
## Deterministic model
#Matrix Calculations

deer.array <- matrix(0,nrow=12, ncol=length(Year))
deer.matrix<-matrix(0,12,12)

Fecundity <- c(0, 0.65, 0.76, 0.76, 0.76, 0.76)
N.0<-c(2000,2000,2000,1000,1000,1000,2000,2000,2000,1000,1000,1000)

# Females
deer.matrix[c(1,7),1:6]<-Fecundity
deer.matrix[2,1] <- 0.52
deer.matrix[3,2] <- 0.93
deer.matrix[4,3] <- 0.84
deer.matrix[5,4] <- 0.84
deer.matrix[6,5] <- 0.84
deer.matrix[6,6] <- 0.84
# Males
deer.matrix[8,7] <- 0.52
deer.matrix[9,8] <- 0.82
deer.matrix[10,9] <- 0.63
deer.matrix[11,10] <- 0.53
deer.matrix[12,11] <- 0.44
deer.matrix[12,12] <- 0.49

deer.array[,1] <- N.0

for (y in 2:length(Year)){
  deer.array[,y]<-deer.matrix%*%deer.array[,y-1] #Make sure to multiply matrix x vector (not vice versa)
}
N.median<-apply(deer.array,2,sum)

results<-data.frame(Year,N.median)

#Plot population with
ggplot(results) +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

# Calculate lambda
lambda <- numeric()
for(y in 1:nrow(results)){
  lambda[y] <- results$N.median[y+1]/results$N.median[y]
}
mean(lambda, na.rm=TRUE)




################################################################################################
##Stochastistic model

Sims<-1000

#Matrix Calculations
deer.array<-array(0,dim=c(6,length(Year),Sims))
deer.matrix<-matrix(0,6,6)

Fecundity<-c(0, 1.14, 1.13, 0.95, 0.79, 0.88)
N.0<-c(4000,3000,2000,1000,1000,1000)

deer.matrix[1,]<-Fecundity
deer.array[,1,]<-N.0

for(i in 1:Sims){

  deer.matrix[2,1] <- rnorm(1, 0.60, 0.05)
  deer.matrix[3,2] <- rnorm(1, 0.82, 0.05)
  deer.matrix[4,3] <- rnorm(1, 0.63, 0.05)
  deer.matrix[5,4] <- rnorm(1, 0.53, 0.05)
  deer.matrix[6,5] <- rnorm(1, 0.44, 0.05)
  deer.matrix[6,6] <- rnorm(1, 0.49, 0.05)
  
  for (y in 2:length(Year)){
    deer.array[,y,i]<-deer.matrix%*%deer.array[,y-1,i] #Make sure to multiply matrix x vector (not vice versa)
  }
}
N.median<-apply(apply(deer.array,c(2,3),sum),1,median)
N.20pct<-apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.20)
N.80pct<-apply(apply(deer.array,c(2,3),sum),1,quantile, probs = 0.80)

results<-data.frame(Year,N.median,N.20pct,N.80pct)

#Plot population projection
ggplot(results) +
  geom_ribbon(aes(x=Year,ymin=N.20pct, ymax=N.80pct), alpha=.2,fill="purple") +
  geom_line(aes(x=Year, y=N.median),colour="purple",alpha=1,linewidth=1.25) +
  scale_y_continuous(name="Abundance (females)")+
  theme_bw() +
  theme(text = element_text(size=16),panel.border = element_blank(), axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
