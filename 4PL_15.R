### R Function for Person Guessing and Person Carelessness Estimation
irt4plpf <- function(theta, c, d, k) {
I<-dim(theta)[1]
G<-rep(NA, I)
G<-exp(-c+(theta*a))/(1+exp(-c+(theta*a))) #Person based ability-based guessing
C<-rep(NA, I)
C<-exp(-d+(theta*a))/(1+exp(-d+(theta*a))) #Person upper asymptote
D<-rep(NA, I)
D<-1-C #Person based ability based pseudo-carelessness
RG<-cor(G, theta) #Correlation between theta and person guessing
RC<-cor(C, theta) #Correlation between theta and person upper asymptote
RD<-cor(D, theta) #Correlation between theta and person carelessness
Results<-list(G[,1], C[,1], D[,1])
mpg=mean(G[,1]) #Mean of Person Guessing
lG<-quantile(G[,1], c(0.025)) #lower 2.5% of the Guessing Distribution
hG<-quantile(G[,1], c(0.975)) #upper 2.5% of the Guessing Distribution
mpgsd=sd(G[,1]) #SD of Person Guessing
mpC=mean(C[,1]) #Mean of Person Upper Asymptote
lC<-quantile(C[,1], c(0.025)) #lower 2.5% of the Person Upper Asymptote Distribution
hC<-quantile(C[,1], c(0.975)) #upper 2.5% of the Person Upper Asymptote Distribution
mpcsd=sd(C[,1]) #SD of Person Upper Asymptote
mc=mean(D[,1])  #Mean of Person Carelessness
lD<-quantile(D[,1], c(0.025)) #lower 2.5% of the Carelessness Distribution
hD<-quantile(D[,1], c(0.975)) #upper 2.5% of the Carelessness Distribution
mcsd=sd(D[,1])  #SD of Person Carelessness
n=nrow(G)   #Samples will be of size equal to the number of rows of guessing
nbs<-10000  #Number of 10000 bootstrap samples
set.seed(12345)
boot.resG<-numeric(nbs) #beginning of bootstrapping Guessing Distribution
for (i in 1:nbs) {
  boot.G<-sample(G[,1],n,replace=TRUE)
  boot.resG[i]=mean(boot.G)}
  meanb.G=mean(boot.G)
  meanb.Gsd=sd(boot.G)
  lowG<-quantile(boot.G, c(0.025))
  highG<-quantile(boot.G, c(0.975))
set.seed(12345)
boot.resC<-numeric(nbs) #beginning of bootstrapping Guessing Distribution
for (i in 1:nbs) {
  boot.C<-sample(C[,1],n,replace=TRUE)
  boot.resC[i]=mean(boot.C)}
  meanb.C=mean(boot.C)
  meanb.Csd=sd(boot.C)
  lowC<-quantile(boot.C, c(0.025))
  highC<-quantile(boot.C, c(0.975))
set.seed(12345)
boot.resD<-numeric(nbs) #beginning of bootstrapping Carelessness Distribution
  for (i in 1:nbs) {
  boot.D<-sample(D[,1],n,replace=TRUE)
  boot.resD[i]=mean(boot.D)}
  meanb.D=mean(boot.D)
  meanb.Dsd=sd(boot.D)
  lowD<-quantile(boot.D, c(0.025))
  highD<-quantile(boot.D, c(0.975))  
cat("\n ######## Results of 4PLPF Function for Estimating Person Parameters ########", "\n") 
cat("The Mean of Person Guessing =", mpg, "\n")
cat("The SD of Person Guessing =", mpgsd, "\n")
cat("The 95% CIs of the Person Guessing Distribution are =", lG, " and =", hG, "\n")
cat("The Mean of Bootstrapped Person Guessing Distribution =", meanb.G, "\n")
cat("The Bootstrapped SD of Person Guessing Distribution =", meanb.Gsd, "\n")
cat("The 95% CIs of the Bootstrapped Person Guessing Distribution are =", lowG, " and =", highG, "\n")
cat("The Mean of Person Upper Asymptote Distribution =", mpC, "\n")
cat("The SD of Person Upper Asymptote Distribution =", mpcsd, "\n")
cat("The 95% CIs of the Person Upper Asymptote Distribution are =", lC, " and =", hC, "\n")
cat("The Mean of Bootstrapped Person Upper Asymptote Distribution =", meanb.C, "\n")
cat("The Bootstrapped SD of Person Upper Asymptote Distribution =", meanb.Csd, "\n")
cat("The 95% CIs of the Bootstrapped Person Asymptote Distribution are =", lowC, " and =", highC, "\n")
cat("The Mean of Person Carelessness =", mc, "\n")
cat("The SD of Person Carelessness =", mcsd, "\n")
cat("The 95% CIs of the Person Carelessness Distribution are =", lD, " and =", hD, "\n")
cat("The Mean of Bootstrapped Person Carelessness Distribution =", meanb.D, "\n")
cat("The Bootstrapped SD of Person Carelessness Distribution =", meanb.Dsd, "\n")
cat("The 95% CIs of the Bootstrapped Person Carelessness Distribution are =", lowD, " and =", highD, "\n")
cat("Pearson's Correlation for Person Theta and Pseuo-Guessing =", RG, "\n")
cat("Pearson's Correlation for Person Theta and Person Upper Asymptote =", RC, "\n")
cat("Pearson's Correlation for Person Theta and Pseudo-Carelessness =", RD, "\n")
par(mfrow=c(4,2))
a<-hist(G[,1], col="yellow", main="Histogram of Person Guessing", xlim=c(0, .5))
   abline(v=c(lG,hG,lowG,highG), col=c("black", "black", "red", "red"), lty=c(2), lwd=c(2))
b<-boxplot(G[,1], main="Boxplot of Person Guessing")
c<-hist (C[,1], col="grey", main="Histogram of Person Upper Asymptote", xlim=c(.5, 1))
   abline(v=c(lC,hC,lowC,highC), col=c("black", "black", "red", "red"), lty=c(2), lwd=c(2))
d<-boxplot(C[,1],  main="Boxplot of Person Upper Asymptote")
e<-hist (D[,1], col="lightblue", main="Histogram of Person Carelessness", xlim=c(0, .5))
   abline(v=c(lD,hD,lowD,highD), col=c("black", "black", "red", "red"), lty=c(2), lwd=c(2))
f<-boxplot(D[,1],  main="Boxplot of Person Carelessness")
g<-plot(G[,1],theta[,1], main="Scatterplot of Person Guessing and Theta", 
  col="blueviolet", xlab="Person Guessing", ylab="Person Theta")
  abline(lm(theta[,1]~G[,1]), col="red", lwd=3, lty=2)
h<-plot(D[,1],theta[,1], main="Scatterplot of Person Carelessness and Theta", 
  col="lightsalmon2", xlab="Person Carelessness", ylab="Person Theta")
  abline(lm(theta[,1]~D[,1]), col="black", lwd=3, lty=2)
names(Results)<-c("Person Guessing", "Person Upper Asymptote", "1-Person Upper Asymptote")
options("scipen=999", digits=4, max.print = 9999)
print(Results)
}