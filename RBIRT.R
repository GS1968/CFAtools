RBIRT <- function(model,Nopt,respv,ip,rbe,est) {
J <- length(ip[,2])
I<-dim(respv)[1]
Theta<-rep(NA,I)
se<-rep(NA,I)
p.mean<-0
p.sd<-1
ip<- as.matrix(ip)
rbe<-as.matrix(rbe)
for (i in 1:I){
x <- sum(respv[i,])
if (x == 0) {
Theta[i] <- -log(2 * J) }
if (x == J) {
Theta[i] <- log(2 * J) }
if (x == 0 | x == J) {
d2 <- 0.0
for (j in 1:J) {
z<-ip[j,1] * (Theta[i] - ip[j,2])
ps1 <- (rbe[i,1])*(1 / (1 + exp(-z)))
ps0<- 1-ps1
d2 <- d2 - ip[j,1]^2 * ps1 * ps0}
se[i] <- 1 / sqrt(-d2)}
if (x != 0 & x != J) {
Theta[i] <- log(x / (J - x))
S <- 10
ccrit <- 0.0001
for (s in 1:S) {
d1 <- 0.0
d2 <- 0.0
for (j in 1:J) {
z<-ip[j,1] * (Theta[i] - ip[j,2])
ps1 <- (rbe[i,1])*(1 / (1 + exp(-z)))
ps0<- 1-ps1
if (model == 1 ) {
d1 <- d1 + ip[j,1] * (respv[i,j] - ps1)
d2 <- d2 - ip[j,1]^2 * ps1 * ps0
}
if (model == 2 ) {
d1 <- d1 + ip[j,1] * (respv[i,j] - ps1)-(p.mean)/p.sd^2
d2 <- d2 - ip[j,1]^2 * ps1 * ps0-1/p.sd^2
}
}
Theta[i] <- Theta[i] - d1/d2
if (!is.na(abs(Theta[i]-d1/d2)) < ccrit | s == S) {
se[i] <- 1 / sqrt(-d2)
if (ps1==0 & est==1) {(Theta[i]=log((1/Nopt)/(1-(1/Nopt))))}
if (ps1==0 & est==2) {
Theta[i]=(th[i,]*ip[j,3])+log((1/Nopt)/(1-(1/Nopt)))
break
}
}
}
}
}
Results<-list(Theta,se)
return(Results)
}