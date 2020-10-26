OmnibusF<-function(Model,chi1,DF1,N1,chi2,DF2,N2) {
  if (Model==1) {
  Probchi1=1-pchisq(chi1,DF1) #probability of chi-square statistic
  ncp1=chi1-DF1 #noncentrality parameter
  {if (ncp1<=0) {(ncp1=0)}} #noncentrality parameter constrained to be zero if negative
  {if (Model<=0|chi1<=0|DF1<=0|N1<=0)stop("Zeros and negative numbers are not parameters of this function")}
  CVchi1=qchisq(0.95,DF1) # Critical value of chi-square statistic for given DF
  powerchi1=1-pchisq(CVchi1,DF1,ncp1) #power of chi-square test
  F_Ratio1=chi1/DF1 #F-Ratio of model based on chi-square estimate
  CVf1<-qf(.95, DF1, N1-1) #critical value of F-dist
  pfvalue1<-pf(F_Ratio1, DF1, N1-1, lower.tail=FALSE) #probability value of F-statistic
  powerf1=1-pf(CVf1, DF1, N1-1, ncp=ncp1) # Power of F-test
  NdfRatio1=N1/DF1 # N to DF Ratio on what constitutes small sized sample data
  F1<-rf(N1, DF1, N1-1) #simulated data on F distribution based on DF1 and DF2
  DensF1<-df(F1, DF1, N1-1) #Density of simulated data
  X1<-rchisq(N1, DF1) #simulated data of chi-square distribution
  Denschi1<-dchisq(X1, DF1, log=FALSE) #Density of chi-square distribution
  cat("\n  OmnibusF Function (Version 1.0)", "\n")
  cat("\nThe OmnibusF Function was developed to provide an index for the evaluation of model fit", "\n") 
  cat("in Confirmatory Factor Analysis (CFA) using small sample sizes as it is well known that the", "\n") 
  cat("Maximum Likelihood Discrepancy measure (Chi-square test) is not distributed as a chi-square", "\n") 
  cat("statistic but rather as an F-Test (McNeish, 2018; Steiger, 2004; Thomas & Rao, 1987; Yuan & Bentler, 1999).", "\n") 
  cat("\nThe R-Function is as Follows:", "\n") 
  cat("\nOmnibusF<-function(Model,chi1,DF1,N1,chi2,DF2,N2)","\n")
  cat("\n         ############ Results ############", "\n")  
  cat("The N to DF Ratio is =", NdfRatio1, "If Ratio<=3 you are strongly advised to use the F-test", "\n")
  cat("The Chi-square Statistic of model =", chi1, ", and the p-value of the Chi-square test =", Probchi1, "\n")
  cat("The power of the Chi-Square Statistic =", powerchi1, "\n")
  cat("The F-Ratio statistic of the model =", F_Ratio1,", and the p-value of the F-Ratio Statistic =", pfvalue1, "\n")
  cat("The power of the F-test =", powerf1, "\n")
  cat(" ", "\n") 
  par(mfrow=c(1,2))
  figTL<-plot(X1,Denschi1, main="Density of Chi-square Test", 
            xlab="Chi-Square Estimates", ylab="Density of Chi-Square", pch = 21, cex=2, col="black", bg="grey", lwd=2)
            abline(v=c(CVchi1, chi1), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
            text(c(CVchi1, chi1), max(Denschi1), c("Crit-Chi","Obs-Chi"), pos=2, srt=90)         
  figTR<-plot(F1,DensF1, main="Density of F-Test",
            xlab="F-Test Estimates", ylab="Density of F-test", pch = 21, cex=2, col="black", bg="grey", lwd=2) 
            abline(v=c(CVf1, F_Ratio1), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
            text(c(CVf1, F_Ratio1), max(DensF1), c("Crit-F","Obs-F"), pos=2, srt=90)
  options(max.print = 999999)
}
if (Model==2) {
  Probchi1=1-pchisq(chi1,DF1) #probability of chi-square statistic 1
  ncp1=chi1-DF1 #noncentrality parameter
  {if (ncp1<=0) {(ncp1=0)}}
  {if (N1!=N2) {N1=(N1+N2)}}
  {if (DF1==DF2) stop("You cannot compare models with equal DF. Use Model 2 for nested models only.")}
  {if (Model<=0|chi1<=0|DF1<=0|N1<=0|chi2<=0|DF2<=0|N2<=0)stop("Zeros and negative numbers are not parameters of this function")}
  CVchi1=qchisq(0.95,DF1) # Critical value of chi-square statistic for given DF
  powerchi1=1-pchisq(CVchi1,DF1,ncp1)
  F_Ratio1=chi1/DF1 #F-Ratio of model 1
  CVf1<-qf(.95, DF1, N1-1) #critical value of F-dist
  pfvalue1<-pf(F_Ratio1, DF1, N1-1, lower.tail=FALSE)
  powerf1=1-pf(CVf1, DF1, N1-1, ncp1) # Power of F-test of Model 1
  NdfRatio1=N1/DF1 # N to DF Ratio on what constitutes small sized sample data
  F1<-rf(N1, DF1, N1-1)
  DensF1<-df(F1, DF1, N1-1)
  X1<-rchisq(N1, DF1)
  Denschi1<-dchisq(X1, DF1, log=FALSE)
  Probchi2=1-pchisq(chi2,DF2) #probability of chi-square statistic 2
  ncp2=chi2-DF2 #noncentrality parameter
  {if (ncp2<=0) {(ncp2=0)}}
  CVchi2=qchisq(0.95,DF2) # Critical value of chi-square statistic for given DF
  powerchi2=1-pchisq(CVchi2,DF2,ncp2)
  F_Ratio2=chi2/DF2 #F-Ratio model 2
  CVf2<-qf(.95, DF2, N2-1) #critical value of F-dist  
  pfvalue2<-pf(F_Ratio2, DF2, N2-1, lower.tail=FALSE)
  powerf2=1-pf(CVf2, DF2, N2-1, ncp2)
  NdfRatio2=N2/DF2 # N to DF Ratio on what constitutes small sized sample data
  F2<-rf(N2, DF2, N2-1)
  DensF2<-df(F2, DF2, N2-1)
  X2<-rchisq(N2, DF2)
  Denschi2<-dchisq(X2, DF2, log=FALSE)
  chiDIFF=abs(chi1-chi2)
  DFDIFF=abs(DF1-DF2)
  CVchi3=qchisq(0.95,DFDIFF) # Critical value of chi-square for DIFF DF
  COMPchi=1-pchisq(chiDIFF,DFDIFF) #probability of chi-square statistic 1
  CVf3<-qf(.95, DFDIFF, N1-1) #critical value of F-dist  
  f_DIFF=(chiDIFF/DFDIFF)
  COMPf<-pf(f_DIFF, DFDIFF, N1-1, lower.tail=FALSE)
  F3<-rf(N1, DFDIFF, N1-1)
  DensF3<-df(F3, DFDIFF, N1-1)
  X3<-rchisq(N1, DFDIFF)
  Denschi3<-dchisq(X3, DFDIFF)
  cat("\n  OmnibusF Function (Version 1.0)", "\n")
  cat("\nThe OmnibusF Function was developed to provide an index for the evaluation of model fit", "\n") 
  cat("in Confirmatory Factor Analysis (CFA) using small sample sizes as it is well known that the", "\n") 
  cat("Maximum Likelihood Discrepancy measure (Chi-square test) is not distributed as a chi-square", "\n") 
  cat("statistic but rather as an F-Test (McNeish, 2018; Steiger, 2004; Thomas & Rao, 1987; Yuan & Bentler, 1999).", "\n") 
  cat("\nThe R-Function is as Follows:", "\n") 
  cat("\nOmnibusF<-function(chi1,DF1,N1,chi2,DF2,N2)","\n")
  cat("\n         ############ Results ############", "\n")  
  cat("The N to DF Ratio of Model 1 =", NdfRatio1, "If Ratio<=3 you are strongly advised to use the F-test", "\n")
  cat("The Chi-square Statistic of model 1 =", chi1, ", and the p-value of the Chi-square test =", Probchi1, "\n")
  cat("The power of the Chi-Square Statistic of model 1 =", powerchi1, "\n")
  cat("The F-Ratio statistic of model 1 =", F_Ratio1,", and the p-value of the F-Ratio Statistic =", pfvalue1, "\n")
  cat("The N to DF Ratio of Model 2 =", NdfRatio2, "If Ratio<=3 you are strongly advised to use the F-test", "\n")
  cat("The Chi-square Statistic of model 2=", chi2, ", and the p-value of the Chi-square test=", Probchi2, "\n")
  cat("The power of the Chi-Square Statistic of model 2 =", powerchi2, "\n")
  cat("The F-Ratio statistic of model 2 =", F_Ratio2,", and the p-value of the F-Ratio Statistic =", pfvalue2, "\n")
  cat("The difference Chi-square Statistic =", chiDIFF, ", with Difference DF =", DFDIFF,", has a p-value =", COMPchi, "\n")  
  cat("The difference F-test =", f_DIFF, ", with Difference DF =", DFDIFF,", has a p-value =", COMPf, "\n")  
  cat(" ", "\n") 
  par(mfrow=c(3,2))
  figTL<-plot(F1,DensF1, main="Density of F-Test of Model 1",
              xlab="F-Test Estimates Model 1", ylab="Density of F-test", pch = 21, cex=2, col="black", 
              bg="grey", lwd=2) 
  abline(v=c(CVf1, F_Ratio1), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
  text(c(CVf1, F_Ratio1), max(DensF1), c("Crit-F","Obs-F"), pos=2, srt=90)
  figTR<-plot(X1,Denschi1, main="Density of Chi-square Test of Model 1", 
              xlab="Chi-Square Estimates Model 1", ylab="Density of Chi-Square", pch = 21, cex=2, col="black", 
              bg="grey", lwd=2)
  abline(v=c(CVchi1, chi1), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
  text(c(CVchi1, chi1), max(Denschi1), c("Crit-Chi","Obs-Chi"), pos=2, srt=90)         
  figML<-plot(F2,DensF2, main="Density of F-Test of Model 2", 
              xlab="F-Test Estimates Model 2", ylab="Density of F-test", pch = 21, cex=2, col="black", 
              bg="grey", lwd=2) 
  abline(v=c(CVf2, F_Ratio2), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
  text(c(CVf2, F_Ratio2), max(DensF2), c("Crit-F","Obs-F"), pos=2, srt=90)      
  figMR<-plot(X2,Denschi2, main="Density of Chi-square Test of Model 2", 
              xlab="Chi-Square Estimates Model 2", ylab="Density of Chi-Square", pch = 21, cex=2, col="black", 
              bg="grey", lwd=2)
  abline(v=c(CVchi2, chi2), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
  text(c(CVchi2, chi2), max(Denschi2), c("Crit-Chi","Obs-Chi"), pos=2, srt=90)            
  figBL<-plot(F3,DensF3, main="Density of F-Test of Model Comparison", 
              xlab="F-Test Estimates Model Comparison", ylab="Density of F-test", pch = 21, cex=2, col="black", 
              bg="grey", lwd=2) 
  abline(v=c(CVf3, f_DIFF), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
  text(c(CVf3, f_DIFF), max(DensF3), c("Crit-F","Obs-F"), pos=2, srt=90)
  figBR<-plot(X3,Denschi3, main="Density of Chi-square Test of Model Comparison", 
              xlab="Chi-Square Estimates of Model Comparison", ylab="Density of Chi-Square", pch = 21, cex=2, 
              col="black", bg="grey", lwd=2)
  abline(v=c(CVchi3, chiDIFF), col = c('black', 'black'), lty=c(2,1), lwd=c(2))
  text(c(CVchi3, chiDIFF), max(Denschi3), c("Crit-Chi","Obs-Chi"), pos=2, srt=90)  
  options(max.print = 999999)
}
}
