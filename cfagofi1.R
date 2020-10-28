cfagofi1<-function(latent, indicators, N, chi.square, DF, chi.squarenull, DFnull, scf, ML.t, BSp) { 
  N.minus.1 = N - 1
  if(latent<=0)stop("The number of latent variables cannot be equal to or smaller than zero")
  if(indicators<=0)stop("The number of indicators cannot be equal to or smaller than zero")
  if(N<=20)stop("The number of participants cannot be less than 20")
  if(chi.square<=0)stop("The Chi-square statistic cannot be equal to or smaller than zero")
  if(DF<=0)stop("Model based DF can not be equal to or smaller than zero")
  if(chi.squarenull<=0)stop("The null model chi-square cannot be equal to or smaller than zero")
  if(DFnull<=0)stop("Null Model DF cannot be equal to or smaller than zero")
  if(scf == FALSE | ML.t == FALSE | BSp == FALSE) (chi.square = chi.square) else {print("Warning: one of the scf, ML.t and/or BSp parameters is missing, regular ML Chi-square will be used instead")}
  if(scf>0 & ML.t>0) stop("Check your syntax, you have defined more than one of the last 3 parameters, where only one is needed")
  if(scf>0 & BSp>0) stop("Check your syntax, you have defined more than one of the last 3 parameters, where only one is needed")
  if(ML.t>0 & BSp>0) stop("Check your syntax, you have defined more than one of the last 3 parameters, where only one is needed")
  if(scf>0 & ML.t>0 & BSp>0) stop("Check your syntax, you have defined more than one of the last 3 parameters, where only one is needed")
  if(scf>0) {(chi.square = chi.square / scf)}
  if(scf<0)stop("Scaling correction factor cannot be less than zero")
  if(ML.t>=0) {(ML.t = chi.square)}
  if(ML.t<0)stop("ML.t cannot be negative")
  p = (indicators*(indicators + 1))/2 #Number of non-redundant elements of S
  q = p - DF #Number of parameters
  BSx2 = qchisq(1 - BSp,DF) #Bollen Stine chi-square
  BSscf = chi.square / BSx2 #Bollen Stine correction factor
  if(BSp>0) {(chi.square = chi.square/BSscf)}
  if(BSp<0)stop("A p-value cannot be negative")
  ncp = chi.square-DF #noncentrality parameter
  Ne1 = (N) - (2.381 + 0.361 * indicators + 0.003 * q) #Yuan's Ne1 sample size correction
  Ne2 = (N) - (2.299 + 0.365 * indicators + 0.038 * latent) #Yuan's Ne2 sample size correction
  Ne3 = (N) - (2.262 + 0.369 * indicators + 0.052 * latent - 0.002 * q) #Yuan's Ne3 sample size correction
  RMSEA = sqrt(max((chi.square - DF)/(DF * N.minus.1),0)) #RMSEA
  Bartletc1 = 1 - ((4 * latent + 2 * indicators + 5) / (6 * (N.minus.1))) #Bartlett's correction
  chi.squareb = chi.square * Bartletc1 #Bartlett's corrected Chi-square 
  Yuanc1 = 1 - ((2 * indicators + 2 * latent + 7) / (6 * (N.minus.1))) #Yuan's correction
  Yuanc1a = 1 - ((2 * latent + 2 * indicators + 7) / (6 * Ne1)) #Yuan's correction for Ne1
  Yuanc2a = 1 - ((2 * latent + 2 * indicators + 7) / (6 * Ne2)) #Yuan's correction for Ne2
  Yuanc3a = 1 - ((2 * latent + 2 * indicators + 7) / (6 * Ne3)) #Yuan's correction for Ne3
  chi.squarey = chi.square * Yuanc1 #Yuan's corrected chi-square
  chi.squareyNe1 = chi.square * Yuanc1a #Yuan's corrected chi-square for Ne1
  chi.squareyNe2 = chi.square * Yuanc2a #Yuan's corrected chi-square for Ne2
  chi.squareyNe3 = chi.square * Yuanc3a #Yuan's corrected chi-square for Ne3
  q1 = (sqrt(1 + 4 * indicators * (indicators + 1) - 8 * DF) - 1) / 2 #Swain's q correction
  Swainc1 = 1 - ((indicators * (2 * indicators^2 + 3 * indicators - 1) - q1 * (2 * q1^2 + 3 * q1 - 1)) / (12 * N * DF)) #Swain's correction
  chi.squaresw = chi.square * Swainc1 #Swain's corrected Chi-square
  RMSEAb = sqrt(max((chi.squareb - DF) / (DF * N.minus.1),0)) #RMSEA Bartlett
  RMSEAy = sqrt(max((chi.squarey - DF) / (DF * N.minus.1),0)) #RMSEA Yuan
  RMSEAyNe1 = sqrt(max((chi.squareyNe1 - DF) / (DF * Ne1),0)) #RMSEA Yuan Ne1
  RMSEAyNe2 = sqrt(max((chi.squareyNe2 - DF) / (DF * Ne2),0)) #RMSEA Yuan Ne2
  RMSEAyNe3 = sqrt(max((chi.squareyNe3 - DF) / (DF * Ne3),0)) #RMSEA Yuan Ne3
  RMSEAsw = sqrt(max((chi.squaresw - DF) / (DF * N.minus.1),0)) #RMSEA Swain
  FF  = (chi.square / N.minus.1) #Fitting Function Uncorrected
  FFnull = (chi.squarenull / N.minus.1) #Fitting Function of Null Model
  FFb = (chi.squareb / N.minus.1) #Fitting Function Bartlett's correction
  FFy = (chi.squarey / N.minus.1) #Fitting Function, Yuan's correction
  FFyNe1 = (chi.squareyNe1 / N.minus.1) #Fitting Function, Yuan's correction Ne1
  FFyNe2 = (chi.squareyNe2 / N.minus.1) #Fitting Function, Yuan's correction Ne2
  FFyNe3 = (chi.squareyNe3 / N.minus.1) #Fitting Function, Yuan's correction Ne3
  FFsw = (chi.squaresw / N.minus.1) #Fitting Function Swain's correction
  LHR = exp(chi.square / (-2 * (N.minus.1))) #LHR index Uncorrected
  LHRb = exp(chi.squareb / (-2 * (N.minus.1))) #LHR index Bartlett's correction
  LHRy = exp(chi.squarey / (-2 * (N.minus.1))) #LHR index Yuan's correction
  LHRyNe1 = exp(chi.squareyNe1 / (-2 * (N.minus.1))) #LHR index Yuan's correction Ne1
  LHRyNe2 = exp(chi.squareyNe2 / (-2 * (N.minus.1))) #LHR index Yuan's correction Ne2
  LHRyNe3 = exp(chi.squareyNe3 / (-2 * (N.minus.1))) #LHR index Yuan's correction Ne3
  LHRsw = exp(chi.squaresw / (-2 * (N.minus.1))) #LHR index Swain's correction
  d = max(chi.square - DF,0)/N.minus.1  #scaled non-centrality parameter uncorrected chi-square
  NCI = exp((-d) / 2) #MacDonald's Non-Central Fit index uncorrected  
  db = max(chi.squareb - DF,0) / N.minus.1  #scaled non-centrality parameter Bartlett's correction on chi-square
  NCIb = exp((-db) / 2) #MacDonald's Non-Central Fit index Bartlett's correction  
  dy = max(chi.squarey - DF,0) / N.minus.1  #scaled non-centrality parameter Yuan's correction on chi-square
  dyNe1 = max(chi.squareyNe1 - DF,0) / N.minus.1  #scaled non-centrality parameter Yuan's correction on chi-square Ne1
  dyNe2 = max(chi.squareyNe2 - DF,0) / N.minus.1  #scaled non-centrality parameter Yuan's correction on chi-square Ne2
  dyNe3 = max(chi.squareyNe3 - DF,0) / N.minus.1  #scaled non-centrality parameter Yuan's correction on chi-square ne3
  NCIy = exp((-dy) / 2) #MacDonald's Non-Central Fit index Yuan's correction  
  NCIyNe1 = exp((-dyNe1) / 2) #MacDonald's Non-Central Fit index Yuan's correction Ne1 
  NCIyNe2 = exp((-dyNe2) / 2) #MacDonald's Non-Central Fit index Yuan's correction Ne2 
  NCIyNe3 = exp((-dyNe3) / 2) #MacDonald's Non-Central Fit index Yuan's correction Ne3 
  dsw = max(chi.squaresw - DF,0) / N.minus.1  #scaled non-centrality parameter SWain's correction on chi-square
  NCIsw = exp((-dsw) / 2) #MacDonald's Non-Central Fit index Swain's correction  
  PRatio = DF/DFnull #PRatio
  NFI = 1 - (chi.square / chi.squarenull) #NFI Uncorrected
  NFIb = 1 - (chi.squareb / chi.squarenull) #NFI Bartlett
  NFIy = 1 - (chi.squarey / chi.squarenull) #NFI Yuan
  NFIyNe1 = 1 - (chi.squareyNe1 / chi.squarenull) #NFI Yuan Ne1
  NFIyNe2 = 1 - (chi.squareyNe2 / chi.squarenull) #NFI Yuan Ne2
  NFIyNe3 = 1 - (chi.squareyNe3 / chi.squarenull) #NFI Yuan Ne3
  NFIsw = 1 - (chi.squaresw / chi.squarenull) #NFI Swain
  PFI = PRatio * NFI #James-Mulaik-Brett Parsimonious Fit Index
  PFIb = PRatio * NFIb #James-Mulaik-Brett Parsimonious Fit Index Bartlett
  PFIy = PRatio * NFIy #James-Mulaik-Brett Parsimonious Fit Index Yuan
  PFIyNe1 = PRatio * NFIyNe1 #James-Mulaik-Brett Parsimonious Fit Index Yuan Ne1
  PFIyNe2 = PRatio * NFIyNe2 #James-Mulaik-Brett Parsimonious Fit Index Yuan Ne2
  PFIyNe3 = PRatio * NFIyNe3 #James-Mulaik-Brett Parsimonious Fit Index Yuan Ne3
  PFIsw = PRatio * NFIsw #James-Mulaik-Brett Parsimonious Fit Index Swain
  NNFI = ((chi.squarenull / DFnull) - (chi.square / DF)) / ((chi.squarenull / DFnull)-1) #NNFI Uncorrected
  NNFIb = ((chi.squarenull / DFnull) - (chi.squareb/DF)) / ((chi.squarenull / DFnull)-1) #NNFI Bartlett
  NNFIy = ((chi.squarenull / DFnull) - (chi.squarey/DF)) / ((chi.squarenull/DFnull)-1) #NNFI Yuan
  NNFIyNe1 = ((chi.squarenull / DFnull) - (chi.squareyNe1/DF)) / ((chi.squarenull / DFnull)-1) #NNFI Yuan ne1
  NNFIyNe2 = ((chi.squarenull / DFnull) - (chi.squareyNe2/DF)) / ((chi.squarenull / DFnull)-1) #NNFI Yuan Ne2
  NNFIyNe3 = ((chi.squarenull / DFnull) - (chi.squareyNe3/DF)) / ((chi.squarenull / DFnull)-1) #NNFI Yuan ne3
  NNFIsw = ((chi.squarenull / DFnull) - (chi.squaresw/DF)) / ((chi.squarenull/DFnull)-1) #NNFI Swain
  RFI = ((chi.squarenull / DFnull) - (chi.square/DF)) / (chi.squarenull / DFnull) #RFI Uncorrected
  RFIb = ((chi.squarenull / DFnull) - (chi.squareb/DF)) / (chi.squarenull / DFnull) #RFI Bartlett
  RFIy = ((chi.squarenull / DFnull) - (chi.squarey/DF)) / (chi.squarenull / DFnull) #RFI Yuan
  RFIyNe1 = ((chi.squarenull / DFnull) - (chi.squareyNe1 / DF))/(chi.squarenull / DFnull) #RFI Yuan Ne1
  RFIyNe2 = ((chi.squarenull / DFnull) - (chi.squareyNe2 / DF))/(chi.squarenull / DFnull) #RFI Yuan Ne2
  RFIyNe3 = ((chi.squarenull / DFnull) - (chi.squareyNe3 / DF))/(chi.squarenull / DFnull) #RFI Yuan Ne3
  RFIsw = ((chi.squarenull / DFnull) - (chi.squaresw / DF)) / (chi.squarenull / DFnull) #RFI Swain
  CFI = max(((chi.squarenull - DFnull) - (chi.square - DF)) / (chi.squarenull - DFnull),0) #CFI Uncorrected
  CFIb = max(((chi.squarenull - DFnull) - (chi.squareb - DF)) / (chi.squarenull - DFnull),0) #CFI Bartlett
  CFIy = max(((chi.squarenull - DFnull) - (chi.squarey - DF)) / (chi.squarenull - DFnull),0) #CFI Yuan
  CFIyNe1 = max(((chi.squarenull - DFnull) - (chi.squareyNe1 - DF))/(chi.squarenull - DFnull),0) #CFI Yuan Ne1
  CFIyNe2 = max(((chi.squarenull - DFnull) - (chi.squareyNe2 - DF))/(chi.squarenull - DFnull),0) #CFI Yuan Ne2
  CFIyNe3 = max(((chi.squarenull - DFnull) - (chi.squareyNe3 - DF))/(chi.squarenull - DFnull),0) #CFI Yuan Ne3
  CFIsw = max(((chi.squarenull - DFnull) - (chi.squaresw - DF)) / (chi.squarenull - DFnull),0) #CFI Swain
  Rho = ((chi.squarenull / DFnull) - (chi.square / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Uncorrected
  Rhob = ((chi.squarenull / DFnull) - (chi.squareb / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Bartlett
  Rhoy = ((chi.squarenull / DFnull) - (chi.squarey / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Yuan
  RhoyNe1 = ((chi.squarenull / DFnull) - (chi.squareyNe1 / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Yuan Ne1
  RhoyNe2 = ((chi.squarenull / DFnull) - (chi.squareyNe2 / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Yuan Ne2
  RhoyNe3 = ((chi.squarenull / DFnull) - (chi.squareyNe3 / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Yuan Ne3
  Rhosw = ((chi.squarenull / DFnull) - (chi.squaresw / DF)) / (chi.squarenull / DFnull) #Bollen's Rho Swain
  Delta = (chi.squarenull - chi.square) / (chi.squarenull - (DF / N.minus.1)) #Bollen's Delta
  Deltab = (chi.squarenull - chi.squareb) / (chi.squarenull - (DF / N.minus.1)) #Bollen's Delta Bartlett
  Deltay = (chi.squarenull - chi.squarey) / (chi.squarenull - (DF / N.minus.1)) #Bollen's Delta Yuan
  DeltayNe1 = (chi.squarenull - chi.squareyNe1) / (chi.squarenull - (DF / Ne1)) #Bollen's Delta Yuan Ne1
  DeltayNe2 = (chi.squarenull - chi.squareyNe2) / (chi.squarenull - (DF / Ne2)) #Bollen's Delta Yuan Ne2
  DeltayNe3 = (chi.squarenull - chi.squareyNe3) / (chi.squarenull - (DF / Ne3)) #Bollen's Delta Yuan Ne3
  Deltasw = (chi.squarenull - chi.squaresw) / (chi.squarenull - (DF / N.minus.1)) #Bollen's Delta Swain
  IFI = (chi.squarenull - chi.square) / (chi.squarenull - DF) #IFI Uncorrected
  IFIb = (chi.squarenull - chi.squareb) / (chi.squarenull - DF) #IFI Bartlett
  IFIy = (chi.squarenull - chi.squarey) / (chi.squarenull - DF) #IFI Yuan
  IFIyNe1 = (chi.squarenull - chi.squareyNe1) / (chi.squarenull - DF) #IFI Yuan Ne1
  IFIyNe2 = (chi.squarenull - chi.squareyNe2) / (chi.squarenull - DF) #IFI Yuan Ne2
  IFIyNe3 = (chi.squarenull - chi.squareyNe3) / (chi.squarenull - DF) #IFI Yuan Ne3
  IFIsw = (chi.squarenull - chi.squaresw) / (chi.squarenull - DF) #IFI Swain
  RNI = 1 - ((chi.square - DF) / (chi.squarenull - DFnull)) #RNI Uncorrected
  RNIb = 1 - ((chi.squareb - DF) / (chi.squarenull - DFnull)) #RNI Bartlett
  RNIy = 1 - ((chi.squarey - DF) / (chi.squarenull - DFnull)) #RNI Yuan
  RNIyNe1 = 1 - ((chi.squareyNe1 - DF) / (chi.squarenull - DFnull)) #RNI Yuan Ne1
  RNIyNe2 = 1 - ((chi.squareyNe2 - DF) / (chi.squarenull - DFnull)) #RNI Yuan Ne2
  RNIyNe3 = 1 - ((chi.squareyNe3 - DF) / (chi.squarenull - DFnull)) #RNI Yuan Ne3
  RNIsw = 1 - ((chi.squaresw - DF) / (chi.squarenull - DFnull)) #RNI Swain
  PRNI = PRatio * RNI #PRNI Uncorrected
  PRNIb = PRatio * RNIb #PRNI Bartlett
  PRNIy = PRatio * RNIy #PRNI Yuan
  PRNIyNe1 = PRatio * RNIyNe1 #PRNI Yuan Ne1
  PRNIyNe2 = PRatio * RNIyNe2 #PRNI Yuan Ne2
  PRNIyNe3 = PRatio * RNIyNe3 #PRNI Yuan Ne3
  PRNIsw = PRatio * RNIsw #PRNI Swain
  PRFI = PRatio * RFI #PRFI Uncorrected
  PRFIb = PRatio * RFIb #PRFI Bartlett
  PRFIy = PRatio * RFIy #PRFI Yuan
  PRFIyNe1 = PRatio * RFIyNe1 #PRFI Yuan Ne1
  PRFIyNe2 = PRatio * RFIyNe2 #PRFI Yuan Ne2
  PRFIyNe3 = PRatio * RFIyNe3 #PRFI Yuan Ne3
  PRFIsw = PRatio * RFIsw #PRFI Swain
  PIFI = PRatio * IFI #PIFI Uncorrected
  PIFIb = PRatio * IFIb #PIFI Bartlett
  PIFIy = PRatio * IFIy #PIFI Yuan
  PIFIyNe1 = PRatio * IFIyNe1 #PIFI Yuan Ne1
  PIFIyNe2 = PRatio * IFIyNe2 #PIFI Yuan Ne2
  PIFIyNe3 = PRatio * IFIyNe3 #PIFI Yuan Ne3
  PIFIsw = PRatio * IFIsw #PIFI Swain
  pchi = 1 - pchisq(chi.square, DF) #probability of uncorrected chi-square
  pchib = 1 - pchisq(chi.squareb, DF) #probability of chi-square using Bartlett's correction
  pchiy = 1 - pchisq(chi.squarey, DF) #probability of chi-square using Yuan's correction
  pchiyNe1 = 1 - pchisq(chi.squareyNe1, DF) #probability of chi-square using Yuan's correction Ne1
  pchiyNe2 = 1 - pchisq(chi.squareyNe2, DF) #probability of chi-square using Yuan's correction Ne2
  pchiyNe3 = 1 - pchisq(chi.squareyNe3, DF) #probability of chi-square using Yuan's correction Ne3
  pchisw = 1 - pchisq(chi.squaresw,DF) #probability of chi-square using Swain's correction
  GFI = indicators/(2 * (chi.square - DF) / (N.minus.1) + indicators) #Unbiased GFI Steiger uncorrected
  GFIb = indicators/(2 * (chi.squareb - DF) / (N.minus.1) + indicators) #Unbiased GFI Steiger Bartlett's correction
  GFIy = indicators/(2 * (chi.squarey - DF) / (N.minus.1) + indicators) #Unbiased GFI Steiger Yuan's correction
  GFIyNe1 = indicators / (indicators + 2 * chi.squareyNe1 / N.minus.1) #Unbiased GFI Steiger Yuan's correction Ne1
  GFIyNe2 = indicators / (indicators + 2 * chi.squareyNe2 / N.minus.1) #Unbiased GFI Steiger Yuan's correction Ne2
  GFIyNe3 = indicators / (indicators + 2 * chi.squareyNe3 / N.minus.1) #Unbiased GFI Steiger Yuan's correction Ne3
  GFIsw = indicators / (indicators + 2 * chi.squaresw / N.minus.1) #GFI Swain's correction
  AGFI = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFI) #AGFI uncorrected
  AGFIb = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFIb) #AGFI corrected Bartlett
  AGFIy = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFIy) #AGFI corrected Yuan
  AGFIyNe1 = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFIyNe1) #AGFI corrected Yuan Ne1
  AGFIyNe2 = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFIyNe2) #AGFI corrected Yuan Ne2
  AGFIyNe3 = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFIyNe3) #AGFI corrected Yuan Ne3
  AGFIsw = 1 - (indicators * (indicators + 1) / (2 * DF)) * (1 - GFIsw) #AGFI corrected Swain
  Gamma = indicators / (indicators + 2 * (chi.square - DF) / N.minus.1) #Gamma population index uncorrected
  Gammab = indicators / (indicators + 2 * (chi.squareb - DF) / N.minus.1) #Gamma population index, Bartlett's correction
  Gammay = indicators / (indicators + 2 * (chi.squarey - DF) / N.minus.1) #Gamma population index, Yuan's correction
  GammayNe1 = indicators / (indicators + 2 * (chi.squareyNe1 - DF) / Ne1) #Gamma population index, Yuan's correction Ne1
  GammayNe2 = indicators / (indicators + 2 * (chi.squareyNe2 - DF) / Ne2) #Gamma population index, Yuan's correction Ne2
  GammayNe3 = indicators / (indicators + 2 * (chi.squareyNe3 - DF) / Ne3) #Gamma population index, Yuan's correction Ne3
  Gammasw = indicators / (indicators + 2 * (chi.squaresw - DF) / N.minus.1) #Gamma population index, Swain's correction
  PCFI = PRatio * CFI #Parsimonious CFI uncorrected
  PCFIb = PRatio * CFIb #Parsimonious CFI corrected Bartlett
  PCFIy = PRatio * CFIy #Parsimonious CFI corrected Yuan
  PCFIyNe1 = PRatio * CFIyNe1 #Parsimonious CFI corrected Yuan Ne1
  PCFIyNe2 = PRatio * CFIyNe2 #Parsimonious CFI corrected Yuan Ne2
  PCFIyNe3 = PRatio * CFIyNe3 #Parsimonious CFI corrected Yuan Ne3
  PCFIsw = PRatio * CFIsw #Parsimonious CFI corrected Swain
  PGFI = PRatio * GFI #Parsimonious GFI Uncorrected
  PGFIb = PRatio * GFIb #Parsimonious GFI corrected Bartlett
  PGFIy = PRatio * GFIy #Parsimonious GFI corrected Yuan
  PGFIyNe1 = PRatio * GFIyNe1 #Parsimonious GFI corrected Yuan Ne1
  PGFIyNe2 = PRatio * GFIyNe2 #Parsimonious GFI corrected Yuan Ne2
  PGFIyNe3 = PRatio * GFIyNe3 #Parsimonious GFI corrected Yuan Ne3
  PGFIsw = PRatio * GFIsw #Parsimonious GFI corrected Swain
  out1 = chi.square #Chi-square value uncorrected
  out2 = pchi #Probability of Chi-square test, uncorrected test
  out3 = RMSEA #RMSEA Uncorrected
  out4 = FF #Fitting Function
  out5 = LHR #LHR
  out6 = chi.square / DF #Chi.square/DF Uncorrected
  out7 = NFI #NFI Uncorrected
  out8 = PFI #James-Mulaik-Brett Parsimonious Fit Index
  out9 = NNFI #NNFI Uncorrected
  out10 = CFI #CFI Uncorrected
  out11 = PCFI #Parsimonious CFI uncorrected
  out12 = Rho #Bollen's Rho Uncorrected
  out13 = Delta #Bollen's Delta
  out14 = RFI #RFI Uncorrected
  out15 = PRFI #PRFI Uncorrected
  out16 = IFI #IFI Uncorrected
  out17 = PIFI #Parsimonious IFI Uncorrected
  out18 = GFI #GFI Uncorrected
  out19 = AGFI #AGFI Uncorrected
  out20 = PGFI #Parsimonious GFI Uncorrected
  out21 = Gamma #Gamma population index
  out22 = RNI #RNI Uncorrected
  out23 = PRNI #PRNI Uncorrected
  out24 = NCI #NCI uncorrected
  out25 = chi.squareb #Chi-square value Bartlett's correction
  out26 = pchib #Probability of Chi-square test, Bartlett's correction
  out27 = RMSEAb #RMSEA Bartlett
  out28 = FFb #Fitting Function Bartlett correction
  out29 = LHRb #LHR Bartlett correcction
  out30 = chi.squareb/DF #Chi.square Bartlett/DF
  out31 = NFIb #NFI Bartlet
  out32 = PFIb #James-Mulaik-Brett Parsimonious Fit Index Bartlett
  out33 = NNFIb #NNFI Bartlett
  out34 = CFIb #CFI Bartlett
  out35 = PCFIb #Parsimonious CFI corrected Bartlett
  out36 = Rhob #Bollen's Rho Bartlett
  out37 = Deltab #Bollen's Delta Bartlett
  out38 = RFIb #RFI Bartlett
  out39 = PRFIb #PRFI corrected Bartlett
  out40 = IFIb #IFI Bartlett
  out41 = PIFIb #Parsimonious IFI corrected Bartlett
  out42 = GFIb #GFI Bartlett
  out43 = AGFIb #AGFI corrected Bartlett
  out44 = PGFIb #Parsimonious GFI corrected Bartlett
  out45 = Gammab #Gamma population index, Bartlett's correction
  out46 = RNIb #RNI Bartlett
  out47 = PRNIb #PRNI Bartlett
  out48 = NCIb #NCI Bartlett
  out49 = chi.squarey #Chi-square value Yuan's correction
  out50 = pchiy #Probability of Chi-square test, Yuan's correction
  out51 = RMSEAy #RMSEA Yuan
  out52 = FFy #Fitting Function Yuan's correction
  out53 = LHRy #LHR Yuan's correction
  out54 = chi.squarey/DF #Chi.square Yuan/DF
  out55 = NFIy #NFI Yuan
  out56 = PFIy #James-Mulaik-Brett Parsimonious Fit Index Yuan
  out57 = NNFIy #NNFI Yuan
  out58 = CFIy #CFI Yuan
  out59 = PCFIy #Parsimonious CFI corrected Yuan
  out60 = Rhoy #Bollen's Rho Yuan's Correction
  out61 = Deltay #Bollen's Delta Yuan's Correction
  out62 = RFIy #RFI Yuan
  out63 = PRFIy #PRFI corrected Yuan
  out64 = IFIy #IFI Yuan
  out65 = PIFIy #Parsimonious IFI corrected Yuan
  out66 = GFIy #GFI Yuan
  out67 = AGFIy #AGFI corrected Yuan
  out68 = PGFIy #Parsimonious GFI corrected Yuan
  out69 = Gammay #Gamma population index, Yuan's correction
  out70 = RNIy #RNI Yuan
  out71 = PRNIy #PRNI Yuan
  out72 = NCIy #NCI Yuan
  out73 = chi.squaresw #Chi-square value Swain's correction
  out74 = pchisw #Probability of Chi-square test Swain's correction
  out75 = RMSEAsw #RMSEA Swain
  out76 = FFsw #Fitting Function Swain's correction
  out77 = LHRsw #LHR Swain's correction
  out78 = chi.squaresw/DF #Chi.square Swain/DF
  out79 = NFIsw #NFI Swain
  out80 = PFIsw #James-Mulaik-Brett Parsimonious Fit Index Swain
  out81 = NNFIsw #NNFI Swain
  out82 = CFIsw #CFI Swain
  out83 = PCFIsw #Parsimonious CFI corrected Swain
  out84 = Rhosw #Bollen's Rho Swain
  out85 = Deltasw #Bollen's Delta Swain
  out86 = RFIsw #RFI Swain
  out87 = PRFIsw #PRFI corrected Swain
  out88 = IFIsw #IFI Swain
  out89 = PIFIsw #Parsimonious IFI corrected Swain
  out90 = GFIsw #GFI Swain
  out91 = AGFIsw #AGFI corrected Swain
  out92 = PGFIsw #Parsimonious GFI corrected Swain
  out93 = Gammasw #Gamma population index Swain's correction
  out94 = RNIsw #RNI Swain
  out95 = PRNIsw #PRNI Swain
  out96 = NCIsw #NCI Swain
  out97 = chi.squareyNe1 #Chi-square value Yuan's Ne1 correction
  out98 = pchiyNe1 #Probability of Chi-square test, Yuan's Ne1 correction
  out99 = RMSEAyNe1 #RMSEA Yuan Ne1 correction
  out100 = FFyNe1 #Fitting Function Yuan's Ne1 correction
  out101 = LHRyNe1 #LHR Yuan's Ne1 correction
  out102 = chi.squareyNe1/DF #Chi.square Yuan/DF Ne1 correction
  out103 = NFIyNe1 #NFI Yuan Ne1 correction
  out104 = PFIyNe1 #James-Mulaik-Brett Parsimonious Fit Index Yuan Ne1 correction
  out105 = NNFIyNe1 #NNFI Yuan Ne1 correction
  out106 = CFIyNe1 #CFI Yuan Ne1 correction
  out107 = PCFIyNe1 #Parsimonious CFI corrected Yuan Ne1 correction
  out108 = RhoyNe1 #Bollen's Rho Yuan's Correction Ne1 correction
  out109 = DeltayNe1 #Bollen's Delta Yuan's Correction Ne1 correction
  out110 = RFIyNe1 #RFI Yuan Ne1 correction
  out111 = PRFIyNe1 #PRFI corrected Yuan Ne1 correction
  out112 = IFIyNe1 #IFI Yuan Ne1 correction
  out113 = PIFIyNe1 #Parsimonious IFI corrected Yuan Ne1 correction
  out114 = GFIyNe1 #GFI Yuan Ne1 correction
  out115 = AGFIyNe1 #AGFI corrected Yuan Ne1 correction
  out116 = PGFIyNe1 #Parsimonious GFI corrected Yuan Ne1 correction
  out117 = GammayNe1 #Gamma population index, Yuan's Ne1 correction
  out118 = RNIyNe1 #RNI Yuan Ne1 correction
  out119 = PRNIyNe1 #PRNI Yuan Ne1 correction
  out120 = NCIyNe1 #NCI Yuan Ne1 correction
  out121 = chi.squareyNe2 #Chi-square value Yuan's Ne1 correction
  out122 = pchiyNe2 #Probability of Chi-square test, Yuan's Ne1 correction
  out123 = RMSEAyNe2 #RMSEA Yuan Ne1 correction
  out124 = FFyNe2 #Fitting Function Yuan's Ne1 correction
  out125 = LHRyNe2 #LHR Yuan's Ne1 correction
  out126 = chi.squareyNe2/DF #Chi.square Yuan/DF Ne2 correction
  out127 = NFIyNe2 #NFI Yuan Ne2 correction
  out128 = PFIyNe2 #James-Mulaik-Brett Parsimonious Fit Index Yuan Ne2 correction
  out129 = NNFIyNe2 #NNFI Yuan Ne2 correction
  out130 = CFIyNe2 #CFI Yuan Ne2 correction
  out131 = PCFIyNe2 #Parsimonious CFI corrected Yuan Ne2 correction
  out132 = RhoyNe2 #Bollen's Rho Yuan's Correction Ne2 correction
  out133 = DeltayNe2 #Bollen's Delta Yuan's Correction Ne2 correction
  out134 = RFIyNe2 #RFI Yuan Ne2 correction
  out135 = PRFIyNe2 #PRFI corrected Yuan Ne2 correction
  out136 = IFIyNe2 #IFI Yuan Ne2 correction
  out137 = PIFIyNe2 #Parsimonious IFI corrected Yuan Ne2 correction
  out138 = GFIyNe2 #GFI Yuan Ne2 correction
  out139 = AGFIyNe2 #AGFI corrected Yuan Ne2 correction
  out140 = PGFIyNe2 #Parsimonious GFI corrected Yuan Ne2 correction
  out141 = GammayNe2 #Gamma population index, Yuan's Ne2 correction
  out142 = RNIyNe2 #RNI Yuan Ne2 correction
  out143 = PRNIyNe2 #PRNI Yuan Ne2 correction
  out144 = NCIyNe2 #NCI Yuan Ne2 correction
  out145 = chi.squareyNe3 #Chi-square value Yuan's Ne3 correction
  out146 = pchiyNe3 #Probability of Chi-square test, Yuan's Ne3 correction
  out147 = RMSEAyNe3 #RMSEA Yuan Ne3 correction
  out148 = FFyNe3 #Fitting Function Yuan's Ne3 correction
  out149 = LHRyNe3 #LHR Yuan's Ne3 correction
  out150 = chi.squareyNe3/DF #Chi.square Yuan/DF Ne3 correction
  out151 = NFIyNe3 #NFI Yuan Ne3 correction
  out152 = PFIyNe3 #James-Mulaik-Brett Parsimonious Fit Index Yuan Ne3 correction
  out153 = NNFIyNe3 #NNFI Yuan Ne3 correction
  out154 = CFIyNe3 #CFI Yuan Ne3 correction
  out155 = PCFIyNe3 #Parsimonious CFI corrected Yuan Ne3 correction
  out156 = RhoyNe3 #Bollen's Rho Yuan's Correction Ne3 correction
  out157 = DeltayNe3 #Bollen's Delta Yuan's Correction Ne3 correction
  out158 = RFIyNe3 #RFI Yuan Ne3 correction
  out159 = PRFIyNe3 #PRFI corrected Yuan Ne3 correction
  out160 = IFIyNe3 #IFI Yuan Ne3 correction
  out161 = PIFIyNe3 #Parsimonious IFI corrected Yuan Ne3 correction
  out162 = GFIyNe3 #GFI Yuan Ne3 correction
  out163 = AGFIyNe3 #AGFI corrected Yuan Ne3 correction
  out164 = PGFIyNe3 #Parsimonious GFI corrected Yuan Ne3 correction
  out165 = GammayNe3 #Gamma population index, Yuan's Ne3 correction
  out166 = RNIyNe3 #RNI Yuan Ne3 correction
  out167 = PRNIyNe3 #PRNI Yuan Ne3 correction
  out168 = NCIyNe3 #NCI Yuan Ne3 correction
  cat("\n    CFAGOFI (Version 1.0)  ------- G.D. Sideridis & F. Al-Jaffari (2020)", "\n")
  cat("\nThe cfagofi Function was designed to correct the Chi-square statistic in CFA as it forms the", "\n") 
  cat("basis for the estimation of several descriptive Fit-Indices such as the CFI, TLI, AGFI, RMSEA,", "\n") 
  cat("and other. It is provided without warranty and for personal use only. Some theory behind the ", "\n") 
  cat("proposed corrections. Bartlett's correction was developed for Principal Components Analysis with", "\n") 
  cat("the goal of providing proper estimates for small sample sizes. His approach involves estimating a", "\n") 
  cat("corrective factor and multiplying it with the Likelihood Ratio Test. The correction involves the", "\n") 
  cat("number of latent variables k the number of observed variables p and the sample size n and is", "\n") 
  cat("estimated as follows: 1-(4k+2p+5)/6n. It is the most prominent among corrective procedures. The", "\n") 
  cat("Yuan correction is based on the proposition that exploratory and confirmatory factor models are", "\n") 
  cat("identical when having one latent variable. Thus in that instance both corrections would be identical.", "\n") 
  cat("However with more complex models Yuan's correction came to reduce the number of estimated parameters", "\n") 
  cat("as in comparison to a CFA model an EFA model estimates factor loadings for an item on all factors.", "\n") 
  cat("It is estimated as follows: Y=1-2k+2p+7/6n. Last Swain introduced several corrective procedures", "\n") 
  cat("with the most popular being estimated as a function of a Q parameter: q=sqrt((1+4p(p+1)-8d-1)/2", "\n") 
  cat("which is applied to the following corrective factor: 1-p(2p^+3p-1)-q(2q^+3q-1)/12dn.", "\n") 
  cat("\nThe CFAGOFI Function is as Follows:", "\n") 
  cat("\ncfagofi1<-function(latent,indicators,N,chi.square,DF,chi.squarenull,DFnull,scf,ML.t,BSp)", "\n")
  cat("\nDescription of cfagofi1 Function's Components:", "\n") 
  cat("Number of Latent Variables (latent), N =", latent, "\n") 
  cat("Number of Measured Variables(indicators), N =", indicators, "\n") 
  cat("Sample size, N-1 (N.minus.1) =", N.minus.1, "\n") 
  cat("Chi-Square Value of Tested Model (chi.Square) =", chi.square, "\n") 
  cat("Degrees of Freedom of Tested Model (DF) =", DF, "\n") 
  cat("Chi-Square Value of Independence Model (chi-squarenull) =", chi.squarenull, "\n") 
  cat("Degrees of Freedom of Independence Model (DFnull) =", DFnull, "\n")
  cat("Bartlett's corrective factor (Bartletc1) =", Bartletc1, "\n")
  cat("Swain's q1 factor (q1) =", q1, "\n")
  cat("Swain's corrective factor (Swainc1) =", Swainc1, "\n")
  cat("Yuan's correction (Yuanc1) =", Yuanc1, "\n")
  cat("The S-B scaling correction factor (scf) =", scf, "\n")
  cat("The ML, t-based chi-square (ML.t) =", ML.t, "\n")
  cat("The Bollen-Stine Chi-square =", BSx2, "\n")
  cat("The Bollen-Stine scaling correction factor =", BSscf, "\n")
  cat("\n                          ################# Results #################", "\n")  
  v1<-matrix(c(out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,out13,out14,out15,out16,out17,out18,
     out19,out20,out21,out22,out23,out24,out25,out26,out27,out28,out29,out30,out31,out32,out33,out34,out35,
     out36,out37,out38,out39,out40,out41,out42,out43,out44,out45,out46,out47,out48,out49,out50,out51,out52,
     out53,out54,out55,out56,out57,out58,out59,out60,out61,out62,out63,out64,out65,out66,out67,out68,out69,
     out70,out71,out72,out73,out74,out75,out76,out77,out78,out79,out80,out81,out82,out83,out84,out86,out87,
     out87,out88,out89,out90,out91,out92,out93,out94,out95,out96,out97,out98,out99,out100,out101,out102,
     out103,out104,out105,out106,out107,out108,out109,out110,out111,out112,out113,out114,out115,out116,
     out117,out118,out119,out120,out121,out122,out123,out124,out125,out126,out127,out128,out129,out130,
     out131,out132,out133,out134,out135,out136,out137,out138,out139,out140,out141,out142,out143,out144,
     out145,out146,out147,out148,out149,out150,out151,out152,out153,out154,out155,out156,out157,out158,out159,
     out160,out161,out162,out163,out164,out165,out166,out167,out168),ncol=24,nrow=7, byrow=TRUE)
  rownames(v1)<-c("1.Uncorrected Estimates","2.Bartlett's Correction","3.Yuan's Correction",
                  "4.Swain's Correction","5.Yuan's Ne1 Correction","6.Yuan's Ne2 Correction","7.Yuan's Ne3 Correction")
  colnames(v1)<-c("Chi-Sq.","Prob.","R.M.S.E.A","FF","LHR","Chi./DF","NFI","PNFI","TLI","CFI","PCFI",
                  "Rho","Delta","RFI","PRFI","IFI","PIFI","GFI","AGFI","PGFI","GAMMA-H","RNI","PRNI","NCI")
  options(scipen=999)
  options(digits=3)
  colors<-c("blue", "orange", "red", "black", "blueviolet", "chocolate1", "black")
  par(mfrow=c(3,4), mar=c(4.1, 4.1, 4.1, 2.1), las=3)
  a<-barplot(as.matrix(v1)[c(1:7)], col=colors, beside = TRUE, 
             main = "Chi-Square", 
             border = TRUE, legend.text = TRUE,
             ylab = "Estimates", axis.lty=1,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 10, 30, 90))
             abline(h=DF, col="Red", lty=5, lwd=2)
  b<-barplot(as.matrix(v1)[c(15:21)], col=colors, beside = TRUE, 
             main = "RMSEA", ylab = "Estimates",
             axis.lty=1,
             border = TRUE, legend.text = TRUE,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.05, col="Red", lty=5, lwd=2)
  c<-barplot(as.matrix(v1)[c(36:42)], col=colors, beside = TRUE, 
             main = "Chi-square/D.F.", ylab = "Estimates",
             axis.lty=1,
             border = TRUE, legend.text = TRUE,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=2, col="Red", lty=5, lwd=2)
  d<-barplot(as.matrix(v1)[c(43:49)], col=colors, beside = TRUE, 
             main = "NFI", ylab = "Estimates",
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  e<-barplot(as.matrix(v1)[c(57:63)], col=colors, beside = TRUE, 
             main = "TLI", ylab = "Estimates",
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  f<-barplot(as.matrix(v1)[c(64:70)], col=colors, beside = TRUE, 
             main = "CFI", ylab = "Estimates",
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  g<-barplot(as.matrix(v1)[c(106:112)], col=colors, beside = TRUE, 
             main = "IFI", ylab = "Estimates",
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE,
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  h<-barplot(as.matrix(v1)[c(120:126)], col=colors, beside = TRUE, 
             main = "GFI", ylab = "Estimates", 
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE, 
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  i<-barplot(as.matrix(v1)[c(127:133)], col=colors, beside = TRUE, 
             main = "AGFI", ylab = "Estimates", 
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE, 
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  j<-barplot(as.matrix(v1)[c(141:147)], col=colors, beside = TRUE, 
             main = "Gamma-H", ylab = "Estimates", 
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE, 
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)
  k<-barplot(as.matrix(v1)[c(148:154)], col=colors, beside = TRUE, 
             main = "RNI", ylab = "Estimates", 
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE, 
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)  
  l<-barplot(as.matrix(v1)[c(162:168)], col=colors, beside = TRUE, 
             main = "NCI", ylab = "Estimates", 
             ylim = c(0,1.1), axis.lty=1,
             border = TRUE, legend.text = TRUE, 
             names.arg = c("Base","Bartlett","Yuan","Swain", "Yuan-Ne1", "Yuan-Ne2", "Yuan-Ne3"),
             density=c(90, 70, 40, 20, 50, 60, 90))
  abline(h=.9, col="Red", lty=5, lwd=2)  
  v1
}

