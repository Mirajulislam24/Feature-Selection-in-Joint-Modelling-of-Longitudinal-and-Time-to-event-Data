
##Model part

jagsmodel<- "model{
   for (i in 1:N) {
    for (j in offset[i]:(offset[i + 1] - 1)) {
      muy1[j] <- inprod(betas1[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 
                                                                  1:ncZ1], Z1[j, 1:ncZ1])
      y1[j] ~ dnorm(muy1[j], tau1)
    }
    for (j in offset[i]:(offset[i + 1] - 1)) {
      muy2[j] <- inprod(betas2[1:ncX2], X2[j, 1:ncX2]) + 
        inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], Z2[j, 1:ncZ2])
      y2[j] ~ dnorm(muy2[j], tau2)
    }
    
    for (j in offset[i]:(offset[i + 1] - 1)) {
      muy3[j] <- inprod(betas3[1:ncX3], X3[j, 1:ncX3]) + 
        inprod(b[i, (ncZ1+ncZ2 + 1):(ncZ1 + ncZ2+ ncZ3)], Z3[j, 1:ncZ3])
      y3[j] ~ dnorm(muy3[j], tau3)
    }
    
    etaBaselineD[i] <- inprod(gammasD[1:(ncWD)], WD[i, 1:ncWD])
    log.h0.TD[i] <- inprod(Bs.gammasD[1:(ncW2D)], W2D[i, 
                                                      1:ncW2D])
    f.T1[i] <- inprod(betas1[1:ncX1], Xtime1[i, 1:ncX1]) + inprod(b[i, 
                                                                    1:ncZ1], Ztime1[i, 1:ncZ1])
    f.T.derivY1[i] <- inprod(betas1[2:3], Xtime.derivY1[i, 
                                                        1:ncX.derivY1]) + inprod(b[i, 2:3], Ztime.derivY1[i, 
                                                                                                          1:ncZ.derivY1])
    fA.T.derivY1[i] <- inprod(betas1[1:3], XAtime.derivY1[i, 
                                                          1:ncXA.derivY1]) + inprod(b[i, 1:3], ZAtime.derivY1[i, 
                                                                                                              1:ncZA.derivY1])
    
    
    f.I1[i]<-step(f.T1[i]-threshold[1])
    
    f.T2[i] <- inprod(betas2[1:ncX2], Xtime2[i, 1:ncX2]) + 
      inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], Ztime2[i, 1:ncZ2])
    f.T.derivY2[i] <- inprod(betas2[2:3], Xtime.derivY2[i, 
                                                        1:ncX.derivY2]) + inprod(b[i, 5:6], Ztime.derivY2[i, 
                                                                                                          1:ncZ.derivY2])
    fA.T.derivY2[i] <- inprod(betas2[1:3], XAtime.derivY2[i, 
                                                          1:ncXA.derivY2]) + inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], 
                                                                                    ZAtime.derivY2[i, 1:ncZA.derivY2])
    
    f.I2[i]<-step(f.T2[i]-threshold[2])
    
    f.T3[i] <- inprod(betas3[1:ncX3], Xtime3[i, 1:ncX3]) + 
      inprod(b[i, (ncZ1 +ncZ2+ 1):(ncZ1 + ncZ2+ncZ3)], 
             Ztime3[i, 1:ncZ3])
    f.T.derivY3[i] <- inprod(betas3[2:3], Xtime.derivY3[i, 
                                                        1:ncX.derivY3]) + inprod(b[i, 8:9], Ztime.derivY3[i, 
                                                                                                          1:ncZ.derivY3])
    fA.T.derivY3[i] <- inprod(betas3[1:3], XAtime.derivY3[i, 
                                                          1:ncXA.derivY3]) + inprod(b[i, (ncZ1+ncZ2 + 1):(ncZ1 + ncZ2+ncZ3)], 
                                                                                    ZAtime.derivY3[i, 1:ncZA.derivY3])
    
    
    f.I3[i]<-step(f.T3[i]-threshold[3])
    
    log.hazardD[i] <- log.h0.TD[i] + etaBaselineD[i] + alphasD1 * 
      ((f.T1[i] - muV1)/stdV1) + DalphasD1 * ((f.T.derivY1[i] - 
                                                 muS1)/stdS1) + DAalphasD1 * ((fA.T.derivY1[i] - muA1)/stdA1) +DIalphasD1 * f.I1[i]+
      alphasD2 * ((f.T2[i] - muV2)/stdV2) + DalphasD2 * 
      ((f.T.derivY2[i] - muS2)/stdS2) + DAalphasD2 * ((fA.T.derivY2[i] - 
                                                         muA2)/stdA2) +DIalphasD2 * f.I2[i]+ 
      alphasD3 * ((f.T3[i] - muV3)/stdV3) + DalphasD3 * 
      ((f.T.derivY3[i] - muS3)/stdS3) + DAalphasD3 * ((fA.T.derivY3[i] - 
                                                         muA3)/stdA3) +DIalphasD3 * f.I3[i]
    
    #########################################################################################################################   
    for (k in 1:K) {
      log.h0.sD[i, k] <- inprod(Bs.gammasD[1:(ncW2D)], 
                                W2sD[K * (i - 1) + k, 1:ncW2D])
      f.s1[i, k] <- inprod(betas1[1:ncX1], Xs1[K * (i - 1) + 
                                                 k, 1:ncX1]) + inprod(b[i, 1:ncZ1], Zs1[K * (i - 
                                                                                               1) + k, 1:ncZ1])
      f.s.derivY1[i, k] <- inprod(betas1[2:3], Xs.derivY1[K * 
                                                            (i - 1) + k, 1:ncX.derivY1]) + inprod(b[i, 2:3], 
                                                                                                  Zs.derivY1[K * (i - 1) + k, 1:ncZ.derivY1])
      
      fA.s.derivY1[i, k] <- inprod(betas1[1:3], XAs.derivY1[K * 
                                                              (i - 1) + k, 1:ncXA.derivY1]) + inprod(b[i, 1:3],
                                                                                                     ZAs.derivY1[K * (i - 1) + k, 1:ncZA.derivY1])
      f.sI1[i,k]<-step(f.s1[i, k]-threshold[1])
      
      f.s2[i, k] <- inprod(betas2[1:ncX2], Xs2[K * (i - 
                                                      1) + k, 1:ncX2]) + inprod(b[i, 4:6], Zs2[K * (i - 1) + k, 
                                                                                               1:ncZ2])
      f.s.derivY2[i, k] <- inprod(betas2[2:3], Xs.derivY2[K * 
                                                            (i - 1) + k, 1:ncX.derivY2]) + inprod(b[i, 5:6], 
                                                                                                  Zs.derivY2[K * (i - 1) + k, 1:ncZ.derivY2])
      
      
      fA.s.derivY2[i, k] <- inprod(betas2[1:3], XAs.derivY2[K * 
                                                              (i - 1) + k, 1:ncXA.derivY2]) + inprod(b[i, 4:6], ZAs.derivY2[K * (i - 1) + k, 
                                                                                                                            1:ncZA.derivY2])    
      
      
      f.sI2[i,k]<-step(f.s2[i, k]-threshold[2])
      f.s3[i, k] <- inprod(betas3[1:ncX3], Xs3[K * (i - 
                                                      1) + k, 1:ncX3]) + inprod(b[i, 7:9], Zs3[K * (i - 1) + k, 
                                                                                               1:ncZ3])
      f.s.derivY3[i, k] <- inprod(betas3[2:3], Xs.derivY3[K * 
                                                            (i - 1) + k, 1:ncX.derivY3]) + inprod(b[i, 8:9], 
                                                                                                  Zs.derivY3[K * (i - 1) + k, 1:ncZ.derivY3])
      
      
      fA.s.derivY3[i, k] <- inprod(betas3[1:3], XAs.derivY3[K * 
                                                              (i - 1) + k, 1:ncXA.derivY3]) + inprod(b[i, 7:9], ZAs.derivY3[K * (i - 1) + k, 
                                                                                                                            1:ncZA.derivY3])                                                            
      
      
      f.sI3[i,k]<-step(f.s3[i, k]-threshold[3])
      SurvLongD[i, k] <-  wk[k] * exp(log.h0.sD[i, k] + 
                                        alphasD1 * ((f.s1[i, k] - muV1)/stdV1) + DalphasD1 * 
                                        ((f.s.derivY1[i, k] - muS1)/stdS1) + DAalphasD1 * 
                                        ((fA.s.derivY1[i, k] - muA1)/stdA1)+DIalphasD1 * f.sI1[i,k] + 
                                        alphasD2 *((f.s2[i, k] - muV2)/stdV2) + DalphasD2 * ((f.s.derivY2[i, 
                                                                                                          k] - muS2)/stdS2) + DAalphasD2 * ((fA.s.derivY2[i, 
                                                                                                                                                          k] - muA2)/stdA2)+DIalphasD2 * f.sI2[i,k] + 
                                        alphasD3 *((f.s3[i, k] - muV3)/stdV3) + DalphasD3 * ((f.s.derivY3[i, 
                                                                                                          k] - muS3)/stdS3) + DAalphasD3 * ((fA.s.derivY3[i, 
                                                                                                                                                          k] - muA3)/stdA3)+DIalphasD3 * f.sI3[i,k]) 
    }
    #zeros-ones trick to get likelihood
    log.survivalD[i] <- -exp(etaBaselineD[i]) * P[i] * sum(SurvLongD[i, 
    ])
    phi[i] <- C - ((eventD[i] * log.hazardD[i])) - (log.survivalD[i])
    zeros[i] ~ dpois(phi[i])
    
    
    b[i,1]<-v[i,1]
    b[i,2]<-v[i,4]
    b[i,3]<-v[i,7]
    b[i,4]<-v[i,2]
    b[i,5]<-v[i,5]
    b[i,6]<-v[i,8]
    b[i,7]<-v[i,3]
    b[i,8]<-v[i,6]
    b[i,9]<-v[i,9]
    
    v[i, 1:3] ~ dmnorm(mu01[], inv.D1[, ])
    v[i, 4:6] ~ dmnorm(mu02[], inv.D2[, ])
    v[i, 7:9] ~ dmnorm(mu03[], inv.D3[, ])
  }
  
  #Priors
  
  betas1[1:ncX1] ~ dmnorm(priorMean.betas1[], priorTau.betas1[, 
  ])
  betas2[1:ncX2] ~ dmnorm(priorMean.betas2[], priorTau.betas2[, 
  ])
  betas3[1:ncX3] ~ dmnorm(priorMean.betas3[], priorTau.betas3[, 
  ])
  tau1 ~ dgamma(priorA.tau1, priorB.tau1)
  tau2 ~ dgamma(priorA.tau2, priorB.tau2)
  tau3 ~ dgamma(priorA.tau3, priorB.tau3)
  
  sigma1<-1/tau1
  sigma2<-1/tau2
  sigma3<-1/tau3
  
  gammasD~ dmnorm(priorMean.gammas[], priorTau.gammas[,])
   ########## SPIKE and SLAB prior specification#####################
  alphasD1<-tau11*d11
  DalphasD1<-tau12*d12
  DAalphasD1<-tau13*d13
  DIalphasD1<-tau14*d14
  alphasD2<-tau21*d21
  DalphasD2<-tau22*d22
  DAalphasD2<-tau23*d23
  DIalphasD2<-tau24*d24
  alphasD3<-tau31*d31
  DalphasD3<-tau32*d32
  DAalphasD3<-tau33*d33
  DIalphasD3<-tau34*d34
  
  ######## KM method for spike(point mass) and slab #######
  invs2~dgamma(1,1)
  d11<-h1[1]*p[1]
  d12<-h1[2]*p[1]
  d13<-h1[3]*p[1]
  d14<-h1[4]*p[1]
  d21<-h2[1]*p[2]
  d22<-h2[2]*p[2]
  d23<-h2[3]*p[2]
  d24<-h2[4]*p[2]
  d31<-h3[1]*p[3]
  d32<-h3[2]*p[3]
  d33<-h3[3]*p[3]
  d34<-h3[4]*p[3]
  
  tau11<-t1[1]*p1[1]
  tau12<-t1[2]*p1[2]
  tau13<-t1[3]*p1[3]
  tau14<-t1[4]*p1[4]
  tau21<-t2[1]*p2[1]
  tau22<-t2[2]*p2[2]
  tau23<-t2[3]*p2[3]
  tau24<-t2[4]*p2[4]
  tau31<-t3[1]*p3[1]
  tau32<-t3[2]*p3[2]
  tau33<-t3[3]*p3[3]
  tau34<-t3[4]*p3[4]
  
  
  for (l in 1:3)#for risk factors
  {
    pi[l]~dbeta(1,1)# equal prob
    p[l]~dbern(pi[l])
  }
  
  for (m in 1:4)# for features
  {
    p1[m]~dbern(pi1[m])#more 0
    p2[m]~dbern(pi2[m])
    p3[m]~dbern(pi3[m])
    t1[m]~dnorm(0, invs2) T(0,)
    t2[m]~dnorm(0, invs2) T(0,)
    t3[m]~dnorm(0, invs2) T(0,)
    h1[m]~dnorm(0,1)
    h2[m]~dnorm(0,1)
    h3[m]~dnorm(0,1)
  }
  
  
# assuming combinations: {1000,0100,0010,0001,1100,1010,1001,0110,0101,0011,1110,1101,1011,0111,1111}
  # a>b, b>c, c>1, d>0
  
  q1[1:15]~ddirch(a_11)#for 1st risk factor
  a_11<-c(a11,a11,a11,a11,a12,a12,a12,a12,a12,a12,a13,a13,a13,a13,a14)
  a14~dt(0,1,1) T(0,)#all 4 features
  a13~dt(0,1,1) T(a14,)# 3 features
  a12~dt(0,1,1) T(a13,)# 2 features 
  a11~dt(0,1,1) T(a12,)# 1 feature
  
  A1=4*a11+6*a12+4*a13+a14
  A11=4*a11/A1
  A12=6*a12/A1
  A13=4*a13/A1
  A14=a14/A1
  pi1[1]=sum(q1[1],q1[5:7],q1[11:13],q1[15])
  pi1[2]=sum(q1[2],q1[5],q1[8:9],q1[11:12],q1[14:15])
  pi1[3]=sum(q1[3],q1[6],q1[8],q1[10:11],q1[13:15])
  pi1[4]=sum(q1[4],q1[7],q1[9:10],q1[12:15])
  
  q2[1:15]~ddirch(a_22)#for 2nd risk factor
  a_22<-c(a21,a21,a21,a21,a22,a22,a22,a22,a22,a22,a23,a23,a23,a23,a24)
  
  a24~dt(0,1,1) T(0,)#all 4 features
  a23~dt(0,1,1) T(a24,)# 3 features
  a22~dt(0,1,1) T(a23,)# 2 features 
  a21~dt(0,1,1) T(a22,)# 1 feature
  
  A2=4*a21+6*a22+4*a23+a24
  A21=4*a21/A2
  A22=6*a22/A2
  A23=4*a23/A2
  A24=a24/A2
  pi2[1]=sum(q2[1],q2[5:7],q2[11:13],q2[15])
  pi2[2]=sum(q2[2],q2[5],q2[8:9],q2[11:12],q2[14:15])
  pi2[3]=sum(q2[3],q2[6],q2[8],q2[10:11],q2[13:15])
  pi2[4]=sum(q2[4],q2[7],q2[9:10],q2[12:15])
  
  q3[1:15]~ddirch(a_33)#for 3rd risk factor
  a_33<-c(a31,a31,a31,a31,a32,a32,a32,a32,a32,a32,a33,a33,a33,a33,a34)
  
  a34~dt(0,1,1) T(0,)#all 4 features
  a33~dt(0,1,1) T(a34,)# 3 features
  a32~dt(0,1,1) T(a33,)# 2 features 
  a31~dt(0,1,1) T(a32,)# 1 feature
  
  A3=4*a31+6*a32+4*a33+a34
  A31=4*a31/A3
  A32=6*a32/A3
  A33=4*a33/A3
  A34=a34/A3
  pi3[1]=sum(q3[1],q3[5:7],q3[11:13],q3[15])
  pi3[2]=sum(q3[2],q3[5],q3[8:9],q3[11:12],q3[14:15])
  pi3[3]=sum(q3[3],q3[6],q3[8],q3[10:11],q3[13:15])
  pi3[4]=sum(q3[4],q3[7],q3[9:10],q3[12:15])
  
  
  p11=p[1]*p1[1]
  p12=p[1]*p1[2]
  p13=p[1]*p1[3]
  p14=p[1]*p1[4]
  p21=p[2]*p2[1]
  p22=p[2]*p2[2]
  p23=p[2]*p2[3]
  p24=p[2]*p2[4]
  p31=p[3]*p3[1]
  p32=p[3]*p3[2]
  p33=p[3]*p3[3]
  p34=p[3]*p3[4]
  
  
  ###############################################################################
  
  Bs.gammasD[1:(ncW2D)] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, 
  ])
  inv.D1[1:3, 1:3] ~ dwish(priorR.D1[, ], priorK.D1)
  inv.D2[1:3, 1:3] ~ dwish(priorR.D2[, ], priorK.D2)
  inv.D3[1:3, 1:3] ~ dwish(priorR.D3[, ], priorK.D3)
  for(i in 1:3){
    Tbetas1[i]=betas1[i]*SD1+mean1[i]
    Tbetas2[i]=betas2[i]*SD2+mean2[i]
    Tbetas3[i]=betas3[i]*SD3+mean3[i]
  }
  Tsigma1=sigma1*SD1[1]^2
  Tsigma2=sigma2*SD2[1]^2
  Tsigma3=sigma3*SD3[1]^2
  TalphasD1=alphasD1/(stdV1*SD1)
  TDalphasD1=DalphasD1/(stdS1*SD1)
  TDAalphasD1=DAalphasD1/(stdA1*SD1)
  TalphasD2=alphasD2/(stdV2*SD2)
  TDalphasD2=DalphasD2/(stdS2*SD2)
  TDAalphasD2=DAalphasD2/(stdA2*SD2)
  TalphasD3=alphasD3/(stdV3*SD3)
  TDalphasD3=DalphasD3/(stdS3*SD3)
  TDAalphasD3=DAalphasD3/(stdA3*SD3)
}"


##Running codes in jags
load("SimData_1R2F.R")
param=c("Tbetas1","Tbetas2","Tbetas3","Tsigma1","Tsigma2","Tsigma3", "TalphasD1",
        "TDalphasD1","TDAalphasD1","DIalphasD1","TalphasD2","TDalphasD2","TDAalphasD2","DIalphasD2","TalphasD3",
        "TDalphasD3","TDAalphasD3","DIalphasD3","invs2","p")

init.rng1<-list(".RNG.seed" = 3333, ".RNG.name" = "base::Super-Duper")
init.rng2<-list(".RNG.seed" = 4444, ".RNG.name" = "base::Super-Duper")
cl <- makeCluster(2)
fit <- run.jags(jagsmodel, n.chains = 2,
                monitor=param,burnin = 10000,sample = 1000,thin=5,
                adapt =100,data=Data,method="parallel",inits = list(init.rng1, init.rng2),cl=cl)

library('coda')
coeff <- as.mcmc.list(fit, vars=param)
library(MCMCvis)  
f<-function(x){(1-sum(x==0)/length(x))*100}
summary<-MCMCsummary(coeff,params=param,func=f, round = 4)
summary

##Sensitivity and specificity
P=summary[13:24,8]#percentage of selection for the features 
P_70=ifelse(P>=70,1,0)#indicator
#observed important among true features/total true features
Sen_70=sum(P_70[1]+P_70[3])/2
#observed unimportant among true unimportant features/total true unimportant features
Spec_70=(10-(sum(P_70)-sum(P_70[1]+P_70[3])))/10

P_80=ifelse(P>=80,1,0)
P_85=ifelse(P>=85,1,0)
P_90=ifelse(P>=90,1,0)
P_95=ifelse(P>=95,1,0)

Sen_80=sum(P_80[1]+P_80[3])/2;Spec_80=1-(sum(P_80)-sum(P_80[1]+P_80[3]))/10
Sen_85=sum(P_85[1]+P_85[3])/2;Spec_85=1-(sum(P_85)-sum(P_85[1]+P_85[3]))/10
Sen_90=sum(P_90[1]+P_90[3])/2;Spec_90=1-(sum(P_90)-sum(P_90[1]+P_90[3]))/10
Sen_95=sum(P_95[1]+P_95[3])/2;Spec_95=1-(sum(P_95)-sum(P_95[1]+P_95[3]))/10
