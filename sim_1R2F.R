set.seed(1234)

x=c("JMbayes","splines","parallel","runjags","coda")
#install.packages(x)
lapply(x, library, character.only = TRUE)

T=50 #maximum time
legP2<-function(t){
  #T<-max(t)
  cbind((2*t/T-1),(-0.5+1.5*(2*t/T-1)^2))
}
dlegP2<-function(t){
  #T<-max(t)
  cbind(2/T,(6/T)*(2*t/T-1))
}
IlegP2<-function(t){
  #T<-max(t)
  cbind(t^2/T-t, (-0.5*t+T/4*((2*t/T-1)^3+1)))
}
#set.seed(4537)
SimJoint <- function (alpha11, alpha12,alpha13,alpha21,alpha22,alpha23,alpha31,alpha32,alpha33, n, upp_Cens) {
  #n=30;upp_Cens=12;set.seed(sample(1000:9999,1))
  K <- 5  # number of planned repeated measurements per subject, per outcome
  t.max <- upp_Cens # maximum follow-up time
  ################################################
  # parameters for the linear mixed effects models
  
  #BMI
  betas1 <- c("Intercept" =  27.897, "Time1" =-0.124, "Time2" =  -2.57 )
  sigma.y1 <- sqrt(1.622) # measurement error standard deviation
  #SBP
  betas2 <- c("Intercept" = 129.239, "Time1" = 12.694, "Time2" = -4.342 )
  sigma.y2 <- sqrt(134.816) # measurement error standard deviation
  #DBP
  betas3 <- c("Intercept" = 190.624 , "Time1" = -38.326  , "Time2" = -9.121)
  sigma.y3 <- sqrt( 540.374) # measurement error standard deviation
  
  
  
  D1=D2=D3<-matrix(0,3,3)
  D1[lower.tri(D1, TRUE)]<-c( 28.16, 7.92, -14.31,186.21,-.108,903.07)
  D1 <- D1 + t(D1)
  diag(D1) <- diag(D1) * 0.5
  
  D2[lower.tri(D2, TRUE)]<-c(14.43, 18.64, 17.29,204.07,156.14,1402.61)
  D2 <- D2 + t(D2)
  diag(D2) <- diag(D2) * 0.5
  
  D3[lower.tri(D3, TRUE)]<-c( 9.62,9.39,24.20,105.01,145.96,225.63)
  D3 <- D3 + t(D3)
  diag(D3) <- diag(D3) * 0.5
  
  
  D<-matrix(0,9,9)
  D[1:3,1:3]<-D1
  D[4:6,4:6]<-D2
  D[7:9,7:9]<-D3
  
  
  
  # design matrices for the longitudinal measurement model
  #In ARIC, 1visit 6%, 2visit 9%, 3visit 11%, 4visit 36%, 5visit 38%
  r<-sample(0:4,n,replace = TRUE,prob=c(.06,.09,.11,.36,.38))
  times<-list()
  
  time0<-c(0,sample(0:19,n-2,replace=T),19)#max(baseAge)=19, max(follow-up)=12, so max(age at end)=31
  for (i in 1:n){
    interval<-c(time0[i]+3, time0[i]+6,time0[i]+9,time0[i]+12)
    times[[i]] <-if(r[i]==0) time0[i]
    else c(time0[i],interval[1:r[i]])
  }
  times1<-unlist(times)
  r=r+1#maximum 5 visits with baseline
  group <- rep(0:1, each = n/2)
  #r1<-r+1
  RACE<-list()
  for (i in 1:n){
    RACE[[i]] <-factor(rep(group[i],each=r[i]))
  }
  RACE=factor(unlist(RACE))
  DF <- data.frame(time = times1,RACE=RACE)
  
  
  X1 <- model.matrix(~legP2(time), data = DF)
  Z1 <- model.matrix(~ legP2(time), data = DF)
  
  X2 <- model.matrix(~ legP2(time), data = DF)
  Z2 <- model.matrix(~ legP2(time), data = DF)
  
  X3 <- model.matrix(~ legP2(time), data = DF)
  Z3 <- model.matrix(~ legP2(time), data = DF)
  
  
  # design matrix for the survival model
  W <- cbind( "RACE" = group)
  
  ################################################
  
  # simulate random effects
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  
  
  # simulate longitudinal responses
  ID<-list()
  for (i in 1:n){
    ID[[i]] <- rep(i, r[i])
  }
  id<-unlist(ID)
  eta.y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b[id,c(1,4,7) ]))
  y1 <- rnorm(length(id), eta.y1, sigma.y1)
  
  eta.y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b[id,c(2,5,8) ]))
  y2 <- rnorm(length(id), eta.y2, sigma.y2)
  
  eta.y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b[id,c(3,6,9) ]))
  y3 <- rnorm(length(id), eta.y3, sigma.y3)
  
  # parameters for the survival model
  gammas <- c("RACE" = 0.15)
  phi <- 1.75
  # simulate event times
  eta.t <- as.vector(W %*% gammas)
  
  #alpha11 =-0.5;alpha13=0.009;alpha12=0; alpha22=alpha21=alpha23=alpha31=alpha32=alpha33=0
  invS <- function (t, u, i) {
    h <- function (s) {
      P2 <-legP2(s)
      dP2 <- dlegP2(s)
      IP2<-IlegP2(s)
      
      
      XX1=XX2=XX3 <- cbind(1,P2[,1],P2[,2])
      ZZ1=ZZ2=ZZ3 <- cbind(1,P2[,1],P2[,2])
      XXd1=XXd2=XXd3 <- cbind(dP2[,1],dP2[,2])
      ZZd1=ZZd2=ZZd3 <- cbind(dP2[,1],dP2[,2])
      XXA1=XXA2=XXA3 <- cbind(s,IP2[,1],IP2[,2])
      ZZA1=ZZA2=ZZA3 <- cbind(s,IP2[,1],IP2[,2])
      
      
      fV1 <- as.vector(XX1 %*% betas1 + rowSums(ZZ1 * b[rep(i, nrow(ZZ1)),c(1,4,7)]))
      fS1 <- as.vector(XXd1 %*% betas1[2:3] + rowSums(ZZd1 * b[rep(i, nrow(ZZd1)),c(4,7)]))
      fA1 <- as.vector(XXA1 %*% betas1 + rowSums(ZZA1 * b[rep(i,nrow(ZZA1)),c(1,4,7)]))
      
      
      fV2 <- as.vector(XX2 %*% betas2 + rowSums(ZZ2 * b[rep(i, nrow(ZZ2)),c(2,5,8)]))
      fS2 <- as.vector(XXd2 %*% betas2[2:3] + rowSums(ZZd2 * b[rep(i, nrow(ZZd2)),c(5,8)]))
      fA2 <- as.vector(XXA2 %*% betas2 + rowSums(ZZA2 * b[rep(i,nrow(ZZA2)),c(2,5,8)]))
      
      
      
      fV3 <- as.vector(XX3 %*% betas3 + rowSums(ZZ3 * b[rep(i, nrow(ZZ3)),c(3,6,9)]))
      fS3 <- as.vector(XXd3 %*% betas3[2:3] + rowSums(ZZd3 * b[rep(i, nrow(ZZd3)),c(6,9)]))
      fA3 <- as.vector(XXA3 %*% betas3 + rowSums(ZZA3 * b[rep(i,nrow(ZZA3)),c(3,6,9)]))
      
      
      exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + fV1 * alpha11 +fS1 * alpha12+fA1 * alpha13+
            fV2 * alpha21 +fS2 * alpha22+fA2 * alpha23+ fV3 * alpha31+fS3* alpha32+fA3 * alpha33)
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  
  u <- runif(n)
  (Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[1], i = 1)$root, TRUE))
  last_visit<-sapply(times,tail,1)
  trueTimes<- numeric(n)
  for (i in 1:n) {
    Up <- 40
    tries <- 5
    Root <- try(uniroot(invS, interval = c(last_visit[i]+0.0001, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 40
      Root <- try(uniroot(invS, interval = c(last_visit[i]+0.0001, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  
  na.ind <- !is.na(trueTimes)
  trueTimes <- trueTimes[na.ind]
  W <- W[na.ind, , drop = FALSE]
  
  long.na.ind<-list()
  for (i in 1:n){
    long.na.ind[[i]] <- rep(na.ind[i], each =r[i] )
  }
  long.na.ind <- unlist(long.na.ind)
  y1 <- y1[long.na.ind]
  X1 <- X1[long.na.ind, , drop = FALSE]
  Z1 <- Z1[long.na.ind, , drop = FALSE]
  y2 <- y2[long.na.ind]
  X2 <- X2[long.na.ind, , drop = FALSE]
  Z2 <- Z2[long.na.ind, , drop = FALSE]
  y3 <- y3[long.na.ind]
  X3 <- X3[long.na.ind, , drop = FALSE]
  Z3 <- Z3[long.na.ind, , drop = FALSE]
  DF <- DF[long.na.ind, , drop = FALSE]
  n <- length(trueTimes)
  
  Ctimes <-  time0+t.max
  Ctimes <- Ctimes[na.ind]
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  
  ################################################
  
  # keep the nonmissing cases, i.e., drop the longitudinal measurements
  # that were taken after the observed event time for each subject.
  
  Time2<-list()
  r1<-r[na.ind]
  for (i in 1:n){
    Time2[[i]]<-rep(Time[i],r1[i])
  }
  
  
  ind<-times1[long.na.ind] <= unlist(Time2)
  y1 <- y1[ind]
  X1 <- X1[ind, , drop = FALSE]
  Z1 <- Z1[ind, , drop = FALSE]
  y2 <- y2[ind]
  X2 <- X2[ind, , drop = FALSE]
  Z2 <- Z2[ind, , drop = FALSE]
  y3 <- y3[ind]
  X3 <- X3[ind, , drop = FALSE]
  Z3 <- Z3[ind, , drop = FALSE]
  id <- id[long.na.ind][ind]
  id <- match(id, unique(id))
  
  dat <- DF[ind,]
  dat$id <- id
  dat$y1 <- y1
  dat$y2 <- y2
  dat$y3 <- y3
  dat$Time <- Time[id]
  dat$event <- event[id]
  dat <- dat[c("id", "y1","y2","y3","RACE", "time", "Time", "event")]
  dat.id <- data.frame(id = unique(id), Time = Time, event = event, RACE = W[, 1])
  trueValues <- list(betas1 = betas1, sigmas1 = sigma.y1,betas2 = betas2, sigmas2 = sigma.y2, betas3 = betas3, sigmas3 = sigma.y3,gammas = gammas,
                     alphas11 = alpha11,alphas12 = alpha12, alphas21 = alpha21,alphas22 = alpha22, alphas31=alpha31,alphas32 = alpha32, sigma.t = phi,
                     D = D, b = b)
  
  # return list
  list(DF = dat, DF.id = dat.id, trueValues = trueValues)
}
simdata <- SimJoint(alpha11 =-0.55, alpha12 =0,alpha13=0.01,alpha21=0,alpha22=0,alpha23=0,alpha31=0,alpha32=0,alpha33=0, n= 800, upp_Cens = 31)
table(simdata$DF.id$event);table(simdata$DF$event)
save(simdata,file="data_1R2F.R")