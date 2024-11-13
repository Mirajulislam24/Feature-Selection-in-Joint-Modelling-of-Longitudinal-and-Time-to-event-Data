load("data_1R2F.R")
data = simdata$DF
data.id=simdata$DF.id

##Standardizing variables
data$z1<-scale(data$y1)
data$z2<-scale(data$y2)
data$z3<-scale(data$y3)
means<-c(mean(data$y1),mean(data$y2),mean(data$y3))
SDs<-c(sd(data$y1),sd(data$y2),sd(data$y3))
threshold<-c(30,120,230)
threshold<-(threshold-means)/SDs

# fit the longitudinal outcomes using mixed-effetcs models
lmeObject1<- lme(z1~ legP2(time), data = data,
                 random=list(id=pdDiag(form=~legP2(time))))
lmeObject2<- lme(z2~ legP2(time), data = data,
                 random=list(id=pdDiag(form=~legP2(time))))
lmeObject3<- lme(z3~ legP2(time), data = data,
                 random=list(id=pdDiag(form=~legP2(time))))

### Set
timeVar <- "time"
lag <- 0
survMod <- "spline-PH"
Time <- data.id$Time###the survival time

# for the continuous longitudinal outcome create the design matrices
id <- data$id
offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

# 1st longitudinal outcome
formYx1 <- formula(lmeObject1)
TermsX1 <- lmeObject1$terms
mfX1 <- model.frame(TermsX1, data = data)
# 2nd longitudinal outcome
formYx2 <- formula(lmeObject2)
TermsX2 <- lmeObject2$terms
mfX2 <- model.frame(TermsX2, data = data)
# 3rd longitudinal outcome
formYx3 <- formula(lmeObject3)
TermsX3 <- lmeObject3$terms
mfX3 <- model.frame(TermsX3, data = data)

# Longitudinal Response variables
y.long1 <- model.response(mfX1, "numeric")
y.long2 <- model.response(mfX2, "numeric")
y.long3 <- model.response(mfX3, "numeric")

y <- list(y1 = y.long1, y2 = y.long2,y3 = y.long3)


#Design matrices for longitudinal and hazard function
X1 <- model.matrix(formYx1, mfX1)
formYz1 <- formula(lmeObject1$modelStruct$reStruct[[1]])
mfZ1 <- model.frame(terms(formYz1), data = data)
TermsZ1 <- attr(mfZ1, "terms")
Z1 <- model.matrix(formYz1, mfZ1)

data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)

mfX.id1 <- model.frame(TermsX1, data = data.id) 
mfZ.id1 <- model.frame(TermsZ1, data = data.id)  
Xtime1 <- model.matrix(formYx1, mfX.id1)
Ztime1 <- model.matrix(formYz1, mfZ.id1)

# 2nd and 3rd longitudinal outcome
X2=X3=X1;Z2=Z3=Z1
Xtime2=Xtime3=Xtime1;Ztime2=Ztime3=Ztime1


#################################
# survival submodel
# design matrices for the survival submodel
WD <- model.matrix(~-1+RACE, data=data.id)
eventD <- data.id$event
nT <- length(Time)
zeros <- numeric(nT)

x <- list(X1 = X1, Z1 = Z1, WD = if (survMod == "weibull-PH") {
  if (is.null(WD)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                            WD)
} else {
  if (is.null(WD)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(WD) == 1) cbind(WD, rep(0, nT)) else WD
  }
})

###################################
# for the longitudinal outcomes - design matrices for the 15-point Gauss-Kronrod quadrature rule approximation
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)
P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))


# 1st longitudinal outcome
mfX1 <- model.frame(TermsX1, data = data.id2)  
mfZ1 <- model.frame(TermsZ1, data = data.id2)    
Xs1 <- model.matrix(formYx1, mfX1)
Zs1 <- model.matrix(formYz1, mfZ1)

# 2nd and 3rd longitudinal outcome
Xs2=Xs3=Xs1;Zs2=Zs3=Zs1
#################################
# set MCMC details
con <- list( K = 100,C = 5000, knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 4)

# design matrices for the baseline hazard (a B-splines baseline hazard function is asssumed)
kn <- if (is.null(con$knots)) {
  pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
  pp <- tail(head(pp, -1), -1)
  tt <- if (con$ObsTimes.knots) {
    Time
  } else {Time[event == 1]    }
  quantile(tt, pp, names = FALSE)
} else {
  con$knots
}
kn <- kn[kn < max(Time)]
rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
con$knots <- rr

W2D <- splineDesign(rr, Time, ord = con$ordSpline)
if (any(colSums(W2D) == 0))
  stop("\nsome of the knots of the B-splines basis are set outside the range",
       "\n   of the observed event times for one of the strata; refit the model",
       "\n   setting the control argument 'equal.strata.knots' to FALSE.")

# design matrices for the baseline hazard for the 15-point Gauss-Kronrod quadrature rule approximation

W2sD <- splineDesign(rr, c(t(st)), ord = con$ordSpline)

x <- c(x, list(W2D = W2D, W2sD = W2sD))

#################################
ncX1 <- ncol(X1)
ncZ1 <- ncol(Z1)
ncWD <- ncol(x$WD)
ncW2D <- ncol(x$W2D)
ncX2 <- ncol(X2)
ncZ2 <- ncol(Z2)
ncX3 <- ncol(X3)
ncZ3 <- ncol(Z3)


C <- con$C
nb <- ncZ1 + ncZ2+ncZ3
b <- cbind(data.matrix(ranef(lmeObject1)),data.matrix(ranef(lmeObject2)),data.matrix(ranef(lmeObject3)))

nY <- nrow(b)
sigma21 <- lmeObject1$sigma^2#residual SD^2
sigma22 <- lmeObject2$sigma^2
sigma23 <- lmeObject3$sigma^2



#################################
# priors/hyperpriors
mu01=mu02=mu03<- rep(0, 3)

betas1=betas2=betas3 <- rep(0, ncX1)
var.betas1=var.betas2=var.betas3 <- rep(con$K, ncX1)

gammasD <- rep(0,ncWD); var.gammasD <- rep(con$K, (ncWD))
Bs.gammasD <- rep(0, (ncW2D)); var.Bs.gammasD <- rep(con$K/10, (ncW2D))


#################################SLOPE PARAMETERIZATION###############
# design matrices for the slope of the 1st longitudinal outcome
extraFormY1 <- list(fixed=~-1+dlegP2(time),indFixed=2:3,random=~-1+dlegP2(time),indRandom=2:3)
mfX.derivY1 <- model.frame(terms(extraFormY1$fixed), data = data)
TermsX.derivY1 <- attr(mfX.derivY1, "terms")
mfZ.derivY1 <- model.frame(terms(extraFormY1$random), data = data)
TermsZ.derivY1 <- attr(mfZ.derivY1, "terms")
mfX.deriv.idY1 <- model.frame(TermsX.derivY1, data = data.id)
mfZ.deriv.idY1 <- model.frame(TermsZ.derivY1, data = data.id)
Xtime.derivY1 <- model.matrix(extraFormY1$fixed, mfX.deriv.idY1)
Ztime.derivY1 <- model.matrix(extraFormY1$random, mfZ.deriv.idY1)
XderivY1 <- model.matrix(extraFormY1$fixed, mfX.derivY1)
ZderivY1 <- model.matrix(extraFormY1$random, mfZ.derivY1)


mfX.derivY1 <- model.frame(TermsX.derivY1, data = data.id2)
mfZ.derivY1 <- model.frame(TermsZ.derivY1, data = data.id2)
Xs.derivY1 <- model.matrix(extraFormY1$fixed, mfX.derivY1)
Zs.derivY1 <- model.matrix(extraFormY1$random, mfZ.derivY1)

# Slope features of 2nd and 3rd longitudinal outcome
Xtime.derivY2=Xtime.derivY3=Xtime.derivY1;Ztime.derivY2=Ztime.derivY3=Ztime.derivY1
Xs.derivY2=Xs.derivY3=Xs.derivY1;Zs.derivY2=Zs.derivY3=Zs.derivY1


# design matrices for the area of the 1st longitudinal outcome
AextraFormY1 <- list(fixed=~-1+time+IlegP2(time),indFixed=1:3,random=~-1+time+IlegP2(time),indRandom=1:3)
mfXA.derivY1 <- model.frame(terms(AextraFormY1$fixed), data = data)
TermsXA.derivY1 <- attr(mfXA.derivY1, "terms")
mfZA.derivY1 <- model.frame(terms(AextraFormY1$random), data = data)
TermsZA.derivY1 <- attr(mfZA.derivY1, "terms")
mfXA.deriv.idY1 <- model.frame(TermsXA.derivY1, data = data.id)
mfZA.deriv.idY1 <- model.frame(TermsZA.derivY1, data = data.id)
XAtime.derivY1 <- model.matrix(AextraFormY1$fixed, mfXA.deriv.idY1)
ZAtime.derivY1 <- model.matrix(AextraFormY1$random, mfZA.deriv.idY1)
XAderivY1 <- model.matrix(AextraFormY1$fixed, mfXA.derivY1)
ZAderivY1 <- model.matrix(AextraFormY1$random, mfZA.derivY1)


mfXA.derivY1 <- model.frame(TermsXA.derivY1, data = data.id2)
mfZA.derivY1 <- model.frame(TermsZA.derivY1, data = data.id2)
XAs.derivY1 <- model.matrix(AextraFormY1$fixed, mfXA.derivY1)
ZAs.derivY1 <- model.matrix(AextraFormY1$random, mfZA.derivY1)

# Area features of 2nd and 3rd longitudinal outcome
XAtime.derivY2=XAtime.derivY3=XAtime.derivY1;ZAtime.derivY2=ZAtime.derivY3=ZAtime.derivY1
XAs.derivY2=XAs.derivY3=XAs.derivY1;ZAs.derivY2=ZAs.derivY3=ZAs.derivY1



x <- c(x, list(Xs.derivY1 = Xs.derivY1, Zs.derivY1 = Zs.derivY1, XAs.derivY1 = XAs.derivY1, ZAs.derivY1 = ZAs.derivY1,
               Xs.derivY2 = Xs.derivY2, Zs.derivY2 = Zs.derivY2, XAs.derivY2 = XAs.derivY2, ZAs.derivY2 = ZAs.derivY2,
               Xs.derivY3 = Xs.derivY3, Zs.derivY3 = Zs.derivY3, XAs.derivY3 = XAs.derivY3, ZAs.derivY3 = ZAs.derivY3))

# 2 stage approach - calculation of the mean and sd of the longitudinal outcomes
# 1st longitudinal outcome
muV1 <- mean(Xtime1%*%fixef(lmeObject1) + rowSums(Ztime1*ranef(lmeObject1)))
stdV1 <- sd(Xtime1%*%fixef(lmeObject1) + rowSums(Ztime1*ranef(lmeObject1)))
muS1 <- mean(Xtime.derivY1%*%fixef(lmeObject1)[c(2:3)] + rowSums(Ztime.derivY1*ranef(lmeObject1)[c(2:3)]))
stdS1 <- sd(Xtime.derivY1%*%fixef(lmeObject1)[c(2:3)] + rowSums(Ztime.derivY1*ranef(lmeObject1)[c(2:3)]))
muA1 <- mean(XAtime.derivY1%*%fixef(lmeObject1) +  rowSums(ZAtime.derivY1*ranef(lmeObject1)))
stdA1 <- sd(XAtime.derivY1%*%fixef(lmeObject1) +  rowSums(ZAtime.derivY1*ranef(lmeObject1)))

# 2nd longitudinal outcome
muV2 <- mean(Xtime2%*%fixef(lmeObject2) + rowSums(Ztime2*ranef(lmeObject2)))
stdV2 <- sd(Xtime2%*%fixef(lmeObject2) + rowSums(Ztime2*ranef(lmeObject2)))
muS2 <- mean(Xtime.derivY2%*%fixef(lmeObject2)[c(2:3)] + rowSums(Ztime.derivY2*ranef(lmeObject2)[c(2:3)]))
stdS2 <- sd(Xtime.derivY2%*%fixef(lmeObject2)[c(2:3)] + rowSums(Ztime.derivY2*ranef(lmeObject2)[c(2:3)]))
muA2 <- mean(XAtime.derivY2%*%fixef(lmeObject2) +  rowSums(ZAtime.derivY2*ranef(lmeObject2)))
stdA2 <- sd(XAtime.derivY2%*%fixef(lmeObject2) +  rowSums(ZAtime.derivY2*ranef(lmeObject2)))

# 3rd longitudinal outcome
muV3 <- mean(Xtime3%*%fixef(lmeObject3) + rowSums(Ztime3*ranef(lmeObject3)))
stdV3 <- sd(Xtime3%*%fixef(lmeObject3) + rowSums(Ztime3*ranef(lmeObject3)))
muS3 <- mean(Xtime.derivY3%*%fixef(lmeObject3)[c(2:3)] + rowSums(Ztime.derivY3*ranef(lmeObject3)[c(2:3)]))
stdS3 <- sd(Xtime.derivY3%*%fixef(lmeObject3)[c(2:3)] + rowSums(Ztime.derivY3*ranef(lmeObject3)[c(2:3)]))
muA3 <- mean(XAtime.derivY3%*%fixef(lmeObject3) +  rowSums(ZAtime.derivY3*ranef(lmeObject3)))
stdA3 <- sd(XAtime.derivY3%*%fixef(lmeObject3) +  rowSums(ZAtime.derivY3*ranef(lmeObject3)))



#################################
Data <- list(N = nY, K = K, offset = offset, X1 = X1, Xtime1 = Xtime1, 
             Xtime.derivY1 = Xtime.derivY1, XAtime.derivY1 = XAtime.derivY1,
             y1 = y$y1, 
             Xs1 = Xs1, Xs.derivY1 = Xs.derivY1, XAs.derivY1 = XAs.derivY1,
             Z1 = Z1, Ztime1 = Ztime1,  Ztime.derivY1 = Ztime.derivY1, ZAtime.derivY1 = ZAtime.derivY1,
             Zs1 = Zs1, Zs.derivY1 = Zs.derivY1,  ZAs.derivY1 = ZAs.derivY1, 
             
             X2 = X2, Xtime2 = Xtime2, 
             Xtime.derivY2 = Xtime.derivY2, XAtime.derivY2 = XAtime.derivY2,
             y2 = y$y2, 
             Xs2 = Xs2, Xs.derivY2 = Xs.derivY2, XAs.derivY2 = XAs.derivY2,
             Z2 = Z2, Ztime2 = Ztime2,  Ztime.derivY2 = Ztime.derivY2, ZAtime.derivY2 = ZAtime.derivY2,
             Zs2 = Zs2, Zs.derivY2 = Zs.derivY2,  ZAs.derivY2 = ZAs.derivY2, 
             
             X3 = X3, Xtime3 = Xtime3, 
             Xtime.derivY3 = Xtime.derivY3, XAtime.derivY3 = XAtime.derivY3,
             y3 = y$y3, 
             Xs3 = Xs3, Xs.derivY3 = Xs.derivY3, XAs.derivY3 = XAs.derivY3,
             Z3 = Z3, Ztime3 = Ztime3,  Ztime.derivY3 = Ztime.derivY3, ZAtime.derivY3 = ZAtime.derivY3,
             Zs3 = Zs3, Zs.derivY3 = Zs.derivY3,  ZAs.derivY3 = ZAs.derivY3, 
             
             
             
             eventD = eventD, zeros = zeros, 
             WD = x$WD, ncZ1 = ncol(Z1), 
             ncX1 = ncol(X1), 
             ncZ2 = ncol(Z2), 
             ncX2 = ncol(X2),
             ncZ3 = ncol(Z3), 
             ncX3 = ncol(X3),
             
             
             ncWD = ncol(x$WD),  
             ncX.derivY1 = ncol(XderivY1), ncZ.derivY1 = ncol(ZderivY1),
             ncXA.derivY1 = ncol(XAderivY1), ncZA.derivY1 = ncol(ZAderivY1),
             
             ncX.derivY2 = ncol(XderivY1), ncZ.derivY2 = ncol(ZderivY1),
             ncXA.derivY2 = ncol(XAderivY1), ncZA.derivY2 = ncol(ZAderivY1),
             
             ncX.derivY3 = ncol(XderivY1), ncZ.derivY3 = ncol(ZderivY1),
             ncXA.derivY3 = ncol(XAderivY1), ncZA.derivY3 = ncol(ZAderivY1),
             threshold=threshold,
             
             W2D = W2D, 
             W2sD = W2sD, ncW2D = ncol(x$W2D), C = C, P = P,
             wk = wk,
             mu01 = mu01,mu02 = mu02, mu03 = mu03,
             priorMean.betas1 = betas1, 
             priorTau.betas1 = diag(1,3),
             
             priorMean.betas2 = betas2, 
             priorTau.betas2 = diag(1,3),
             
             priorMean.betas3 = betas3, 
             priorTau.betas3 = diag(1,3),
             
             priorA.tau1 = (1/sigma21)^2/10,
             priorB.tau1 = (1/sigma21)/10, 
             
             priorA.tau2 = (1/sigma22)^2/10,
             priorB.tau2 = (1/sigma22)/10, 
             
             priorA.tau3 = (1/sigma23)^2/10,
             priorB.tau3 = (1/sigma23)/10, 
             
             
             priorMean.gammas = gammasD,
             priorTau.gammas = diag(1/var.gammasD),
             
             priorMean.Bs.gammas = Bs.gammasD,
             priorTau.Bs.gammas = diag(1/var.Bs.gammasD),
             priorR.D1 = diag(1,3), priorK.D1 = 3,
             priorR.D2 = diag(1,3), priorK.D2 = 3,
             priorR.D3 = diag(1,3), priorK.D3 = 3,
             muV1 = muV1, stdV1 = stdV1, muS1 = muS1, stdS1 = stdS1, muA1 = muA1, stdA1 = stdA1,             
             muV2 = muV2, stdV2 = stdV2, muS2 = muS2, stdS2 = stdS2, muA2 = muA2, stdA2 = stdA2,
             muV3 = muV3, stdV3 = stdV3, muS3 = muS3, stdS3 = stdS3, muA3 = muA3, stdA3 = stdA3,
             mean1=c(means[1],rep(0,2)),SD1=SDs[1],
             mean2=c(means[2],rep(0,2)),SD2=SDs[2],
             mean3=c(means[3],rep(0,2)),SD3=SDs[3])

save(Data,file="SimData_1R2F.R")
