
library(rjags)
library(daymetr)
library(zoo)
library(ecoforecastR)
library(dplyr)
jags.out <- readRDS("data/cal.jags.out.rds")





ourmodel <- "  model{
  
#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)

#### Random Effects
# tau_alpha~dgamma(0.1,0.1)
# for(i in 1:n){                  
# alpha[i]~dnorm(0,tau_alpha)
# }

#### Fixed Effects
beta_IC~dnorm(0,0.001)
betaIntercept~dnorm(0,0.001)
betaTmin~dnorm(0,0.001)
betaPrecip~dgamma(0.01,0.01)
muTmin~dnorm(0,0.001)
muPrecip~dnorm(0,0.001)
tauTmin~dgamma(0.01,0.01)
tauPrecip~dgamma(0.01,0.01)


#### Data Model
for(t in 1:n){
OBS[t] ~ dnorm(x[t],tau_obs)
Z[t,2] ~ dnorm(muTmin,tauTmin)
Z[t,3] ~ dnorm(muPrecip,tauPrecip)
}

#### Process Model
for(t in 2:n){
mu[t] <- beta_IC*x[t-1]  + betaIntercept*Z[t,1] + betaTmin*Z[t,2] + betaPrecip*Z[t,3] 
# + alpha[t[i]]
x[t]~dnorm(mu[t],tau_add)
}

}"

ppt <- readRDS("data/climate_ensembles.rds") # created using load_ERA5 function - only runs on SCC
# view temperature ensemble
temp <- ppt[[1]]
plot(1:length(temp[1,]),temp[1,],type='n',ylim=range(temp[1,]),xlab="time",ylab="temperature (c)")
for(i in 1:10){
lines(1:length(temp[1,]),temp[i,],lwd=0.5,col="grey")
}

#view precipitation ensemble
precip <- ppt[[2]]
plot(1:length(precip[1,]),precip[1,],type='n',ylim=range(precip[1,]),xlab="time",ylab="precip (cm)")
for(i in 1:10){
lines(1:length(precip[1,]),precip[i,],lwd=0.5,col="grey")
}


real_data <- readRDS("data/calibration_model_data.rds")
STER <- real_data[real_data$siteID=="STER",]
STER <- STER[STER$dateID > "2013-05" & STER$dateID < "2014-12",]

STER$date <- as.Date(as.yearmon(STER$dateID))
STER <- pad(STER)
STER$times <- seq(1:18)
y <- STER$ratio



# settings
NT = 48
s <- 1             ## Focal site for forward simulation
Nmc = 5000         ## set number of Monte Carlo draws
ylim = c(4,30) ## set Y range on plot
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy
time = 1:63   ## total time
time1 = 1:18      ## calibration period
time2 = 19:63   ## forecast period
NS <- 10
dates <- seq(as.Date("2013-06-01"), by = "month", length.out = 63)
dates2 <- seq(as.Date("2014-12-01"), by = "month", length.out = 45)

driver_ensemble <- readRDS("data/climate_ensembles.rds")
ppt_ensemble <- driver_ensemble[[2]]
ppt.mean <- matrix(apply(ppt_ensemble,2,mean),1,NT) ## driver
temp_ensemble <- driver_ensemble[[1]]
temp.mean <- matrix(apply(temp_ensemble,2,mean),1,NT) ## driver


# our function
forecastN <- function(IC,betaIntercept,betaPrecip,betaTmin,beta_IC,ppt,tmin,Q=0,n=Nmc) {
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  N[,1] <- Nprev
  for(t in 2:NT){
    mu = beta_IC*Nprev + betaIntercept*1 + betaPrecip*ppt[,t] + betaTmin*tmin[,t]
    N[,t] <- rnorm(n,mu,Q)  
    Nprev <- N[,t]        ## update IC
  }
  return(N)
}

plot.run <- function(ylim=c(4,31)){
  sel = seq(1,18)
  plot(dates,dates,type='n',ylim=ylim, ylab="N")
  ecoforecastR::ciEnvelope(observed$date,ci[1,sel],ci[3,sel],col=col.alpha("lightBlue",0.6))
  #lines(time1,ci[2,sel],col="blue")
  points(observed$date,log(observed$ratio), pch=18)
}

#### HISTORICAL TIME-SERIES AND DETERMINISTIC FORECAST
## parameters
out <- fitModel(input.data = STER)
params <- as.matrix(out$params)
param.mean <- apply(params,2,mean)
## initial conditions
IC.orig <- as.matrix(out$predict)
NT <- 45
NS <- 10
s <- 5
ci <- apply(as.matrix(out$predict),2,quantile,c(0.025,0.5,0.975))


#observed$times <- seq(1:16)
#observed <- STER[!is.na(STER[,"ratio"]),]
N.det <- forecastN(IC=mean(IC.orig[,"x[10]"]),
                   betaIntercept=param.mean["betaIntercept"],
                   betaPrecip=param.mean["betaPrecip"],
                   betaTmin=param.mean["betaTmin"],
                   beta_IC = param.mean["beta_IC"],
                   ppt=ppt.mean,
                   tmin=temp.mean,
                   Q=0,  ## process error off
                   n=1)

observed <- STER
plot.run()
lines(dates2,N.det,col="purple",lwd=3)



val <- readRDS("data/valibration_model_data.rds")
val_STER <- val[val$siteID=="STER",]
val_date1 <- val_STER[1:16,]

STER_1 <- rbind.fill(STER, val_date1)


out <- fitModel(input.data = STER_1)
params <- as.matrix(out$params)
param.mean <- apply(params,2,mean)
## initial conditions
IC.orig <- as.matrix(out$predict)
NT <- 45
NS <- 10
s <- 5
ci <- apply(as.matrix(out$predict),2,quantile,c(0.025,0.5,0.975))


N.fcast <- forecastN(IC=mean(IC.orig[,"x[34]"]),
                   betaIntercept=param.mean["betaIntercept"],
                   betaPrecip=param.mean["betaPrecip"],
                   betaTmin=param.mean["betaTmin"],
                   beta_IC = param.mean["beta_IC"],
                   ppt=ppt.mean,
                   tmin=temp.mean,
                   Q=0,  ## process error off
                   n=1)

observed <- STER_1
plot.run()
lines(dates2,N.fcast,col="purple",lwd=3)



prow = sample.int(nrow(params),Nmc,replace=TRUE)

N.I <- forecastN(IC=IC.orig[prow,"x[10]"],  ## sample IC
                 betaIntercept=param.mean["betaIntercept"],
                 betaPrecip=param.mean["betaPrecip"],
                 betaTmin=param.mean["betaTmin"],
                 beta_IC = param.mean["beta_IC"],
                 ppt=ppt.mean,
                 tmin = temp.mean,
                 Q=0,
                 n=Nmc)

# create plot
plot.run()
N.I.ci = apply(N.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)





# Forecast Step:
#   P(Xt+1) ~ N(mXt ,m2τa2 + q)
# ✤ Analysis Step
# P(Xt+1 |Yt+1) ~ N(µa,τa2)
# ✤ 1/pa = n/pf + 1/r
# ✤ µa = (µf /pf +nY/r)·pa

#mXt = mean of current X
mXt <- log(observed$ratio)[16]

# process error
Q = 0
Q = mean(params[,"tau_add"]) #1.079387

# observation error
a = mean(params[,"tau_obs"]) # 1.680847

# current covariates?
tmin = observed$min_temp.C_avg[16]
prec = observed$precip.mm_avg[16]



params$tau_add

## log transform data
Y   = log10(observed$ratio)

## options for process model 
alpha = 0       ## assume no spatial flux
adj <- matrix(1,1,1)
#alpha = 0.05    ## assume a large spatial flux
M = adj*alpha + diag(1-alpha*apply(adj,1,sum))  ## random walk with flux

## options for process error covariance
#Q = tau_proc            ## full covariance matrix
Q = mean(params[,"tau_add"])
#mean(params[,"tau_add"])
params$tau_add
  matrix(0,1,1)
#Q = diag(diag(Q))       ## diagonal covariance matrix

## observation error covariance (assumed independent)  
R = mean(params[,"tau_obs"])
  #diag(tau_obs,nstates) 

## prior on first step, initialize with long-term mean and covariance
# > mean(Y, na.rm=TRUE)
# [1] 8.150998
mu0 = matrix(8.15,1,1)#apply(Y,1,mean,na.rm=TRUE)
y.na <- Y[!is.na(Y)]
Y[is.na(Y)] <- 8

x <- log10(observed$ratio)
x[is.na(x)] <- 8
Y <- matrix(x,x,x)
P0 = cov(t(Y),use="pairwise.complete.obs")

## Run Kalman Filter
KF00 = KalmanFilter(M,mu0,P0,Q,R,Y)
