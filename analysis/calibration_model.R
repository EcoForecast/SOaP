### SOaP HISTORICAL TIME-SERIES FIT - looping over five NEON sites

library(ggplot2)
library(rjags)
#library(rnoaa)
library(daymetr)
library(zoo)
library("ecoforecastR")

df <- readRDS("data/calibration_model_data.rds")

### read in dataframe with response variable and covariates
df <- readRDS("data/calibration_model_data.rds")
sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")
par(mfrow=c(5,1), mar=c(1.7,3.1,1.5,1.1), mgp=c(3,.6,0))


### JAGS code adapted from the output of ecoforecastR::fit_dlm()

ourmodel <- "model{

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)

#### Random Effects
tau_alpha~dgamma(0.1,0.1)
for(i in 1:n){
alpha[i]~dnorm(0,tau_alpha)
}

#### Fixed Effects
beta_IC~dnorm(0,0.001)
betaIntercept~dnorm(0,0.001)
betaTmin~dnorm(0,0.001)
betaPrecip~dnorm(0,0.001)
betapH~dnorm(0,0.001)
betaLitter~dnorm(0,0.001)
muTmin~dnorm(0,0.001)
muPrecip~dnorm(0,0.001)
mupH~dnorm(0,0.001)
mulitter~dnorm(0,0.001) 
tauTmin~dgamma(0.01,0.01)
tauPrecip~dgamma(0.01,0.01)
taupH~dgamma(0.01,0.01)
taulitter~dgamma(0.01,0.01)

#### Data Model
for(t in 1:n){
y[t] ~ dnorm(x[t],tau_obs)
# missing data model:
Z[t,2] ~ dnorm(muTmin,tauTmin)
Z[t,3] ~ dnorm(muPrecip,tauPrecip)
Z[t,4] ~ dnorm(mupH,taupH)
Z[t,5] ~ dnorm(mulitter,taulitter)
}

#### Process Model
for(t in 2:n){
mu[t] <- beta_IC*x[t-1]  + betaIntercept*Z[t,1] + betaTmin*Z[t,2] + betaPrecip*Z[t,3] + betapH*Z[t,4] + betaLitter*Z[t,5] #+ alpha[t[i]]
x[t]~dnorm(mu[t],tau_add)
}

}"

# loop through all 5 sites
for(s in 1:length(sites)) {
    
  site.data <- df[df$siteID==sites[s],]
  y <- site.data$ratio
  
  # convert dateID to date
  time <- site.data$dateID
  time <- as.yearmon(time, format="%Y-%m")
  time <- as.Date(time, format="%Y-%m")
  
  # set up data object
  z <- cbind(rep(1,length(y)), site.data$min_temp.C_avg, site.data$precip.mm_avg,
             site.data$pH, site.data$litterDepth)
  colnames(z) <- c("(Intercept)", "Tmin", "Precip", "pH", "Litter")
  data1 <- list(y=log(y),n=length(y), x_ic = 0,tau_ic = 0.00001,a_obs=0.1,
                r_obs=0.1,a_add=0.1,r_add=0.1)
  data1[["Z"]] <- z
  
  # # set initial conditions
  # nchain = 3
  # init <- list()
  # for(i in 1:nchain){
  #   y.samp = sample(y,length(y),replace=TRUE)
  #   init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=5/var(log(y.samp)))
  # }
  
  j.model   <- jags.model (file = textConnection(ourmodel),
                           data = data1,
                          # inits = init,
                           n.chains = 3)
  
  # hide diagnostic plots for now
  # jags.out   <- coda.samples (model = j.model,
  #                             variable.names = c("tau_add","tau_obs"),
  #                             n.iter = 1000)
  # plot(jags.out)
  
  jags.out   <- coda.samples (model = j.model,
                              variable.names = c("x","tau_add",
                              "tau_obs", "beta_IC", "betaIntercept", 
                              "betaTmin", "betaPrecip", "betapH", "betaLitter", 
                              "tau_alpha","alpha"),
                              n.iter = 10000)
  
  # plot fit with confidence interval
  time.rng = c(1,length(time)) ## adjust to zoom in and out
  out <- as.matrix(jags.out)
  x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
  ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
  
  plot(time, ci[2,], type='n', ylim=range(y,na.rm=TRUE), ylab="", 
       log='y', xlim=time[time.rng], xaxt="n", 
       main=paste0("Time-series for ",sites[s]))
  mtext(text = "Bacteria+Archaea : Fungi", side=2, line=2,cex=0.6)
  
  ## adjust x-axis labels
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
  ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
  points(time,y,pch="+",cex=0.5)
}

