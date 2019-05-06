rm(list=ls())

library(animation)
library(rjags)
library(daymetr)
library(zoo)
library(ecoforecastR)
library(dplyr)
library(padr)

source("functions/fitModel.R")

# read in driver data
all_driver_ensemble <- readRDS("data/climate_ensembles.rds") # ensemble data for forecasting: 2014-2018

# read in calibration data
cal <- readRDS("data/calibration_model_data.rds")
cal$date <- as.Date(as.yearmon(cal$dateID))
cal <- cal[cal$dateID > "2013-05" & cal$dateID <= "2014-12",] # make sure we're not including val data

# read in validation data
val <- readRDS("data/valibration_model_data.rds")
val$date <- as.Date(as.yearmon(val$dateID))

# global settings (not forecast-specific)
Nmc = 50000         ## set number of Monte Carlo draws
trans <- 0.8       ## set transparancy
N.cols <- c("black","red","green","blue","orange") ## set colors

# forecasting function
forecastN <- function(IC,betaIntercept,betaPrecip,beta_IC,ppt,Q=0,n=Nmc) {
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  N[,1] <- Nprev
  for(t in 2:NT){
    mu = beta_IC*Nprev + betaIntercept*1 + betaPrecip*ppt[,t]
    N[,t] <- rnorm(n,mu,Q)  
    Nprev <- N[,t]        ## update IC
  }
  return(N)
}

# plotting function
plot.run <- function(observed = STER, ts = 10, ylim=c(0,20)){
  ts.col <- rownames(observed[which(observed$timestep == ts),])[[1]]
  ts.col <- as.numeric(ts.col)
  in.ts <- observed[1:ts.col,]
  sel = seq(1,length(dates_fit))
  plot(dates_all,dates_all,type='n',ylim=ylim, ylab = "", xlab="")
  ecoforecastR::ciEnvelope(in.ts$date,ci[1,sel],ci[3,sel],col=col.alpha("lightBlue",0.6))
  lines(dates_fit,ci[2,sel],col="blue")
  points(in.ts$date,log(in.ts$ratio), pch=18)
  axis(2, at=seq(from = 0, to = 20, by=5), las=2)
  title(paste0("Forecast for ", site_name), xlab="time", ylab="log(Bacteria:Fungi)")
}

n.iter <- 30000
#n.iter <- 10000 # sufficient for STER, but other sites require many more iterations
#n.iter <- 170000 
#n.iter <- 300000

# subset to our site
#site_name = "DSNY"
site_name = "STER"
# site_name = "OSBS"
# site_name = "HARV"
# site_name = "CPER"

STER <- cal[cal$siteID==site_name,]

## remove leading, empty rows from our sites
if (site_name == "DSNY" | site_name == "HARV" | site_name=="CPER"){ # remove 1 leading empty row 
  STER <- STER[-c(1),]
} else if (site_name == "OSBS"){ # remove 3 leading empty rows
  STER <- STER[-c(1:3),]
}

# this section (messy, i know) creates "timestep" columns for every new observation added. this is because we have so much missing data, and want to start each forecast *from* the most recent timestep.
    # create timestep columns for calibration data
    to_align <- STER[!is.na(STER$ratio), c("ratio", "log_BF_ratio")]
    to_align$timestep <- seq(1:nrow(to_align))
    STER$timestep <- NA
    STER <- merge(STER,to_align, by=c("ratio", "log_BF_ratio"), all=T)
    STER$timestep.x <- NULL
    STER <- STER[order(STER$dateID),]
    STER <- pad(STER)
    STER$timestep <- na.locf(STER$timestep.y)
    end_2014 <- nrow(STER[STER$timestep == max(STER$timestep),]) - 1 # count the number of empty months at the end of 2014 (AKA where forecast will start)
    NT = 48 + end_2014 # months we start out forecasting (4 years = 48 months, plus some of 2014)
    STER_fit <- STER[-c((nrow(STER)-end_2014+1):nrow(STER)),]
    
    # create timestep columns for validation data
    val_STER <- val[val$siteID==site_name,]
    to_align <- val_STER[!is.na(val_STER$ratio), c("ratio", "log_BF_ratio")]
    to_align$timestep <- seq(from = max(STER$timestep)+1, to = max(STER$timestep) + nrow(to_align))
    val_STER$timestep <- NA
    val_STER <- merge(val_STER,to_align, by=c("ratio", "log_BF_ratio"), all=T)
    val_STER$timestep.x <- NULL
    val_STER <- val_STER[order(val_STER$dateID),]
    val_STER <- pad(val_STER)
    val_STER$timestep <- na.locf(val_STER$timestep.y,na.rm = F)

STER_1 <- plyr::rbind.fill(STER, val_STER)
STER_1$timestep <- na.locf(STER_1$timestep.y,na.rm = F) 
dates_all <- seq(as.Date(STER$date[1]), by = "month", length.out = (nrow(STER_1))) # all dates

# subset driver ensembles, which include all of 2014, even though we have *some* of 2014 for each site
#substract.n <- 13
driver_ensemble <- all_driver_ensemble[[which(names(all_driver_ensemble)==site_name)]]
ppt_ensemble <- driver_ensemble[[2]][,(13-end_2014):60] # subset to after last 2014 date
ppt.mean <- matrix(apply(ppt_ensemble,2,mean),1,ncol(ppt_ensemble)) ## take means per month

# list of timesteps to forecast (starting after our calibration period)
ts_loop <- (length(unique(na.omit(STER_fit$ratio)))):length(unique(na.omit(STER_1$ratio)))

ani.options(interval=.25) # speed of animation
saveGIF({
for (t in 1:length(ts_loop)){
  ts <- ts_loop[t] 
ts.col <- rownames(STER_1[which(STER_1$timestep == ts),])[[1]] # get row of new data point
ts.col <- as.numeric(ts.col)
NT <- length(dates_all) - ts.col # get number of forecast rows
ts_min1.col <- as.numeric(rownames(STER_1[which(STER_1$timestep == ts-1),])[[1]]) 
new.driver.start <- (ts.col-ts_min1.col)+1 # adjust ensemble size to match the new forecast period
driver.subset <- new.driver.start:ncol(ppt.mean)
if(t == 1){
  driver.subset <-  c(1:ncol(ppt.mean)) # at the first timestep, we already have the right size ensemble
  } else driver.subset <- new.driver.start:ncol(ppt.mean)
dates_fcast <- seq(as.Date(STER_1[(which(STER_1$timestep.y==ts)+1),]$date), by = "month", length.out = NT) # forecast period dates
dates_fit <- seq(as.Date(STER$date[1]), by = "month", length.out = ts.col) # model fit dates


out <- fitModel(input.data = STER_1[1:ts.col,], n.iter=n.iter, climate = "era5")
params <- as.matrix(out$params)
param.mean <- apply(params,2,mean)
## initial conditions
IC.orig <- as.matrix(out$predict)
ci <- apply(as.matrix(out$predict),2,quantile,c(0.025,0.5,0.975))
np <- ncol(IC.orig)

prow = sample.int(nrow(params),Nmc,replace=TRUE)
Qmc <- 1/sqrt(params[prow,"tau_add"])  ## convert from precision to standard deviation
drow = sample.int(nrow(ppt_ensemble),Nmc,replace=TRUE)

N.I <- forecastN(IC=IC.orig[prow,paste0("x[",ts.col,"]")],  ## Sample IC
                 betaIntercept=param.mean["betaIntercept"],
                 betaPrecip=param.mean["betaPrecip"],
                 beta_IC = param.mean["beta_IC"],
                 ppt=ppt.mean,
                 Q=0,
                 n=Nmc)
N.IP <- forecastN(IC=IC.orig[prow,paste0("x[",ts.col,"]")],  ## Sample IC
                  betaIntercept=params[prow,"betaIntercept"], ## Sample parameters
                  betaPrecip=params[prow,"betaPrecip"],
                  beta_IC = params[prow,"beta_IC"],
                  ppt=ppt.mean,
                  Q=0,
                  n=Nmc)


N.IPD.p <- forecastN(IC=IC.orig[prow,paste0("x[",ts.col,"]")],  ## Sample IC
                     betaIntercept=params[prow,"betaIntercept"], ## Sample parameters
                     betaPrecip=params[prow,"betaPrecip"],
                     beta_IC = params[prow,"beta_IC"],
                     ppt=ppt_ensemble[drow,],   ## Sample drivers
                     Q=0,
                     n=Nmc)


N.IPDE <- forecastN(IC=IC.orig[prow,paste0("x[",ts.col,"]")],  ## Sample IC
                    betaIntercept=params[prow,"betaIntercept"],
                    betaPrecip=params[prow,"betaPrecip"],
                    beta_IC = params[prow,"beta_IC"],
                    ppt=ppt_ensemble[drow,],   ## Sample drivers
                    Q=Qmc, ## Sample process error
                    n=Nmc)
N.I.ci = apply(N.I,2,quantile,c(0.025,0.5,0.975))
N.IP.ci = apply(N.IP,2,quantile,c(0.025,0.5,0.975), na.rm=T)
N.IPD.p.ci = apply(N.IPD.p,2,quantile,c(0.025,0.5,0.975))
N.IPDE.ci = apply(N.IPDE,2,quantile,c(0.025,0.5,0.975), na.rm=TRUE)

plot.run(observed = STER_1[1:ts.col,], ts = ts)
ecoforecastR::ciEnvelope(dates_fcast,N.IPDE.ci[1,],N.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(dates_fcast,N.IPD.p.ci[1,],N.IPD.p.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(dates_fcast,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(dates_fcast,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(dates_fcast,N.I.ci[2,],lwd=0.5)
}
}, movie.name = paste0(site_name, "_forecast.gif"))
