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
site_name = "STER"
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

#  for (t in 1:length(ts_loop)){
  t <- 1
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
    

## create variance partitioning figure, only for STER
png("rel_variance_STER.png", width     = 5.25,
    height    = 4.25,
    units     = "in",
    res       = 200,
    pointsize = 8)
### calculation of variances
varI     <- apply(N.I,2,var, na.rm=T)
varIP    <- apply(N.IP,2,var, na.rm=T)
varIPD.p   <- apply(N.IPD.p,2,var, na.rm=T)
#varIPD.t   <- apply(N.IPD.t,2,var, na.rm=T)
varIPDE  <- apply(N.IPDE,2,var, na.rm=T)
varMat   <- rbind(varI,varIP, varIPD.p,varIPDE)

## in-sample stacked area plot
V.pred.rel.in <- apply(varMat,2,function(x) {x/max(x)})
plot(dates_fcast,V.pred.rel.in[1,],ylim=c(0,1),type='n',main="Relative Variance",ylab="Proportion of Variance",xlab="time")
ciEnvelope(dates_fcast,rep(0,ncol(V.pred.rel.in)),V.pred.rel.in[1,],col=N.cols[1]) # IC error = black
ciEnvelope(dates_fcast,V.pred.rel.in[1,],V.pred.rel.in[2,],col=N.cols[2]) # parameter error = red
ciEnvelope(dates_fcast,V.pred.rel.in[2,],V.pred.rel.in[3,],col=N.cols[3]) # driver error = green = precipitation
ciEnvelope(dates_fcast,V.pred.rel.in[3,],V.pred.rel.in[4,],col=N.cols[4]) # process error = blue
legend("topright",legend=c("Process","Driver - Precip","Parameter","InitCond"),col=N.cols[4:1],lty=1,lwd=5)
dev.off()