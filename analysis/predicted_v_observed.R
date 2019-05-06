## data exploration
# 
# 
# # read in calibration data
# cal <- readRDS("data/calibration_model_data.rds")
# cal$date <- as.Date(as.yearmon(cal$dateID))
# cal <- cal[cal$dateID > "2013-05" & cal$dateID <= "2014-12",]
# 
# # read in validation data
# val <- readRDS("data/valibration_model_data.rds")
# val$date <- as.Date(as.yearmon(val$dateID))
# 
# raw_abun <- readRDS("data/raw_abundances.rds")
# 
# abun <- raw_abun[raw_abun$ratio >= 1,]
# abun$collectDate <- as.Date(abun$collectDate)
# 


library(animation)
library(rjags)
library(daymetr)
library(zoo)
library(ecoforecastR)
library(dplyr)
library(padr)

source("functions/fitModel.R")

# read in calibration data
cal <- readRDS("data/calibration_model_data.rds")
cal$date <- as.Date(as.yearmon(cal$dateID))
cal <- cal[cal$dateID > "2013-05" & cal$dateID <= "2014-12",]

# read in validation data
val <- readRDS("data/valibration_model_data.rds")
val$date <- as.Date(as.yearmon(val$dateID))

# read in driver data
ensembles <- readRDS("data/climate_ensembles2.rds")
all_driver_ensemble <- readRDS("data/climate_ensembles.rds")

pred_obs_allsites <- list()
sites <- c("DSNY", "STER", "OSBS", "HARV", "CPER")
for (s in 1:5){
  site_name <- sites[s]
# subset to our site
#site_name = "DSNY"
#site_name = "STER"
# site_name = "OSBS"
# site_name = "HARV"
# site_name = "CPER"

STER <- cal[cal$siteID==site_name,]


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

#n.iter <- 30000
n.iter <- 170000
#n.iter <- 300000

if (site_name == "DSNY" | site_name == "HARV" | site_name=="CPER"){ # remove leading empty row from DSNY
  STER <- STER[-c(1),]
} else if (site_name == "OSBS"){ # remove leading empty row from DSNY
  STER <- STER[-c(1:3),]
}

# create timestep columns
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

STER_fcast <- STER[-c((nrow(STER)-end_2014+1):nrow(STER)),]


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

substract.n <- 13
driver_ensemble <- all_driver_ensemble[[which(names(all_driver_ensemble)==site_name)]]
ppt_ensemble <- driver_ensemble[[2]][,(substract.n-end_2014):60] # subset to after last 2014 date
ppt.mean <- matrix(apply(ppt_ensemble,2,mean),1,ncol(ppt_ensemble)) ## take means per month

ts_loop <- (length(unique(na.omit(STER_fcast$ratio)))):length(unique(na.omit(STER_1$ratio)))

  t <- length(ts_loop)
  t <- 2
    ts <- ts_loop[t] 
    ts.col <- rownames(STER_1[which(STER_1$timestep == ts),])[[1]]
    ts.col <- as.numeric(ts.col)
    NT <- 67-ts.col
    NT <- length(dates_all) - ts.col
    ts_min1.col <- as.numeric(rownames(STER_1[which(STER_1$timestep == ts-1),])[[1]])
    new.driver.start <- (ts.col-ts_min1.col)+1
    driver.subset <- new.driver.start:ncol(ppt.mean)
    if(t == 1){
      driver.subset <-  c(1:ncol(ppt.mean)) 
    } else driver.subset <- new.driver.start:ncol(ppt.mean)
    dates_fcast <- seq(as.Date(STER_1[(which(STER_1$timestep.y==ts)+1),]$date), by = "month", length.out = NT) # forecast period dates
    dates_fit <- seq(as.Date(STER$date[1]), by = "month", length.out = ts.col) # model fit dates
    

    out <- fitModel(input.data = STER_1[1:ts.col,], n.iter=n.iter, climate = "era5")
    params <- as.matrix(out$params)
    param.mean <- apply(params,2,mean)
    IC.orig <- as.matrix(out$predict)
    ci <- apply(as.matrix(out$predict),2,quantile,c(0.025,0.5,0.975))
    np <- ncol(IC.orig)
    
    prow = sample.int(nrow(params),Nmc,replace=TRUE)
    Qmc <- 1/sqrt(params[prow,"tau_add"])  ## convert from precision to standard deviation
    drow = sample.int(nrow(ppt_ensemble),Nmc,replace=TRUE)
    
    N.IPDE <- forecastN(IC=IC.orig[prow,paste0("x[",ts.col,"]")],  ## sample IC
                        betaIntercept=params[prow,"betaIntercept"],
                        betaPrecip=params[prow,"betaPrecip"],
                        beta_IC = params[prow,"beta_IC"],
                        ppt=ppt_ensemble[drow,],   ## Sample drivers
                        Q=Qmc,
                        n=Nmc)

    N.IPDE.ci = apply(N.IPDE,2,quantile,c(0.025,0.5,0.975), na.rm=TRUE)
    
    STER_pred <- STER_1[(ts.col+1):nrow(STER_1),]
    pred <- N.IPDE.ci[2,]
    pred <- N.IPDE.ci
    pred <- pred[,which(!is.na(STER_pred$log_BF_ratio))]
    obs <- unique(na.omit(STER_pred$log_BF_ratio))
    plot.df <- data.frame(siteID = rep(site_name, length(obs)), obs = obs, pred = pred[2,])
    pred_obs_allsites[[s]] <- plot.df
}

pred_obs_allsites

pred_obs <- plyr::rbind.fill(pred_obs_allsites)

png("predicted_v_observed.png", width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 200,
    pointsize = 7)
plot(pred_obs$pred, pred_obs$obs, pch = 16, cex=.7, xlim = c(3,15), ylim=c(3,15), xlab = "", ylab="")
fit <- lm(pred_obs$obs ~ pred_obs$pred)
abline(fit, col="red")
abline(1,1, lwd=3)
title(paste0("observed ~ predicted, all sites"), xlab="predicted", ylab="observed")
dev.off()
