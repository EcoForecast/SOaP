# This script downloads and aggregates ERA5 daily estimates into monthly means (10 ensembles)
# It only runs on the SCC, and contains code written by Hamze Dokoohaki and Kathryn Wheeler

load_ERA5 <- function(lat,long,years) {
  source('/projectnb/dietzelab/hamzed/ERA5/nc_extractor_Ens.R')
  out <- ERA5_extract_ENS(lat=lat, long=long,years=years,var=c("t2m", "tp"))
  precip_output <- numeric()
  temp_output <- numeric()
  for(e in 1:10){
    t <- out[[e]]
    tDaily <- as.data.frame(xts::apply.daily(t,mean))
    tMonthly <- as.data.frame(xts::apply.monthly(t,mean))
    tMonthly <- tMonthly %>% mutate(Date=rownames(.))
    precip <- tMonthly[,2] *1000 #convert m to mm
    temp <- tMonthly[,1]-273 #convert kelvin to celsius
    precip_output <- rbind(precip_output, precip)
    temp_output <- rbind(temp_output, temp)
  }
  out <- list(temp_output, precip_output)
  return(out)
}

site = c('HARV', 'DSNY', 'OSBS', 'STER', 'CPER')
coords <- data.frame(x=c(-72.17266, -81.4362, -81.99343, -103.0293, -104.7456 ), y=c(42.5369, 28.12504, 29.68927, 40.4619, 40.81553))
coords$site <- site

# run function for 6 years (for creating calibration/observation means)
clim <- list()
for (i in 1:5){
  clim[[i]] <- load_ERA5(lat=coords$y[[i]], long=coords$x[[i]],years=c(2013, 2014, 2015, 2016, 2017, 2018)) #download using above function
}
names(clim) <- site
saveRDS(clim, "SOaP/data/climate_ensembles2.rds")

# run function for 5 years (for forecasting with full ensembles)
tMonthly <- list()
for (i in 1:5){
  tMonthly[[i]] <- load_ERA5(lat=coords$y[[i]], long=coords$x[[i]],years=c(2014, 2015, 2016, 2017, 2018)) #download using above function
saveRDS(tMonthly, "climate_ensembles.rds")