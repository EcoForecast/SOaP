##' Calculates weather params (daily for all 10 ensembles)
##' Initially written by Kathryn Wheeler, sources code written by Hamze Dokoohaki,
##' Further adapted by Zoey Werbin
##' 
##' @param lat The site latitude
##' @param long The site longitude
##' @param years The desired years to download
##' @import xts
##' @export
load_ERA5 <- function(lat,long,years) {
  source('/projectnb/dietzelab/hamzed/ERA5/nc_extractor_Ens.R')
  ##Could probably save this in a csv file so we don't have to process it all of the time
  out <- ERA5_extract_ENS(lat=lat, long=long,years=years,var=c("t2m", "tp"))
  #TairsOutput <- matrix(ncol=10,nrow=366*length(years))
  precip_output <- numeric()
  temp_output <- numeric()
  for(e in 1:10){
    t <- out[[e]]
    tDaily <- as.data.frame(xts::apply.daily(t,mean))
    tDaily <- tDaily %>% mutate(Date=rownames(.))
    print("length(tdaily) all")
    print(length(tDaily[,1]))
    tDaily <- as.data.frame(tDaily)
    precip <- tDaily[,2]*1000 #convert m to cm
    temp <- tDaily[,1]-273 #convert kelvin to celsius
    precip_output <- rbind(precip_output, precip)
    temp_output <- rbind(temp_output, temp)
  }
  return(list(temp_output, precip_output))
}
