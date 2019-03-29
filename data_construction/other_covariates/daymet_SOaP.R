
## call relevant packages
library(daymetr)
library(dplyr)

## I saw the code from Exercise_06 and expanded it to grab the data from all of our sites...
#LIST OF sites
#lat and longs come from: https://www.neonscience.org/field-sites/field-sites-map/list
locs <- data.frame(site = c('HARV', 'DSNY', 'OSBS', 'STER', 'CPER'),
                   lat = c(42.5369, 28.12504, 29.68927, 40.4619, 40.81553),
                   lon = c(-72.17266, -81.4362, -81.99343, -103.0293, -104.7456 ))


#WRITE DATA TO CSV FILE
write.table(locs, paste0(tempdir(),"/locations.csv"),
            sep=",",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

#This line of code makes sure you get the most recent data since daymet is updated once a year. 
end.year  <- as.numeric(substr(Sys.Date(),1,4)) - 1


#download the data here. 
daymet_data <- download_daymet_batch(file_location = paste0(tempdir(),
                                                            "/locations.csv"),
                                     start=2010,
                                     end=end.year,
                                     internal=TRUE,
                                     silent=TRUE)
#HARV can call each site 1-5
HARV <- as.matrix(daymet_data[[1]]$data)
harv_wx <- as.data.frame(HARV)

days.harv<-format(strptime(harv_wx$yday, format="%j"), format="%m-%d")
dayyear <- paste(days.harv, harv_wx$year, sep="-")
harv_wx$dayyear <-dayyear
HARVdta <- as.Date(as.character(dayyear), format="%m-%d-%Y")

head(harv_wx)
#DSNY
DSNY <- as.matrix(daymet_data[[2]]$data)
dsny_wx <- as.data.frame(DSNY)

days.dsny <- format(strptime(dsny_wx$yday, format="%j"), format="%m-%d")
dayyear <- paste(days.dsny, harv_wx$year, sep="-")
dsny_wx$dayyear <-dayyear
DSNYdta <- as.Date(as.character(dayyear), format="%m-%d-%Y")

#OSBS
OSBS <- as.matrix(daymet_data[[3]]$data)
osbs_wx <- as.data.frame(OSBS)

days.osbs <- format(strptime(osbs_wx$yday, format="%j"), format="%m-%d")
dayyear <- paste(days.osbs, osbs_wx$year, sep="-")
osbs_wx$dayyear <-dayyear
OSBSdta <- as.Date(as.character(dayyear), format="%m-%d-%Y")

#STER
STER <- as.matrix(daymet_data[[4]]$data)
ster_wx <- as.data.frame(STER)

days.ster <- format(strptime(ster_wx$yday, format="%j"), format="%m-%d")
dayyear <- paste(days.ster, ster_wx$year, sep="-")
ster_wx$dayyear <-dayyear
STERdta <- as.Date(as.character(dayyear), format="%m-%d-%Y")

#CPER
CPER <- as.matrix(daymet_data[[5]]$data)
cper_wx <- as.data.frame(CPER)

days.cper <- format(strptime(cper_wx$yday, format="%j"), format="%m-%d")
dayyear <- paste(days.ster, cper_wx$year, sep="-")
cper_wx$dayyear <-dayyear
CPERdta <- as.Date(as.character(dayyear), format="%m-%d-%Y")

#plot of min temps for each site
plot(HARVdta, harv_wx$tmin..deg.c., type='l', ylim=c(-30,25))
points(DSNYdta, dsny_wx$tmin..deg.c., type='l', col=rgb(red=1,green=0,blue=0, alpha=0.5))
points(OSBSdta, osbs_wx$tmin..deg.c., type='l', col=rgb(red=1,green=1,blue=0, alpha=0.3))
points(STERdta, ster_wx$tmin..deg.c., type='l', col=rgb(red=0,green=1,blue=1, alpha=0.5))
points(CPERdta, cper_wx$tmin..deg.c., type='l', col=rgb(red=0,green=0,blue=1, alpha=0.5))


#plot of max temps for each site
plot(HARVdta, harv_wx$tmax..deg.c., type='l', ylim=c(-20,40))
points(DSNYdta, dsny_wx$tmax..deg.c., type='l', col=rgb(red=1,green=0,blue=0, alpha=0.5))
points(OSBSdta, osbs_wx$tmax..deg.c., type='l', col=rgb(red=1,green=1,blue=0, alpha=0.3))
points(STERdta, ster_wx$tmax..deg.c., type='l', col=rgb(red=0,green=1,blue=1, alpha=0.5))
points(CPERdta, cper_wx$tmax..deg.c., type='l', col=rgb(red=0,green=0,blue=1, alpha=0.5))

#plot of precip for each site
plot(HARVdta, harv_wx$prcp..mm.day., type='l', ylim=c(0,200))
points(DSNYdta, dsny_wx$prcp..mm.day, type='l', col=rgb(red=1,green=0,blue=0, alpha=0.5))
points(OSBSdta, osbs_wx$prcp..mm.day, type='l', col=rgb(red=1,green=1,blue=0, alpha=0.3))
points(STERdta, ster_wx$prcp..mm.day, type='l', col=rgb(red=0,green=1,blue=1, alpha=0.5))
points(CPERdta, cper_wx$prcp..mm.day, type='l', col=rgb(red=0,green=0,blue=1, alpha=0.5))


harv_wx$site <- rep('HARV', length(harv_wx$dayyear))
dsny_wx$site <- rep('DSNY', length(dsny_wx$dayyear))
osbs_wx$site <- rep('OSBS', length(osbs_wx$dayyear))
ster_wx$site <- rep('STER', length(ster_wx$dayyear))
cper_wx$site <- rep('CPER', length(cper_wx$dayyear))


all_daymet <- rbind(harv_wx, dsny_wx, osbs_wx, ster_wx, cper_wx)

precip_mintemp <- cbind.data.frame(all_daymet$site, all_daymet$dayyear, all_daymet$prcp..mm.day., all_daymet$tmin..deg.c.)
colnames(precip_mintemp) <- c('siteID', 'Date', 'Precipitation.mm.day', 'Average.minimum.temperature.deg.C')

#precip_mintemp$Date <- as.Date(precip_mintemp$Date, format = "%m-%d-%Y")
precip_mintemp$dateID <- substr(precip_mintemp$Date, 1, 7)
monthly_avg <- precip_mintemp %>% 
  group_by(siteID, dateID) %>%                           
  summarise(min_temp.C_avg = mean(Average.minimum.temperature.deg.C), 
            min_temp.C_sd = sd(Average.minimum.temperature.deg.C),
            precip.mm_avg = mean(Precipitation.mm.day),
            precip.mm_sd = sd(Precipitation.mm.day))

saveRDS(monthly_avg, "data/daymet_monthly.rds")
<<<<<<< HEAD
=======

>>>>>>> 6d999e4144362ad9ee17492069d8711ef1078d50
