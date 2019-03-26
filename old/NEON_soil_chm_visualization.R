#Script completed by Zoey and Ryan
library("ggplot2")
library("scales")
soil_chm_merge <- readRDS("NEON_soil_chm_merge.rds")

for (i in 1:length(sites)){
  soil_chm_site <- soil_chm_merge[which(soil_chm_merge$siteID==sites[i]),]
  chm.time <-as.Date(soil_chm_site$collectDate, format = "%Y-%m-%d")
  p <- ggplot(soil_chm_site, aes(x = chm.time, y = soil_chm_site$organicCPercent)) + 
    geom_point(aes(colour = soil_chm_site$plotID, fill = soil_chm_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i], "   % Organic C")))
}

for (i in 1:length(sites)){
  soil_chm_site <- soil_chm_merge[which(soil_chm_merge$siteID==sites[i]),]
  chm.time <-as.Date(soil_chm_site$collectDate, format = "%Y-%m-%d")
  p <- ggplot(soil_chm_site, aes(x = chm.time, y = soil_chm_site$nitrogenPercent)) + 
    geom_point(aes(colour = soil_chm_site$plotID, fill = soil_chm_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i], "   % N")))
}

for (i in 1:length(sites)){
  soil_chm_site <- soil_chm_merge[which(soil_chm_merge$siteID==sites[i]),]
  chm.time <-as.Date(soil_chm_site$collectDate, format = "%Y-%m-%d")
  p <- ggplot(soil_chm_site, aes(x = chm.time, y = soil_chm_site$CNratio)) + 
    geom_point(aes(colour = soil_chm_site$plotID, fill = soil_chm_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i], "   C:N Ratio")))
}

