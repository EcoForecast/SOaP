#read in data
soil <- readRDS("data/soil_chm_phys.rds")
colnames(soil)
mean.ph.water <- aggregate(list(soil$waterpHRatio), by = list(soil$siteID, soil$newdate), "mean")
mean.ph.cacl <- aggregate(list(soil$caclpHRatio), by = list(soil$siteID, soil$newdate), "mean")
mean.litter.depth <- aggregate(list(soil$litterDepth), by = list(soil$siteID, soil$newdate), "mean")          
mean.soil.temp <- aggregate(list(soil$soilTemp), by = list(soil$siteID, soil$newdate), "mean")          
mean.water.depth <- aggregate(list(soil$standingWaterDepth), by = list(soil$siteID, soil$newdate), "mean")          
mean.soil.moisture <- aggregate(list(soil$soilMoisture), by = list(soil$siteID, soil$newdate), "mean")          

for(i in names(soil)){
  
}