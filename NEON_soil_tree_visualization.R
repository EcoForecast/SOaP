library("ggplot2")
library("scales")
soil_tree_merge <- readRDS("data/NEON_soil_tree_merge.rds")

for (i in 1:length(sites)){
  soil_tree_site <- soil_tree_merge[which(soil_tree_merge$siteID==sites[i]),]
  tree.time <-as.Date(soil_tree_site$date, format = "%Y-%m-%d")
  p <- ggplot(soil_tree_site, aes(x = tree.time, y = soil_tree_site$stemDistance)) + 
    geom_point(aes(colour = soil_tree_site$plotID, fill = soil_tree_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i])))
}

for (i in 1:length(sites)){
  soil_tree_site <- soil_tree_merge[which(soil_tree_merge$siteID==sites[i]),]
  tree.time <-as.Date(soil_tree_site$date, format = "%Y-%m-%d")
  p <- ggplot(soil_tree_site, aes(x = tree.time, y = soil_tree_site$stemDiameter)) + 
    geom_point(aes(colour = soil_tree_site$plotID, fill = soil_tree_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i])))
}

for (i in 1:length(sites)){
  soil_tree_site <- soil_tree_merge[which(soil_tree_merge$siteID==sites[i]),]
  tree.time <-as.Date(soil_tree_site$date, format = "%Y-%m-%d")
  p <- ggplot(soil_tree_site, aes(x = tree.time, y = soil_tree_site$height)) + 
    geom_point(aes(colour = soil_tree_site$plotID, fill = soil_tree_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i])))
}

for (i in 1:length(sites)){
  soil_tree_site <- soil_tree_merge[which(soil_tree_merge$siteID==sites[i]),]
  tree.time <-as.Date(soil_tree_site$date, format = "%Y-%m-%d")
  p <- ggplot(soil_tree_site, aes(x = tree.time, y = soil_tree_site$basalStemDiameterMsrmntHeight)) + 
    geom_point(aes(colour = soil_tree_site$plotID, fill = soil_tree_site$plotID)) + 
    scale_x_date("", labels = date_format("%b %Y"))
  print(p + ggtitle(paste0(sites[i])))
}

