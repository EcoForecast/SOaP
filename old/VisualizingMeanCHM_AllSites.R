library(rgdal)
library(raster)
library(neonUtilities)

##Visualizing mean CHM of CPER plots

a=read.csv("/usr3/graduate/wangytj/GE585-Share/CPER_Soilcore_CorrespondingMeanCHM_OfEachPlot_2013.csv",header = TRUE,sep=",");

hist(a[,9],xlab="mean CHM",ylab="frequency",main="Corresponding mean CHM of all plots in CPER in 2013");




##Visualizing mean CHM of DSNY plots

a=read.csv("/usr3/graduate/wangytj/GE585-Share/DSNY_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv",header = TRUE,sep=",");

hist(a[,9],xlab="mean CHM",ylab="frequency",main="Corresponding mean CHM of all plots in DSNY in 2014");




##Visualizing mean CHM of HARV plots

a=read.csv("/usr3/graduate/wangytj/GE585-Share/HARV_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv",header = TRUE,sep=",");

hist(a[,9],xlab="mean CHM",ylab="frequency",main="Corresponding mean CHM of all plots in HARV in 2014");




##Visualizing mean CHM of OSBS plots

a=read.csv("/usr3/graduate/wangytj/GE585-Share/OSBS_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv",header = TRUE,sep=",");

hist(a[,9],xlab="mean CHM",ylab="frequency",main="Corresponding mean CHM of all plots in OSBS in 2014");




##Visualizing mean CHM of STER plots

a=read.csv("/usr3/graduate/wangytj/GE585-Share/STER_Soilcore_CorrespondingMeanCHM_OfEachPlot_2013.csv",header = TRUE,sep=",");

hist(a[,9],xlab="mean CHM",ylab="frequency",main="Corresponding mean CHM of all plots in STER in 2013")
