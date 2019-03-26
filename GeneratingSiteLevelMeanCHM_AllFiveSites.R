## Calculating mean CHM for each of the five sites

library(neonUtilities)
library(sp)
library(rgdal)

## calculating mean CHM of all 1km by 1km areas with soil cores for each of the five sites
D1 <- readRDS("PlotLevelCovariates_CPER.rds")
D2 <- readRDS("PlotLevelCovariates_DSNY.rds")
D3 <- readRDS("PlotLevelCovariates_HARV.rds")
D4 <- readRDS("PlotLevelCovariates_OSBS.rds")
D5 <- readRDS("PlotLevelCovariates_STER.rds")

name=t(c("siteID","Mean CHM","dateID"));
result=matrix(0,5,3);
result[,3]=c("2013-06","2014-05","2014-06","2014-05","2013-06");
result[,2]=c(mean(D1[,3]),mean(D2[,3]),mean(D3[,3]),mean(D4[,3]),mean(D5[,3]));
result[,1]=c("CPER","DSNY","HARV","OSBS","STER");
colnames(result) <- name

write.csv(result,file="MeanCHM_FiveSites_OnlyIncludingAreasWithSoilCores.csv");
saveRDS(result, "MeanCHM_FiveSites_OnlyIncludingAreasWithSoilCores.rds");



## calculating mean CHM of all 1km by 1km areas for each of the five sites
D1 <- read.csv("/usr3/graduate/wangytj/GE585-Share/CPER_CHM_2013.csv")
D2 <- read.csv("/usr3/graduate/wangytj/GE585-Share/DSNY_CHM_2014.csv")
D3 <- read.csv("/usr3/graduate/wangytj/GE585-Share/HARV_CHM_2014.csv")
D4 <- read.csv("/usr3/graduate/wangytj/GE585-Share/OSBS_CHM_2014.csv")
D5 <- read.csv("/usr3/graduate/wangytj/GE585-Share/STER_CHM_2013.csv")

name=t(c("siteID","Mean CHM","dateID"));
result=matrix(0,5,3);
result[,3]=c("2013-06","2014-05","2014-06","2014-05","2013-06");
result[,2]=c(mean(as.numeric(as.character(D1[2:dim(D1)[1],2]))),mean(as.numeric(as.character(D2[2:dim(D2)[1],2]))),mean(as.numeric(as.character(D3[2:dim(D3)[1],2]))),mean(as.numeric(as.character(D4[2:dim(D4)[1],2]))),mean(as.numeric(as.character(D5[2:dim(D5)[1],2]))));
result[,1]=c("CPER","DSNY","HARV","OSBS","STER");
colnames(result) <- name

write.csv(result,file="MeanCHM_FiveSites_AllAreas.csv");
saveRDS(result, "MeanCHM_FiveSites_AllAreas.rds")
