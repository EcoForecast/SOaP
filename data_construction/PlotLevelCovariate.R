## Placing CHM data (as one covariate data), date and corresponding plots together

library(neonUtilities)
library(sp)
library(rgdal)

D=read.csv("/usr3/graduate/wangytj/GE585-Share/CPER_Soilcore_CorrespondingMeanCHM_OfEachPlot_2013.csv",header=TRUE, sep=",");
D1<-D[!duplicated(D$Plot.ID),];
siteID=matrix("CPER",dim(D1)[1],1);
dateID=matrix("2013-06",dim(D1)[1],1);
D2=cbind(siteID,D1[,c(2,9)],dateID);

write.csv(D2, file="PlotLevelCovariates_CPER.csv")
saveRDS(D2, "PlotLevelCovariates_CPER.rds")
##D3 <- readRDS("PlotLevelCovariates_CPER.rds")
##identical(D2,D3)



D=read.csv("/usr3/graduate/wangytj/GE585-Share/DSNY_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv",header=TRUE, sep=",");
D1<-D[!duplicated(D$Plot.ID),];
siteID=matrix("DSNY",dim(D1)[1],1);
dateID=matrix("2014-05",dim(D1)[1],1);
D2=cbind(siteID,D1[,c(2,9)],dateID);

write.csv(D2, file="PlotLevelCovariates_DSNY.csv")
saveRDS(D2, "PlotLevelCovariates_DSNY.rds")



D=read.csv("/usr3/graduate/wangytj/GE585-Share/HARV_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv",header=TRUE, sep=",");
D1<-D[!duplicated(D$Plot.ID),];
siteID=matrix("HARV",dim(D1)[1],1);
dateID=matrix("2014-06",dim(D1)[1],1);
D2=cbind(siteID,D1[,c(2,9)],dateID);

write.csv(D2, file="PlotLevelCovariates_HARV.csv")
saveRDS(D2, "PlotLevelCovariates_HARV.rds")




D=read.csv("/usr3/graduate/wangytj/GE585-Share/OSBS_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv",header=TRUE, sep=",");
D1<-D[!duplicated(D$Plot.ID),];
siteID=matrix("OSBS",dim(D1)[1],1);
dateID=matrix("2014-05",dim(D1)[1],1);
D2=cbind(siteID,D1[,c(2,9)],dateID);

write.csv(D2, file="PlotLevelCovariates_OSBS.csv")
saveRDS(D2, "PlotLevelCovariates_OSBS.rds")




D=read.csv("/usr3/graduate/wangytj/GE585-Share/STER_Soilcore_CorrespondingMeanCHM_OfEachPlot_2013.csv",header=TRUE, sep=",");
D1<-D[!duplicated(D$Plot.ID),];
siteID=matrix("STER",dim(D1)[1],1);
dateID=matrix("2013-06",dim(D1)[1],1);
D2=cbind(siteID,D1[,c(2,9)],dateID);

write.csv(D2, file="PlotLevelCovariates_STER.csv")
saveRDS(D2, "PlotLevelCovariates_STER.rds")





D1 <- readRDS("data/PlotLevelCovariates_CPER.rds")
D2 <- readRDS("data/PlotLevelCovariates_DSNY.rds")
D3 <- readRDS("data/PlotLevelCovariates_HARV.rds")
D4 <- readRDS("data/PlotLevelCovariates_OSBS.rds")
D5 <- readRDS("data/PlotLevelCovariates_STER.rds")

write.csv(rbind(D1,D2,D3,D4,D5), file="PlotLevelCovariates.csv")
saveRDS(rbind(D1,D2,D3,D4,D5), "PlotLevelCovariates.rds")




