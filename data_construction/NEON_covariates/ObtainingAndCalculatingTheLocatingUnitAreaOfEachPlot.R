## When running the script, please run the line starting with "zipsByProduct" function first, type "y" 
## in the console, and then run the following codes for step 1 to 5.

## The UTM northing and easting of the lower left corner of the 1km by 1km unit areas where plots of each site
## are located will be stored.

## Each plot is followed by the coordinates of the UTM northing and easting of the lower left corner of the unit area 
## it is in and the corresponding mean CHM. 

library(neonUtilities)
library(sp)
library(rgdal)

## Coordinate conversion from latitude-longitude to UTM
## Run this function before running following commands
LatitudeLongitudeToUTM<-function(x,y,ID,Zone){
  xy <- data.frame(ID = 1:length(y), X=x,Y=y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## default projection is considered to be WGS84
  
  Res <- spTransform(xy, CRS(paste("+proj=utm +zone=",Zone,"ellps=WGS84",sep='')))
  return(as.data.frame(Res));
}


## 1. Getting the spatial coordinates of plots in CPER 
zipsByProduct(dpID="DP1.10086.001", site="CPER", package="expanded", check.size = T)
stackByTable("filesToStack10086", folder = T)
soilcore <- read.csv("filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
dim(soilcore) 
# grab lat/lon, coordinate uncertainty, elevation, standingWaterDepth, litterDepth, etc...
soilcore.df <- soilcore[,c(1:6,8:9,11:15,17,19,23:26,31:33,43)]
Head=head(soilcore.df);
coord_CPER=soilcore.df[,9:10];
UTM_CPER_Soilcore <- LatitudeLongitudeToUTM(coord_CPER[,2],coord_CPER[,1],1:6,13);  ##CPER is in the UTM zone of 13

## Calculating corresponding 1km by 1km unit area
UTM_CPER_CorrespondingUnitArea <- cbind(floor(UTM_CPER_Soilcore[,2]/1000)*1000,floor(UTM_CPER_Soilcore[,3]/1000)*1000);

## Finding the corresponding mean CHM value
D=read.csv("/usr3/graduate/wangytj/GE585-Share/CPER_CHM_2013.csv",header=TRUE, sep=",");
CHM_CPER=D[2:dim(D)[1],2:dim(D)[2]];
CHM_Value=matrix(0,dim(UTM_CPER_CorrespondingUnitArea)[1],1);
for (i in 1:dim(UTM_CPER_CorrespondingUnitArea)[1])
{
  m=intersect(which(CHM_CPER[,2]==UTM_CPER_CorrespondingUnitArea[i,1]),which(CHM_CPER[,3]==UTM_CPER_CorrespondingUnitArea[i,2]));
  CHM_Value[i]=as.numeric(as.character(CHM_CPER[m,1]));
}
name=c("Plot ID","Long","Lat",c("UTM east","UTM north"),c("UTM east of lower left","UTM north of lower left"),"Mean CHM");
Output=cbind(soilcore.df$plotID,coord_CPER[,2],coord_CPER[,1],UTM_CPER_Soilcore[,2:3],UTM_CPER_CorrespondingUnitArea,CHM_Value);
colnames(Output) <- name
write.csv(Output,file="/usr3/graduate/wangytj/GE585-Share/CPER_Soilcore_CorrespondingMeanCHM_OfEachPlot_2013.csv");






## 2. Getting the spatial coordinates of plots in DSNY
zipsByProduct(dpID="DP1.10086.001", site="DSNY", package="expanded", check.size = T)
stackByTable("filesToStack10086", folder = T)
soilcore <- read.csv("filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
dim(soilcore) 
# grab lat/lon, coordinate uncertainty, elevation, standingWaterDepth, litterDepth, etc...
soilcore.df <- soilcore[,c(1:6,8:9,11:15,17,19,23:26,31:33,43)]
Head=head(soilcore.df);
coord_DSNY=soilcore.df[,9:10];
UTM_DSNY_Soilcore <- LatitudeLongitudeToUTM(coord_DSNY[,2],coord_DSNY[,1],1:6,17);  ##DSNY is in the UTM zone of 17

## Calculating corresponding 1km by 1km unit area
UTM_DSNY_CorrespondingUnitArea <- cbind(floor(UTM_DSNY_Soilcore[,2]/1000)*1000,floor(UTM_DSNY_Soilcore[,3]/1000)*1000);

## Finding the corresponding mean CHM value
D=read.csv("/usr3/graduate/wangytj/GE585-Share/DSNY_CHM_2014.csv",header=TRUE, sep=",");
CHM_DSNY=D[2:dim(D)[1],2:dim(D)[2]];
CHM_Value=matrix(0,dim(UTM_DSNY_CorrespondingUnitArea)[1],1);
for (i in 1:dim(UTM_DSNY_CorrespondingUnitArea)[1])
{
  m=intersect(which(CHM_DSNY[,2]==UTM_DSNY_CorrespondingUnitArea[i,1]),which(CHM_DSNY[,3]==UTM_DSNY_CorrespondingUnitArea[i,2]));
  CHM_Value[i]=as.numeric(as.character(CHM_DSNY[m,1]));
}
name=c("Plot ID","Long","Lat",c("UTM east","UTM north"),c("UTM east of lower left","UTM north of lower left"),"Mean CHM");
Output=cbind(soilcore.df$plotID,coord_DSNY[,2],coord_DSNY[,1],UTM_DSNY_Soilcore[,2:3],UTM_DSNY_CorrespondingUnitArea,CHM_Value);
colnames(Output) <- name
write.csv(Output,file="/usr3/graduate/wangytj/GE585-Share/DSNY_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv")






## 3. Getting the spatial coordinates of plots in HARV
zipsByProduct(dpID="DP1.10086.001", site="HARV", package="expanded", check.size = T)
stackByTable("filesToStack10086", folder = T)
soilcore <- read.csv("filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
dim(soilcore) 
# grab lat/lon, coordinate uncertainty, elevation, standingWaterDepth, litterDepth, etc...
soilcore.df <- soilcore[,c(1:6,8:9,11:15,17,19,23:26,31:33,43)]
Head=head(soilcore.df);
coord_HARV=soilcore.df[,9:10];
UTM_HARV_Soilcore <- LatitudeLongitudeToUTM(coord_HARV[,2],coord_HARV[,1],1:6,18);  ##HARV is in the UTM zone of 18

## Calculating corresponding 1km by 1km unit area
UTM_HARV_CorrespondingUnitArea <- cbind(floor(UTM_HARV_Soilcore[,2]/1000)*1000,floor(UTM_HARV_Soilcore[,3]/1000)*1000);

## Finding the corresponding mean CHM value
D=read.csv("/usr3/graduate/wangytj/GE585-Share/HARV_CHM_2014.csv",header=TRUE, sep=",");
CHM_HARV=D[2:dim(D)[1],2:dim(D)[2]];
CHM_Value=matrix(0,dim(UTM_HARV_CorrespondingUnitArea)[1],1);
for (i in 1:dim(UTM_HARV_CorrespondingUnitArea)[1])
{
  m=intersect(which(CHM_HARV[,2]==UTM_HARV_CorrespondingUnitArea[i,1]),which(CHM_HARV[,3]==UTM_HARV_CorrespondingUnitArea[i,2]));
  CHM_Value[i]=as.numeric(as.character(CHM_HARV[m,1]));
}
name=c("Plot ID","Long","Lat","UTM east","UTM north","UTM east of lower left","UTM north of lower left","Mean CHM");
Output=cbind(soilcore.df$plotID,coord_HARV[,2],coord_HARV[,1],UTM_HARV_Soilcore[,2:3],UTM_HARV_CorrespondingUnitArea,CHM_Value);
colnames(Output) <- name
write.csv(Output,file="/usr3/graduate/wangytj/GE585-Share/HARV_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv")




## 4. Getting the spatial coordinates of plots in OSBS
zipsByProduct(dpID="DP1.10086.001", site="OSBS", package="expanded", check.size = T)
stackByTable("filesToStack10086", folder = T)
soilcore <- read.csv("filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
dim(soilcore) 
# grab lat/lon, coordinate uncertainty, elevation, standingWaterDepth, litterDepth, etc...
soilcore.df <- soilcore[,c(1:6,8:9,11:15,17,19,23:26,31:33,43)]
Head=head(soilcore.df);
coord_OSBS=soilcore.df[,9:10];
UTM_OSBS_Soilcore <- LatitudeLongitudeToUTM(coord_OSBS[,2],coord_OSBS[,1],1:6,17);  ##OSBS is in the UTM zone of 17

## Calculating corresponding 1km by 1km unit area
UTM_OSBS_CorrespondingUnitArea <- cbind(floor(UTM_OSBS_Soilcore[,2]/1000)*1000,floor(UTM_OSBS_Soilcore[,3]/1000)*1000);

## Finding the corresponding mean CHM value
D=read.csv("/usr3/graduate/wangytj/GE585-Share/OSBS_CHM_2014.csv",header=TRUE, sep=",");
CHM_OSBS=D[2:dim(D)[1],2:dim(D)[2]];
CHM_Value=matrix(0,dim(UTM_OSBS_CorrespondingUnitArea)[1],1);
for (i in 1:dim(UTM_OSBS_CorrespondingUnitArea)[1])
{
  m=intersect(which(CHM_OSBS[,2]==UTM_OSBS_CorrespondingUnitArea[i,1]),which(CHM_OSBS[,3]==UTM_OSBS_CorrespondingUnitArea[i,2]));
  CHM_Value[i]=as.numeric(as.character(CHM_OSBS[m,1]));
}
name=c("Plot ID","Long","Lat","UTM east","UTM north","UTM east of lower left","UTM north of lower left","Mean CHM");
Output=cbind(soilcore.df$plotID,coord_OSBS[,2],coord_OSBS[,1],UTM_OSBS_Soilcore[,2:3],UTM_OSBS_CorrespondingUnitArea,CHM_Value);
colnames(Output) <- name
write.csv(Output,file="/usr3/graduate/wangytj/GE585-Share/OSBS_Soilcore_CorrespondingMeanCHM_OfEachPlot_2014.csv")





## 5. Getting the spatial coordinates of plots in STER
zipsByProduct(dpID="DP1.10086.001", site="STER", package="expanded", check.size = T)
stackByTable("filesToStack10086", folder = T)
soilcore <- read.csv("filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
dim(soilcore) 
# grab lat/lon, coordinate uncertainty, elevation, standingWaterDepth, litterDepth, etc...
soilcore.df <- soilcore[,c(1:6,8:9,11:15,17,19,23:26,31:33,43)]
Head=head(soilcore.df);
coord_STER=soilcore.df[,9:10];
UTM_STER_Soilcore <- LatitudeLongitudeToUTM(coord_STER[,2],coord_STER[,1],1:6,13);  ##STER is in the UTM zone of 13

## Calculating corresponding 1km by 1km unit area
UTM_STER_CorrespondingUnitArea <- cbind(floor(UTM_STER_Soilcore[,2]/1000)*1000,floor(UTM_STER_Soilcore[,3]/1000)*1000);

## Finding the corresponding mean CHM value
D=read.csv("/usr3/graduate/wangytj/GE585-Share/STER_CHM_2013.csv",header=TRUE, sep=",");
CHM_STER=D[2:dim(D)[1],2:dim(D)[2]];
CHM_Value=matrix(0,dim(UTM_STER_CorrespondingUnitArea)[1],1);
for (i in 1:dim(UTM_STER_CorrespondingUnitArea)[1])
{
  m=intersect(which(CHM_STER[,2]==UTM_STER_CorrespondingUnitArea[i,1]),which(CHM_STER[,3]==UTM_STER_CorrespondingUnitArea[i,2]));
  CHM_Value[i]=as.numeric(as.character(CHM_STER[m,1]));
}
name=c("Plot ID","Long","Lat","UTM east","UTM north","UTM east of lower left","UTM north of lower left","Mean CHM");
Output=cbind(soilcore.df$plotID,coord_STER[,2],coord_STER[,1],UTM_STER_Soilcore[,2:3],UTM_STER_CorrespondingUnitArea,CHM_Value);
colnames(Output) <- name
write.csv(Output,file="/usr3/graduate/wangytj/GE585-Share/STER_Soilcore_CorrespondingMeanCHM_OfEachPlot_2013.csv")



## To visualize the histograms of mean CHM for each site, open and run "VisualizingMeanCHM_AllSites.R"


