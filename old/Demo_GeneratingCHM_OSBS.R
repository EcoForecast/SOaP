library(rgdal)
library(raster)
library(neonUtilities)

## range of coordinates of 1km by 1km unit areas
cd1=seq(350000,450000,by=1000)
cd2=seq(3000000,4000000,by=1000)

head=c("EcoStructure","UTM east lower left","UTM north lower left");

## downloading the data (Hint: the time needed ford download can be hours. Users are not encouraged to use the NEON
## package to download remote sensing data. It can be 3 times faster to mannually download the remote sensing data.)

byFileAOP("DP3.30015.001", site="OSBS", year="2014", check.size=T);

## 2014 data are available
for (time in 2014:2014) 
{
  for (i in cd1)
  {
    for (j in cd2)
    {
      inputpath=sprintf("DP3.30015.001/%i/FullSite/D03/%i_OSBS_1/L3/DiscreteLidar/CanopyHeightModelGtif/%i_OSBS_1_%i_%i_CHM.tif",time,time,time,i,j);
      
      # check if file exists       
      if(file.exists(inputpath))
      {
        img=raster(inputpath);
        
        value=as.matrix(img);
        
        # excluding NA pixels        
        p=which(!is.na(value));
        
        # calculate and record the mean CHM of the 1km by 1km unit area         
        if(length(p)>0){
          head=rbind(head,c(mean(value[p]),i,j));}}
    }    
  }  
}
write.csv(head,file=sprintf("/usr3/graduate/wangytj/GE585-Share/OSBS_CHM_%i.CSV",time))