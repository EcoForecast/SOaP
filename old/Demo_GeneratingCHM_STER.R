library(rgdal)
library(raster)
library(neonUtilities)

## range of coordinates of 1km by 1km unit areas
cd1=seq(650000,750000,by=1000)
cd2=seq(4000000,5000000,by=1000)

head=c("EcoStructure","UTM east lower left","UTM north lower left");

## downloading the data (Hint: the time needed ford download can be hours. Users are not encouraged to use the NEON
## package to download remote sensing data. It can be 3 times faster to mannually download the remote sensing data.)

byFileAOP("DP3.30015.001", site="STER", year="2013", check.size=T);

## 2013 data are available
for (time in 2013:2013) 
{
  for (i in cd1)
  {
    for (j in cd2)
    {
      inputpath=sprintf("DP3.30015.001/%i/FullSite/D10/%i_STER_1/L3/DiscreteLidar/CanopyHeightModelGtif/%i_STER_1_%i_%i_pit_free_CHM.tif",time,time,time,i,j);
      
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
write.csv(head,file=sprintf("/usr3/graduate/wangytj/GE585-Share/STER_CHM_%i.CSV",time))