#Worldclim
#call packages
library(raster)
library(sp)

#Define sites, variables of interest, coordinates of sites
#lat and longs come from: https://www.neonscience.org/field-sites/field-sites-map/list
site = c('HARV', 'DSNY', 'OSBS', 'STER', 'CPER')
r <- getData('worldclim', var='bio',res=10)
r <- r[[c(1,12)]]
names(r) <- c("MAT","MAP")

coords <- data.frame(x=c(-72.17266, -81.4362, -81.99343, -103.0293, -104.7456 ), y=c(42.5369, 28.12504, 29.68927, 40.4619, 40.81553))

#Grap the data and put it in 'values'
points <- SpatialPoints(coords, proj4string=r@crs)
values <- extract(r,points)

#combine climate data with the site names. 
site_climate_values <- as.data.frame( cbind(values, site))
colnames(site_climate_values) <- c('MAT.degree.C.multiplied.by.ten', 'MAP.mm', 'Site')

saveRDS(site_climate_values, "data/site_climate_values.rds")