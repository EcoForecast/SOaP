#Script completed by Zoey and Ryan

# script expects that you keep the raw data within its filesToStack10086 folder, 
# and everything saves within /data. 
# will throw errors if you delete the raw files but keep the merged output.

# script to download soil physical data, distributed periodic, dp.10086.001, from NEON
# clear environment, load packages and set file paths
rm(list=ls())
library(zoo)
library(neonUtilities)

sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")

# this script:
# 1. checks if the output file already exists
# 2. if it does, we look online and see if there's anything new and combine that with the old.
# 3. if no output file exists, download everything.
# 4. once files have been downloaded, we merge them all together.

# 1. check if output file exists
if(file.exists("data/NEON_soil_phys_merge.rds")) {
  
  # 2. if it exists, check for new files to download
  
  # read in file, if this script has been run before
  downloaded_data <- readRDS("data/NEON_soil_phys_merge.rds")
  downloaded_dates <- downloaded_data[,c("siteID", "dateID")]
  
  #connect to NEON API for DP1.10086.00.----
  req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10086.001")
  req.text <- httr::content(req, as="text")
  avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)
  
  #grab a vector of the urls to data. One per unique site-date combination.-----
  urls <- unlist(avail$data$siteCodes$availableDataUrls)
  
  all_site_dates <- substr(urls, nchar(urls)-11, nchar(urls))
  all_site_dates <- data.frame(siteID=substr(all_site_dates, 1,4),
                               dateID=substr(all_site_dates, 6,12))
  # subset to our 5 sites
  site_dates_5 <- all_site_dates[which(all_site_dates$siteID %in% sites),]
  
  # check if there are sites online that aren't in our previous download
  to_download <- site_dates_5[!do.call(paste0, site_dates_5) %in% do.call(paste0, downloaded_dates),] # have to paste the columns together to check
  
  if (nrow(to_download) > 0){ # if there ARE files in to_download
    
    dir.create("data/new_soil_phys_files")
    
    #if there are dates in to_download, let's use getPackage() on those dates 
    for (d in nrow(to_download)){
      getPackage(dpID = "DP1.10086.001", site_code = to_download[d,1], 
                 year_month = to_download[d,2], package="expanded", 
                 savepath = "data/new_soil_phys_files/")
    }
    stackByTable("data/new_soil_phys_files", folder=T)
    
    # read in new files
    soil_pH_new <- read.csv("data/new_soil_phys_files/stackedFiles/sls_soilpH.csv")
    soil_moisture_new <- read.csv("data/new_soil_phys_files/stackedFiles/sls_soilMoisture.csv")
    soil_core_new <- read.csv("data/new_soil_phys_files/stackedFiles/sls_soilCoreCollection.csv")
    
    # read in older files
    soil_pH_old <- read.csv("data/filesToStack10086/stackedFiles/sls_soilpH.csv")
    soil_moisture_old <- read.csv("data/filesToStack10086/stackedFiles/sls_soilMoisture.csv")
    soil_core_old <- read.csv("data/filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
    
    # combine them 
    if(exists("soil_pH_new")){
      soil_pH <- rbind(soil_pH_old, soil_pH_new)
    } else {
      soil_pH <- soil_pH_old
    }
    if(exists("soil_moisture_new")){
      soil_moisture <- rbind(soil_moisture_old, soil_moisture_new)
    } else {
      soil_moisture <- soil_moisture_old
    }
    if(exists("soil_core_new")){
      soil_core <- rbind(soil_core_old, soil_core_new)
    } else {
      soil_core <- soil_core_old
    }
    
    # remove the directory we created
    unlink("data/new_soil_phys_files", recursive=T)
    
  } else { # if output file exists but to_download is empty, read in older files
    
    # read in older files
    soil_pH <- read.csv("data/filesToStack10086/stackedFiles/sls_soilpH.csv")
    soil_moisture <- read.csv("data/filesToStack10086/stackedFiles/sls_soilMoisture.csv")
    soil_core <- read.csv("data/filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
    
  }
} else { #3. if no output file exists, download everything
  
  # loop through 5 sites to download data from each
  for (s in 1:length(sites)){
    zipsByProduct(dpID="DP1.10086.001", sites[s], package="expanded", check.size = FALSE, savepath = "data")
  }
  # combine them all into fewer files
  stackByTable("data/filesToStack10086", folder = T)
  
  # read in older files
  soil_pH <- read.csv("data/filesToStack10086/stackedFiles/sls_soilpH.csv")
  soil_moisture <- read.csv("data/filesToStack10086/stackedFiles/sls_soilMoisture.csv")
  soil_core <- read.csv("data/filesToStack10086/stackedFiles/sls_soilCoreCollection.csv")
  
}

###run this if there are any files that have not previously been downloaded (incl. if the script has never been run)

# 4. merge downloaded files

soil_pH.merge <- soil_pH[,!colnames(soil_pH) %in% colnames(soil_core)] # reduce soil_pH to columns that AREN'T in soil_moisture 
dim(soil_pH.merge) 
soil_pH.merge$sampleID <- soil_pH$sampleID # add back in the identifier of x (since y has indiv id, we need to add it back into x because it was taken one line up)
soil_pH_core.merge <- merge(soil_pH.merge, soil_core, by = "sampleID", all = TRUE) # merge everything and keep all columns (all unique data from x, merged into y)
dim(soil_pH_core.merge)#row and col volume
#subset soil moisture
soil_moisture.merge <- soil_moisture[,!colnames(soil_moisture) %in% colnames(soil_pH_core.merge)] # reduce a to columns that AREN'T in b (Take out everything thats in y)
dim(soil_moisture.merge)
soil_moisture.merge$sampleID <- soil_moisture$sampleID
soil_phys_merge <- merge(soil_moisture.merge, soil_pH_core.merge,  by = "sampleID", all = TRUE)
dim(soil_phys_merge)
soil_phys_merge$newdate <- as.Date(as.yearmon(soil_phys_merge$collectDate)) #make a new column with the collection date, titled "new date"
soil_phys_merge$dateID <- substr(soil_phys_merge$newdate, 1, 7)
saveRDS(soil_phys_merge, "data/NEON_soil_phys_merge.rds") # save data as a .rds file (smaller than a .csv, easier to read into R) 

#below is for chopping certain dates off (validation)
#just_calibration <- soil_phys_merge[which(soil_phys_merge$newdate < "2015-01-01"),] #get rid of everything before 01/2015
#saveRDS(just_calibration, "NEON_soil_phys_merge_dates.rds") #save as RDS 

