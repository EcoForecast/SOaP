#Script completed by Zoey and Ryan
# script to download soil chemical data, distributed periodic, dp.10078.001, from NEON

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
if(file.exists("data/NEON_soil_chm_merge.rds")) {
  
  # 2. if it exists, check for new files to download
  
  # read in file, if this script has been run before
  downloaded_data <- readRDS("data/NEON_soil_chm_merge.rds")
  downloaded_dates <- downloaded_data[,c("siteID", "dateID")]
  
  #connect to NEON API for DP1.10078.00.----
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
    
    dir.create("data/new_soil_chm_files")
    
    #if there are dates in to_download, let's use getPackage() on those dates 
    for (d in nrow(to_download)){
      getPackage(dpID = "DP1.10078.001", site_code = to_download[d,1], 
                 year_month = to_download[d,2], package="expanded", 
                 savepath = "data/new_soil_chm_files")
    }
    stackByTable("data/new_soil_chm_files", folder=T)
    
    # read in new files
    soil_chm_new <- read.csv("data/new_soil_chm_files/stackedFiles/sls_soilChemistry.csv")
    
    # read in older files
    soil_chm_old <- read.csv("data/filesToStack10078/stackedFiles/sls_soilChemistry.csv")
    
    # combine them 
    if(exists("data/soil_chm_new")){
      soil_chm <- rbind(soil_chm_old, soil_chm_new)
    } else {
      soil_chm <- soil_chm_old
    }
    
    # remove the directory we created
    unlink("data/new_soil_chm_files", recursive=T)
    
  } else { # if output file exists but to_download is empty, read in older files
    
    # read in older files
    soil_chm <- read.csv("data/filesToStack10086/stackedFiles/sls_soilChemistry.csv")
    
  }
} else { #3. if no output file exists, download everything
  
  # loop through 5 sites to download data from each
  for (s in 1:length(sites)){
    zipsByProduct(dpID="DP1.10078.001", sites[s], package="expanded", 
                  savepath = "data", check.size = FALSE)
  }
  # combine them all into fewer files
  stackByTable("data/filesToStack10078", folder = T)
  
  # read in older files
  soil_chm <- read.csv("data/filesToStack10078/stackedFiles/sls_soilChemistry.csv")
}

###run this if there are any files that have not previously been downloaded (incl. if the script has never been run)

# 4. merge downloaded files

dim(soil_chm) 
soil_chm$sampleID <- soil_chm$sampleID # add back in the identifier of x (since y has indiv id, we need to add it back into x because it was taken one line up)
soil_chm$newdate <- as.Date(as.yearmon(soil_chm$collectDate)) #make a new column with the collection date, titled "new date"
soil_chm$dateID <- substr(soil_chm$newdate, 1, 7)
saveRDS(soil_chm, "data/NEON_soil_chm_merge.rds") # save data as a .rds file (smaller than a .csv, easier to read into R) 

#below is for chopping certain dates off (validation)
#just_calibration.chm <- soil_chm[which(soil_chm$newdate < "2015-01-01"),] #get rid of everything before 01/2015
#saveRDS(just_calibration.chm, "NEON_soil_chm_merge_dates.rds") #save as RDS 



