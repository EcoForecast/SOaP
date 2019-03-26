rm(list=ls())
setwed("")
dir.create("data")

# script to download soil physical data, distributed periodic, dp.10098.001, from NEON

# clear environment, load packages and set file paths
library(zoo)
library(neonUtilities)

sites <- c("DSNY", "HARV", "OSBS", "CPER")

# this script:
# 1. checks if the output file already exists
# 2. if it does, we look online and see if there's anything new and combine that with the old.
# 3. if no output file exists, download everything.
# 4. once files have been downloaded, we merge them all together.

# 1. check if output file exists
if(file.exists("data/NEON_soil_tree_merge.rds")) {
  
  # 2. if it exists, check for new files to download
  
  # read in file, if this script has been run before
  downloaded_data <- readRDS("data/NEON_soil_tree_merge.rds")
  downloaded_dates <- downloaded_data[,c("siteID", "dateID")]
  
  #connect to NEON API for DP1.10086.00.----
  req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10098.001")
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
    
    dir.create("data/new_soil_tree_files")
    
    #if there are dates in to_download, let's use getPackage() on those dates 
    for (d in nrow(to_download)){
      getPackage(dpID = "DP1.10098.001", site_code = to_download[d,1], year_month = to_download[d,2], package="expanded", savepath = "data/new_soil_tree_files/")
    }
    stackByTable("data/new_soil_tree_files", folder=T)
    
    # read in new files
    soil_maptag_new <- read.csv("data/new_soil_tree_files/stackedFiles/vst_mappingandtagging.csv")
    soil_appindiv_new <- read.csv("data/new_soil_tree_files/stackedFiles/vst_apparentindividual.csv")
    
    # read in older files
    soil_maptag_old <- read.csv("data/filesToStack10098/stackedFiles/vst_mappingandtagging.csv")
    soil_appindiv_old <- read.csv("data/filesToStack10098/stackedFiles/vst_apparentindividual.csv")
    
    # combine them 
    if(exists("soil_maptag_new")){
      soil_pH <- rbind(soil_maptag_old, soil_maptag_new)
    } else {
      soil_maptag <- soil_maptag_old
    }
    if(exists("soil_appindiv_new")){
      soil_appindiv <- rbind(soil_appindiv_old, soil_appindiv_new)
    } else {
      soil_appindiv <- soil_appindiv_old
    }
    
    # remove the directory we created
    unlink("data/new_soil_tree_files", recursive=T)
    
  } else { # if output file exists but to_download is empty, read in older files
    
    # read in older files
    soil_maptag <- read.csv("data/filesToStack10098/stackedFiles/vst_mappingandtagging.csv")
    soil_appindiv <- read.csv("data/filesToStack10098/stackedFiles/vst_apparentindividual.csv")
    
  }
} else { #3. if no output file exists, download everything
  
  # loop through 5 sites to download data from each
  for (s in 1:length(sites)){
    zipsByProduct(dpID="DP1.10098.001", sites[s], package="expanded", check.size = T, savepath = "data")
  }
  # combine them all into fewer files
  stackByTable("data/filesToStack10098", folder = T)
  
  # read in older files
  soil_maptag <- read.csv("data/filesToStack10098/stackedFiles/vst_mappingandtagging.csv")
  soil_appindiv <- read.csv("data/filesToStack10098/stackedFiles/vst_apparentindividual.csv")
  
}

###run this if there are any files that have not previously been downloaded (incl. if the script has never been run)

# 4. merge downloaded files
dim(soil_maptag)
dim(soil_appindiv)
soil_maptag.merge <- soil_maptag[,!colnames(soil_maptag) %in% colnames(soil_appindiv)] # reduce soil_pH to columns that AREN'T in soil_moisture 
dim(soil_maptag.merge) 
soil_maptag.merge$individualID <- soil_maptag$individualID # add back in the identifier of x (since y has indiv id, we need to add it back into x because it was taken one line up)
soil_tree.merge <- merge(soil_maptag.merge, soil_appindiv, by = "individualID", all = TRUE) # merge everything and keep all columns (all unique data from x, merged into y)
head(soil_tree.merge)
soil_tree.merge$newdate <- as.Date(as.yearmon(soil_tree.merge$date)) #make a new column with the collection date, titled "new date"
soil_tree.merge$dateID <- substr(soil_tree.merge$newdate, 1, 7)
saveRDS(soil_tree.merge, "data/NEON_soil_tree_merge.rds") # save data as a .rds file (smaller than a .csv, easier to read into R) 

#below is for chopping certain dates off (validation)
#just_calibration.tree <- soil_tree.merge[which(soil_tree.merge$newdate < "2015-01-01"),] #get rid of everything before 01/2015
#saveRDS(just_calibration.tree, "data/NEON_soil_tree_merge_dates.rds") #save as RDS 


