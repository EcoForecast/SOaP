# Script to download metadata of NEON ITS sequences (DP1.10108.001) 
# (next script) we'll use that metadata to access fungal sequence data from MG-RAST

# load packages and set file paths
rm(list=ls())
library(zoo)
library(neonUtilities)

# change to local file path
foldername <- "/usr3/graduate/zrwerbin/NEON_data/"

# set working directory to above file path
#setwd(foldername)
# list of our 5 sites
sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")

# first product: core-level sequence metadata. includes the following 5 products:
# mmg_soilDnaExtraction, mmg_soilMarkerGeneSequencing_16S, mmg_soilMarkerGeneSequencing_ITS, 
# mmg_soilPcrAmplification_16S, mmg_soilPcrAmplification_ITS

# check if data has already been downloaded 
if (!file.exists("filesToStack10108")) {
# loop through 5 sites to download data from each
for (s in 1:length(sites)){
zipsByProduct(dpID="DP1.10108.001", site=sites[s], package="expanded", check.size = T)
}
# combine them all into fewer files
stackByTable("filesToStack10108", folder = T)
}

# now let's format all this metadata!
pcr <- read.csv("filesToStack10108/stackedFiles/mmg_soilPcrAmplification_ITS.csv")
dna.data <- read.csv("filesToStack10108/stackedFiles/mmg_soilDnaExtraction.csv")

pcrITS <- pcr[which(pcr$targetGene=="ITS"),]

# merge the DNA data into the PCR data.
dna.merge <- dna.data[,!(colnames(dna.data) %in% colnames(pcrITS))]
dna.merge$dnaSampleID <- dna.data$dnaSampleID
dna.metadata <- merge(pcrITS,dna.merge, by="dnaSampleID", all = T)
dim(dna.metadata) # check dimensions

# remove some extra characters
dna.metadata$dateID <- substr(dna.metadata$collectDate,1,7)
dna.metadata$geneticSampleID <- as.character(dna.metadata$geneticSampleID)
dna.metadata$geneticSampleID <- substr(dna.metadata$geneticSampleID,1,nchar(dna.metadata$geneticSampleID)-4)

# save metadata
saveRDS(dna.metadata,"ITS_metadata.rds")


#get nested lsited of sites and dates sampled.
sites <- unique(dna.metadata$siteID)
sites <- sites[!is.na(sites)]
site_dates <- list()
for(i in 1:length(sites)){
  dates <- unique(dna.metadata[dna.metadata$siteID == sites[i],]$dateID)
  dates <- dates[!is.na(dates)]
  site_dates[[i]] <- dates
}
names(site_dates) <- sites

#save site_dates output.
saveRDS(site_dates,ITS_site_dates.path)
