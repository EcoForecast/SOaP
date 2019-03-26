#install libraries necessary to run script.
library(readr)
library(data.table)
#load data 
phys <- readRDS("data/NEON_soil_phys_merge.rds")
chm <- readRDS("data/NEON_soil_chm_merge.rds")
chm <- chm[!duplicated(chm$sampleID),] #remove duplicate samples 
length(intersect(phys$sampleID, chm$sampleID)) #should be 202 if all duplicate samples were dropped 
length(unique(chm$sampleID)) 

#merge chm and phys into one file 
chm.2 <- chm[,!colnames(chm) %in% colnames(phys)] # reduce chm to columns that AREN'T in phys 
dim(chm.2) #reduced columns by 12   
chm.2$sampleID <- chm$sampleID # add back in the identifier of x (since y has indiv id, we need to add it back into x because it was taken one line up)
neon <- merge(chm.2, phys, by = "sampleID", all = TRUE) # merge everything and keep all columns (all unique data from x, merged into y)
nrow(neon) == nrow(phys) # should be TRUE
ncol(neon) == (ncol(phys) + ncol(chm.2) -1) # should be TRUE
#want to keep 
neon.2 <- neon[,c("sampleID", "newdate", "plotID", "siteID", "geneticSampleID", "sampleID", "soilTemp", "litterDepth", "standingWaterDepth", "soilMoisture", "waterpHRatio", "caclpHRatio")]

