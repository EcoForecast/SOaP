# previous script downloaded metadata for our sequences.
# now, let's match this metadata with the sequences that Lee Stanish sent via Sharefile.
# in April, data should be all accessible online somehow, and this script will be updated to reflect that

rm(list=ls())
library(neonUtilities)
library(zoo)

# read in metadata
metadata <- readRDS("ITS_metadata.rds")

# 2014-2017 data. one set of forward and reverse reads (C3D5Y); another set of only forward reads (C3D5Y). 
C3D5Y_fwd <- list.files("/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_raw_2014-2017/hpc/home/minardsmitha/NEON/16S_ITS_2018/C3D5Y_ITS_Oct_5_2018/RAW_Upload_to_BOX/R1/")
C3D5Y_rev <- list.files("/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_raw_2014-2017/hpc/home/minardsmitha/NEON/16S_ITS_2018/C3D5Y_ITS_Oct_5_2018/RAW_Upload_to_BOX/R2/")
#C3CNF_fwd <- list.files("/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_raw_2014-2017/hpc/home/minardsmitha/NEON/16S_ITS_2018/C3CNF_ITS_Sep_21_2018/RAW_Upload_to_BOX/R1/")
#all_newer_samples <- c(C3D5Y_fwd, C3CNF_fwd) # can't use the fwd-only reads right now
newer_samples <- substr(C3D5Y_fwd,1,nchar(C3D5Y_fwd)-13) # remove extra characters

# bring in the 2013-2014 
seqs_150225 <- list.files("/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_raw_ITS_fastq/ITS_run150225/r1_per_sample_demux")
seqs_150922 <- list.files("/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_raw_ITS_fastq/ITS_run150922/r1_per_sample_demux")
older_samples <- c(seqs_150225, seqs_150922)
older_samples <- substr(older_samples,1,nchar(older_samples)-10) # remove extra characters

# combine 'em
samples <- c(older_samples, newer_samples) 

# how many site-date combos samples do we have fwd/reverse reads for?
we_have <- metadata[which(metadata$internalLabID %in% samples |
                            metadata$geneticSampleID %in% samples),]
#unique(we_have[,c("siteID", "dateID")]); dim(we_have) # check it out

# # adjust dates. not sure what timezone NEON reports in; below, dates are converted to EST/EDT
# we_have$date_day <- as.POSIXct(we_have$collectDate, format = "%Y-%m-%dT%H:%MZ")
# 
# # subset to calibration/validation metadata dates
# calibration_metadata <- we_have[which(we_have$date_day < "2015-01-01 12:00:00 EST"),]
# cal <- unique(calibration_metadata[,c("siteID", "dateID")])
# 
# validation_metadata <- we_have[which(we_have$date_day > "2015-01-01 12:00:00 EST"),]
# val <- unique(validation_metadata[,c("siteID", "dateID")])
# 
# # check which samples exist in the metadata but not in our files
# list_all <-list(
#   CPER = c("2013-07", "2013-08", "2013-09",  "2013-10", "2013-11", "2013-12", "2014-01", "2014-03", "2014-04", "2014-05", "2014-06", "2014-07", "2014-08", "2014-09", "2016-04", "2016-07", "2016-11", "2017-04", "2017-07"), 
#   OSBS = c("2013-09", "2013-11", "2013-12", "2014-02", "2014-03", "2014-05", "2014-06", "2014-07", "2014-08", "2014-09", "2016-04", "2016-08", "2016-10", "2017-04", "2017-07", "2017-10"), 
#   HARV = c("2013-07", "2013-08", "2013-09", "2013-10", "2013-11", "2014-05", "2014-06", "2014-07", "2014-08", "2016-03", "2016-08", "2016-12", "2017-04", "2017-07", "2017-08", "2017-10"), 
#   STER = c("2013-06", "2013-07", "2013-08", "2013-09", "2013-10", "2013-11", "2013-12", "2014-01", "2014-03", "2014-04", "2014-05", "2014-06", "2014-07", "2014-08", "2014-09", "2016-04", "2016-08", "2016-11", "2017-04", "2017-06", "2017-10"),            
#   DSNY = c("2013-07", "2013-11", "2013-12", "2014-01", "2014-02","2014-03", "2014-04", "2014-05", "2014-06", "2014-07", "2014-08",  "2014-09", "2014-10", "2016-05", "2016-07", "2016-11", "2017-03", "2017-08"))
# list_all <- stack(list_all)
# list_all <- list_all[,c(2,1)]
# 
# missing_sample_check$source <- "we have"
# status <- merge(list_all, missing_sample_check, by.x = c("ind","values"), by.y=c("siteID", "dateID"), all=T)
# need <- status[which(is.na(status$source)),]
# need$source <- NULL
# need <- need[-c(2:3,22:23),]
# rownames(need) <- NULL
# need # print list of sampling dates that we appear to be missing