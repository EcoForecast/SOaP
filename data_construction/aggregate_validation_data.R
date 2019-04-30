# combine all validation data for initial models
library(dplyr)
library(zoo)

# daymet data
if(!file.exists("data/daymet_monthly.rds")){
  source("data_construction/other_covariates/daymet_SOaP.R")
}
daymet <- readRDS("data/daymet_monthly.rds")

# ITS:16S ratios
if(!file.exists("data/validation_abundances.rds")){
  source("data_construction/microbial_data/01_download_abundance_data.R")
}
microbes <- readRDS("data/validation_abundances.rds")


df2 <- merge(daymet, microbes, all=T)
df2$log_BF_ratio <- log(df2$ratio)
df2 <- df2[df2$dateID > "2014-12",]

saveRDS(df2, "data/calibration_model_data.rds")
