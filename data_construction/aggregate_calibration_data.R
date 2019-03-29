# combine all calibration data for initial models
library(dplyr)
library(zoo)

# daymet data
source("data_construction/other_covariates/daymet_SOaP.R")
daymet <- readRDS("data/daymet_monthly.rds")

# ITS:16S ratios
source("data_construction/sequence_processing/01_download_abundance_data.R")
microbes <- readRDS("data/calibration_abundances.rds")

# soil phys
source("data_construction/NEON_covariates/NEON_soil_phys_iterate.R")
soil_phys <- readRDS("data/NEON_soil_phys_merge.rds")

# soil chem
source("data_construction/NEON_covariates/NEON_soil_chem_iterate.R")
soil_chem <- readRDS("data/NEON_soil_chm_merge.rds")

#### spatial covariates - not used for initial time-series ####

# climate data from WorldClim
source("data_construction/other_covariates/worldClim_SOaP.R")
worldclim <- readRDS("data/site_climate_values.rds")
worldclim$siteID <- worldclim$Site
worldclim$Site <- NULL

# #CHM - can only generate on SCC
# source("data_function/NEON_covariates/PlotLevelCovariate.R")
if(file.exists("data/MeanCHM_FiveSites_AllAreas.rds")){
  CHM <- as.data.frame(readRDS("data/MeanCHM_FiveSites_AllAreas.rds"))
  CHM <- as.data.frame(readRDS("data/MeanCHM_FiveSites_AllAreas.rds"))
  CHM$CHM <- CHM$`Mean CHM`
  CHM$`Mean CHM` <- NULL
  CHM$dateID <- NULL
}

# aggregate data by month and site

# soil physical properties
soil_phys <- soil_phys %>% 
  group_by(siteID, dateID) %>% 
  summarise(pH = mean(soilInCaClpH),
            pH_sd = sd(soilInCaClpH),
            standingWaterDepth = mean(standingWaterDepth),
            standingWaterDepth_sd = sd(standingWaterDepth),
            soilTemp = mean(soilTemp),
            soilTemp_sd = sd(soilTemp),
            litterDepth = mean(litterDepth),
            litterDepth_sd = sd(litterDepth))

# soil chemical properties
soil_chem <- soil_chem %>% 
  group_by(siteID, dateID) %>% 
  summarise(percentC = mean(organicCPercent),
            percentC_sd = sd(organicCPercent),
            CNratio = mean(CNratio),
            CNratio_sd = sd(CNratio))

df1 <- merge(soil_phys, soil_chem, all=T)
df2 <- merge(daymet, microbes) # since we don't have all=T, we will drop any dates not in the calibration abundances
df3 <- merge(df1, df2)

if(exists("CHM")){
  df4 <- merge(worldclim, CHM)
  master.df <- merge(df3, df4, all=T)
} else {
  master.df <- merge(df3, worldclim, all=T)
}

master.df$log_BF_ratio <- log(master.df$ratio)

saveRDS(master.df, "data/calibration_model_data.rds")
