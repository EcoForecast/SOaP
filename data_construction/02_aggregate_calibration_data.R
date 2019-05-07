# combine all calibration data for initial models
library(dplyr)
library(zoo)

# daymet data
if(!file.exists("data/daymet_monthly.rds")){
source("data_construction/other_covariates/daymet_SOaP.R")
}
daymet <- readRDS("data/daymet_monthly.rds")

# ensemble data
ensembles <- readRDS("data/climate_ensembles2.rds")


# ITS:16S ratios
if(!file.exists("data/calibration_abundances.rds")){
  source("data_construction/microbial_data/01_download_abundance_data.R")
}
microbes <- readRDS("data/calibration_abundances.rds")

# soil phys
if(!file.exists("data/NEON_soil_phys_merge.rds")){
  source("data_construction/NEON_covariates/NEON_soil_phys_iterate.R")
}
soil_phys_raw <- readRDS("data/NEON_soil_phys_merge.rds")

# soil chem
if(!file.exists("data/NEON_soil_chm_merge.rds")){
  source("data_construction/NEON_covariates/NEON_soil_chem_iterate.R")
}
soil_chem_raw <- readRDS("data/NEON_soil_chm_merge.rds")

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
soil_phys <- soil_phys_raw[,colnames(soil_phys_raw) %in% c("siteID", "soilInCaClpH", "litterDepth", "soilTemp", "dateID")]

soil_phys <- soil_phys %>% 
  dplyr::group_by(siteID, dateID) %>% 
  dplyr::summarise_all(funs(mean, sd), na.rm = TRUE)

soil_chem <- soil_chem_raw[,colnames(soil_chem_raw) %in% c("siteID", "organicCPercent", "CNratio", "dateID")]
soil_chem <- soil_chem %>% 
  dplyr::group_by(siteID, dateID) %>% 
  dplyr::summarise_all(funs(mean, sd), na.rm = TRUE)

df1 <- merge(soil_phys, soil_chem, all=T)
df2 <- merge(daymet, microbes, all=T) # since we don't have all=T, we will drop any dates not in the calibration abundances
df3 <- merge(df1, df2, all = T)

if(exists("CHM")){
  df4 <- merge(worldclim, CHM)
  master.df <- merge(df3, df4, all=T)
} else {
  master.df <- merge(df3, worldclim, all=T)
}

master.df$log_BF_ratio <- log(master.df$ratio)

sites = c('HARV', 'DSNY', 'OSBS', 'STER', 'CPER')
site.ensemble.summary <- list()
for (i in 1:length(sites)){
  site.name <- sites[[i]]
  site.ensembles <- ensembles[[which(names(ensembles)==site.name)]]
site.temp <- site.ensembles[[1]]
site.precip <- site.ensembles[[2]]
dateIDs <- substr(x = as.character(seq(as.Date("2013-01-01"), by = "month", length.out = 72)),  1,  7)
colnames(site.temp) <- dateIDs
colnames(site.precip) <- dateIDs
site.temp.df <- t(site.temp)
temp.mean <- rowMeans(site.temp.df)
temp.sd <- apply(site.temp.df,1, sd)

site.precip.df <- t(site.precip)
precip.mean <- rowMeans(site.precip.df)
precip.sd <- apply(site.precip.df,1, sd)

climate_summary <- data.frame(siteID= rep(site.name, 72), dateID = rownames(site.temp.df), temp.mean, temp.sd, precip.mean, precip.sd)
site.ensemble.summary[[i]] <- climate_summary
}

site.ensembles <- plyr::rbind.fill(site.ensemble.summary)

master.df <- merge(site.ensembles, master.df, by = c("siteID", "dateID"))
master.df$dateID <- as.character(master.df$dateID)

saveRDS(master.df, "data/calibration_model_data.rds")
