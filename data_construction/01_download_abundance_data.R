# download and process microbial group abundances from 2013-2017
# outputs the site-level average values of each group 

library(neonUtilities)
library(plyr)
library(dplyr)
library(data.table)
library(zoo)
library(ggplot2)

# do you need to download the data?
needDownload <- F

if (needDownload == T) {
  # download data
  sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")
  
  for (s in 1:5){
    zipsByProduct(dpID="DP1.10109.001", site=sites[s], package="expanded", check.size = F, savepath = "data")
  }
  # combine data
  stackByTable(filepath = "data/filesToStack10109", folder=T)
  
}

# read in from csv
df_raw <- read.csv("data/filesToStack10109/stackedFiles/mga_soilGroupAbundances.csv")

# subset to columns of interest
df <- df_raw[,c("plotID", "siteID", "collectDate", "targetTaxonGroup", "meanCopyNumber", 
                "copyNumberStandardDeviation", "nucleicAcidConcentration", "geneticSampleID","processedDate", "dnaSampleID", "dataQF")]

df <- df[which(df$dataQF!="<LLOQ"),]
# first lets get the bacteria/archaea columns fixed.
# sometimes they were processed together, other times separately.

# JUST the samples where bac and arc were analyzed separately - combine values but keep SD columns
df_16S_raw <- df[df$targetTaxonGroup== "archaea" | df$targetTaxonGroup=="bacteria",]

df_16S <- df_16S_raw %>% 
  dplyr::group_by(dnaSampleID) %>%                           
  dplyr::summarise(bacteriaAndArchaea = sum(meanCopyNumber))
arc_sd <- df[df$targetTaxonGroup== "archaea",c("dnaSampleID", "copyNumberStandardDeviation")]
bac_sd <- df[df$targetTaxonGroup== "bacteria",c("dnaSampleID", "copyNumberStandardDeviation")]
setnames(arc_sd, "copyNumberStandardDeviation", "archaea_sd")
setnames(bac_sd, "copyNumberStandardDeviation", "bacteria_sd")
df_16S <- merge(merge(arc_sd, bac_sd), df_16S)

# now the samples where bacteria and archaea were analyzed together
df_16S_together <- df[which(df$targetTaxonGroup== "bacteria and archaea"),
                      c("dnaSampleID", "meanCopyNumber", "copyNumberStandardDeviation")]
names(df_16S_together) <- c("dnaSampleID", "bacteriaAndArchaea", "bacteriaArchaea_sd")  
df_16S_all_16S_cases <- rbind.fill(df_16S, df_16S_together)

# extract ITS rows
df_ITS <- df[df$targetTaxonGroup== "fungi",]
setnames(df_ITS, c("meanCopyNumber", "copyNumberStandardDeviation"), c("fungi", "fungi_sd"))

# merge ITS with 16S
df_merged <- merge(df_ITS, df_16S_all_16S_cases)
df_merged$targetTaxonGroup <- NULL

# add dateID column
df_merged$dateID <- substr(df_merged$collectDate, 1, 7)
unique(df_merged$dateID)

df_merged$fungi_old <- df_merged$fungi
df_merged$bacteriaAndArchaea_old <- df_merged$bacteriaAndArchaea

# add one so that none of our division causes infinite values
df_merged$fungi <- df_merged$fungi + 1
df_merged$bacteriaAndArchaea <- df_merged$bacteriaAndArchaea + 1

# remove any rows with no copy numbers for either group
df_merged <- df_merged[-which(df_merged$fungi == 1 | df_merged$bacteriaAndArchaea == 1),]

# calculate ratio of bacteria and archaea to fungi
df_merged$ratio <-  df_merged$bacteriaAndArchaea / df_merged$fungi
#df_mergeddateID<−as.yearmon(dfmergeddateID)
rownames(df_merged) <- df_merged$dnaSampleID
df_merged <- df_merged[, !names(df_merged) %in% c("nucleicAcidConcentration","dnaSampleID","geneticSampleID")]
saveRDS(df_merged, "data/raw_abundances.rds")



##### METHOD 1 - average ratios from each sample/date, to the site level
# df_method1 <- aggregate(list(ratio = df_merged$ratio), 
#                         list(siteID = df_merged$siteID, dateID = df_merged$dateID), 
#                         mean, na.action = na.pass)

df_aggregate <- df_merged %>% 
  dplyr::group_by(siteID, dateID) %>%                           
  dplyr::summarise(ratio = mean(ratio))

df_aggregate$date <- as.Date(as.yearmon(df_aggregate$dateID))

# visualize data
ggplot(data = df_aggregate, aes(x = date, y = log(ratio), color = siteID)) +       
  geom_line(aes(group = siteID)) + geom_point() + ggtitle("Calibration period") 


# SUBSET to calibration data
calibration <- df_aggregate[which(df_aggregate$date < "2015-01-01"),]
validation <- df_aggregate[which(df_aggregate$date >= "2015-01-01"),]

ggplot(data = calibration, aes(x = date, y = log(ratio), color = siteID)) +       
  geom_line(aes(group = siteID)) + geom_point() + ggtitle("Calibration period") 

saveRDS(calibration, "data/calibration_abundances.rds")
saveRDS(validation, "data/validation_abundances.rds")

## CREATE spline figure

library("ggformula")

df_merged$collectDate <- as.Date(df_merged$collectDate)
df_merged$y <- log(df_merged$ratio) 
df_merged$x <- as.numeric(df_merged$collectDate)
df_merged$x <-df_merged$collectDate
abun <- df_merged
abun$collectDate <- as.Date(abun$collectDate)
abun <- abun[abun$ratio >= 1,]

cals <- abun[abun$collectDate < "2015-01-01",]
vals <- abun[abun$collectDate > "2015-01-01",]
# p <- ggplot(abun, aes(x = as.numeric(collectDate), y = log(ratio)))+ geom_point() + ggtitle("Calibration period") 
#  
#  #p + geom_spline(aes(x,y),  colour="blue", spar=0.85, lwd=3)
#  p + geom_spline(data = cals, aes(x,y),  colour="blue", spar=.85, lwd=3) + geom_spline(data = vals, aes(x,y),  colour="blue", spar=.85, lwd=3) 
#  
#  
 p <- ggplot(abun, aes(x = collectDate, y = log(ratio)))+ geom_point() + ggtitle("Log B:F ratios from all samples") 
 
 #p + geom_spline(aes(x,y),  colour="blue", spar=0.85, lwd=3)
 p + geom_spline(data = cals, aes(x,y),  colour="blue", spar=.85, lwd=3) + geom_spline(data = vals, aes(x,y),  colour="blue", spar=.85, lwd=2) + scale_x_date(date_labels = "%b-%Y")
 ggsave("spline.png")
# ##### METHOD 2 - actually this is the same thing
# 
# df_method2 <- df_merged
# df_method2Total<−dfmethod2fungi + df_method2$bacteriaAndArchaea
# 
# df_method2relfungi<−dfmethod2fungi/df_method2$Total
# df_method2relbacarc<−dfmethod2bacteriaAndArchaea/df_method2$Total
# 
# 
# df_method2 <- aggregate(list(rel_fungi = df_method2relfungi,relbacarc=dfmethod2rel_bac_arc), 
#                         list(siteID = df_mergedsiteID,dateID=dfmergeddateID), 
#                         mean, na.action = na.pass)
# 
# # SUBSET to calibration data
# calibration <- df_method1[which(df_method1$date < "2015-01"),]
# validation <- df_method1[which(df_method1$date > "2015-01"),]
# 
# ggplot(data = calibration, aes(x = factor(date), y = log(ratio), color = siteID)) +       
#   geom_line(aes(group = siteID)) + geom_point()
# 
# ggplot(data = calibration, aes(x = factor(date), y = log(ratio), color = siteID)) +       
#   geom_line(aes(group = siteID)) + geom_point()