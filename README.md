# SOaP
GE585 Spring 2019: "Soil Organism chAnging Predictably" fungal functional group forecast

Name: Yetianjian Wang
Email: wangytj@bu.edu

Name: Steve Gougherty
Email: gougher@bu.edu

Name: Ryan Quinn 
Email: rkq@bu.edu 

Name: Zoey Werbin
Email: zrwerbin@gmail.com

## Updated 2019/03/29

## How to run our calibration model

The entire model can currently be run using the "analysis/Calibration_model_milestone.Rmd" script.
It calls a script called "data_construction/aggregate_calibration_data.R", which, in turn, calls the scripts for downloading a variety of NEON, Daymet, and WorldClim products.

Our remote-sensing product, Canopy Height Model, should only be downloaded via the SCC, so we recommend reading in the hard-coded paths to Yetianjian Wang's SCC directory. If you want to download them yourself (please don't), here's what to do (from the data_construction/NEON_covariates/):

1. To generate the CHM data for each site, run the 5 r scripts starting with "Demo", then run the r script "ObtainingAndCalculatingTheLocatingUnitAreaOfEachPlot.R" and csv files "XXXX_Soilcore_CorrespondingMeanCHM_OfEachPlot_YYYY.csv" (XXXX=CPER/DSNY/HARV/OSBS/STER, YYYY=2013/2014) will generate.
2. To visualize the histograms of correponding mean Canopy Height Model (CHM) of each sitem run the r script "VisualizingMeanCHM_AllSites.R".

Note: 
1. CPER and STER only have 2013 data for CHM, and DSNY, HARV, OSBS only have 2014 data for CHM.
2. Running through 5 "Demo" scripts can take hours. 5 output csv files from r scripts starting with "Demo" are in advance generated and stored in the same folder with "XXXX_Soilcore_CorrespondingMeanCHM_OfEachPlot_YYYY.csv" that    is to be generated. One can directly run "ObtainingAndCalculatingTheLocatingUnitAreaOfEachPlot.R" to generate the 
"XXXX_Soilcore_CorrespondingMeanCHM_OfEachPlot_YYYY.csv" files.
