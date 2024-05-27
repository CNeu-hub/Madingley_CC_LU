###--------------------------------------------------------------------------------------------------------------------------------------------------------------------###
### TITLE:       Averaged Climate Data.                                                                                                                                ###
### DESCRIPTION: Calculates monthly averages over last 30 years time frame for climate variables as input for MadingleyR & export them                                 ###
###              to monthly raster bands compatible with direct MadingleyR input                                                                                       ###
### DATE: 08.09.2023                                                                                                                                                   ###
###--------------------------------------------------------------------------------------------------------------------------------------------------------------------###

source("Pre-processing/Averaged_Climate_Data_Function.R")

#create path for easier file selection 
climatepath <- "Path/to/climate/data"
outpath <- "Path/to/store/output"

#Load nc raster data for historical climate scenario 
near_surface_temperature_ <- rast(climatepath %+% "Historical/tas_FINAL.nc")
diurnal_temperature_range_ <- rast(climatepath %+% "Historical/diurnal.nc")
precipitation_ <- rast(climatepath %+% "Historical/pr_FINAL.nc")
terrestrial_net_primary_productivity_ <- rast(climatepath %+% "Historical/npp_FINAL.nc")

#calculate monthly averages & save as separate tiff files, historical has different time frame 
tas <- averages(near_surface_temperature_, timeframe = 1621:1980, output_directory = outpath%+%"Historical/")
diurnal <- averages(diurnal_temperature_range_, timeframe = 1621:1980, output_directory = outpath%+%"Historical/")
pr <- averages(precipitation_, timeframe = 1621:1980, output_directory = outpath%+%"Historical/")
npp <- averages(terrestrial_net_primary_productivity_, timeframe = 1621:1980, output_directory = outpath%+%"Historical/")

#Load nc raster data for ssp126 climate scenario 
near_surface_temperature_ <- rast(climatepath %+% "SSP2.6/tas_FINAL.nc")
diurnal_temperature_range_ <- rast(climatepath %+% "SSP2.6/diurnal.nc")
precipitation_ <- rast(climatepath %+% "SSP2.6/pr_FINAL.nc")
terrestrial_net_primary_productivity_ <- rast(climatepath %+% "SSP2.6/npp_FINAL.nc")

#calculate monthly averages & save as separate tiff files 
tas <- averages(near_surface_temperature_,output_directory = outpath%+%"SSP126/")
diurnal <- averages(diurnal_temperature_range_,output_directory = outpath%+%"SSP126/")
pr <- averages(precipitation_,output_directory = outpath%+%"SSP126/")
npp <- averages(terrestrial_net_primary_productivity_,output_directory = outpath%+%"SSP126/")

#Load nc raster data for historical ssp585 scenario 
near_surface_temperature_ <- rast(climatepath %+% "SSP8.5/tas_FINAL.nc")
diurnal_temperature_range_ <- rast(climatepath %+% "SSP8.5/diurnal.nc")
precipitation_ <- rast(climatepath %+% "SSP8.5/pr_FINAL.nc")
terrestrial_net_primary_productivity_ <- rast(climatepath %+% "SSP8.5/npp_FINAL.nc")

#calculate monthly averages & save as separate tiff files 
tas <- averages(near_surface_temperature_,output_directory = outpath%+%"SSP585/")
diurnal <- averages(diurnal_temperature_range_,output_directory = outpath%+%"SSP585/")
pr <- averages(precipitation_,output_directory = outpath%+%"SSP585/")
npp <- averages(terrestrial_net_primary_productivity_,output_directory = outpath%+%"SSP585/")


