###--------------------------------------------------------------------------------------------------------------------------------------------------------------------###
### TITLE:       MadingleyR Simulations                                                                                                                                ###
### DESCRIPTION: Runs MadingleyR for different regions / climate datasets of:                                                                                          ###                              
###              "Model-based Impact Analysis of Climate Change and Land-use Intensification on Trophic Networks."                                                     ###
### NOTES:       - Has to be applied for each region separately (spatial extents are below in the code).                                                               ###                                              
###              - Creates climate scenarios + HANPP (Current land-use intensity) + max land-use intensity scenarios for historical, SSP1-2.6, and SSP5-8.5 climate    ###
###              - Output is saved as rdata list (stage 1 of list = climate scenario, entry 2 = climate scenario under current land use, last entry of list = climate  ###  
###              scenario under max. land use); length of list depends on time until only 10 % vegetation = left in region                                             ###
###             - Output = rdata list, is saved in input folder within repository.                                                                                     ###
###             - Climate data can be obtained at the sources provided in the repository description, or SI, has to be pre-processed before usage in model.            ###
### DATE: 27.05.2024.                                                                                                                                                  ###
###--------------------------------------------------------------------------------------------------------------------------------------------------------------------###

###Create environment 
library(MadingleyR)
library(terra)
library(tidyverse)

#Source: Hoeks, S. (2022): 
#Function to crop spatial raster input of MadingleyR using spatial window
source("https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/crop_spatial_rasters_to_window.r")

#Insert climate data that matches MadingleyR requirements
source("InsertClimateDataFunction.R") 

#define spatial extent of region, meeting requirements of MadingleyR (min. longitude, max. longitude, min. latitude, max. latitude)
#spatial_window <- c(25.5,28.5,61,69) #96 grid cells, Finland
#spatial_window <- c(16,21,-22,-17) #90 grid cells, Namibia 
#spatial_window <- c(-69,-61,-3,0) #96 grid cells, Brazil 
spatial_window <- c(-1,6,46,49) #84 grid cells, France

#define region
region <- "France"

#create output path for regional output 
regionpath <- paste0(getwd(), sep = "/", "Input/", region, sep="/")

#check if region extent is chosen correctly 
plot_spatialwindow(spatial_window)

#define no. of replicates (there are some stochastic processes inside the model, so we use replicates & average across them to generate more robust analysis of output)
reps <- 10

###----------------------------------------------###
###         1.) HISTORICAL CLIMATE SCENARIO      ###
###----------------------------------------------###

for(i in 1:reps) {

  cat("Replicate", i, "historical climate is running\n")
  
#Create model inputs
sp_inputs_histo = insert_climate_data("/Path/to/historical/climate/data") #here it is possible to insert any climate data matching the default madingleyR requirements in 0.5degree resolution
chrt_def = madingley_inputs("cohort definition")
stck_def = madingley_inputs("stock definition")
mdl_prms = madingley_inputs("model parameters")

# crop the raster to spatial window extent
# helps to speed up loading times and makes working with the inputs easier
sp_inputs_histo = crop_spatial_rasters_to_window(inputs = sp_inputs_histo, spatial_window = spatial_window)

###------------------###
### INITIALIZE MODEL ###
###------------------###

historical_2014 = madingley_init(spatial_window = spatial_window, 
                                 cohort_def = chrt_def,
                                 stock_def = stck_def,
                                 spatial_inputs = sp_inputs_histo, 
                                 max_cohort = 1000)

###---------------###
### RUN THE MODEL ###
###---------------###
###-------------------------------------------------------###
### Simulation Experiment 1: Climate (Spin-up)            ###
###-------------------------------------------------------###
###create a list file, where all outputs will be stored
historical_2014_list = list()

historical_2014_list[[1]] = madingley_run(out_dir = regionpath,
                                          madingley_data = historical_2014,
                                          years = 200, 
                                          cohort_def = chrt_def,
                                          stock_def = stck_def,
                                          spatial_inputs = sp_inputs_histo,
                                          model_parameters = mdl_prms,
                                          max_cohort = 1000) 

###-------------------------------------------------------###
### Simulation Experiment 2: Current Land Use (HANPP)     ###
###-------------------------------------------------------###

#run model for additional 200 years with applying hanpp as input raster with gC/m-2/year
historical_2014_list[[2]] = madingley_run(out_dir = regionpath, 
                                          madingley_data = historical_2014_list[[1]],
                                          years = 200, 
                                          cohort_def = chrt_def,
                                          stock_def = stck_def,
                                          spatial_inputs = sp_inputs_histo,
                                          model_parameters = mdl_prms,
                                          max_cohort = 1000, 
                                          apply_hanpp = 2)  #apply_hanpp = 2 is to apply hanpp with gC/m-2/year as input

###-------------------------------------------------------###
### Simulation Experiment 3: Maximum Land Use (maxHANPP)  ###
###-------------------------------------------------------###

#calculation of fractional hanpp raster
sp_inputs_histo$hanpp[] = sp_inputs_histo$hanpp[] + abs(min(sp_inputs_histo$hanpp[], na.rm = TRUE)) 
sp_inputs_histo$hanpp[] = sp_inputs_histo$hanpp[] / max(sp_inputs_histo$hanpp[], na.rm = TRUE)  
sp_inputs_histo$hanpp[] = 1-sp_inputs_histo$hanpp[]  # The subtraction 1-HANPP leads to swapped values, this means the lower values are the values with the highest HANPP and the higher values are the values with lowest HANPP --> So, we can assume that a reduction in these values is a reduction in biomass, and not a reduction in HANPP.

#remove NA's
sp_inputs_histo$hanpp[is.na(sp_inputs_histo$hanpp[])] = 0.001

###VEGETATION REDUCTION###
#Create HANPP backup to compare with reduced HANPP later
hanpp_backup = sp_inputs_histo$hanpp

#Reduce vegetation using a while loop, stops when needed (no more values above 0.1 in hanpp raster)
while(max(sp_inputs_histo$hanpp[])>0.1) {
  
  #use ifelse to reduce vegetation by 0.1 per 5 years (if its > 0.1), but not below the threshold 0.1
  sp_inputs_histo$hanpp[] = ifelse(sp_inputs_histo$hanpp[] > 0.1,
                                        ifelse(sp_inputs_histo$hanpp[] - 0.1 < 0.1,0.1,sp_inputs_histo$hanpp[] - 0.1),
                                        sp_inputs_histo$hanpp[])
  plot(sp_inputs_histo$hanpp)
  
  #run the model for 10 years with HANPP reduced by 0.1 (10 %)
  historical_2014_list[[length(historical_2014_list)+1]] = madingley_run(out_dir = regionpath,
                                                                         years = 10, 
                                                                         madingley_data = historical_2014_list[[length(historical_2014_list)]],
                                                                         spatial_inputs = sp_inputs_histo,
                                                                         silenced = TRUE,
                                                                         max_cohort = 1000, 
                                                                         apply_hanpp = 1) 
}

#Compare start HANPP with end HANPP 
par(mfrow=c(1,2))
plot(hanpp_backup,main="Starting HANPP", xlab = "Longitude", "Latitude")
plot(sp_inputs_histo$hanpp,main="HANPP Reduced to Maximum of 0.1", xlab = "Longitude", ylab = "Latitude")
par(mfrow=c(1,1))

###POST VEGETATION REDUCTION###
historical_2014_list[[length(historical_2014_list)+1]]  = madingley_run(out_dir = regionpath,
                                                                        madingley_data = historical_2014_list[[length(historical_2014_list)]], 
                                                                        years = 200, #at least 200
                                                                        model_parameters = mdl_prms, 
                                                                        spatial_inputs = sp_inputs_histo,
                                                                        cohort_def = chrt_def,
                                                                        max_cohort = 1000, #at least 500
                                                                        apply_hanpp = 1)

#save image for further use in next scripts
save.image(paste0(regionpath %+% i %+% "Historical_Run.Rdata"))

#clear all data frames from environment for next scenario
rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "data.frame"])
#clear all lists from environment for next scenario
rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "list"])

}

###----------------------------------------------###
###         2.) SSP1-2.6 CLIMATE SCENARIO        ###
###----------------------------------------------###

for(i in 1:reps) {
  
  cat("Replicate", i, "SSP1-2.6 climate is running\n")

  #Create model inputs  
  sp_inputs_SSP126 = insert_climate_data("/Path/to/SSP1-2.6/climate/data") #here it is possible to insert any climate data matching the default madingleyR requirements in 0.5degree resolution
  chrt_def = madingley_inputs("cohort definition")
  stck_def = madingley_inputs("stock definition")
  mdl_prms = madingley_inputs("model parameters")
  
  # crop the raster to spatial window extent
  # helps to speed up loading times and makes working with the inputs easier
  sp_inputs_SSP126 = crop_spatial_rasters_to_window(inputs = sp_inputs_SSP126, spatial_window = spatial_window)
  
  ###------------------###
  ### INITIALIZE MODEL ###
  ###------------------###
  
  historical_2014 = madingley_init(spatial_window = spatial_window, 
                                   cohort_def = chrt_def,
                                   stock_def = stck_def,
                                   spatial_inputs = sp_inputs_SSP126, 
                                   max_cohort = 1000)
  
  ###---------------###
  ### RUN THE MODEL ###
  ###---------------###
  ###-------------------------------------------------------###
  ### Simulation Experiment 1: Climate (Spin-up)            ###
  ###-------------------------------------------------------###
  ###create a list file, where all outputs will be stored
  historical_2014_list = list()
  
  historical_2014_list[[1]] = madingley_run(out_dir = regionpath,
                                            madingley_data = historical_2014,
                                            years = 200, 
                                            cohort_def = chrt_def,
                                            stock_def = stck_def,
                                            spatial_inputs = sp_inputs_SSP126,
                                            model_parameters = mdl_prms,
                                            max_cohort = 1000) 
  
  ###-------------------------------------------------------###
  ### Simulation Experiment 2: Current Land Use (HANPP)     ###
  ###-------------------------------------------------------###
  
  #run model for additional 200 years with applying hanpp as input raster with gC/m-2/year
  historical_2014_list[[2]] = madingley_run(out_dir = regionpath, 
                                            madingley_data = historical_2014_list[[1]],
                                            years = 200, 
                                            cohort_def = chrt_def,
                                            stock_def = stck_def,
                                            spatial_inputs = sp_inputs_SSP126,
                                            model_parameters = mdl_prms,
                                            max_cohort = 1000, 
                                            apply_hanpp = 2)  #apply_hanpp = 2 is to apply hanpp with gC/m-2/year as input
  
  ###-------------------------------------------------------###
  ### Simulation Experiment 3: Maximum Land Use (maxHANPP)  ###
  ###-------------------------------------------------------###
  
  #calculation of fractional hanpp raster
  sp_inputs_SSP126$hanpp[] = sp_inputs_SSP126$hanpp[] + abs(min(sp_inputs_SSP126$hanpp[], na.rm = TRUE)) 
  sp_inputs_SSP126$hanpp[] = sp_inputs_SSP126$hanpp[] / max(sp_inputs_SSP126$hanpp[], na.rm = TRUE)  
  sp_inputs_SSP126$hanpp[] = 1-sp_inputs_SSP126$hanpp[]  # The subtraction 1-HANPP leads to swapped values, this means the lower values are the values with the highest HANPP and the higher values are the values with lowest HANPP --> So, we can assume that a reduction in these values is a reduction in biomass, and not a reduction in HANPP.
  
  #remove NA's
  sp_inputs_SSP126$hanpp[is.na(sp_inputs_SSP126$hanpp[])] = 0.001
  
  ###VEGETATION REDUCTION###
  #Create HANPP backup to compare with reduced HANPP later
  hanpp_backup = sp_inputs_SSP126$hanpp
  
  #Reduce vegetation using a while loop, stops when needed (no more values above 0.1 in hanpp raster)
  while(max(sp_inputs_SSP126$hanpp[])>0.1) {
    
    #use ifelse to reduce vegetation by 0.1 per 5 years (if its > 0.1), but not below the threshold 0.1
    sp_inputs_SSP126$hanpp[] = ifelse(sp_inputs_SSP126$hanpp[] > 0.1,
                                          ifelse(sp_inputs_SSP126$hanpp[] - 0.1 < 0.1,0.1,sp_inputs_SSP126$hanpp[] - 0.1),
                                          sp_inputs_SSP126$hanpp[])
    plot(sp_inputs_SSP126$hanpp)
    
    #run the model for 10 years with HANPP reduced by 0.1 (10 %)
    historical_2014_list[[length(historical_2014_list)+1]] = madingley_run(out_dir = regionpath,
                                                                           years = 10, 
                                                                           madingley_data = historical_2014_list[[length(historical_2014_list)]],
                                                                           spatial_inputs = sp_inputs_SSP126,
                                                                           silenced = TRUE,
                                                                           max_cohort = 1000, 
                                                                           apply_hanpp = 1) 
  }
  
  #Compare start HANPP with end HANPP 
  par(mfrow=c(1,2))
  plot(hanpp_backup,main="Starting HANPP", xlab = "Longitude", "Latitude")
  plot(sp_inputs_SSP126$hanpp,main="HANPP Reduced to Maximum of 0.1", xlab = "Longitude", ylab = "Latitude")
  par(mfrow=c(1,1))
  
  ###POST VEGETATION REDUCTION###
  historical_2014_list[[length(historical_2014_list)+1]]  = madingley_run(out_dir = regionpath,
                                                                          madingley_data = historical_2014_list[[length(historical_2014_list)]], 
                                                                          years = 200, #at least 200
                                                                          model_parameters = mdl_prms, 
                                                                          spatial_inputs = sp_inputs_SSP126,
                                                                          cohort_def = chrt_def,
                                                                          max_cohort = 1000, #at least 500
                                                                          apply_hanpp = 1)
  
  #save image for further use in next scripts
  save.image(paste0(regionpath %+% i %+% "SSP126_Run.Rdata"))
  
  #clear all data frames from environment for next scenario
  rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "data.frame"])
  #clear all lists from environment for next scenario
  rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "list"])
  
}

###----------------------------------------------###
###         2.) SSP5-8.5 CLIMATE SCENARIO        ###
###----------------------------------------------###

for(i in 1:reps) {
  
  cat("Replicate", i, "SSP5-8.5 climate is running\n")

  #Create model inputs
  sp_inputs_SSP585 = insert_climate_data("/Path/to/SSP1-2.6/climate/data") #here it is possible to insert any climate data matching the default madingleyR requirements in 0.5degree resolution
  chrt_def = madingley_inputs("cohort definition")
  stck_def = madingley_inputs("stock definition")
  mdl_prms = madingley_inputs("model parameters")
  
  # crop the raster to spatial window extent
  # helps to speed up loading times and makes working with the inputs easier
  sp_inputs_SSP585 = crop_spatial_rasters_to_window(inputs = sp_inputs_SSP585, spatial_window = spatial_window)
  
  ###------------------###
  ### INITIALIZE MODEL ###
  ###------------------###
  
  historical_2014 = madingley_init(spatial_window = spatial_window, 
                                   cohort_def = chrt_def,
                                   stock_def = stck_def,
                                   spatial_inputs = sp_inputs_SSP585, 
                                   max_cohort = 1000)
  
  ###---------------###
  ### RUN THE MODEL ###
  ###---------------###
  ###-------------------------------------------------------###
  ### Simulation Experiment 1: Climate (Spin-up)            ###
  ###-------------------------------------------------------###
  ###create a list file, where all outputs will be stored
  historical_2014_list = list()
  
  historical_2014_list[[1]] = madingley_run(out_dir = regionpath,
                                            madingley_data = historical_2014,
                                            years = 200, 
                                            cohort_def = chrt_def,
                                            stock_def = stck_def,
                                            spatial_inputs = sp_inputs_SSP585,
                                            model_parameters = mdl_prms,
                                            max_cohort = 1000) 
  
  ###-------------------------------------------------------###
  ### Simulation Experiment 2: Current Land Use (HANPP)     ###
  ###-------------------------------------------------------###
  
  #run model for additional 200 years with applying hanpp as input raster with gC/m-2/year
  historical_2014_list[[2]] = madingley_run(out_dir = regionpath, 
                                            madingley_data = historical_2014_list[[1]],
                                            years = 200, 
                                            cohort_def = chrt_def,
                                            stock_def = stck_def,
                                            spatial_inputs = sp_inputs_SSP585,
                                            model_parameters = mdl_prms,
                                            max_cohort = 1000, 
                                            apply_hanpp = 2)  #apply_hanpp = 2 is to apply hanpp with gC/m-2/year as input
  
  ###-------------------------------------------------------###
  ### Simulation Experiment 3: Maximum Land Use (maxHANPP)  ###
  ###-------------------------------------------------------###
  
  #calculation of fractional hanpp raster
  sp_inputs_SSP585$hanpp[] = sp_inputs_SSP585$hanpp[] + abs(min(sp_inputs_SSP585$hanpp[], na.rm = TRUE)) 
  sp_inputs_SSP585$hanpp[] = sp_inputs_SSP585$hanpp[] / max(sp_inputs_SSP585$hanpp[], na.rm = TRUE)  
  sp_inputs_SSP585$hanpp[] = 1-sp_inputs_SSP585$hanpp[]  # The subtraction 1-HANPP leads to swapped values, this means the lower values are the values with the highest HANPP and the higher values are the values with lowest HANPP --> So, we can assume that a reduction in these values is a reduction in biomass, and not a reduction in HANPP.
  
  #remove NA's
  sp_inputs_SSP585$hanpp[is.na(sp_inputs_SSP585$hanpp[])] = 0.001
  
  ###VEGETATION REDUCTION###
  #Create HANPP backup to compare with reduced HANPP later
  hanpp_backup = sp_inputs_SSP585$hanpp
  
  #Reduce vegetation using a while loop, stops when needed (no more values above 0.1 in hanpp raster)
  while(max(sp_inputs_SSP585$hanpp[])>0.1) {
    
    #use ifelse to reduce vegetation by 0.1 per 5 years (if its > 0.1), but not below the threshold 0.1
    sp_inputs_SSP585$hanpp[] = ifelse(sp_inputs_SSP585$hanpp[] > 0.1,
                                      ifelse(sp_inputs_SSP585$hanpp[] - 0.1 < 0.1,0.1,sp_inputs_SSP585$hanpp[] - 0.1),
                                      sp_inputs_SSP585$hanpp[])
    plot(sp_inputs_SSP585$hanpp)
    
    #run the model for 10 years with HANPP reduced by 0.1 (10 %)
    historical_2014_list[[length(historical_2014_list)+1]] = madingley_run(out_dir = regionpath,
                                                                           years = 10, 
                                                                           madingley_data = historical_2014_list[[length(historical_2014_list)]],
                                                                           spatial_inputs = sp_inputs_SSP585,
                                                                           silenced = TRUE,
                                                                           max_cohort = 1000, 
                                                                           apply_hanpp = 1) 
  }
  
  #Compare start HANPP with end HANPP 
  par(mfrow=c(1,2))
  plot(hanpp_backup,main="Starting HANPP", xlab = "Longitude", "Latitude")
  plot(sp_inputs_SSP585$hanpp,main="HANPP Reduced to Maximum of 0.1", xlab = "Longitude", ylab = "Latitude")
  par(mfrow=c(1,1))
  
  ###POST VEGETATION REDUCTION###
  historical_2014_list[[length(historical_2014_list)+1]]  = madingley_run(out_dir = regionpath,
                                                                          madingley_data = historical_2014_list[[length(historical_2014_list)]], 
                                                                          years = 200, #at least 200
                                                                          model_parameters = mdl_prms, 
                                                                          spatial_inputs = sp_inputs_SSP585,
                                                                          cohort_def = chrt_def,
                                                                          max_cohort = 1000, #at least 500
                                                                          apply_hanpp = 1)
  
  #save image for further use in next scripts
  save.image(paste0(regionpath %+% i %+% "SSP585_Run.Rdata"))
  
  #clear all data frames from environment for next scenario
  rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "data.frame"])
  #clear all lists from environment for next scenario
  rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "list"])
  
}


