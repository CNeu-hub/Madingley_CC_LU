###-----------------------------------------------------------------------------------------###
###Insert climate data function                                                             ###
### DESCRIPTION: Function to load climate data of:                                          ###
###                                                                                         ### 
###              "Model-based Impact Analysis of Climate Change and Land-use Intensification###
###              on Trophic Networks."                                                      ###
###-----------------------------------------------------------------------------------------###
###Author: Hoeks, S. (2021)                                                                 ###

madingley_inputs()
# I just erased the part which is checking the 0.5 degree data, and downloading 
# it if it's not available, and modified the path to the rasters
# This script can be used to load any raster data out of a folder, which is 
# structured like madingley default inputs (layer, layer names, monthly data of 
# climate variables) can be loaded as function like the "loadhalfdegreeinputdata" 
# function of Selwyn to the Environment to insert the climate data into Madingley sptl inputs

# From Christian: This is the source of the modifie CLimate Data
# I've got it from the official climate data store of WCRP for CMIP6: 
# https://esgf-data.dkrz.de/search/cmip6-dkrz/ .
# There I searched for "CNRM-CM6-1-HR" which is a french high resolution model, 
# the SSP, monthly frequency and the variable.This is the doi for 
# SSP585 https://cera-www.dkrz.de/WDCC/ui/cerasearch/cmip6?input=CMIP6.ScenarioMIP.CNRM-CERFACS.CNRM-CM6-1-HR.ssp585 , 
# and this for SSP126 
# https://cera-www.dkrz.de/WDCC/ui/cerasearch/cmip6?input=CMIP6.ScenarioMIP.CNRM-CERFACS.CNRM-CM6-1-HR.ssp126 . 

insert_climate_data = function(wd) 
{
  if(substr(wd,(nchar(wd)+1)-1,nchar(wd))=='/')  wd=substr(wd,1,nchar(wd)-1)
  if(substr(wd,(nchar(wd)+1)-1,nchar(wd))=='\\') wd=substr(wd,1,nchar(wd)-1)
  
  rasters_path = paste0(wd)
  
  input = list()
  
  spatial_path = paste0(rasters_path)
  file_names = list.files(spatial_path, full.names = T, recursive = T)
  list_names = gsub("\\..*", "", list.files(spatial_path, full.names = F, recursive = T))
  FILES_sp = c("realm_classification", "land_mask", "hanpp", 
               "available_water_capacity") #, "Ecto_max", "Endo_C_max", "Endo_H_max", "Endo_O_max")
  FILES_sp_temp = c("terrestrial_net_primary_productivity", 
                    "near-surface_temperature", "precipitation", "ground_frost_frequency", 
                    "diurnal_temperature_range")
  for (i in FILES_sp) {
    if (length(grep(i, file_names, value = T)) != 1) {
      stop("Could not find raster: ", i, ".tif \n")
    }
  }
  for (i in FILES_sp_temp) {
    if (length(grep(i, file_names, value = T)) != 12) {
      stop("Could not find raster all 12 monthly rasters containing data on: ", 
           i, "\n")
    }
  }
  cat("Reading default input rasters from: ", spatial_path)
  for (i in FILES_sp) {
    file_name = grep(i, file_names, value = T)
    cat(".")
    input[[i]] = raster(file_name)
  }
  for (i in FILES_sp_temp) {
    file_name = grep(i, file_names, value = T)
    if (length(grep("_1.tif", file_name, value = T)) == 
        0) {
      file_name_sort = file_name
    }
    else {
      if (substr(spatial_path, nchar(spatial_path), 
                 nchar(spatial_path)) == "/") {
        file_name_sort = paste0(spatial_path, i, "_", 
                                1:12)
      }
      else {
        file_name_sort = paste0(spatial_path, "/", 
                                i, "_", 1:12, ".tif")
      }
    }
    if (length(file_name_sort) == 12) {
      input[[i]] = brick(lapply(file_name_sort, raster))
    }
    cat(".")
  }
  cat("\n")
  
  # no max body mass raster available yet, setting single value
  input$Ecto_max = input$Endo_C_max = input$Endo_H_max = input$Endo_O_max = input$land_mask
  input$Ecto_max[] = 50000
  input$Endo_C_max[] = 200000
  input$Endo_H_max[] = 1000000
  input$Endo_O_max[] = 500000
  
  return(input)
}
