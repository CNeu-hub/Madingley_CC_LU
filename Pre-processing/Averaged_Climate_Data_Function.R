###--------------------------------------------------------------------------------------------------------------------------------------------------------------------###
### TITLE:       Averaged Climate Data Function                                                                                                                        ###
### DESCRIPTION: Function to calculate monthly averages of netcdf climate data & save result as monthly tiff raster output                                             ###
### DATE: 08.09.2023                                                                                                                                                   ###
###--------------------------------------------------------------------------------------------------------------------------------------------------------------------###

# x = spatial raster dataset 
# timeframe = timeframe to be considered (default = last 30 years) 
# output_directory = where do you want to save data??

averages <- function(x, timeframe = 673:1032, output_directory = "output") {
  
  #load libraries
  require(terra)
  require(tidyverse)
  
  #Create the output directory if it doesn't exist
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  #Subset to the specified timeframe (e.g., 673:1032 for 2071 - 2100)
  subset <- x[[timeframe]]
  
  #Initialize a list to store the monthly averages
  averaged <- list()
  
  #Loop through each month (1 to 12) and calculate the average
  for (i in 1:12) {
    #Extract all bands corresponding to the current month
    month_bands <- subset[[seq(i, 360, 12)]]
    
    #Calculate the average for the current month
    averaged[[i]] <- mean(month_bands, na.rm = TRUE)
    
    #Define the output file path for the current month using the input dataset name
    input_name <- deparse(substitute(x))
    
    month_str <- ifelse(input_name == "diurnal_temperature_range_",i,sprintf("%02d", i))
    #month_str <- sprintf("%02d", i)  # Ensure two digits for month (01, 02, ..., 12)
    
    output_file <- file.path(output_directory, paste0(ifelse(input_name == "near_surface_temperature_", "near-surface_temperature_",input_name), month_str, ".tif"))
    
    #Save the monthly average as a raster TIFF file
    writeRaster(averaged[[i]], filename = output_file, overwrite = TRUE)
    
    cat("Saved:", output_file, "\n")
  }
  
  #Create a SDS (spatial data stack) from the list of monthly averages
  averaged_sds <- sds(averaged)
  
  #Create a 3x4 layout for 12 plots
  layout(matrix(1:12, nrow = 3, ncol = 4))
  
  #Plot all 12 months
  for (i in 1:12) {
    plot(averaged_sds[[i]], main = month.abb[i])
  }
  
  # Return the SDS
  return(averaged_sds)
}
