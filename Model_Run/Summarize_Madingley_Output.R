#Script to analyse Madingley Outputs for further Analysis


library(dplyr)

source("Settings.R")

#define region
#define region
region <- "Brazil"
#NEW REGIONS (27.10.2023):
#spatial_window = c(25.5,28.5,61,69) #96 grid cells, works brilliant, Finland
#spatial_window = c(16,21,-22,-17) #90 grid cells, works brilliant, Namibia 
#spatial_window = c(-69,-61,-3,0) #96 grid cells, works brilliant, Brazil 
#spatial_window = c(-1,6,46,49) #84 grid cells, works brilliant, France

#create new path for region
regionpath <- paste0("/Volumes/scratch_60_days/neumanch/Outpath_Madingley/",region,sep="/")

#create vectors for filenames of scenarios 
filenames_historical <- list.files(path=regionpath, pattern="Historical", full.names=TRUE)
filenames_ssp126 <- list.files(path=regionpath, pattern="SSP126", full.names = TRUE)
filenames_ssp585 <- list.files(path=regionpath, pattern="SSP585", full.names = TRUE)

#load data for all scenarios with filenames
replicates_data_historical <- lapply(filenames_historical, function(x) mget(load(x)))
replicates_data_ssp126 <- lapply(filenames_ssp126, function(x) mget(load(x)))
replicates_data_ssp585 <- lapply(filenames_ssp585, function(x) mget(load(x)))

##subset to only needed output data.frames (output list objects)
replicates_data_historical <- lapply(replicates_data_historical[1:10], function(x) x$historical_2014_list)
replicates_data_ssp126 <- lapply(replicates_data_ssp126[1:10], function(x) x$SSP126_2100_list)
replicates_data_ssp585 <- lapply(replicates_data_ssp585[1:10], function(x) x$SSP585_2100_list)

#define simulation stage to loop over stages 
simulation_stage <- length(replicates_data_historical[[1]])

#create list object, for collecting aggregated results 
output <- vector("list", length = simulation_stage)

#create function that summarizes the data across replicates 
summarize_replicates <- function(replicates_data, i, variable) {
  
  #create object containing all replicates as rbind
  result <- do.call(rbind, lapply(replicates_data, function(replicate) replicate[[i]][[variable]]))
  
  #calculate mean across replicates, depending on condition across grid cells/functional groups, or across month
  if(variable == "cohorts") {
    result <- aggregate(result, list(result$GridcellIndex,result$FunctionalGroupIndex), mean)
    return(result)
  } else {
    result <- aggregate(result, list(result$Month), mean)
    return(result)
  }
}

#loop across all simulation stages & calculate replicates mean for each simulation stage 
#for historical scenarios 
for (i in 1:simulation_stage) {
  output[[i]]$time_line_cohorts <- summarize_replicates(replicates_data_historical, i, "time_line_cohorts")
  output[[i]]$time_line_stocks <- summarize_replicates(replicates_data_historical, i, "time_line_stocks")
  output[[i]]$cohorts <- summarize_replicates(replicates_data_historical, i, "cohorts")
  historical_2014_list <- output
}

#for SSP126 scenarios 
for (i in 1:simulation_stage) {
  output[[i]]$time_line_cohorts <- summarize_replicates(replicates_data_ssp126, i, "time_line_cohorts")
  output[[i]]$time_line_stocks <- summarize_replicates(replicates_data_ssp126, i, "time_line_stocks")
  output[[i]]$cohorts <- summarize_replicates(replicates_data_ssp126, i, "cohorts")
  SSP126_2100_list <- output
}

#for SSP585 scenarios 
for (i in 1:simulation_stage) {
  output[[i]]$time_line_cohorts <- summarize_replicates(replicates_data_ssp585, i, "time_line_cohorts")
  output[[i]]$time_line_stocks <- summarize_replicates(replicates_data_ssp585, i, "time_line_stocks")
  output[[i]]$cohorts <- summarize_replicates(replicates_data_ssp585, i, "cohorts")
  SSP585_2100_list <- output
}

save(historical_2014_list,SSP126_2100_list,SSP585_2100_list, file = outpath %+% region %+% "/Output_Summary.Rdata")


