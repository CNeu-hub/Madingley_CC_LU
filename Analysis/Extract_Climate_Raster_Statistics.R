###-----------------------------------------------------------------------------------------###
### TITLE:       EXTRACT RASTER STATISTICS                                                  ###
### DESCRIPTION: Extracts raster statistics for climate data of:                            ###
###              "Model-based Impact Analysis of Climate Change and Land-use Intensification###
###              on Trophic Networks."                                                      ###                                                                 ###
### PROCEDURE:   Crops climate data input to extent of study regions and extracts summary   ###
###              statistics for variables: HANPP, NPP, TAS, PR -> Base for climate bar plots###
### DATE:        21.09.2022                                                                 ###
###-----------------------------------------------------------------------------------------###    

library(ggplot2)
library(dplyr)
library(psych)
library(raster)
library(MadingleyR)

###---------------------------------------------------------------------------###
###                          CREATE SETTINGS / ENVIRONMENT                    ###
###---------------------------------------------------------------------------###
source("InsertClimateDataFunction.R") 

#create datapath
datapath <- paste(getwd(), "Output", sep = "/")
  
#load spatial input madingley 
sp_inputs = insert_climate_data("/Path/to/climate/data/")

#manual definition of climate scenario
scenario <- "Historical"
#scenario <- "SSP1-2.6"
#scenario <- "SSP5-8.5"

#create folder to store data
ifelse(!dir.exists(file.path(datapath)),dir.create(file.path(datapath)),print("This folder already exists"))

#generate yearly climate rasters for each variable 
tas <- sum(sp_inputs$`near-surface_temperature`)/12
pr <- sum(sp_inputs$precipitation)/12
NPP <- sum(sp_inputs$terrestrial_net_primary_productivity)/12
HANPP <- sp_inputs$hanpp

###---------------------------------------------------------------------------###
###                              CREATE STUDY AREAS                           ###
###---------------------------------------------------------------------------###
#create spatial window extent for each region  
spatial_window_finland = c(25.5,28.5,61,69) #96 grid cells, works brilliant, Finland
spatial_window_namibia = c(16,21,-22,-17) #90 grid cells, works brilliant, Namibia 
spatial_window_brazil = c(-69,-61,-3,0) #96 grid cells, works brilliant, Brazil 
spatial_window_france = c(-1,6,46,49) #84 grid cells, works brilliant, France


#create bounding box for cropping out of spatial_window for each region 
#Brazil
ext_br <- as(extent(spatial_window_brazil), 'SpatialPolygons')
crs(ext_br) <- "+proj=longlat +datum=WGS84 +no_defs"
#Namibia
ext_nam <- as(extent(spatial_window_namibia), 'SpatialPolygons')
crs(ext_nam) <- "+proj=longlat +datum=WGS84 +no_defs"
#France
ext_fr <- as(extent(spatial_window_france), 'SpatialPolygons')
crs(ext_fr) <- "+proj=longlat +datum=WGS84 +no_defs"
#Finland 
ext_fi <- as(extent(spatial_window_finland), 'SpatialPolygons')
crs(ext_fi) <- "+proj=longlat +datum=WGS84 +no_defs"

#Chris: crop rasters to spatial_window // region of concern
#Brazil
brazil_tas <- crop(tas,ext_br)
brazil_pr <- crop(pr,ext_br)
brazil_NPP <- crop(NPP,ext_br)
brazil_HANPP <- crop(HANPP,ext_br)
#Namibia
namibia_tas <- crop(tas,ext_nam)
namibia_pr <- crop(pr,ext_nam)
namibia_NPP <- crop(NPP,ext_nam)
namibia_HANPP <- crop(HANPP,ext_nam)
#France
france_tas <- crop(tas,ext_fr)
france_pr <- crop(pr,ext_fr)
france_NPP <- crop(NPP,ext_fr)
france_HANPP <- crop(HANPP,ext_fr)
#Finland
finland_tas <- crop(tas,ext_fi)
finland_pr <- crop(pr,ext_fi)
finland_NPP <- crop(NPP,ext_fi)
finland_HANPP <- crop(HANPP,ext_fi)

###---------------------------------------------------------------------------###
###                           EXTRAXT CLIMATE DATA                            ###
###---------------------------------------------------------------------------###
###    
#1. HANPP 
#analysis
hanpp_summary <- data.frame("Region" = c("Brazil","Namibia","France","Finland"),
                            "Mean" = c(cellStats(brazil_HANPP$hanpp_2005,stat=mean,na.rm=T),cellStats(namibia_HANPP$hanpp_2005,stat=mean,na.rm=T),cellStats(france_HANPP$hanpp_2005,stat=mean,na.rm=T),cellStats(finland_HANPP$hanpp_2005,stat=mean,na.rm=T)),
                            "SD" = c(cellStats(brazil_HANPP$hanpp_2005,stat=sd,na.rm=T),cellStats(namibia_HANPP$hanpp_2005,stat=sd,na.rm=T),cellStats(france_HANPP$hanpp_2005,stat=sd,na.rm=T),cellStats(finland_HANPP$hanpp_2005,stat=sd,na.rm=T)),
                            "Max" =c(cellStats(brazil_HANPP$hanpp_2005,stat=max,na.rm=T),cellStats(namibia_HANPP$hanpp_2005,stat=max,na.rm=T),cellStats(france_HANPP$hanpp_2005,stat=max,na.rm=T),cellStats(finland_HANPP$hanpp_2005,max,na.rm=T)),
                            "Min" =c(cellStats(brazil_HANPP$hanpp_2005,stat=min,na.rm=T),cellStats(namibia_HANPP$hanpp_2005,stat=min,na.rm=T),cellStats(france_HANPP$hanpp_2005,stat=min,na.rm=T),cellStats(finland_HANPP$hanpp_2005,min,na.rm=T)),
                            "Scenario" = paste0(scenario),
                            "Variable" = "Human Appropriation of net Primary Productivity",
                            "Unit" = "gC/m2/yr",
                            "Extent" = c("-69, -61, -3, 0", "16, 21, -22, -17", "-1, 6, 46, 49", "25.5, 28.5, 61, 69"))
write.csv(hanpp_summary,datapath %+% paste0("HANPP_summary","_",scenario,".csv"))

#2. NPP 
#analysis
npp_summary <- data.frame("Region" = c("Brazil","Namibia","France","Finland"),
                          "Mean" = c(cellStats(brazil_NPP$layer,stat=mean,na.rm=T),cellStats(namibia_NPP$layer,stat=mean,na.rm=T),cellStats(france_NPP$layer,stat=mean,na.rm=T),cellStats(finland_NPP$layer,stat=mean,na.rm=T)),
                          "SD" = c(cellStats(brazil_NPP$layer,stat=sd,na.rm=T),cellStats(namibia_NPP$layer,stat=sd,na.rm=T),cellStats(france_NPP$layer,stat=sd,na.rm=T),cellStats(finland_NPP$layer,stat=sd,na.rm=T)),
                          "Max" =c(cellStats(brazil_NPP$layer,stat=max,na.rm=T),cellStats(namibia_NPP$layer,stat=max,na.rm=T),cellStats(france_NPP$layer,stat=max,na.rm=T),cellStats(finland_NPP$layer,max,na.rm=T)),
                          "Min" =c(cellStats(brazil_NPP$layer,stat=min,na.rm=T),cellStats(namibia_NPP$layer,stat=min,na.rm=T),cellStats(france_NPP$layer,stat=min,na.rm=T),cellStats(finland_NPP$layer,min,na.rm=T)),
                          "Scenario" = paste0(scenario),
                          "Variable" = "Net Primary Productivity",
                          "Unit" = "gC/m2/d",
                          "Extent" = c("-69, -61, -3, 0", "16, 21, -22, -17", "-1, 6, 46, 49", "25.5, 28.5, 61, 69"))
write.csv(npp_summary,datapath %+% paste0("NPP_summary","_",scenario,".csv"))

#3. TAS 
#analysis
tas_summary <- data.frame("Region" = c("Brazil","Namibia","France","Finland"),
                          "Mean" = c(cellStats(brazil_tas$layer,stat=mean,na.rm=T),cellStats(namibia_tas$layer,stat=mean,na.rm=T),cellStats(france_tas$layer,stat=mean,na.rm=T),cellStats(finland_tas$layer,stat=mean,na.rm=T)),
                          "SD" = c(cellStats(brazil_tas$layer,stat=sd,na.rm=T),cellStats(namibia_tas$layer,stat=sd,na.rm=T),cellStats(france_tas$layer,stat=sd,na.rm=T),cellStats(finland_tas$layer,stat=sd,na.rm=T)),
                          "Max" =c(cellStats(brazil_tas$layer,stat=max,na.rm=T),cellStats(namibia_tas$layer,stat=max,na.rm=T),cellStats(france_tas$layer,stat=max,na.rm=T),cellStats(finland_tas$layer,max,na.rm=T)),
                          "Min" =c(cellStats(brazil_tas$layer,stat=min,na.rm=T),cellStats(namibia_tas$layer,stat=min,na.rm=T),cellStats(france_tas$layer,stat=min,na.rm=T),cellStats(finland_tas$layer,min,na.rm=T)),
                          "Scenario" = paste0(scenario),
                          "Variable" = "Near-Surface Temperature",
                          "Unit" = "°C",
                          "Extent" = c("-69, -61, -3, 0", "16, 21, -22, -17", "-1, 6, 46, 49", "25.5, 28.5, 61, 69"))
write.csv(tas_summary,datapath %+% paste0("tas_summary","_",scenario,".csv"))


#4. PR 
#analysis
pr_summary <- data.frame("Region" = c("Brazil","Namibia","France","Finland"),
                         "Mean" = c(cellStats(brazil_pr$layer,stat=mean,na.rm=T),cellStats(namibia_pr$layer,stat=mean,na.rm=T),cellStats(france_pr$layer,stat=mean,na.rm=T),cellStats(finland_pr$layer,stat=mean,na.rm=T)),
                         "SD" = c(cellStats(brazil_pr$layer,stat=sd,na.rm=T),cellStats(namibia_pr$layer,stat=sd,na.rm=T),cellStats(france_pr$layer,stat=sd,na.rm=T),cellStats(finland_pr$layer,stat=sd,na.rm=T)),
                         "Max" =c(cellStats(brazil_pr$layer,stat=max,na.rm=T),cellStats(namibia_pr$layer,stat=max,na.rm=T),cellStats(france_pr$layer,stat=max,na.rm=T),cellStats(finland_pr$layer,max,na.rm=T)),
                         "Min" =c(cellStats(brazil_pr$layer,stat=min,na.rm=T),cellStats(namibia_pr$layer,stat=min,na.rm=T),cellStats(france_pr$layer,stat=min,na.rm=T),cellStats(finland_pr$layer,min,na.rm=T)),
                         "Scenario" = paste0(scenario),
                         "Variable" = "Precipitation",
                         "Unit" = "mm/month",
                         "Extent" = c("-69, -61, -3, 0", "16, 21, -22, -17", "-1, 6, 46, 49", "25.5, 28.5, 61, 69"))
write.csv(pr_summary,datapath %+% paste0("pr_summary","_",scenario,".csv"))

all <- rbind(hanpp_summary, npp_summary, tas_summary, pr_summary)
all[2:5] <- round(all[2:5], digits = 2)

write.csv(all,datapath %+% paste0("Climate_Summary","_",scenario,".csv"))

