###LOAD LIBRARIES
library(ggplot2)
library(dplyr)
library(psych)
library(RColorBrewer)
library(ggpubr)
library(ggpmisc)
library(MASS)
library(scales)
library(viridis)
library(pheatmap)
library(stringr)
library(ggh4x) #allows for in depth facet grid customization 

#define input path (Output_Data in repository)
myoutpath <- "/Users/neumanch/Desktop/Thesis/bpm-madingley/Christians Projects/MadingleyR_Sep_2023/Outpath"
outpath <- myoutpath
#define output path (Output_CSV in repository)
myfigpath <- "/Users/neumanch/Desktop/Thesis/bpm-madingley/Christians Projects/MadingleyR_Sep_2023/Figpath"
figpath <-myfigpath

#Activate aggregated FGs, or deactivate (Results for all FGs will be plotted)
#chose either "aggregated", or "all" 
FGs <- "all"


###--------###
### Finland ###
###--------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","Finland",sep="/","Output_Summary.Rdata"))
outpath <- myoutpath
figpath <- myfigpath

regionpath <- paste0(figpath,sep="/",dir.create(file.path(figpath, "Finland"), showWarnings = FALSE),sep="/")

#get output directory path
out_dir <- outpath 

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#just load data again 
#list[1] = control; list[2] = HANPP; last list entry [individual on region] = After vegreduction
#climate
cohorts_histo <- data.frame(historical_2014_list[[1]]$cohorts) 
cohorts_SSP126 <- data.frame(SSP126_2100_list[[1]]$cohorts)
cohorts_SSP585 <- data.frame(SSP585_2100_list[[1]]$cohorts) 

#HANPP
cohorts_histo_HANPP <- data.frame(historical_2014_list[[2]]$cohorts) 
cohorts_SSP126_HANPP <- data.frame(SSP126_2100_list[[2]]$cohorts)
cohorts_SSP585_HANPP <- data.frame(SSP585_2100_list[[2]]$cohorts) 

#Max. HANPP
cohorts_histo_maxHANPP <- data.frame(historical_2014_list[[10]]$cohorts) 
cohorts_SSP126_maxHANPP <- data.frame(SSP126_2100_list[[10]]$cohorts)
cohorts_SSP585_maxHANPP <- data.frame(SSP585_2100_list[[10]]$cohorts) 

#subset data 
#climate
cohorts_histo <- cohorts_histo[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126 <- cohorts_SSP126[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585 <- cohorts_SSP585[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#HANPP
cohorts_histo_HANPP <- cohorts_histo_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_HANPP <- cohorts_SSP126_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_HANPP <- cohorts_SSP585_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#maxHANPP
cohorts_histo_maxHANPP <- cohorts_histo_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_maxHANPP <- cohorts_SSP126_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_maxHANPP <- cohorts_SSP585_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#create new grouping variable for 
#climate
cohorts_histo$group = "Historical"
cohorts_SSP126$group = "SSP1-2.6"
cohorts_SSP585$group = "SSP5-8.5"
cohorts_histo$group2 = "Climate"
cohorts_SSP126$group2 = "Climate"
cohorts_SSP585$group2 = "Climate"

#HANPP
cohorts_histo_HANPP$group = "Historical"
cohorts_SSP126_HANPP$group = "SSP1-2.6"
cohorts_SSP585_HANPP$group = "SSP5-8.5"
cohorts_histo_HANPP$group2 = "Current Land Use"
cohorts_SSP126_HANPP$group2 = "Current Land Use"
cohorts_SSP585_HANPP$group2 = "Current Land Use"

#maxHANPP
cohorts_histo_maxHANPP$group = "Historical"
cohorts_SSP126_maxHANPP$group = "SSP1-2.6"
cohorts_SSP585_maxHANPP$group = "SSP5-8.5"
cohorts_histo_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP126_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP585_maxHANPP$group2 = "Maximum Land Use"

#combine dataframes for plotting dataset
test <- rbind(cohorts_histo,cohorts_SSP126,cohorts_SSP585,
              cohorts_histo_HANPP,cohorts_SSP126_HANPP,cohorts_SSP585_HANPP,
              cohorts_histo_maxHANPP,cohorts_SSP126_maxHANPP,cohorts_SSP585_maxHANPP)

#convert to kg 
test[3] <- test[3]/1000 

#change naming
names(test) <- c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance","group","group2")

#FG naming
if(FGs == "all") {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores (i.)", '1' = "Endotherm Carnivores (i.)", '2' =  "Endotherm Omnivores (i.)",
                                      '3'="Ectotherm Herbivores (s.)",'4'="Ectotherm Carnivores (s.)",'5'="Ectotherm Omnivores (s.)","6" = "Ectotherm Herbivores (i.)",
                                      '7' = "Ectotherm Carnivores (i.)",'8'="Ectotherm Omnivores (i.)" )
} else {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores", '1' = "Endotherm Carnivores", '2' =  "Endotherm Omnivores",
                                      '3'="Ectotherm Herbivores",'4'="Ectotherm Carnivores",'5'="Ectotherm Omnivores","6" = "Ectotherm Herbivores",
                                      '7' = "Ectotherm Carnivores",'8'="Ectotherm Omnivores" )
  
}

#manual model fitting
fitted_models_finland <- test %>%
  group_by(FunctionalGroupIndex,group, group2) %>%
  do(broom::tidy(lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .),conf.int=T)) %>%
  ungroup() %>%
  mutate(Region = "Finland")

models_finland <- fitted_models_finland %>%
  filter(term == "log10(IndividualBodyMass)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate,std.error,statistic,conf.high,conf.low,p.value) %>%
  mutate(estimate = round(estimate,digits = 2),
         std.error = round(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         conf.high = round(conf.high, digits = 2),
         conf.low = round(conf.low, digits = 2),
         p.value = format.pval(p.value, digits = 2)) %>%
  rename(Slope = estimate) 

#add intercept 
intercepts <- fitted_models_finland %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate) %>%
  rename(Intercept = estimate)

models_finland <- left_join(models_finland, intercepts)

###--------###
### Brazil ###
###--------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","Brazil",sep="/","Output_Summary.Rdata"))
outpath <- myoutpath
figpath <- myfigpath

regionpath <- paste0(figpath,sep="/",dir.create(file.path(figpath, "Brazil"), showWarnings = FALSE),sep="/")

#get output directory path
out_dir <- outpath 

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#just load data again 
#list[1] = control; list[2] = HANPP; last list entry [individual on region] = After vegreduction
#climate
cohorts_histo <- data.frame(historical_2014_list[[1]]$cohorts) 
cohorts_SSP126 <- data.frame(SSP126_2100_list[[1]]$cohorts)
cohorts_SSP585 <- data.frame(SSP585_2100_list[[1]]$cohorts) 

#HANPP
cohorts_histo_HANPP <- data.frame(historical_2014_list[[2]]$cohorts) 
cohorts_SSP126_HANPP <- data.frame(SSP126_2100_list[[2]]$cohorts)
cohorts_SSP585_HANPP <- data.frame(SSP585_2100_list[[2]]$cohorts) 

#Max. HANPP
cohorts_histo_maxHANPP <- data.frame(historical_2014_list[[11]]$cohorts) 
cohorts_SSP126_maxHANPP <- data.frame(SSP126_2100_list[[11]]$cohorts)
cohorts_SSP585_maxHANPP <- data.frame(SSP585_2100_list[[11]]$cohorts) 

#subset data 
#climate
cohorts_histo <- cohorts_histo[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126 <- cohorts_SSP126[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585 <- cohorts_SSP585[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#HANPP
cohorts_histo_HANPP <- cohorts_histo_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_HANPP <- cohorts_SSP126_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_HANPP <- cohorts_SSP585_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#maxHANPP
cohorts_histo_maxHANPP <- cohorts_histo_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_maxHANPP <- cohorts_SSP126_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_maxHANPP <- cohorts_SSP585_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#create new grouping variable for 
#climate
cohorts_histo$group = "Historical"
cohorts_SSP126$group = "SSP1-2.6"
cohorts_SSP585$group = "SSP5-8.5"
cohorts_histo$group2 = "Climate"
cohorts_SSP126$group2 = "Climate"
cohorts_SSP585$group2 = "Climate"

#HANPP
cohorts_histo_HANPP$group = "Historical"
cohorts_SSP126_HANPP$group = "SSP1-2.6"
cohorts_SSP585_HANPP$group = "SSP5-8.5"
cohorts_histo_HANPP$group2 = "Current Land Use"
cohorts_SSP126_HANPP$group2 = "Current Land Use"
cohorts_SSP585_HANPP$group2 = "Current Land Use"

#maxHANPP
cohorts_histo_maxHANPP$group = "Historical"
cohorts_SSP126_maxHANPP$group = "SSP1-2.6"
cohorts_SSP585_maxHANPP$group = "SSP5-8.5"
cohorts_histo_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP126_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP585_maxHANPP$group2 = "Maximum Land Use"

#combine dataframes for plotting dataset
test <- rbind(cohorts_histo,cohorts_SSP126,cohorts_SSP585,
              cohorts_histo_HANPP,cohorts_SSP126_HANPP,cohorts_SSP585_HANPP,
              cohorts_histo_maxHANPP,cohorts_SSP126_maxHANPP,cohorts_SSP585_maxHANPP)

#convert to kg 
test[3] <- test[3]/1000 

#change naming
names(test) <- c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance","group","group2")

#FG naming
if(FGs == "all") {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores (i.)", '1' = "Endotherm Carnivores (i.)", '2' =  "Endotherm Omnivores (i.)",
                                      '3'="Ectotherm Herbivores (s.)",'4'="Ectotherm Carnivores (s.)",'5'="Ectotherm Omnivores (s.)","6" = "Ectotherm Herbivores (i.)",
                                      '7' = "Ectotherm Carnivores (i.)",'8'="Ectotherm Omnivores (i.)" )
} else {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores", '1' = "Endotherm Carnivores", '2' =  "Endotherm Omnivores",
                                      '3'="Ectotherm Herbivores",'4'="Ectotherm Carnivores",'5'="Ectotherm Omnivores","6" = "Ectotherm Herbivores",
                                      '7' = "Ectotherm Carnivores",'8'="Ectotherm Omnivores" )
  
}

#manual model fitting
fitted_models_Brazil<- test %>%
  group_by(FunctionalGroupIndex,group, group2) %>%
  do(broom::tidy(lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .),conf.int=T)) %>%
  ungroup() %>%
  mutate(Region = "Brazil")

models_brazil <- fitted_models_Brazil %>%
  filter(term == "log10(IndividualBodyMass)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate,std.error,statistic,conf.high,conf.low,p.value) %>%
  mutate(estimate = round(estimate,digits = 2),
         std.error = round(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         conf.high = round(conf.high, digits = 2),
         conf.low = round(conf.low, digits = 2),
         p.value = format.pval(p.value, digits = 2)) %>%
  rename(Slope = estimate) 

#add intercept 
intercepts <- fitted_models_Brazil%>%
  filter(term == "(Intercept)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate) %>%
  rename(Intercept = estimate)

models_brazil <- left_join(models_brazil, intercepts)

###--------###
### Namibia ###
###--------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","Namibia",sep="/","Output_Summary.Rdata"))
outpath <- myoutpath
figpath <- myfigpath

regionpath <- paste0(figpath,sep="/",dir.create(file.path(figpath, "Namibia"), showWarnings = FALSE),sep="/")

#get output directory path
out_dir <- outpath 

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#just load data again 
#list[1] = control; list[2] = HANPP; last list entry [individual on region] = After vegreduction
#climate
cohorts_histo <- data.frame(historical_2014_list[[1]]$cohorts) 
cohorts_SSP126 <- data.frame(SSP126_2100_list[[1]]$cohorts)
cohorts_SSP585 <- data.frame(SSP585_2100_list[[1]]$cohorts) 

#HANPP
cohorts_histo_HANPP <- data.frame(historical_2014_list[[2]]$cohorts) 
cohorts_SSP126_HANPP <- data.frame(SSP126_2100_list[[2]]$cohorts)
cohorts_SSP585_HANPP <- data.frame(SSP585_2100_list[[2]]$cohorts) 

#Max. HANPP
cohorts_histo_maxHANPP <- data.frame(historical_2014_list[[13]]$cohorts) 
cohorts_SSP126_maxHANPP <- data.frame(SSP126_2100_list[[13]]$cohorts)
cohorts_SSP585_maxHANPP <- data.frame(SSP585_2100_list[[13]]$cohorts) 

#subset data 
#climate
cohorts_histo <- cohorts_histo[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126 <- cohorts_SSP126[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585 <- cohorts_SSP585[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#HANPP
cohorts_histo_HANPP <- cohorts_histo_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_HANPP <- cohorts_SSP126_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_HANPP <- cohorts_SSP585_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#maxHANPP
cohorts_histo_maxHANPP <- cohorts_histo_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_maxHANPP <- cohorts_SSP126_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_maxHANPP <- cohorts_SSP585_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#create new grouping variable for 
#climate
cohorts_histo$group = "Historical"
cohorts_SSP126$group = "SSP1-2.6"
cohorts_SSP585$group = "SSP5-8.5"
cohorts_histo$group2 = "Climate"
cohorts_SSP126$group2 = "Climate"
cohorts_SSP585$group2 = "Climate"

#HANPP
cohorts_histo_HANPP$group = "Historical"
cohorts_SSP126_HANPP$group = "SSP1-2.6"
cohorts_SSP585_HANPP$group = "SSP5-8.5"
cohorts_histo_HANPP$group2 = "Current Land Use"
cohorts_SSP126_HANPP$group2 = "Current Land Use"
cohorts_SSP585_HANPP$group2 = "Current Land Use"

#maxHANPP
cohorts_histo_maxHANPP$group = "Historical"
cohorts_SSP126_maxHANPP$group = "SSP1-2.6"
cohorts_SSP585_maxHANPP$group = "SSP5-8.5"
cohorts_histo_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP126_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP585_maxHANPP$group2 = "Maximum Land Use"

#combine dataframes for plotting dataset
test <- rbind(cohorts_histo,cohorts_SSP126,cohorts_SSP585,
              cohorts_histo_HANPP,cohorts_SSP126_HANPP,cohorts_SSP585_HANPP,
              cohorts_histo_maxHANPP,cohorts_SSP126_maxHANPP,cohorts_SSP585_maxHANPP)

#convert to kg 
test[3] <- test[3]/1000 

#change naming
names(test) <- c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance","group","group2")

#FG naming
#FG naming
if(FGs == "all") {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores (i.)", '1' = "Endotherm Carnivores (i.)", '2' =  "Endotherm Omnivores (i.)",
                                      '3'="Ectotherm Herbivores (s.)",'4'="Ectotherm Carnivores (s.)",'5'="Ectotherm Omnivores (s.)","6" = "Ectotherm Herbivores (i.)",
                                      '7' = "Ectotherm Carnivores (i.)",'8'="Ectotherm Omnivores (i.)" )
} else {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores", '1' = "Endotherm Carnivores", '2' =  "Endotherm Omnivores",
                                      '3'="Ectotherm Herbivores",'4'="Ectotherm Carnivores",'5'="Ectotherm Omnivores","6" = "Ectotherm Herbivores",
                                      '7' = "Ectotherm Carnivores",'8'="Ectotherm Omnivores" )
  
}

#manual model fitting
fitted_models_Namibia <- test %>%
  group_by(FunctionalGroupIndex,group, group2) %>%
  do(broom::tidy(lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .),conf.int=T)) %>%
  ungroup() %>%
  mutate(Region = "Namibia")

models_namibia <- fitted_models_Namibia %>%
  filter(term == "log10(IndividualBodyMass)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate,std.error,statistic,conf.high,conf.low,p.value) %>%
  mutate(estimate = round(estimate,digits = 2),
         std.error = round(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         conf.high = round(conf.high, digits = 2),
         conf.low = round(conf.low, digits = 2),
         p.value = format.pval(p.value, digits = 2)) %>%
  rename(Slope = estimate) 

#add intercept 
intercepts <- fitted_models_Namibia %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate) %>%
  rename(Intercept = estimate)

models_namibia <- left_join(models_namibia, intercepts)

###--------###
### France ###
###--------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","France",sep="/","Output_Summary.Rdata"))
outpath <- myoutpath
figpath <- myfigpath

regionpath <- paste0(figpath,sep="/",dir.create(file.path(figpath, "France"), showWarnings = FALSE),sep="/")

#get output directory path
out_dir <- outpath 

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#just load data again 
#list[1] = control; list[2] = HANPP; last list entry [individual on region] = After vegreduction
#climate
cohorts_histo <- data.frame(historical_2014_list[[1]]$cohorts) 
cohorts_SSP126 <- data.frame(SSP126_2100_list[[1]]$cohorts)
cohorts_SSP585 <- data.frame(SSP585_2100_list[[1]]$cohorts) 

#HANPP
cohorts_histo_HANPP <- data.frame(historical_2014_list[[2]]$cohorts) 
cohorts_SSP126_HANPP <- data.frame(SSP126_2100_list[[2]]$cohorts)
cohorts_SSP585_HANPP <- data.frame(SSP585_2100_list[[2]]$cohorts) 

#Max. HANPP
cohorts_histo_maxHANPP <- data.frame(historical_2014_list[[8]]$cohorts) 
cohorts_SSP126_maxHANPP <- data.frame(SSP126_2100_list[[8]]$cohorts)
cohorts_SSP585_maxHANPP <- data.frame(SSP585_2100_list[[8]]$cohorts) 

#subset data 
#climate
cohorts_histo <- cohorts_histo[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126 <- cohorts_SSP126[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585 <- cohorts_SSP585[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#HANPP
cohorts_histo_HANPP <- cohorts_histo_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_HANPP <- cohorts_SSP126_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_HANPP <- cohorts_SSP585_HANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#maxHANPP
cohorts_histo_maxHANPP <- cohorts_histo_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP126_maxHANPP <- cohorts_SSP126_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]
cohorts_SSP585_maxHANPP <- cohorts_SSP585_maxHANPP[ ,c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance")]

#create new grouping variable for 
#climate
cohorts_histo$group = "Historical"
cohorts_SSP126$group = "SSP1-2.6"
cohorts_SSP585$group = "SSP5-8.5"
cohorts_histo$group2 = "Climate"
cohorts_SSP126$group2 = "Climate"
cohorts_SSP585$group2 = "Climate"

#HANPP
cohorts_histo_HANPP$group = "Historical"
cohorts_SSP126_HANPP$group = "SSP1-2.6"
cohorts_SSP585_HANPP$group = "SSP5-8.5"
cohorts_histo_HANPP$group2 = "Current Land Use"
cohorts_SSP126_HANPP$group2 = "Current Land Use"
cohorts_SSP585_HANPP$group2 = "Current Land Use"

#maxHANPP
cohorts_histo_maxHANPP$group = "Historical"
cohorts_SSP126_maxHANPP$group = "SSP1-2.6"
cohorts_SSP585_maxHANPP$group = "SSP5-8.5"
cohorts_histo_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP126_maxHANPP$group2 = "Maximum Land Use"
cohorts_SSP585_maxHANPP$group2 = "Maximum Land Use"

#combine dataframes for plotting dataset
test <- rbind(cohorts_histo,cohorts_SSP126,cohorts_SSP585,
              cohorts_histo_HANPP,cohorts_SSP126_HANPP,cohorts_SSP585_HANPP,
              cohorts_histo_maxHANPP,cohorts_SSP126_maxHANPP,cohorts_SSP585_maxHANPP)

#convert to kg 
test[3] <- test[3]/1000 

#change naming
names(test) <- c("GridcellIndex","FunctionalGroupIndex","IndividualBodyMass","CohortAbundance","group","group2")

#FG naming
#FG naming
if(FGs == "all") {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores (i.)", '1' = "Endotherm Carnivores (i.)", '2' =  "Endotherm Omnivores (i.)",
                                      '3'="Ectotherm Herbivores (s.)",'4'="Ectotherm Carnivores (s.)",'5'="Ectotherm Omnivores (s.)","6" = "Ectotherm Herbivores (i.)",
                                      '7' = "Ectotherm Carnivores (i.)",'8'="Ectotherm Omnivores (i.)" )
} else {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores", '1' = "Endotherm Carnivores", '2' =  "Endotherm Omnivores",
                                      '3'="Ectotherm Herbivores",'4'="Ectotherm Carnivores",'5'="Ectotherm Omnivores","6" = "Ectotherm Herbivores",
                                      '7' = "Ectotherm Carnivores",'8'="Ectotherm Omnivores" )
  
}

#manual model fitting
fitted_models_France <- test %>%
  group_by(FunctionalGroupIndex,group, group2) %>%
  do(broom::tidy(lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .),conf.int=T)) %>%
  ungroup() %>%
  mutate(Region = "France")

models_france <- fitted_models_France %>%
  filter(term == "log10(IndividualBodyMass)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate,std.error,statistic,conf.high,conf.low,p.value) %>%
  mutate(estimate = round(estimate,digits = 2),
         std.error = round(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         conf.high = round(conf.high, digits = 2),
         conf.low = round(conf.low, digits = 2),
         p.value = format.pval(p.value, digits = 2)) %>%
  rename(Slope = estimate) 

#add intercept 
intercepts <- fitted_models_France %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(FunctionalGroupIndex,Region,group,group2,estimate) %>%
  rename(Intercept = estimate)

models_france <- left_join(models_france, intercepts)

models_all <- rbind(models_brazil, models_finland, models_namibia, models_france)

#reorder and rename columns
models_all <- models_all[,c(1,2,3,4,11,5,6,7,8,9,10)]

colnames(models_all) <- c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "F Statistic", "Upper CI", "Lower CI", "P Value")

models_all$Intercept <- round(models_all$Intercept, digits = 2)

#html table output
models_all %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "F Statistic", "Upper CI", "Lower CI", "P Value"))

climate <- subset(models_all, Scenario == "Climate")
climate %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "F Statistic", "Upper CI", "Lower CI", "P Value"))

curr_lu <- subset(models_all, Scenario == "Current Land Use")
curr_lu %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "F Statistic", "Upper CI", "Lower CI", "P Value"))

max_lu <- subset(models_all, Scenario == "Maximum Land Use")
max_lu %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "F Statistic", "Upper CI", "Lower CI", "P Value"))

# Create a matrix with NA values
heatmap_matrix <- matrix(NA, nrow = length(unique(paste(models_all$FunctionalGroupIndex, models_all$Region))), ncol = length(unique(paste(models_all$Climate, models_all$Scenario))))

# Fill the matrix with slope values
for (i in 1:nrow(models_all)) {
  row_index <- which(unique(paste(models_all$FunctionalGroupIndex), models_all$Region) == paste(models_all$FunctionalGroupIndex[i], models_all$Region[i]))
  col_index <- which(colnames(heatmap_matrix) == paste(models_all$Climate[i], models_all$Scenario[i]))
}

# Set row and column names
row_names <- unique(paste(models_all$FunctionalGroupIndex, models_all$Region))
#row_names <- rep(paste(unique(models_all$FunctionalGroupIndex)), 4)
col_names <- unique(paste(models_all$Climate, models_all$Scenario))
colnames(heatmap_matrix) <- col_names
rownames(heatmap_matrix) <- row_names

# Fill the matrix with slope values
for (i in 1:nrow(models_all)) {
  row_index <- which(unique(paste(models_all$FunctionalGroupIndex, models_all$Region)) == paste(models_all$FunctionalGroupIndex[i], models_all$Region[i]))
  col_index <- which(colnames(heatmap_matrix) == paste(models_all$Climate[i], models_all$Scenario[i]))
  print(paste("Row index:", row_index))
  print(paste("Col index:", col_index))
  print(paste("Slope value:", models_all$Slope[i]))
  heatmap_matrix[row_index, col_index] <- models_all$Slope[i]
}

col_order <- c("Historical Climate", "SSP1-2.6 Climate", "SSP5-8.5 Climate",
               "Historical Current Land Use", "SSP1-2.6 Current Land Use", "SSP5-8.5 Current Land Use",
               "Historical Maximum Land Use", "SSP1-2.6 Maximum Land Use", "SSP5-8.5 Maximum Land Use")

heatmap_matrix <- heatmap_matrix[, col_order]

if(FGs == "all") {
  rownames(heatmap_matrix) <- str_replace_all(rownames(heatmap_matrix), "Brazil|Finland|Namibia|France", c(rep("a",9),rep("b", 9),rep("c",9),rep("d",9)))
} else {
  rownames(heatmap_matrix) <- str_replace_all(rownames(heatmap_matrix), "Brazil|Finland|Namibia|France", c(rep("a",6),rep("b", 6),rep("c",6),rep("d",6)))
}
colnames(heatmap_matrix) <- str_replace_all(colnames(heatmap_matrix), "Climate|Current Land Use|Maximum Land Use", c(rep("1",3),rep("2", 3),rep("3", 3)))

if(FGs == "all") {
  annot_row <- data.frame(c(rep("Brazil", 9), rep("Finland", 9), rep("Namibia", 9), rep("France", 9)))
} else {
  annot_row <- data.frame(c(rep("Brazil", 6), rep("Finland", 6), rep("Namibia", 6), rep("France", 6)))
}
colnames(annot_row) <- c("Region")
annot_row$Region <- as.factor(annot_row$Region)
rownames(annot_row) <- rownames(heatmap_matrix)

annot_col <- data.frame(c(rep("Climate", 3), rep("Current Land Use", 3), rep("Maximum Land Use", 3)))
colnames(annot_col) <- c("Scenario")
annot_col$Scenario <- as.factor(annot_col$Scenario)
rownames(annot_col) <- colnames(heatmap_matrix)

if(FGs == "all") {
  row_labels <- rep(c("Ectotherm  Carnivores (it.)", "Ectotherm  Carnivores (s.)", "Ectotherm  Herbivores (it.)", "Ectotherm  Herbivores (s.)", 
                      "Ectotherm  Omnivores (it.)", "Ectotherm  Omnivores (s.)", "Endotherm Carnivores", "Endotherm Herbivores", "Endotherm Omnivores"),4)
} else {
  row_labels <- rep(c("Ectotherm  Carnivores", "Ectotherm  Herbivores", 
                      "Ectotherm  Omnivores", "Endotherm Carnivores", 
                      "Endotherm Herbivores", "Endotherm Omnivores"),4)
}

#row_labels <- rep(c("Ectotherm  Carnivores (it.)", "Ectotherm  Carnivores (s.)", "Ectotherm  Herbivores (it.)", "Ectotherm  Herbivores (s.)", "Ectotherm  Omnivores (it.)", "Ectotherm  Omnivores (s.)", "Endotherm Carnivores", "Endotherm Herbivores", "Endotherm Omnivores"),4)

col_labels <- c(rep(c("Historical", "SSP1-2.6", "SSP5-8.5"),3))

pdf(paste0(figpath,"/",FGs,"_","Heatmap_Slopes.pdf"),width = 20, height = 20, paper="special")

# Plot heatmap
pheatmap(heatmap_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_row = annot_row,
         annotation_col = annot_col,
         labels_row = row_labels,
         labels_col = col_labels,
         show_colnames = T,
        main = "Slope Values",
         annotation_colors = list(Region = c(Brazil = "#DCDCDC", Finland = "#C0C0C0", Namibia = "#808080", France = "#696969"),
                                  Scenario = c(Climate = "#DCDCDC", `Current Land Use` = "#C0C0C0", `Maximum Land Use` = "#808080")),
         border_color = NA,
         drop_levels  = TRUE,
         fontsize = 25,
         #scale = "row",
         number_color = ifelse(heatmap_matrix > 0, "black", "white"),
         display_numbers = T,
         color = inferno(10000))

# Add a legend title manually
#legend_title <- "Z-Score Slope"
#grid.text(legend_title, x = 0.78, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


###ggplot heatmap: 

test3 <- reshape::melt(heatmap_matrix)
test3$Region <- annot_row$Region[match(test3$X1, rownames(annot_row))]
test3$"Simulation Experiment" <- annot_col$Scenario[match(test3$X2, rownames(annot_col))]

#using ggh4x package here, allows for in depth customization of facet grids 
strip_variation <- ggh4x::strip_themed(
  #NULL,
  # Horizontal strips
  background_x = elem_list_rect(fill = c("#DCDCDC", "#C0C0C0", "#808080"), colour = c(rep("white",3))),
  text_x = elem_list_text(face = c("bold")),
  by_layer_x = FALSE,
  # Vertical strips
  background_y = elem_list_rect(fill = c("white"), colour = c(rep("grey",4))),
  text_y = elem_list_text(face = c("bold")),
  by_layer_y = FALSE)

pdf(paste0(figpath,"/",FGs,"_","Heatmap_Slopes_ggplot.pdf"),width = 10, height = 10, paper="special")

#ggplot heatmap with geom_tile & facet_grid
ggplot(test3, aes(x = X2, y = X1, fill = value)) +
  geom_tile() +  
  geom_text(aes(label = value), color = ifelse(test3$value > 0, "black", "white"), size = 4) +
  #coord_fixed()
  scale_fill_viridis(option = "inferno", breaks = seq(-6, 1, 1)) +
  ggh4x::facet_grid2(rows = vars(Region),cols = vars(`Simulation Experiment`), scales = "free", space = "free", switch = "x", strip = strip_variation) +
  labs(fill = "Slope Value", y = "Functional Group", x = "Simulation Experiment") +
  scale_y_discrete(labels = rev(row_labels), limits = rev, expand = c(0,0)) +
  scale_x_discrete(labels = col_labels, expand = c(0,0)) +
  theme(axis.text.y = element_text(size=12,hjust=0))+
  #geom_vline(xintercept=c(3.5,6.5), color='black', size = 1) +
  theme_classic()+
  theme(legend.position = "top") +
  theme(panel.spacing=unit(0,units = "cm")) +
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_rect(fill = "white"),
        axis.title = element_blank())  +
  theme(axis.text.y = element_text(size=12,hjust=0))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=14))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=12,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=12,face="bold")) +
  theme(legend.text = element_text(size=12),legend.title = element_text(size=12,face="bold")) +       
  theme(legend.margin=margin(grid::unit(0, "cm")), legend.key.height=grid::unit(0.3, "cm"), legend.key.width=grid::unit(1, "cm")) 
dev.off()

