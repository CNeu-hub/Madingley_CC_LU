###-----------------------------------------------------------------------------------------###
### TITLE:       Heatmap plot of slopes                                                     ###
### DESCRIPTION: Plots figure 4 and figure S3.1 of:                                         ###
###              "Model-based Impact Analysis of Climate Change and Land-use Intensification###
###              on Trophic Networks."                                                      ###                                                                 
### PROCEDURE:   Performs combined spatial autoregression, autocorrelation tests, neighbour ###
###              creation for spatial regression, and extracts model results to display     ###
###              slopes & significance in a heatmap plot                                    ###
### DATE:        12.08.2024                                                                 ###
###-----------------------------------------------------------------------------------------###

###LOAD LIBRARIES
library(tidyverse)
library(RColorBrewer)
library(spatialreg) #package for performing combined spatial autoregression (SAC)
library(spdep) #package to create neighbour objects
library(ggh4x) #needed for heatmap (labeling of facets)

#define input path (Output_Data in repository)
input <- paste(getwd(), "Input", sep = "/")

#define output path (Output_CSV in repository)
outpath <- paste(getwd(), "Output", sep = "/")

#Activate aggregated FGs, or deactivate (Results for all FGs will be plotted)
#chose either "aggregated", or "all" 
#if aggregated is chosen, mean abundance and body mass for ectotherms will be calculated summarizing iteroparous and semelparous ectotherm functional groups
FGs <- "aggregated"

#Define if impacts for spatial sac model should be calculated (if impact == "Impacts") than total effect of impacts() function is used as slope coefficient (incorporating direct and indirect spatial effects into consideration)
impact <- "none"

#set seed to allow reproducibility of random processes (Morans I Monte Carlo Simulations, etc.)
set.seed(1234)

###--------###
### Finland ###
###--------###
#load workspace images from input 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(input,sep="/","Finland",sep="/","Output_Summary.Rdata"))

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
  test <- test %>%
    group_by(GridcellIndex, FunctionalGroupIndex, group, group2) %>%
    summarise(
      CohortAbundance = mean(CohortAbundance, na.rm = TRUE),
      IndividualBodyMass = mean(IndividualBodyMass, na.rm = TRUE)
    ) %>% ungroup()
}

###create grouped df with underlying data & models information we need to perform loops later on for identifying neighbours!
lm_fit_finland <- test %>%
  group_by(FunctionalGroupIndex, group, group2) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .x))) %>%
  ungroup() %>%
  mutate(Region = "Finland")

#maintain df order for plotting later
#lm_fit_finland <- lm_fit_finland %>%
#  group_by(FunctionalGroupIndex, group, group2)

###-----------------------------------------------------------------###
###  2.                      NEIGHBOUR LISTS                        ###
###-----------------------------------------------------------------###
###Step 1: Identify Neighbours 
#create neigbours object for gridcells based on spatial grid dimensions
#spatial_window <- c(25.5,28.5,61,69) #96 grid cells, dimensions Finland
min.x <- 25.5
max.x <- 28.5
min.y <- 61
max.y <- 69

#define cell size (0.5 degree in our case)
cellsize <- 0.5

#calculate x y dimensions of grid
cells.dim.x <- (max.x - min.x) / cellsize
cells.dim.y <- (max.y - min.y) / cellsize

#calculate cells dimensions (combine x y information)
cells.dim <- c(cells.dim.x, cells.dim.y)

#create neigbours object
neighbours <- spdep::grid2nb(d = cells.dim)

###-----------------------------------------------------------------###
###  3.          AUTOCORRELATION: STATISTICS AND MODELS             ###
###-----------------------------------------------------------------###
###from here i loop over datasets
output <- vector("list", length = length(lm_fit_finland$data))
output <- lapply(output, function(x) list(group = NULL, lm = NULL, MoransI = NULL, MoransI_resid = NULL, Rao_lag = NULL, Rao_error = NULL, Spa_lag = NULL, Spa_error = NULL, Spa_sac = NULL, Spa_sac_effects = NULL, Spa_lac_effects = NULL, AIC = NULL))

for(i in seq_along(lm_fit_finland$data)) {
  
  dat <- lm_fit_finland$data[[i]]
  #dat$IndividualBodyMass <- log10(dat$IndividualBodyMass)
  #dat$CohortAbundance <- log10(dat$CohortAbundance)
  
  #create group output
  output[[i]]$group <- c(lm_fit_finland$FunctionalGroupIndex[[i]], lm_fit_finland$group[[i]], lm_fit_finland$group2[[i]], lm_fit_finland$Region[[i]])
  
  #ifelse condition (some combinations (e.g. ectotherm carnivores france max land use), have too less observations so that they dont have neigbours, thus I exclude all combinations with < 5 cells/observations from modelling)
  if(length(dat$GridcellIndex) < 5) {
    
    cat(lm_fit_france$FunctionalGroupIndex[[i]], "has to less observations, not used for modelling.")
    
  } else{
    
    ###1. Subset for each group in i neigbours list (because sometimes we have missing cells in between)
    #create logical cells vector (dimensions like neigbours list, but with FALSE as default)
    observed_cells <- rep(FALSE, length(neighbours))
    
    #add up +1 to gridcellindex, because it starts at 0 (is not matching neighbours list dimensions because R starts at 1), set all cells that can be found in both to TRUE
    observed_cells[dat$GridcellIndex+1] <- TRUE
    
    #calculate no. of missing cells per i 
    missing_cells <- which(!observed_cells)
    cat("Missing cells:", sum(n_distinct(missing_cells)), "\n")
    
    #subset neighbours list to include only observed cells
    nb <- spdep::subset.nb(neighbours, subset = observed_cells)
    
    #create spatial weights object 
    listw <- spdep::nb2listw(nb, style = "W")
    
    cat("Cell correction done, neigbours created. Calculation of Morans I and spatial regression.", "Step:", i, length(lm_fit_finland$data))
    
    fit.lm <- lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = dat)
    output[[i]]$lm <- fit.lm
    
    ###2. calculate global morans I (for dependent variable, and residuals separately) and raos score for spatial error and spatial lag model
    #test for autocorrelation of dependent variable if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI <- spdep::moran.mc(dat$CohortAbundance, listw, nsim = 999)
    
    #test for autocorrelation of residuals if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI_resid <- (spdep::lm.morantest(fit.lm, listw, zero.policy = T)) #test = significant --> spatial autocorrelation detected
    
    #test if spatial error or spatial lag model more appropriate with Raos Score (Lagrange multiplier test)
    output[[i]]$Rao_error <- spdep::lm.RStests(fit.lm, listw, test="RSerr", zero.policy = T) # test not significant --> spatial lag model more appropriate
    
    output[[i]]$Rao_lag <- spdep::lm.RStests(fit.lm, listw, test="RSlag", zero.policy = T) # test not significant --> spatial error model more appropriate
    
    ###3. Perform spatial regression (spatial lag model)
    output[[i]]$Spa_error <- spatialreg::errorsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                        data = dat, 
                                        listw = listw) 
    
    Spa_lag <- spatialreg::lagsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    Spa_sac <- spatialreg::sacsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    output[[i]]$Spa_lag <- Spa_lag
    output[[i]]$Spa_sac <- Spa_sac
    
    #calculate effects/impacts (indirect, direct, total) for spatial lag model 
    output[[i]]$Spa_lag_effects <- impacts(Spa_lag, listw = listw, R = 999)
    output[[i]]$Spa_sac_effects <- impacts(Spa_sac, listw = listw, R = 999)
    
  }}


AIC_list <- vector("list", length(output))

for(i in seq_along(output)) {
  
  AIC_values <- c(lm = AIC(output[[i]]$lm),
                  Spa_error = AIC(output[[i]]$Spa_error),
                  Spa_lag = AIC(output[[i]]$Spa_lag),
                  Spa_sac = AIC(output[[i]]$Spa_sac))
  AIC_list[[i]] <- AIC_values
}

AIC_df <- do.call(rbind, AIC_list)

AIC_df <- data.frame(Index = seq_along(output), AIC_df)

#create vector to calculate which model has lowest AIC
min_AIC_column <- apply(AIC_df[, -1], 1, function(x) names(x)[which.min(x)])

#store minimum aic values
AIC_df$min_AIC <- min_AIC_column

#count how often for each model AIC was minimum 
min_AIC_count <- table(AIC_df$min_AIC)

#store results in txt.file 
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_AICs_Finland", FGs, ".txt", sep = ""))

# Print the counts
print(min_AIC_count)

sink()

###-----------------------------------------------------------------###
###  4.          AUTOCORRELATION: DIAGNOSTICS OUTPUT.               ###
###-----------------------------------------------------------------###
#calculate how many cases are autocorrelated and how many cases have significant raos score to determine which model to use
determine_Moran_dependent <- vector(length = length(output))
determine_Moran_resid <- vector(length = length(output))
determine_Rao_error <- vector(length = length(output))
determine_Rao_lag <- vector(length = length(output))

for(i in seq_along(output)) {
  
  if(!is.null(output[[i]]$MoransI)) {
    print(output[[i]]$MoransI)
    print(output[[i]]$MoransI_resid)
    print(output[[i]]$Rao_error)
    print(output[[i]]$Rao_lag)
    
    determine_Moran_dependent[i] <- ifelse(output[[i]]$MoransI$p.value < 0.05, 1, 0)
    determine_Moran_resid[i] <- ifelse(output[[i]]$MoransI_resid$p.value < 0.05, 1, 0)
    determine_Rao_error[i] <- ifelse(output[[i]]$Rao_error$RSerr$p.value < 0.05, 1, 0)
    determine_Rao_lag[i] <- ifelse(output[[i]]$Rao_lag$RSlag$p.value < 0.05, 1, 0)
  }
}

sum(determine_Moran_dependent) #59 dependent variables show autocorrelation
sum(determine_Moran_resid) #50 models show autocorrelation in residuals
sum(determine_Rao_error) #39 models show significance (spatial error model is preferred)
sum(determine_Rao_lag) #44 models show significance (spatial lag model is preferred)

###Export autocorrelation diagnostics summary
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_Finland", FGs, ".txt", sep = ""))
print(paste("Morans I Monte Carlo Simulation:", sum(determine_Moran_dependent), "of", length(output), "P-Values are significant"))
print(paste("Morans I for Spatial Autocorrelation in Residuals:", sum(determine_Moran_resid), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Error Model):", sum(determine_Rao_error), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Lag Model):", sum(determine_Rao_lag), "of", length(output), "P-Values are significant"))
sink()

###-----------------------------------------------------------------###
###  5.                    PLOTTING SLOPES                          ###
###-----------------------------------------------------------------###

#create dataframe for plotting spatial lag results
final_data_sp_sac_finl <- data.frame("FunctionalGroupIndex" = character(), "Climate" = character(), "Scenario" = character(), "Region" = character(),"Intercept" = numeric(), "Slope" = numeric(), "Std.Err" = numeric(), "Z.Stat" = numeric(), "P.Val" = numeric(), stringsAsFactors = FALSE)

#extract output for plotting
for(i in seq_along(output)) { 
  
  if(!is.null(output[[i]]$MoransI)) {
    
    if(impact == "Impacts") {
    sums <- summary(output[[i]]$Spa_sac_effects, short = T, zstats = T)
    
    final_data <- data.frame("FunctionalGroupIndex" = c(output[[i]]$group[1]),
                             "Climate" = c(output[[i]]$group[2]),
                             "Scenario" = c(output[[i]]$group[3]),
                             "Region" = c(output[[i]]$group[4]),
                             "Intercept" = c(round(sums$total_sum$statistics[1])),
                             "Slope" = c(round(sums$res$total[[1]], digits = 2)),
                             "Std.Err" = c(round(sums$semat[[3]], digits = 2)),
                             "Z.Stat" = round(sums$zmat[[3]], digits = 2),
                             "P.Val" = format.pval(sums$pzmat[[3]], digits = 2))
    
    final_data_sp_sac_finl <- rbind(final_data_sp_sac_finl, final_data)
    } else {
      sums <- summary(output[[i]]$Spa_sac)
      summary(output[[1]]$Spa_sac)
      final_data <- data.frame(
        "FunctionalGroupIndex" = output[[i]]$group[1],
        "Climate" = output[[i]]$group[2],
        "Scenario" = output[[i]]$group[3],
        "Region" = output[[i]]$group[4],
        "Intercept" = round(sums$Coef[1], digits = 2),
        "Slope" = round(sums$Coef[2], digits = 2),
        "Std.Err" = round(sums$Coef[4], digits = 2),
        "Z.Stat" = round(sums$Coef[6], digits = 2),
        "P.Val" = format.pval(sums$Coef[8], digits = 2)
      )
      final_data_sp_sac_finl <- rbind(final_data_sp_sac_finl, final_data)
    }
  }
}

###--------###
### France ###
###--------###
#load workspace images from input 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(input,sep="/","France",sep="/","Output_Summary.Rdata"))

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
if(FGs == "all") {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores (i.)", '1' = "Endotherm Carnivores (i.)", '2' =  "Endotherm Omnivores (i.)",
                                      '3'="Ectotherm Herbivores (s.)",'4'="Ectotherm Carnivores (s.)",'5'="Ectotherm Omnivores (s.)","6" = "Ectotherm Herbivores (i.)",
                                      '7' = "Ectotherm Carnivores (i.)",'8'="Ectotherm Omnivores (i.)" )
  
  
} else {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores", '1' = "Endotherm Carnivores", '2' =  "Endotherm Omnivores",
                                      '3'="Ectotherm Herbivores",'4'="Ectotherm Carnivores",'5'="Ectotherm Omnivores","6" = "Ectotherm Herbivores",
                                      '7' = "Ectotherm Carnivores",'8'="Ectotherm Omnivores" )
  test <- test %>%
    group_by(GridcellIndex, FunctionalGroupIndex, group, group2) %>%
    summarise(
      CohortAbundance = mean(CohortAbundance, na.rm = TRUE),
      IndividualBodyMass = mean(IndividualBodyMass, na.rm = TRUE)
    ) %>% ungroup()
}

###create grouped df with underlying data & models information we need to perform loops later on for identifying neighbours!
lm_fit_france <- test %>%
  group_by(FunctionalGroupIndex, group, group2) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .x))) %>%
  ungroup() %>%
  mutate(Region = "France")

#maintain df order for plotting later
#lm_fit_france <- lm_fit_france %>%
#  group_by(FunctionalGroupIndex, group, group2)

###-----------------------------------------------------------------###
###  2.                      NEIGHBOUR LISTS                        ###
###-----------------------------------------------------------------###
###Step 1: Identify Neighbours 
#create neigbours object for gridcells based on spatial grid dimensions
#spatial_window <- c(-1,6,46,49) #84 grid cells, France
min.x <- -1
max.x <- 6
min.y <- 46
max.y <- 49

#define cell size (0.5 degree in our case)
cellsize <- 0.5

#calculate x y dimensions of grid
cells.dim.x <- (max.x - min.x) / cellsize
cells.dim.y <- (max.y - min.y) / cellsize

#calculate cells dimensions (combine x y information)
cells.dim <- c(cells.dim.x, cells.dim.y)

#create neigbours object
neighbours <- spdep::grid2nb(d = cells.dim)

###-----------------------------------------------------------------###
###  3.          AUTOCORRELATION: STATISTICS AND MODELS             ###
###-----------------------------------------------------------------###
###from here i loop over datasets
output <- vector("list", length = length(lm_fit_france$data))
output <- lapply(output, function(x) list(group = NULL, lm = NULL, MoransI = NULL, MoransI_resid = NULL, Rao_lag = NULL, Rao_error = NULL, Spa_lag = NULL, Spa_error = NULL))

for(i in seq_along(lm_fit_france$data)) {
 # dat <- lm_fit_france$data[[52]]
  dat <- lm_fit_france$data[[i]]
  #dat$IndividualBodyMass <- log10(dat$IndividualBodyMass)
  #dat$CohortAbundance <- log10(dat$CohortAbundance)
  
  #create group output
  output[[i]]$group <- c(lm_fit_france$FunctionalGroupIndex[[i]], lm_fit_france$group[[i]], lm_fit_france$group2[[i]], lm_fit_france$Region[[i]])
  
  #ifelse condition (some combinations (e.g. ectotherm carnivores france max land use), have too less observations so that they dont have neigbours, thus I exclude all combinations with < 5 cells/observations from modelling)
  if(length(dat$GridcellIndex) < 5) {
    
    cat(lm_fit_france$FunctionalGroupIndex[[i]], "has to less observations, not used for modelling.")
    
  } else{
    
    ###1. Subset for each group in i neigbours list (because sometimes we have missing cells in between)
    #create logical cells vector (dimensions like neigbours list, but with FALSE as default)
    observed_cells <- rep(FALSE, length(neighbours))
    
    #add up +1 to gridcellindex, because it starts at 0 (is not matching neighbours list dimensions because R starts at 1), set all cells that can be found in both to TRUE
    observed_cells[dat$GridcellIndex+1] <- TRUE
    
    #calculate no. of missing cells per i 
    missing_cells <- which(!observed_cells)
    cat("Missing cells:", sum(n_distinct(missing_cells)), "\n")
    
    #subset neighbours list to include only observed cells
    nb <- spdep::subset.nb(neighbours, subset = observed_cells)
    
    #create spatial weights object 
    listw <- spdep::nb2listw(nb, style = "W")
    
    cat("Cell correction done, neigbours created. Calculation of Morans I and spatial regression.", "Step:", i, length(lm_fit_france$data))
    
    fit.lm <- lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = dat)
    output[[i]]$lm <- fit.lm
    
    ###2. calculate global morans I (for dependent variable, and residuals separately) and raos score for spatial error and spatial lag model
    #test for autocorrelation of dependent variable if p < 0.05 then autocorrelation detected
    
    output[[i]]$MoransI <- spdep::moran.mc(dat$CohortAbundance, listw, nsim = 999)
    
    #test for autocorrelation of residuals if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI_resid <- (spdep::lm.morantest(fit.lm, listw, zero.policy = T)) #test = significant --> spatial autocorrelation detected
    
    #test if spatial error or spatial lag model more appropriate with Raos Score (Lagrange multiplier test)
    output[[i]]$Rao_error <- spdep::lm.RStests(fit.lm, listw, test="RSerr", zero.policy = T) # test not significant --> spatial lag model more appropriate
    
    output[[i]]$Rao_lag <- spdep::lm.RStests(fit.lm, listw, test="RSlag", zero.policy = T) # test not significant --> spatial error model more appropriate
    
    ###3. Perform spatial regression (spatial lag model)
    output[[i]]$Spa_error <- spatialreg::errorsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                                    data = dat, 
                                                    listw = listw) 
    
    Spa_lag <- spatialreg::lagsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    Spa_sac <- spatialreg::sacsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    output[[i]]$Spa_lag <- Spa_lag
    output[[i]]$Spa_sac <- Spa_sac

    #calculate effect#calculate effect#calculate effects/impacts (indirect, direct, total) for spatial lag model 
    output[[i]]$Spa_lag_effects <- impacts(Spa_lag, listw = listw, R = 999)
    output[[i]]$Spa_sac_effects <- impacts(Spa_sac, listw = listw, R = 999)
    
  }}


AIC_list <- vector("list", length(output))

for(i in seq_along(output)) {
  
  if(!is.null(output[[i]]$lm)) {
  AIC_values <- c(lm = AIC(output[[i]]$lm),
                  Spa_error = AIC(output[[i]]$Spa_error),
                  Spa_lag = AIC(output[[i]]$Spa_lag),
                  Spa_sac = AIC(output[[i]]$Spa_sac))
  AIC_list[[i]] <- AIC_values 
  }
}

AIC_df <- do.call(rbind, AIC_list)

AIC_df <- data.frame(AIC_df)

#create vector to calculate which model has lowest AIC
min_AIC_column <- apply(AIC_df[, -1], 1, function(x) names(x)[which.min(x)])

#store minimum aic values
AIC_df$min_AIC <- min_AIC_column

#count how often for each model AIC was minimum 
min_AIC_count <- table(AIC_df$min_AIC)

#store results in txt.file 
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_AICs_France", FGs, ".txt", sep = ""))

# Print the counts
print(min_AIC_count)

sink()

###-----------------------------------------------------------------###
###  4.          AUTOCORRELATION: DIAGNOSTICS OUTPUT.               ###
###-----------------------------------------------------------------###
#calculate how many cases are autocorrelated and how many cases have significant raos score to determine which model to use
determine_Moran_dependent <- vector(length = length(output))
determine_Moran_resid <- vector(length = length(output))
determine_Rao_error <- vector(length = length(output))
determine_Rao_lag <- vector(length = length(output))

for(i in seq_along(output)) {
  
  if(!is.null(output[[i]]$MoransI)) {
    print(output[[i]]$MoransI)
    print(output[[i]]$MoransI_resid)
    print(output[[i]]$Rao_error)
    print(output[[i]]$Rao_lag)
    
    determine_Moran_dependent[i] <- ifelse(output[[i]]$MoransI$p.value < 0.05, 1, 0)
    determine_Moran_resid[i] <- ifelse(output[[i]]$MoransI_resid$p.value < 0.05, 1, 0)
    determine_Rao_error[i] <- ifelse(output[[i]]$Rao_error$RSerr$p.value < 0.05, 1, 0)
    determine_Rao_lag[i] <- ifelse(output[[i]]$Rao_lag$RSlag$p.value < 0.05, 1, 0)
  }
}

sum(determine_Moran_dependent) #59 dependent variables show autocorrelation
sum(determine_Moran_resid) #50 models show autocorrelation in residuals
sum(determine_Rao_error) #39 models show significance (spatial error model is preferred)
sum(determine_Rao_lag) #44 models show significance (spatial lag model is preferred)

###Export autocorrelation diagnostics summary
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_France", FGs, ".txt", sep = ""))
print(paste("Morans I Monte Carlo Simulation:", sum(determine_Moran_dependent), "of", length(output), "P-Values are significant"))
print(paste("Morans I for Spatial Autocorrelation in Residuals:", sum(determine_Moran_resid), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Error Model):", sum(determine_Rao_error), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Lag Model):", sum(determine_Rao_lag), "of", length(output), "P-Values are significant"))
sink()

###-----------------------------------------------------------------###
###  5.                    PLOTTING SLOPES                          ###
###-----------------------------------------------------------------###

#create dataframe for plotting spatial lag results
final_data_sp_sac_france <- data.frame("FunctionalGroupIndex" = character(), "Climate" = character(), "Scenario" = character(), "Region" = character(),"Intercept" = numeric(), "Slope" = numeric(), "Std.Err" = numeric(), "Z.Stat" = numeric(), "P.Val" = numeric(), stringsAsFactors = FALSE)

#extract output for plotting
for(i in seq_along(output)) { 
  
  if(!is.null(output[[i]]$MoransI)) {
    
    if(impact == "Impacts") {
      sums <- summary(output[[i]]$Spa_sac_effects, short = T, zstats = T)
      
      final_data <- data.frame("FunctionalGroupIndex" = c(output[[i]]$group[1]),
                               "Climate" = c(output[[i]]$group[2]),
                               "Scenario" = c(output[[i]]$group[3]),
                               "Region" = c(output[[i]]$group[4]),
                               "Intercept" = c(round(sums$total_sum$statistics[1])),
                               "Slope" = c(round(sums$res$total[[1]], digits = 2)),
                               "Std.Err" = c(round(sums$semat[[3]], digits = 2)),
                               "Z.Stat" = round(sums$zmat[[3]], digits = 2),
                               "P.Val" = format.pval(sums$pzmat[[3]], digits = 2))
      
      final_data_sp_sac_france <- rbind(final_data_sp_sac_france, final_data)
    } else {
      sums <- summary(output[[i]]$Spa_sac)
      summary(output[[1]]$Spa_sac)
      final_data <- data.frame(
        "FunctionalGroupIndex" = output[[i]]$group[1],
        "Climate" = output[[i]]$group[2],
        "Scenario" = output[[i]]$group[3],
        "Region" = output[[i]]$group[4],
        "Intercept" = round(sums$Coef[1], digits = 2),
        "Slope" = round(sums$Coef[2], digits = 2),
        "Std.Err" = round(sums$Coef[4], digits = 2),
        "Z.Stat" = round(sums$Coef[6], digits = 2),
        "P.Val" = format.pval(sums$Coef[8], digits = 2)
      )
      final_data_sp_sac_france <- rbind(final_data_sp_sac_france, final_data)
    }
  }}

###---------###
### Namibia ###
###---------###
#load workspace images from input 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 2 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(input,sep="/","Namibia",sep="/","Output_Summary.Rdata"))

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
if(FGs == "all") {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores (i.)", '1' = "Endotherm Carnivores (i.)", '2' =  "Endotherm Omnivores (i.)",
                                      '3'="Ectotherm Herbivores (s.)",'4'="Ectotherm Carnivores (s.)",'5'="Ectotherm Omnivores (s.)","6" = "Ectotherm Herbivores (i.)",
                                      '7' = "Ectotherm Carnivores (i.)",'8'="Ectotherm Omnivores (i.)" )
  
  
} else {
  test$FunctionalGroupIndex <- recode(test$FunctionalGroupIndex, '0' = "Endotherm Herbivores", '1' = "Endotherm Carnivores", '2' =  "Endotherm Omnivores",
                                      '3'="Ectotherm Herbivores",'4'="Ectotherm Carnivores",'5'="Ectotherm Omnivores","6" = "Ectotherm Herbivores",
                                      '7' = "Ectotherm Carnivores",'8'="Ectotherm Omnivores" )
  test <- test %>%
    group_by(GridcellIndex, FunctionalGroupIndex, group, group2) %>%
    summarise(
      CohortAbundance = mean(CohortAbundance, na.rm = TRUE),
      IndividualBodyMass = mean(IndividualBodyMass, na.rm = TRUE)
    ) %>% ungroup()
}

###create grouped df with underlying data & models information we need to perform loops later on for identifying neighbours!
lm_fit_namibia <- test %>%
  group_by(FunctionalGroupIndex, group, group2) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .x))) %>%
  ungroup() %>%
  mutate(Region = "Namibia")

#maintain df order for plotting later
#lm_fit_namibia <- lm_fit_namibia %>%
#  group_by(FunctionalGroupIndex, group, group2)

###-----------------------------------------------------------------###
###  2.                      NEIGHBOUR LISTS                        ###
###-----------------------------------------------------------------###
###Step 1: Identify Neighbours 
#create neigbours object for gridcells based on spatial grid dimensions
#spatial_window <- c(16,21,-22,-17) #90 grid cells, Namibia 
min.x <- 16
max.x <- 21
min.y <- -22
max.y <- -17

#define cell size (0.5 degree in our case)
cellsize <- 0.5

#calculate x y dimensions of grid
cells.dim.x <- (max.x - min.x) / cellsize
cells.dim.y <- (max.y - min.y) / cellsize

#calculate cells dimensions (combine x y information)
cells.dim <- c(cells.dim.x, cells.dim.y)

#create neigbours object
neighbours <- spdep::grid2nb(d = cells.dim)

###-----------------------------------------------------------------###
###  3.          AUTOCORRELATION: STATISTICS AND MODELS             ###
###-----------------------------------------------------------------###
###from here i loop over datasets
output <- vector("list", length = length(lm_fit_namibia$data))
output <- lapply(output, function(x) list(group = NULL, lm = NULL, MoransI = NULL, MoransI_resid = NULL, Rao_lag = NULL, Rao_error = NULL, Spa_lag = NULL, Spa_error = NULL))

for(i in seq_along(lm_fit_namibia$data)) {
  
  dat <- lm_fit_namibia$data[[i]]
  #dat$IndividualBodyMass <- log10(dat$IndividualBodyMass)
  #dat$CohortAbundance <- log10(dat$CohortAbundance)
  
  #create group output
  output[[i]]$group <- c(lm_fit_namibia$FunctionalGroupIndex[[i]], lm_fit_namibia$group[[i]], lm_fit_namibia$group2[[i]], lm_fit_namibia$Region[[i]])
  
  #ifelse condition (some combinations (e.g. ectotherm carnivores namibia max land use), have too less observations so that they dont have neigbours, thus I exclude all combinations with < 5 cells/observations from modelling)
  if(length(dat$GridcellIndex) < 5) {
    
    cat(lm_fit_namibia$FunctionalGroupIndex[[i]], "has to less observations, not used for modelling.")
    
  } else{
    
    ###1. Subset for each group in i neigbours list (because sometimes we have missing cells in between)
    #create logical cells vector (dimensions like neigbours list, but with FALSE as default)
    observed_cells <- rep(FALSE, length(neighbours))
    
    #add up +1 to gridcellindex, because it starts at 0 (is not matching neighbours list dimensions because R starts at 1), set all cells that can be found in both to TRUE
    observed_cells[dat$GridcellIndex+1] <- TRUE
    
    #calculate no. of missing cells per i 
    missing_cells <- which(!observed_cells)
    cat("Missing cells:", sum(n_distinct(missing_cells)), "\n")
    
    #subset neighbours list to include only observed cells
    nb <- spdep::subset.nb(neighbours, subset = observed_cells, zero.policy=TRUE)
    
    
    #remove cells without neighbours
    #nb <- nb[!sapply(nb, function(x) any(x == 0))]    
    
    #class(nb) <- "nb"
    
    #create spatial weights object 
    listw <- spdep::nb2listw(nb, style = "W", zero.policy=TRUE)
    
    cat("Cell correction done, neigbours created. Calculation of Morans I and spatial regression.", "Step:", i, length(lm_fit_namibia$data))
    
    fit.lm <- lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = dat)
    output[[i]]$lm <- fit.lm
    
    ###2. calculate global morans I (for dependent variable, and residuals separately) and raos score for spatial error and spatial lag model
    #test for autocorrelation of dependent variable if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI <- spdep::moran.mc(dat$CohortAbundance, listw, zero.policy = TRUE, nsim = 999)
    
    #test for autocorrelation of residuals if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI_resid <- (spdep::lm.morantest(fit.lm, listw, zero.policy = T)) #test = significant --> spatial autocorrelation detected
    
    #test if spatial error or spatial lag model more appropriate with Raos Score (Lagrange multiplier test)
    output[[i]]$Rao_error <- spdep::lm.RStests(fit.lm, listw, test="RSerr", zero.policy = T) # test not significant --> spatial lag model more appropriate
    
    output[[i]]$Rao_lag <- spdep::lm.RStests(fit.lm, listw, test="RSlag", zero.policy = T) # test not significant --> spatial error model more appropriate
    
    ###3. Perform spatial regression (spatial lag model)
    output[[i]]$Spa_error <- spatialreg::errorsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                                    data = dat, 
                                                    listw = listw, zero.policy = TRUE) 
    
    Spa_lag <- spatialreg::lagsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    Spa_sac <- spatialreg::sacsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    output[[i]]$Spa_lag <- Spa_lag
    output[[i]]$Spa_sac <- Spa_sac
    
    #calculate effects/impacts (indirect, direct, total) for spatial lag model 
    output[[i]]$Spa_lag_effects <- impacts(Spa_lag, listw = listw, R = 999)
    output[[i]]$Spa_sac_effects <- impacts(Spa_sac, listw = listw, R = 999)
    
  }}


AIC_list <- vector("list", length(output))

for(i in seq_along(output)) {
  
  AIC_values <- c(lm = AIC(output[[i]]$lm),
                  Spa_error = AIC(output[[i]]$Spa_error),
                  Spa_lag = AIC(output[[i]]$Spa_lag),
                  Spa_sac = AIC(output[[i]]$Spa_sac))
  AIC_list[[i]] <- AIC_values
}

AIC_df <- do.call(rbind, AIC_list)

AIC_df <- data.frame(Index = seq_along(output), AIC_df)

#create vector to calculate which model has lowest AIC
min_AIC_column <- apply(AIC_df[, -1], 1, function(x) names(x)[which.min(x)])

#store minimum aic values
AIC_df$min_AIC <- min_AIC_column

#count how often for each model AIC was minimum 
min_AIC_count <- table(AIC_df$min_AIC)

#store results in txt.file 
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_AICs_Namibia", FGs, ".txt", sep = ""))

# Print the counts
print(min_AIC_count)

sink()

###-----------------------------------------------------------------###
###  4.          AUTOCORRELATION: DIAGNOSTICS OUTPUT.               ###
###-----------------------------------------------------------------###
#calculate how many cases are autocorrelated and how many cases have significant raos score to determine which model to use
determine_Moran_dependent <- vector(length = length(output))
determine_Moran_resid <- vector(length = length(output))
determine_Rao_error <- vector(length = length(output))
determine_Rao_lag <- vector(length = length(output))

for(i in seq_along(output)) {
  
  if(!is.null(output[[i]]$MoransI)) {
    print(output[[i]]$MoransI)
    print(output[[i]]$MoransI_resid)
    print(output[[i]]$Rao_error)
    print(output[[i]]$Rao_lag)
    
    determine_Moran_dependent[i] <- ifelse(output[[i]]$MoransI$p.value < 0.05, 1, 0)
    determine_Moran_resid[i] <- ifelse(output[[i]]$MoransI_resid$p.value < 0.05, 1, 0)
    determine_Rao_error[i] <- ifelse(output[[i]]$Rao_error$RSerr$p.value < 0.05, 1, 0)
    determine_Rao_lag[i] <- ifelse(output[[i]]$Rao_lag$RSlag$p.value < 0.05, 1, 0)
  }
}

sum(determine_Moran_dependent) #59 dependent variables show autocorrelation
sum(determine_Moran_resid) #50 models show autocorrelation in residuals
sum(determine_Rao_error) #39 models show significance (spatial error model is preferred)
sum(determine_Rao_lag) #44 models show significance (spatial lag model is preferred)

###Export autocorrelation diagnostics summary
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_Namibia", FGs, ".txt", sep = ""))
print(paste("Morans I Monte Carlo Simulation:", sum(determine_Moran_dependent), "of", length(output), "P-Values are significant"))
print(paste("Morans I for Spatial Autocorrelation in Residuals:", sum(determine_Moran_resid), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Error Model):", sum(determine_Rao_error), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Lag Model):", sum(determine_Rao_lag), "of", length(output), "P-Values are significant"))
sink()

###-----------------------------------------------------------------###
###  5.                    PLOTTING SLOPES                          ###
###-----------------------------------------------------------------###

#create dataframe for plotting spatial lag results
final_data_sp_sac_namibia <- data.frame("FunctionalGroupIndex" = character(), "Climate" = character(), "Scenario" = character(), "Region" = character(),"Intercept" = numeric(), "Slope" = numeric(), "Std.Err" = numeric(), "Z.Stat" = numeric(), "P.Val" = numeric(), stringsAsFactors = FALSE)

#extract output for plotting
for(i in seq_along(output)) { 
  
  if(!is.null(output[[i]]$MoransI)) {
    
    if(impact == "Impacts") {
      sums <- summary(output[[i]]$Spa_sac_effects, short = T, zstats = T)
      
      final_data <- data.frame("FunctionalGroupIndex" = c(output[[i]]$group[1]),
                               "Climate" = c(output[[i]]$group[2]),
                               "Scenario" = c(output[[i]]$group[3]),
                               "Region" = c(output[[i]]$group[4]),
                               "Intercept" = c(round(sums$total_sum$statistics[1])),
                               "Slope" = c(round(sums$res$total[[1]], digits = 2)),
                               "Std.Err" = c(round(sums$semat[[3]], digits = 2)),
                               "Z.Stat" = round(sums$zmat[[3]], digits = 2),
                               "P.Val" = format.pval(sums$pzmat[[3]], digits = 2))
      
      final_data_sp_sac_namibia <- rbind(final_data_sp_sac_namibia, final_data)
    } else {
      sums <- summary(output[[i]]$Spa_sac)
      summary(output[[1]]$Spa_sac)
      final_data <- data.frame(
        "FunctionalGroupIndex" = output[[i]]$group[1],
        "Climate" = output[[i]]$group[2],
        "Scenario" = output[[i]]$group[3],
        "Region" = output[[i]]$group[4],
        "Intercept" = round(sums$Coef[1], digits = 2),
        "Slope" = round(sums$Coef[2], digits = 2),
        "Std.Err" = round(sums$Coef[4], digits = 2),
        "Z.Stat" = round(sums$Coef[6], digits = 2),
        "P.Val" = format.pval(sums$Coef[8], digits = 2)
      )
      final_data_sp_sac_namibia <- rbind(final_data_sp_sac_namibia, final_data)
    }
  }}

###---------###
### Brazil  ###
###---------###
#load workspace images from input 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(input,sep="/","Brazil",sep="/","Output_Summary.Rdata"))

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
  test <- test %>%
    group_by(GridcellIndex, FunctionalGroupIndex, group, group2) %>%
    summarise(
      CohortAbundance = mean(CohortAbundance, na.rm = TRUE),
      IndividualBodyMass = mean(IndividualBodyMass, na.rm = TRUE)
    ) %>% ungroup()
}

###create grouped df with underlying data & models information we need to perform loops later on for identifying neighbours!
lm_fit_brazil <- test %>%
  group_by(FunctionalGroupIndex, group, group2) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = .x))) %>%
  ungroup() %>%
  mutate(Region = "Brazil")

#maintain df order for plotting later
#lm_fit_brazil <- lm_fit_brazil %>%
#  group_by(FunctionalGroupIndex, group, group2)

###-----------------------------------------------------------------###
###  2.                      NEIGHBOUR LISTS                        ###
###-----------------------------------------------------------------###
###Step 1: Identify Neighbours 
#create neigbours object for gridcells based on spatial grid dimensions
#spatial_window <- c(-69,-61,-3,0) #96 grid cells, Brazil 
min.x <- -69
max.x <- -61
min.y <- -3
max.y <- 0

#define cell size (0.5 degree in our case)
cellsize <- 0.5

#calculate x y dimensions of grid
cells.dim.x <- (max.x - min.x) / cellsize
cells.dim.y <- (max.y - min.y) / cellsize

#calculate cells dimensions (combine x y information)
cells.dim <- c(cells.dim.x, cells.dim.y)

#create neigbours object
neighbours <- spdep::grid2nb(d = cells.dim)

###-----------------------------------------------------------------###
###  3.          AUTOCORRELATION: STATISTICS AND MODELS             ###
###-----------------------------------------------------------------###
###from here i loop over datasets
output <- vector("list", length = length(lm_fit_brazil$data))
output <- lapply(output, function(x) list(group = NULL, lm = NULL, MoransI = NULL, MoransI_resid = NULL, Rao_lag = NULL, Rao_error = NULL, Spa_lag = NULL, Spa_error = NULL))

for(i in seq_along(lm_fit_brazil$data)) {
  
  dat <- lm_fit_brazil$data[[i]]
  #dat$IndividualBodyMass <- log10(dat$IndividualBodyMass)
  #dat$CohortAbundance <- log10(dat$CohortAbundance)
  
  #create group output
  output[[i]]$group <- c(lm_fit_brazil$FunctionalGroupIndex[[i]], lm_fit_brazil$group[[i]], lm_fit_brazil$group2[[i]], lm_fit_brazil$Region[[i]])
  
  #ifelse condition (some combinations (e.g. ectotherm carnivores namibia max land use), have too less observations so that they dont have neigbours, thus I exclude all combinations with < 5 cells/observations from modelling)
  if(length(dat$GridcellIndex) < 5) {
    
    cat(lm_fit_brazil$FunctionalGroupIndex[[i]], "has to less observations, not used for modelling.")
    
  } else{
    
    ###1. Subset for each group in i neigbours list (because sometimes we have missing cells in between)
    #create logical cells vector (dimensions like neigbours list, but with FALSE as default)
    observed_cells <- rep(FALSE, length(neighbours))
    
    #add up +1 to gridcellindex, because it starts at 0 (is not matching neighbours list dimensions because R starts at 1), set all cells that can be found in both to TRUE
    observed_cells[dat$GridcellIndex+1] <- TRUE
    
    #calculate no. of missing cells per i 
    missing_cells <- which(!observed_cells)
    cat("Missing cells:", sum(n_distinct(missing_cells)), "\n")
    
    #subset neighbours list to include only observed cells
    nb <- spdep::subset.nb(neighbours, subset = observed_cells, zero.policy = T)
    
    #create spatial weights object 
    listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    
    cat("Cell correction done, neigbours created. Calculation of Morans I and spatial regression.", "Step:", i, length(lm_fit_namibia$data))
    
    fit.lm <- lm(log10(CohortAbundance) ~ log10(IndividualBodyMass), data = dat)
    output[[i]]$lm <- fit.lm
    
    ###2. calculate global morans I (for dependent variable, and residuals separately) and raos score for spatial error and spatial lag model
    #test for autocorrelation of dependent variable if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI <- spdep::moran.mc(dat$CohortAbundance, listw, nsim = 999)
    
    #test for autocorrelation of residuals if p < 0.05 then autocorrelation detected
    output[[i]]$MoransI_resid <- (spdep::lm.morantest(fit.lm, listw, zero.policy = T)) #test = significant --> spatial autocorrelation detected
    
    #test if spatial error or spatial lag model more appropriate with Raos Score (Lagrange multiplier test)
    output[[i]]$Rao_error <- spdep::lm.RStests(fit.lm, listw, test="RSerr", zero.policy = T) # test not significant --> spatial lag model more appropriate
    
    output[[i]]$Rao_lag <- spdep::lm.RStests(fit.lm, listw, test="RSlag", zero.policy = T) # test not significant --> spatial error model more appropriate
    
    ###3. Perform spatial regression (spatial lag model)
    output[[i]]$Spa_error <- spatialreg::errorsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                                    data = dat, 
                                                    listw = listw, zero.policy = TRUE) 
    
    Spa_lag <- spatialreg::lagsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    Spa_sac <- spatialreg::sacsarlm(log10(CohortAbundance) ~ log10(IndividualBodyMass),  
                                    data = dat, 
                                    listw = listw, zero.policy = TRUE)
    
    output[[i]]$Spa_lag <- Spa_lag
    output[[i]]$Spa_sac <- Spa_sac
    
    #calculate effects/impacts (indirect, direct, total) for spatial lag model 
    output[[i]]$Spa_lag_effects <- impacts(Spa_lag, listw = listw, R = 999)
    output[[i]]$Spa_sac_effects <- impacts(Spa_sac, listw = listw, R = 999)
    
  }}


AIC_list <- vector("list", length(output))

for(i in seq_along(output)) {
  
  AIC_values <- c(lm = AIC(output[[i]]$lm),
                  Spa_error = AIC(output[[i]]$Spa_error),
                  Spa_lag = AIC(output[[i]]$Spa_lag),
                  Spa_sac = AIC(output[[i]]$Spa_sac))
  AIC_list[[i]] <- AIC_values
}

AIC_df <- do.call(rbind, AIC_list)

AIC_df <- data.frame(Index = seq_along(output), AIC_df)

#create vector to calculate which model has lowest AIC
min_AIC_column <- apply(AIC_df[, -1], 1, function(x) names(x)[which.min(x)])

#store minimum aic values
AIC_df$min_AIC <- min_AIC_column

#count how often for each model AIC was minimum 
min_AIC_count <- table(AIC_df$min_AIC)

#store results in txt.file 
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_AICs_Brazil", FGs, ".txt", sep = ""))

# Print the counts
print(min_AIC_count)

sink()

###-----------------------------------------------------------------###
###  4.          AUTOCORRELATION: DIAGNOSTICS OUTPUT.               ###
###-----------------------------------------------------------------###
#calculate how many cases are autocorrelated and how many cases have significant raos score to determine which model to use
determine_Moran_dependent <- vector(length = length(output))
determine_Moran_resid <- vector(length = length(output))
determine_Rao_error <- vector(length = length(output))
determine_Rao_lag <- vector(length = length(output))

for(i in seq_along(output)) {
  
  if(!is.null(output[[i]]$MoransI)) {
    print(output[[i]]$MoransI)
    print(output[[i]]$MoransI_resid)
    print(output[[i]]$Rao_error)
    print(output[[i]]$Rao_lag)
    
    determine_Moran_dependent[i] <- ifelse(output[[i]]$MoransI$p.value < 0.05, 1, 0)
    determine_Moran_resid[i] <- ifelse(output[[i]]$MoransI_resid$p.value < 0.05, 1, 0)
    determine_Rao_error[i] <- ifelse(output[[i]]$Rao_error$RSerr$p.value < 0.05, 1, 0)
    determine_Rao_lag[i] <- ifelse(output[[i]]$Rao_lag$RSlag$p.value < 0.05, 1, 0)
  }
}

sum(determine_Moran_dependent) #59 dependent variables show autocorrelation
sum(determine_Moran_resid) #50 models show autocorrelation in residuals
sum(determine_Rao_error) #39 models show significance (spatial error model is preferred)
sum(determine_Rao_lag) #44 models show significance (spatial lag model is preferred)

###Export autocorrelation diagnostics summary
sink(paste(outpath, "/Autocorrelation_Diagnostics/Summary_Brazil", FGs, ".txt", sep = ""))
print(paste("Morans I Monte Carlo Simulation:", sum(determine_Moran_dependent), "of", length(output), "P-Values are significant"))
print(paste("Morans I for Spatial Autocorrelation in Residuals:", sum(determine_Moran_resid), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Error Model):", sum(determine_Rao_error), "of", length(output), "P-Values are significant"))
print(paste("Raos Score for Spatial Dependence (Spatial Lag Model):", sum(determine_Rao_lag), "of", length(output), "P-Values are significant"))
sink()

###-----------------------------------------------------------------###
###  5.                    PLOTTING SLOPES                          ###
###-----------------------------------------------------------------###

#create dataframe for plotting spatial lag results
final_data_sp_sac_brazil <- data.frame("FunctionalGroupIndex" = character(), "Climate" = character(), "Scenario" = character(), "Region" = character(),"Intercept" = numeric(), "Slope" = numeric(), "Std.Err" = numeric(), "Z.Stat" = numeric(), "P.Val" = numeric(), stringsAsFactors = FALSE)

#extract output for plotting
for(i in seq_along(output)) { 
  
  if(!is.null(output[[i]]$MoransI)) {
    
    if(impact == "Impacts") {
      sums <- summary(output[[i]]$Spa_sac_effects, short = T, zstats = T)
      
      final_data <- data.frame("FunctionalGroupIndex" = c(output[[i]]$group[1]),
                               "Climate" = c(output[[i]]$group[2]),
                               "Scenario" = c(output[[i]]$group[3]),
                               "Region" = c(output[[i]]$group[4]),
                               "Intercept" = c(round(sums$total_sum$statistics[1])),
                               "Slope" = c(round(sums$res$total[[1]], digits = 2)),
                               "Std.Err" = c(round(sums$semat[[3]], digits = 2)),
                               "Z.Stat" = round(sums$zmat[[3]], digits = 2),
                               "P.Val" = format.pval(sums$pzmat[[3]], digits = 2))
      
      final_data_sp_sac_brazil <- rbind(final_data_sp_sac_brazil, final_data)
    } else {
      sums <- summary(output[[i]]$Spa_sac)
      summary(output[[1]]$Spa_sac)
      final_data <- data.frame(
        "FunctionalGroupIndex" = output[[i]]$group[1],
        "Climate" = output[[i]]$group[2],
        "Scenario" = output[[i]]$group[3],
        "Region" = output[[i]]$group[4],
        "Intercept" = round(sums$Coef[1], digits = 2),
        "Slope" = round(sums$Coef[2], digits = 2),
        "Std.Err" = round(sums$Coef[4], digits = 2),
        "Z.Stat" = round(sums$Coef[6], digits = 2),
        "P.Val" = format.pval(sums$Coef[8], digits = 2)
      )
      final_data_sp_sac_brazil <- rbind(final_data_sp_sac_brazil, final_data)
    }
  }}

#####Heatmapping
#First reorder df for consistent data structure
final_data_sp_sac_finl <- final_data_sp_sac_finl %>%
  arrange(FunctionalGroupIndex, Climate, Scenario)

final_data_sp_sac_france <- final_data_sp_sac_france %>%
  arrange(FunctionalGroupIndex, Climate, Scenario)

final_data_sp_sac_namibia <- final_data_sp_sac_namibia %>%
  arrange(FunctionalGroupIndex, Climate, Scenario)

final_data_sp_sac_brazil <- final_data_sp_sac_brazil %>%
  arrange(FunctionalGroupIndex, Climate, Scenario)

models_all <- rbind(final_data_sp_sac_brazil, final_data_sp_sac_finl, final_data_sp_sac_namibia, final_data_sp_sac_france)

#reorder and rename columns
models_all <- models_all[,c(1,4,2,3,5,6,7,8,9)]

colnames(models_all) <- c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "Z. Stat", "P Value")

#remove < signs to convert p-values to numeric
models_all$`P Value` <- gsub("[^0-9.eE-]", "", models_all$`P Value`)
models_all$`P Value` <- as.numeric(models_all$`P Value`)

#create new significance labelling for more intuitive data plotting
models_all$`Significance Lvl` <- ifelse(models_all$`P Value` < 0.001, "***", #highly significant
                                    ifelse(models_all$`P Value` < 0.01, "**", #very significant
                                           ifelse(models_all$`P Value` < 0.05, "*", ""))) #significant

#html table output
models_all %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "Z. Stat",  "P Value", "Significance Lvl"))

climate <- subset(models_all, Scenario == "Climate")
climate %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "Z. Stat", "P Value", "Significance Lvl"))

curr_lu <- subset(models_all, Scenario == "Current Land Use")
curr_lu %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "Z. Stat", "P Value", "Significance Lvl"))

max_lu <- subset(models_all, Scenario == "Maximum Land Use")
max_lu %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable(header = c("FunctionalGroupIndex", "Region", "Climate", "Scenario", "Intercept", "Slope", "Std. Error", "Z. Stat", "P Value", "Significance Lvl"))

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

###ggplot heatmap: 

test3 <- reshape::melt(heatmap_matrix)
test3$Region <- annot_row$Region[match(test3$X1, rownames(annot_row))]
test3$"Simulation Experiment" <- annot_col$Scenario[match(test3$X2, rownames(annot_col))]

#create a key column in old dataframe and new dataframe to match p-values for each specific combination within dataframe for plotting
models_all$key <- paste(models_all$FunctionalGroupIndex, models_all$Region, models_all$Scenario, models_all$Climate, models_all$Slope)
test3$key <- paste(test3$X1, test3$Region, test3$`Simulation Experiment`, test3$X2, test3$value)

if(FGs == "aggregated") {
  test3$key <- gsub(" Carnivores( a| b| c| d)", " Carnivores", test3$key)
  test3$key <- gsub(" Herbivores( a| b| c| d)", " Herbivores", test3$key)
  test3$key <- gsub(" Omnivores( a| b| c| d)", " Omnivores", test3$key)
  test3$key <- gsub(" Historical (1|2|3)", " Historical", test3$key)
  test3$key <- gsub(" SSP1-2.6 (1|2|3)", " SSP1-2.6", test3$key)
  test3$key <- gsub(" SSP5-8.5 (1|2|3)", " SSP5-8.5", test3$key)
} else {
  #remove specific patterns from key column to match old and new key column 
  test3$key <- gsub(" \\(i\\.\\)( a| b| c| d)", " (i.)", test3$key)
  test3$key <- gsub(" \\(s\\.\\)( a| b| c| d)", " (s.)", test3$key)
  test3$key <- gsub(" Historical (1|2|3)", " Historical", test3$key)
  test3$key <- gsub(" SSP1-2.6 (1|2|3)", " SSP1-2.6", test3$key)
  test3$key <- gsub(" SSP5-8.5 (1|2|3)", " SSP5-8.5", test3$key)
}

#match pvalue with key columns to add up pvalues in dataframe for plotting 
test3$"P.Value" <- models_all$`P Value`[match(test3$key, models_all$key)]

#remove < signs to convert p-values to numeric
test3$P.Value <- gsub("[^0-9.eE-]", "", test3$P.Value)
test3$P.Value <- as.numeric(test3$P.Value)

#create new significance labelling for more intuitive data plotting
test3$significance_levels <- ifelse(test3$P.Value < 0.001, "***", #highly significant
                             ifelse(test3$P.Value < 0.01, "**", #very significant
                                    ifelse(test3$P.Value < 0.05, "*", ""))) #significant

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

if(FGs == "aggregated") {
  pdf(paste(getwd(),"Output","Figure4.pdf", sep = "/"),width = 4.33071, height = 4.33071, paper="special") #7.87402 inches = 20 cm
} else {
  pdf(paste(getwd(),"Output","FigureS3.1.pdf", sep = "/"),width = 4.33071, height = 4.33071, paper="special") ##7.87402 inches = 20 cm
}

col <- rev(RColorBrewer::brewer.pal(9, name = "Blues"))

#ggplot heatmap with geom_tile & facet_grid
ggplot(test3, aes(x = X2, y = X1, fill = value)) +
  geom_tile() +  
  geom_text(aes(label = paste0(value, "", significance_levels)), color = ifelse(test3$value > -4, "black", "white"), size = 5/.pt) + # fontface = ifelse(!is.na(test3$P.Value) & test3$P.Value < 0.05, "bold", "plain")) +
  #coord_fixed()
  #scale_fill_viridis(option = "inferno", breaks = seq(-6, 1, 1)) +
  scale_fill_gradientn(colours = col) + 
  ggh4x::facet_grid2(rows = vars(Region),cols = vars(`Simulation Experiment`), scales = "free", space = "free", switch = "x", strip = strip_variation) +
  labs(fill = "Slope Value", y = "Functional Group", x = "Simulation Experiment") +
  scale_y_discrete(labels = rev(row_labels), limits = rev, expand = c(0,0)) +
  scale_x_discrete(labels = col_labels, expand = c(0,0)) +
  theme(axis.text.y = element_text(size=5,hjust=0))+
  #geom_vline(xintercept=c(3.5,6.5), color='black', size = 1) +
  theme_classic()+
  theme(legend.position = "top") +
  theme(panel.spacing=unit(0,units = "cm")) +
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_rect(fill = "white"),
        axis.title = element_blank())  +
  theme(axis.text.y = element_text(size=5,hjust=0))+
  theme(axis.text.x = element_text(size=5))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=8))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=5,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=5,face="bold")) +
  theme(legend.text = element_text(size=5),legend.title = element_text(size=5,face="bold")) +       
  theme(legend.margin=margin(grid::unit(0, "cm")), legend.key.height=grid::unit(0.3, "cm"), legend.key.width=grid::unit(1, "cm")) +
  theme(strip.text = element_text(size = 5))

dev.off()
