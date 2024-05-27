###-----------------------------------------------------------------------------------------###
### TITLE:       Percentage change for current land use (HANPP) scenarios                   ###
### DESCRIPTION: Calculates percentage change for current land use (HANPP) Phase of:        ###
###              "Model-based Impact Analysis of Climate Change and Land-use Intensification###
###              on Trophic Networks."                                                      ###                                                                 
### PROCEDURE:   Calculates percentage change based on bootstrapped logRR & CIs between     ### 
###              the distinct scenarios for each region and exports results as csv datasets ###
### DATE:        13.12.2023                                                                 ###
###-----------------------------------------------------------------------------------------###

###LOAD LIBRARIES
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(boot)#bootstrap CIs

###---------------------------------------------------------------------------###
###                          CREATE SETTINGS / ENVIRONMENT                    ###
###---------------------------------------------------------------------------###
#define input path (Output_Data in repository)
outpath <- paste(getwd(), "Input", sep = "/")

#define output path (Output_CSV in repository)
figpath <- paste(getwd(), "Output", sep = "/")

###------###
###BRAZIL###
###------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","Brazil",sep="/","Output_Summary.Rdata"))

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#1. Create df for control + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_control <- historical_2014_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_control <- SSP126_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_control <- SSP585_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#2. Create df for HANPP scenarios + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_HANPP <- historical_2014_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_HANPP <- SSP126_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_HANPP <- SSP585_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#3. Reorganize data: to create 1 df per climate scenario --> for logRR calculation 
#creates df with control and hanpp application observations (can be distinguished with group variable)
#stick control & treatment dataframe for historical
historical <- rbind(cohorts_historical_control,cohorts_historical_HANPP)
unique(historical$Group)
#calculation of log10 values
historical[3:9] <- log10(historical[3:9]) 

#stick control & treatment dataframe for SSP2.6
SSP2.6 <- rbind(cohorts_SSP2.6_2100_control,cohorts_SSP2.6_2100_HANPP)
unique(SSP2.6$Group)
#calculation of log10 values 
SSP2.6[3:9] <- log10(SSP2.6[3:9])

#stick control & treatment dataframe for SSP8.5
SSP8.5 <- rbind(cohorts_SSP8.5_2100_control,cohorts_SSP8.5_2100_HANPP)
unique(SSP8.5$Group)
#calculation of log10 values 
SSP8.5[3:9] <- log10(SSP8.5[3:9]) 

###-----------------------------------------------------------------###
###  2.       EFFECT SIZE CALCULATION // logRR                      ###
###-----------------------------------------------------------------###

###Set up function which calculates logRR effect sizes & variances between samples for bootstrapping 
logRR <- function(x,y) {
  
  resample <- x[y]
  control <- resample[1:119]
  treatment <- resample[120:238]
  
  mean_c <- mean(control)
  mean_t <- mean(treatment)
  
  var_c <- var(control)
  var_t <- var(treatment)
  
  variance <- as.numeric(c(var_c,var_t))
  mean <- as.numeric(c(mean_c,mean_t))
  n <- as.numeric(c(length(control),length(treatment)))
  
  V_lRR <- sum(variance/(n*mean^2))
  
  V_lRR <- sqrt(V_lRR)
  
  logRR <- log((mean_t))-log((mean_c))
  
  #CI <- logRR + c(-1, 1) * stats::qnorm(1-(1-0.95)/2) * sqrt(V_lRR)
  
  return(c(logRR,V_lRR))
}

###Historical natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_historical_control$Biomass_FG_0)
n2 <- length(cohorts_historical_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(historical$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(historical$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(historical$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(historical$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(historical$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(historical$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(historical$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_historical <- data.frame(c("Brazil"),
                                c("End. Herbivores","End. Carnivores","End. Omnivores",
                                  "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                                c(1,2,3,4,5,6,7),
                                c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                                c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                                c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                                c("Historical"),
                                c(sum(historical$Biomass_FG_0),sum(historical$Biomass_FG_1),sum(historical$Biomass_FG_2),sum(historical$Biomass_FG_3),sum(historical$Biomass_FG_4),sum(historical$Biomass_FG_5),sum(historical$BiomassSum)))

names(effect_historical) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_historical[4:6] <- (exp(effect_historical[4:6])-1)*100

###SSP1-2.6 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP2.6_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP2.6_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP2.6$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP2.6$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP2.6$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP2.6$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP2.6$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP2.6$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP2.6$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP2.6 <- data.frame(c("Brazil"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP1-2.6"),
                            c(sum(SSP2.6$Biomass_FG_0),sum(SSP2.6$Biomass_FG_1),sum(SSP2.6$Biomass_FG_2),sum(SSP2.6$Biomass_FG_3),sum(SSP2.6$Biomass_FG_4),sum(SSP2.6$Biomass_FG_5),sum(SSP2.6$BiomassSum)))

names(effect_SSP2.6) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP2.6[4:6] <- (exp(effect_SSP2.6[4:6])-1)*100

###SSP5-8.5 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP8.5_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP8.5_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP8.5$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP8.5$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP8.5$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP8.5$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP8.5$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP8.5$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP8.5$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP8.5 <- data.frame(c("Brazil"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP5-8.5"),
                            c(sum(SSP8.5$Biomass_FG_0),sum(SSP8.5$Biomass_FG_1),sum(SSP8.5$Biomass_FG_2),sum(SSP8.5$Biomass_FG_3),sum(SSP8.5$Biomass_FG_4),sum(SSP8.5$Biomass_FG_5),sum(SSP8.5$BiomassSum)))

names(effect_SSP8.5) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP8.5[4:6] <- (exp(effect_SSP8.5[4:6])-1)*100 

#create df for later csv output:
HANPP <- rbind(effect_historical,effect_SSP2.6,effect_SSP8.5)

###-------###
###NAMIBIA###
###-------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","Namibia",sep="/","Output_Summary.Rdata"))

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#1. Create df for control + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_control <- historical_2014_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_control <- SSP126_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_control <- SSP585_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#2. Create df for HANPP scenarios + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_HANPP <- historical_2014_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_HANPP <- SSP126_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_HANPP <- SSP585_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#3. Reorganize data: to create 1 df per climate scenario --> for logRR calculation 
#creates df with control and hanpp application observations (can be distinguished with group variable)
#stick control & treatment dataframe for historical
historical <- rbind(cohorts_historical_control,cohorts_historical_HANPP)
unique(historical$Group)
#calculation of log10 values
historical[3:9] <- log10(historical[3:9]) 

#stick control & treatment dataframe for SSP2.6
SSP2.6 <- rbind(cohorts_SSP2.6_2100_control,cohorts_SSP2.6_2100_HANPP)
unique(SSP2.6$Group)
#calculation of log10 values 
SSP2.6[3:9] <- log10(SSP2.6[3:9])

#stick control & treatment dataframe for SSP8.5
SSP8.5 <- rbind(cohorts_SSP8.5_2100_control,cohorts_SSP8.5_2100_HANPP)
unique(SSP8.5$Group)
#calculation of log10 values 
SSP8.5[3:9] <- log10(SSP8.5[3:9]) 

###-----------------------------------------------------------------###
###  2.       EFFECT SIZE CALCULATION // logRR                      ###
###-----------------------------------------------------------------###

###Historical natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_historical_control$Biomass_FG_0)
n2 <- length(cohorts_historical_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(historical$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(historical$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(historical$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(historical$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(historical$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(historical$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(historical$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_historical <- data.frame(c("Namibia"),
                                c("End. Herbivores","End. Carnivores","End. Omnivores",
                                  "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                                c(1,2,3,4,5,6,7),
                                c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                                c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                                c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                                c("Historical"),
                                c(sum(historical$Biomass_FG_0),sum(historical$Biomass_FG_1),sum(historical$Biomass_FG_2),sum(historical$Biomass_FG_3),sum(historical$Biomass_FG_4),sum(historical$Biomass_FG_5),sum(historical$BiomassSum)))

names(effect_historical) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_historical[4:6] <- (exp(effect_historical[4:6])-1)*100

###SSP1-2.6 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP2.6_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP2.6_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP2.6$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP2.6$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP2.6$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP2.6$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP2.6$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP2.6$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP2.6$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP2.6 <- data.frame(c("Namibia"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP1-2.6"),
                            c(sum(SSP2.6$Biomass_FG_0),sum(SSP2.6$Biomass_FG_1),sum(SSP2.6$Biomass_FG_2),sum(SSP2.6$Biomass_FG_3),sum(SSP2.6$Biomass_FG_4),sum(SSP2.6$Biomass_FG_5),sum(SSP2.6$BiomassSum)))

names(effect_SSP2.6) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP2.6[4:6] <- (exp(effect_SSP2.6[4:6])-1)*100

###SSP5-8.5 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP8.5_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP8.5_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP8.5$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP8.5$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP8.5$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP8.5$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP8.5$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP8.5$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP8.5$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP8.5 <- data.frame(c("Namibia"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP5-8.5"),
                            c(sum(SSP8.5$Biomass_FG_0),sum(SSP8.5$Biomass_FG_1),sum(SSP8.5$Biomass_FG_2),sum(SSP8.5$Biomass_FG_3),sum(SSP8.5$Biomass_FG_4),sum(SSP8.5$Biomass_FG_5),sum(SSP8.5$BiomassSum)))

names(effect_SSP8.5) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP8.5[4:6] <- (exp(effect_SSP8.5[4:6])-1)*100 

#create df for later csv output:
HANPP <- rbind(HANPP,effect_historical,effect_SSP2.6,effect_SSP8.5)

###------###
###FRANCE###
###------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","France",sep="/","Output_Summary.Rdata"))

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#1. Create df for control + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_control <- historical_2014_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_control <- SSP126_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_control <- SSP585_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#2. Create df for HANPP scenarios + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_HANPP <- historical_2014_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_HANPP <- SSP126_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_HANPP <- SSP585_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#3. Reorganize data: to create 1 df per climate scenario --> for logRR calculation 
#creates df with control and hanpp application observations (can be distinguished with group variable)
#stick control & treatment dataframe for historical
historical <- rbind(cohorts_historical_control,cohorts_historical_HANPP)
unique(historical$Group)
#calculation of log10 values
historical[3:9] <- log10(historical[3:9]) 

#stick control & treatment dataframe for SSP2.6
SSP2.6 <- rbind(cohorts_SSP2.6_2100_control,cohorts_SSP2.6_2100_HANPP)
unique(SSP2.6$Group)
#calculation of log10 values 
SSP2.6[3:9] <- log10(SSP2.6[3:9])

#stick control & treatment dataframe for SSP8.5
SSP8.5 <- rbind(cohorts_SSP8.5_2100_control,cohorts_SSP8.5_2100_HANPP)
unique(SSP8.5$Group)
#calculation of log10 values 
SSP8.5[3:9] <- log10(SSP8.5[3:9]) 

###-----------------------------------------------------------------###
###  2.       EFFECT SIZE CALCULATION // logRR                      ###
###-----------------------------------------------------------------###

###Historical natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_historical_control$Biomass_FG_0)
n2 <- length(cohorts_historical_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(historical$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(historical$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(historical$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(historical$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(historical$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(historical$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(historical$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_historical <- data.frame(c("France"),
                                c("End. Herbivores","End. Carnivores","End. Omnivores",
                                  "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                                c(1,2,3,4,5,6,7),
                                c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                                c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                                c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                                c("Historical"),
                                c(sum(historical$Biomass_FG_0),sum(historical$Biomass_FG_1),sum(historical$Biomass_FG_2),sum(historical$Biomass_FG_3),sum(historical$Biomass_FG_4),sum(historical$Biomass_FG_5),sum(historical$BiomassSum)))

names(effect_historical) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_historical[4:6] <- (exp(effect_historical[4:6])-1)*100

###SSP1-2.6 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP2.6_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP2.6_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP2.6$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP2.6$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP2.6$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP2.6$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP2.6$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP2.6$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP2.6$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP2.6 <- data.frame(c("France"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP1-2.6"),
                            c(sum(SSP2.6$Biomass_FG_0),sum(SSP2.6$Biomass_FG_1),sum(SSP2.6$Biomass_FG_2),sum(SSP2.6$Biomass_FG_3),sum(SSP2.6$Biomass_FG_4),sum(SSP2.6$Biomass_FG_5),sum(SSP2.6$BiomassSum)))

names(effect_SSP2.6) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP2.6[4:6] <- (exp(effect_SSP2.6[4:6])-1)*100

###SSP5-8.5 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP8.5_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP8.5_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP8.5$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP8.5$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP8.5$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP8.5$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP8.5$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP8.5$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP8.5$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP8.5 <- data.frame(c("France"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP5-8.5"),
                            c(sum(SSP8.5$Biomass_FG_0),sum(SSP8.5$Biomass_FG_1),sum(SSP8.5$Biomass_FG_2),sum(SSP8.5$Biomass_FG_3),sum(SSP8.5$Biomass_FG_4),sum(SSP8.5$Biomass_FG_5),sum(SSP8.5$BiomassSum)))

names(effect_SSP8.5) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP8.5[4:6] <- (exp(effect_SSP8.5[4:6])-1)*100 

#create df for later csv output:
HANPP <- rbind(HANPP,effect_historical,effect_SSP2.6,effect_SSP8.5)

###-------###
###FINLAND###
###-------###
#load workspace images from outpath 
#list file can be read: 
#entry 1 = spinup, entry 2 = apply hanpp, 
#entry 3 - 12 = vegetation reduction stages, 
#last entry = 13 = post reduction run (with 0.1 vegetation left)
load(paste0(outpath,sep="/","Finland",sep="/","Output_Summary.Rdata"))

###-----------------------------------------------------------------###
###  1.                      DATA PREPARATION                       ###
###-----------------------------------------------------------------###
#1. Create df for control + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_control <- historical_2014_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_control <- SSP126_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_control <- SSP585_2100_list[[1]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "control",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#2. Create df for HANPP scenarios + sum up monthly biomass of each scenario (overall effect), filter = > Month 2280, last 10 years
#historical
cohorts_historical_HANPP <- historical_2014_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="Historical",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()  

#SSP2.6_2100
cohorts_SSP2.6_2100_HANPP <- SSP126_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP1-2.6",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()
#SSP8.5_2100
cohorts_SSP8.5_2100_HANPP <- SSP585_2100_list[[2]]$time_line_cohorts %>% filter(Month > 2280) %>% rowwise() %>% 
  dplyr::mutate(BiomassSum = sum(c(Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,Biomass_FG_6,Biomass_FG_7,Biomass_FG_8)),
                Scenario="SSP5-8.5",Group = "HANPP",
                Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5,BiomassSum,Scenario,Group) %>%
  ungroup() %>%   as.data.frame()

#3. Reorganize data: to create 1 df per climate scenario --> for logRR calculation 
#creates df with control and hanpp application observations (can be distinguished with group variable)
#stick control & treatment dataframe for historical
historical <- rbind(cohorts_historical_control,cohorts_historical_HANPP)
unique(historical$Group)
#calculation of log10 values
historical[3:9] <- log10(historical[3:9]) 

#stick control & treatment dataframe for SSP2.6
SSP2.6 <- rbind(cohorts_SSP2.6_2100_control,cohorts_SSP2.6_2100_HANPP)
unique(SSP2.6$Group)
#calculation of log10 values 
SSP2.6[3:9] <- log10(SSP2.6[3:9])

#stick control & treatment dataframe for SSP8.5
SSP8.5 <- rbind(cohorts_SSP8.5_2100_control,cohorts_SSP8.5_2100_HANPP)
unique(SSP8.5$Group)
#calculation of log10 values 
SSP8.5[3:9] <- log10(SSP8.5[3:9]) 

###-----------------------------------------------------------------###
###  2.       EFFECT SIZE CALCULATION // logRR                      ###
###-----------------------------------------------------------------###

###Historical natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_historical_control$Biomass_FG_0)
n2 <- length(cohorts_historical_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(historical$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(historical$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(historical$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(historical$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(historical$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(historical$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(historical$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_historical <- data.frame(c("Finland"),
                                c("End. Herbivores","End. Carnivores","End. Omnivores",
                                  "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                                c(1,2,3,4,5,6,7),
                                c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                                c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                                c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                                c("Historical"),
                                c(sum(historical$Biomass_FG_0),sum(historical$Biomass_FG_1),sum(historical$Biomass_FG_2),sum(historical$Biomass_FG_3),sum(historical$Biomass_FG_4),sum(historical$Biomass_FG_5),sum(historical$BiomassSum)))

names(effect_historical) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_historical[4:6] <- (exp(effect_historical[4:6])-1)*100

###SSP1-2.6 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP2.6_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP2.6_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP2.6$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP2.6$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP2.6$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP2.6$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP2.6$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP2.6$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP2.6$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP2.6 <- data.frame(c("Finland"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP1-2.6"),
                            c(sum(SSP2.6$Biomass_FG_0),sum(SSP2.6$Biomass_FG_1),sum(SSP2.6$Biomass_FG_2),sum(SSP2.6$Biomass_FG_3),sum(SSP2.6$Biomass_FG_4),sum(SSP2.6$Biomass_FG_5),sum(SSP2.6$BiomassSum)))

names(effect_SSP2.6) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP2.6[4:6] <- (exp(effect_SSP2.6[4:6])-1)*100

###SSP5-8.5 natural ecosystem state vs. HANPP application###
#define length for bootstrapping samples (used for strata)
n1 <- length(cohorts_SSP8.5_2100_control$Biomass_FG_0)
n2 <- length(cohorts_SSP8.5_2100_HANPP$Biomass_FG_0)

#calculate ordinary bootstrap replicates with replacement because: 
FG1 <- boot(SSP8.5$Biomass_FG_0,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG2 <- boot(SSP8.5$Biomass_FG_1,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG3 <- boot(SSP8.5$Biomass_FG_2,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG4 <- boot(SSP8.5$Biomass_FG_3,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG5 <- boot(SSP8.5$Biomass_FG_4,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
FG6 <- boot(SSP8.5$Biomass_FG_5,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")
Overall <- boot(SSP8.5$BiomassSum,logRR, R=10000, strata= rep(1:2,c(n1,n2)),sim="ordinary")

#calculate stud bootstrap cis (bias corrected, accelerated confidence interval), because: 
FG1_CI <- boot.ci(FG1,conf=0.95,type = "stud")
FG2_CI <- boot.ci(FG2,conf=0.95,type = "stud")
FG3_CI <- boot.ci(FG3,conf=0.95,type = "stud")
FG4_CI <- boot.ci(FG4,conf=0.95,type = "stud")
FG5_CI <- boot.ci(FG5,conf=0.95,type = "stud")
FG6_CI <- boot.ci(FG6,conf=0.95,type = "stud")
Overall_CI <- boot.ci(Overall,conf=0.95,type = "stud")

#create a effect size data frame = base for forest plot
effect_SSP8.5 <- data.frame(c("Finland"),
                            c("End. Herbivores","End. Carnivores","End. Omnivores",
                              "Ect. Herbivores","Ect. Carnivores","Ect. Omnivores","Overall Biomass"),
                            c(1,2,3,4,5,6,7),
                            c(FG1$t0[1],FG2$t0[1],FG3$t0[1],FG4$t0[1],FG5$t0[1],FG6$t0[1],Overall$t0[1]),
                            c(FG1_CI$stud[4],FG2_CI$stud[4],FG3_CI$stud[4],FG4_CI$stud[4],FG5_CI$stud[4],FG6_CI$stud[4],Overall_CI$stud[4]),
                            c(FG1_CI$stud[5],FG2_CI$stud[5],FG3_CI$stud[5],FG4_CI$stud[5],FG5_CI$stud[5],FG6_CI$stud[5],Overall_CI$stud[5]),
                            c("SSP5-8.5"),
                            c(sum(SSP8.5$Biomass_FG_0),sum(SSP8.5$Biomass_FG_1),sum(SSP8.5$Biomass_FG_2),sum(SSP8.5$Biomass_FG_3),sum(SSP8.5$Biomass_FG_4),sum(SSP8.5$Biomass_FG_5),sum(SSP8.5$BiomassSum)))

names(effect_SSP8.5) = c("Region","Functional Group","Index","Eff_size","lower","upper","Scenario","Biomass")

#exponentiate / convert to percentage change
effect_SSP8.5[4:6] <- (exp(effect_SSP8.5[4:6])-1)*100 

#create df for later csv output:
HANPP <- rbind(HANPP,effect_historical,effect_SSP2.6,effect_SSP8.5)

#subset df to relevant information for csv output:
HANPP <- HANPP %>% 
  summarize(Region, Scenario,`Functional Group`,Eff_size,lower,upper)
#new column names
names(HANPP)=c("Region","Scenario","Functional Group","Effect Size","Lower CI","Upper CI")
#export df as csv 
write.csv(HANPP,paste0(figpath,sep="/","HANPP_Effect.csv"))

#round values
HANPP[4:6] <- round(HANPP[4:6], digits = 2)

#plot html table
HANPP %>%
  htmlTable::addHtmlTableStyle(css.cell = c("width: 150;","width: 150;","width: 150;","width: 150;","width: 150;","width: 150;")) %>%
  htmlTable::htmlTable() 

