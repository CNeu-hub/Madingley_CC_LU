library(MadingleyR)
library(dplyr)
library(ggplot2)

source("settings.R")
load(paste0(outpath,sep="/","France",sep="/","Output_Summary.Rdata"))

#max_year <- max(historical_2014_list[[1]]$time_line_stocks$Year)
#max_month <- max(historical_2014_list[[1]]$time_line_stocks$Month)

autotroph_biomass <- data.frame()  # Initialize an empty data frame

# Concatenate the data frames and keep track of the total years
for (i in 1:length(historical_2014_list)) {
  x <- historical_2014_list[[i]]$time_line_stocks
  
  autotroph_biomass <- rbind(autotroph_biomass, x)
  
}

autotroph_biomass$Month <- seq(1,7792,1)
autotroph_biomass$climate <- "Historical"

#max_year <- max(historical_2014_list[[1]]$time_line_stocks$Year)
#max_month <- max(historical_2014_list[[1]]$time_line_stocks$Month)

FG_biomass <- data.frame()  # Initialize an empty data frame

# Concatenate the data frames and keep track of the total years
for (i in 1:length(historical_2014_list)) {
  x <- historical_2014_list[[i]]$time_line_cohorts
  
  FG_biomass <- rbind(FG_biomass, x)
  
}

FG_biomass$Month <- seq(1,7792,1)
FG_biomass$climate <- "Historical"

#summarize ectotherms biomass: 
#historical
FG_biomass <- FG_biomass %>% rowwise() %>% 
  dplyr::mutate(Biomass_FG_3 = sum(c(Biomass_FG_3,Biomass_FG_6)),Biomass_FG_4 = sum(c(Biomass_FG_4,Biomass_FG_7)),Biomass_FG_5 = sum(c(Biomass_FG_5, Biomass_FG_8))) %>%
  summarize(Month,Year,Biomass_FG_0,Biomass_FG_1,Biomass_FG_2,Biomass_FG_3,Biomass_FG_4,Biomass_FG_5) %>%
  ungroup() %>%   as.data.frame()

#add autotroph biomass to heterotroph biomass for timelines plot
Biomass <- cbind(FG_biomass,autotroph_biomass$TotalStockBiomass)

pdf("/Users/neumanch/Desktop/Thesis/Paper/Figures/BiomassTimelines.pdf", width = 13, height = 5)
ggplot(data=Biomass, aes(x=Month/12))+
  geom_line(aes(y=log10(autotroph_biomass$TotalStockBiomass),colour="Autotrophs"),linewidth=0.8)+
  geom_line(aes(y=log10(Biomass_FG_0),colour="End. Herbivores"),linewidth=0.8)+
  geom_line(aes(y=log10(Biomass_FG_1),colour="End. Carnivores"),linewidth=0.8)+
  geom_line(aes(y=log10(Biomass_FG_2),colour="End. Omnivores"),linewidth=0.8)+
  geom_line(aes(y=log10(Biomass_FG_3),colour="Ect. Herbivores"),linewidth=0.8)+
  geom_line(aes(y=log10(Biomass_FG_4),colour="Ect. Carnivores"),linewidth=0.8)+
  geom_line(aes(y=log10(Biomass_FG_5),colour="Ect. Omnivores"),linewidth=0.8)+
  scale_colour_manual(name = "Functional Group", values = c("Autotrophs"="#009E73",
                                                            "End. Herbivores" = "#CC79A7",
                                                            "End. Carnivores" = "#E69F00",
                                                            "End. Omnivores" = "#56B4E9",
                                                            "Ect. Herbivores" = "#F0E442",
                                                            "Ect. Carnivores" = "#0072B2",
                                                            "Ect. Omnivores" = "#D55E00"))+
  geom_vline(xintercept = c(200,400,410,420,430,440,450),linetype="dashed",alpha =0.2)+
  scale_x_continuous(breaks=c(0,100,200,300,400,500,600,650),limits = c(0,650), expand = c(0,0,0.01,0))+
  labs(x='Year of Simulation', y = "Log10 Biomass [kg]")+
  annotate("text", label=c("1. Climate","2. Current Land Use","Land-use Intensification", "3. Maximum Land Use"),
           x=c(90,290,450,590),y=12,size=20,fontface=c("bold","bold","italic","bold"))+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", hjust =0.5))+
  theme(axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=14)) +
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=14))+
  theme(legend.text = element_text(size=14),legend.title = element_text(size=14,face="bold")) +
  theme(legend.position = "bottom") 

dev.off()
ggsave(figpath %+% "Biomass_Timelines_France.png",dpi=1200,width=360,height=125,units="mm")

citation(package = "MadingleyR")

?coord_cartesian
