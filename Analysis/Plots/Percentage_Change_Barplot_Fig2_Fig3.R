###-----------------------------------------------------------------------------------------###
### TITLE:       Percentage change barplots                                                 ###
### DESCRIPTION: Plots figure 2 and figure 3 of:                                            ###
###              "Model-based Impact Analysis of Climate Change and Land-use Intensification###
###              on Trophic Networks."                                                      ###                                                                 
### PROCEDURE:   Plots results of percentage change scripts for each scenario, functional   ### 
###              group and region as barplot.                                               ###
### DATE:        04.08.2023                                                                 ###
###-----------------------------------------------------------------------------------------###

library(ggplot2)
library(RColorBrewer)
library(ggpubr)

#define figpath for figure output
figpath <- paste(getwd(), "Output/", sep = "/")

#load & modify data #set here the repository output_csv folder
Climate <- read.csv(paste0(figpath,"Climate_Effect.csv"))
Climate$Climate <- ifelse(Climate$Scenario=="SSP1-2.6","SSP1-2.6","SSP5-8.5")
Climate$Scenario <- "Climate"

Current_LU <- read.csv(paste0(figpath, "HANPP_Effect.csv"))
Current_LU$Climate <- paste(Current_LU$Scenario)
Current_LU$Scenario <- "Current Land Use"

Max_LU <- read.csv(paste0(figpath,"maxHANPP_Effect.csv"))
Max_LU$Climate <- paste(Max_LU$Scenario)
Max_LU$Scenario <- "Maximum Land Use"

All <- rbind(Climate,Current_LU,Max_LU)

#Decapitalize groups for clean plotting
All <- All %>%
  mutate(Functional.Group = stringr::str_replace(Functional.Group, "Omn", "omn")) %>%
  mutate(Functional.Group = stringr::str_replace(Functional.Group, "Herb", "herb")) %>%
  mutate(Functional.Group = stringr::str_replace(Functional.Group, "Carn", "carn")) %>%
  mutate(Functional.Group = stringr::str_replace(Functional.Group, "Biomass", "biomass"))

#create labels for facet wraps 
labels_tot_bio <- c(
  `Climate` = "(a) Climate",
  `Current Land Use` = "(b) Current land use",
  `Maximum Land Use` = "(c) Maximum land use"
)

#check color blind palette #we use Rcolorbrewer Set2 here
display.brewer.all(colorblindFriendly = T)

Brazil <- subset(All,Region=="Brazil")
Brazil$Scenario <- factor(Brazil$Scenario,levels=c("Climate","Current Land Use","Maximum Land Use"))
Brazil$Functional.Group <- factor(Brazil$Functional.Group, levels=rev(c("Ect. omnivores","Ect. carnivores","Ect. herbivores","End. omnivores","End. carnivores","End. herbivores","Overall biomass")))

Brazil_plot <- ggplot(data=Brazil)+
  geom_bar(stat="identity",aes(x=Effect.Size,y=Functional.Group,fill=Climate),position="dodge")+
  geom_errorbar(stat="identity",aes(xmax=Upper.CI,xmin=Lower.CI,y=Functional.Group,group=Climate),position=position_dodge(0.9),linewidth=0.2,width=0.2)+
  facet_wrap(vars(Scenario), labeller = as_labeller(labels_tot_bio))+#,scales = "free_x")+
  scale_x_continuous(breaks=seq(-12,1,2))+
  scale_fill_brewer(palette="Set2",direction = 1)+  
  labs(title="Brazil", x='Effect size (% change)', y = "Functional group (FG)",col="Climate scenario") +
  geom_vline(xintercept=0, color='black', linetype='dashed') +
  theme_classic()+
  theme(axis.text.y = element_text(face = c('bold','plain', 'plain', 'plain', 'plain', 'plain', 'plain'),size=6,hjust=0))+
  theme(axis.text.x = element_text(size=6))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=8))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=6,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=6,face="bold")) +
  theme(legend.text = element_text(size=6),legend.title = element_text(size=6,face="bold")) +
  theme(strip.text = element_text(size = 6))

Brazil_plot

Namibia <- subset(All,Region=="Namibia")
Namibia$Scenario <- factor(Namibia$Scenario,levels=c("Climate","Current Land Use","Maximum Land Use"))
Namibia$Functional.Group <- factor(Namibia$Functional.Group, levels=rev(c("Ect. omnivores","Ect. carnivores","Ect. herbivores","End. omnivores","End. carnivores","End. herbivores","Overall biomass")))

Namibia_plot <- ggplot(data=Namibia)+
  geom_bar(stat="identity",aes(x=Effect.Size,y=Functional.Group,fill=Climate),position="dodge")+
  geom_errorbar(stat="identity",aes(xmax=Upper.CI,xmin=Lower.CI,y=Functional.Group,group=Climate),position=position_dodge(0.9),linewidth=0.2,width=0.2)+
  facet_wrap(vars(Scenario), labeller = as_labeller(labels_tot_bio))+#,scales = "free_x")+
  scale_x_continuous(breaks=seq(-25,1,5))+
  scale_fill_brewer(palette="Set2",direction = 1)+  
  labs(title="Namibia", x='Effect size (% change)', y = "Functional group (FG)",col="Climate scenario") +
  geom_vline(xintercept=0, color='black', linetype='dashed') +
  theme_classic()+
  theme(axis.text.y = element_text(face = c('bold','plain', 'plain', 'plain', 'plain', 'plain', 'plain'),size=6,hjust=0))+
  theme(axis.text.x = element_text(size=6))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=8))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=6,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=6,face="bold")) +
  theme(legend.text = element_text(size=6),legend.title = element_text(size=6,face="bold"))+
  theme(strip.text = element_text(size = 6))

Namibia_plot

France <- subset(All,Region=="France")
France$Scenario <- factor(France$Scenario,levels=c("Climate","Current Land Use","Maximum Land Use"))
France$Functional.Group <- factor(France$Functional.Group, levels=rev(c("Ect. omnivores","Ect. carnivores","Ect. herbivores","End. omnivores","End. carnivores","End. herbivores","Overall biomass")))

France_plot <- ggplot(data=France)+
  geom_bar(stat="identity",aes(x=Effect.Size,y=Functional.Group,fill=Climate),position="dodge")+
  geom_errorbar(stat="identity",aes(xmax=Upper.CI,xmin=Lower.CI,y=Functional.Group,group=Climate),position=position_dodge(0.9),linewidth=0.2,width=0.2)+
  facet_wrap(vars(Scenario), labeller = as_labeller(labels_tot_bio))+#,scales = "free_x")+
  scale_x_continuous(breaks=seq(-20,1,5))+
  scale_fill_brewer(palette="Set2",direction = 1)+  
  labs(title="France", x='Effect size (% change)', y = "Functional group (FG)",col="Climate scenario") +
  geom_vline(xintercept=0, color='black', linetype='dashed') +
  theme_classic()+
  theme(axis.text.y = element_text(face = c('bold','plain', 'plain', 'plain', 'plain', 'plain', 'plain'),size=6,hjust=0))+
  theme(axis.text.x = element_text(size=6))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=8))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=6,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=6,face="bold")) +
  theme(legend.text = element_text(size=6),legend.title = element_text(size=6,face="bold"))+
  theme(strip.text = element_text(size = 6))

France_plot

Finland <- subset(All,Region=="Finland")
Finland$Scenario <- factor(Finland$Scenario,levels=c("Climate","Current Land Use","Maximum Land Use"))
Finland$Functional.Group <- factor(Finland$Functional.Group, levels=rev(c("Ect. omnivores","Ect. carnivores","Ect. herbivores","End. omnivores","End. carnivores","End. herbivores","Overall biomass")))

Finland_plot <- ggplot(data=Finland)+
  geom_bar(stat="identity",aes(x=Effect.Size,y=Functional.Group,fill=Climate),position="dodge")+
  geom_errorbar(stat="identity",aes(xmax=Upper.CI,xmin=Lower.CI,y=Functional.Group,group=Climate),position=position_dodge(0.9),linewidth=0.2,width=0.2)+
  facet_wrap(vars(Scenario), labeller = as_labeller(labels_tot_bio))+#,scales = "free_x")+
  scale_x_continuous(breaks=seq(-15,2,5))+
  scale_fill_brewer(palette="Set2",direction = 1)+  
  labs(title="Finland", x='Effect size (% change)', y = "Functional group (FG)",col="Climate scenario") +
  geom_vline(xintercept=0, color='black', linetype='dashed') +
  theme_classic()+
  theme(axis.text.y = element_text(face = c('bold','plain', 'plain', 'plain', 'plain', 'plain', 'plain'),size=6,hjust=0))+
  theme(axis.text.x = element_text(size=6))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=8))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=6,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=6,face="bold")) +
  theme(legend.text = element_text(size=6),legend.title = element_text(size=6,face="bold"))+
  theme(strip.text = element_text(size = 6))

Finland_plot

###-----------------------------------------------------------------###
###         ARRANGE ALL PLOTS TO DINA4 AND SAVE AS PDF              ###
###-----------------------------------------------------------------###
#Makes a big plot of all 9 ggplot objects (FG's) + delete redundant information 
a <- ggarrange(Brazil_plot+rremove("xlab")+rremove("ylab"),Finland_plot+rremove("ylab")+rremove("xlab"),France_plot+rremove("ylab")+rremove("xlab"),Namibia_plot+rremove("ylab")+rremove("xlab"), 
               ncol = 1, nrow = 4,
               labels = c("(a)", "(b)", "(c)","(d)","(e)","(f)"),
               font.label = list(size=8),
               legend = "bottom",
               common.legend = T,
               align = "h")#dev.off()

#add title + x and y axis captions 
annotate_figure(a,
                bottom = text_grob("Effect size (% change)",hjust=0.2,face="bold",size=6),
                left = text_grob("Functional group",rot=90,hjust=0.5,face="bold",size=6))

#save pdf file in dina4
ggsave(paste0(figpath,sep="_","Fig3.pdf"),width=110, height=180, units="mm")

###-----------------------------------------------------------------###
###                         REGIONS PLOT                            ###
###-----------------------------------------------------------------###

Regions <- subset(All,Functional.Group=="Overall biomass")
All$Scenario <- factor(All$Scenario,levels=c("Climate","Current Land Use","Maximum Land Use"))
All$Region <- factor(All$Region, levels = c("Brazil", "Finland", "France", "Namibia"))

All_Regions <- ggplot(data=Regions)+
  geom_bar(stat="identity",aes(x=Effect.Size,y=Region,fill=Climate),position="dodge")+
  geom_errorbar(stat="identity",aes(xmax=Upper.CI,xmin=Lower.CI,y=Region,group=Climate),position=position_dodge(0.9),linewidth=0.2,width=0.2)+
  facet_wrap(vars(Scenario), ncol = 3, labeller = as_labeller(labels_tot_bio))+#,scales = "free_x")+
  scale_fill_brewer(palette="Set2",direction = 1)+  
  scale_x_continuous(limits=c(-10,1),breaks=seq(-10,1,2))+
  scale_y_discrete(limits = rev(levels(All$Region))) +  # Reverse the y-axis order
  labs(x='Effect size (% change)', y = "Region",col="Climate scenario")+
  theme(axis.text.y = element_text(size=6,hjust=0))+
  geom_vline(xintercept=0, color='black', linetype='dashed') +
  theme_classic()+
  theme(axis.text.y = element_text(size=6,hjust=0))+
  theme(axis.text.x = element_text(size=6,hjust=0))+
  theme(plot.title = element_text(face = "bold", hjust =0.5,size=8))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),size=6,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),size=6,face="bold")) +
  theme(legend.position = "bottom",legend.text = element_text(size=6),legend.title = element_text(size=6,face="bold"))+
  theme(strip.text = element_text(size = 6))

All_Regions

#save pdf file in dina4
ggsave(paste0(figpath,sep="_","Fig2.pdf"),width=110, height=75, units="mm")
