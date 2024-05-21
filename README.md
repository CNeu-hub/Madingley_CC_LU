# Model-based Impact Analysis of Climate Change and Land-use Intensification on Trohpic Networks 

## Overview 
The repository is divided into three parts: 

* The input folder contains the already aggregated simulation output for each region and climate scenario (the original files were too large to upload to a repository). 
* The analysis folder contains the scripts needed to run the analysis. 
* The model run folder contains the scripts needed to run the original model simulations. In order to run this script, additional climate data must be provided by the user (as required by MadingleyR). 
* The output folder contains all output figures and output tables created in the repository. 


The climate data used in this study can be found here: 

* [Historical climate scenario](https://doi.org/http://doi.org/10.22033/ESGF/CMIP6.4067)
* [SSP1-2.6 climate scenario](https://doi.org/http://doi.org/10.22033/ESGF/CMIP6.4067)
* [SSP5-8.5 climate scenario](https://doi.org/http://doi.org/10.22033/ESGF/CMIP6.4067)

Before the data can be used, they need to be pre-processed to meet the requirements of the MadingleyR package, for example as described in the supplementary material of the study (also available in this repository).

To give an insight into our simulations, we provide here an example of how the simulations in our study are performed with the MadingleyR package. To do this, we use 0.5 degree data (to match the region sizes of our original model simulations) from the standard MadingleyR input datasets.

We perform the example here with a

## Create the environment 
First we have to create the environment that is necessary to perform the simulations. 

1.) We need to install "MadingleyR", following the [instructions](https://github.com/MadingleyR/MadingleyR/blob/master/README.md).

2.) We need to load MadingleyR for performing the simulations and some additional functions provided on Github for preparing data input. 

```
library("MadingleyR")
library("terra")

#Source: Hoeks, S. (2022): 
#Function to crop spatial raster input of MadingleyR using spatial window
source("https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/crop_spatial_rasters_to_window.r")

# download and load 0.5 degree inputs
source("https://raw.githubusercontent.com/SHoeks/MadingleyR_0.5degree_inputs/master/DownloadLoadHalfDegreeInputs.R")

```

## Define study region 
For this example we use France as our study region.

```
#define spatial extent of region, meeting requirements of MadingleyR (min. longitude, max. longitude, min. latitude, max. latitude)
#spatial_window = c(25.5,28.5,61,69) #96 grid cells, works brilliant, Finland
#spatial_window = c(16,21,-22,-17) #90 grid cells, works brilliant, Namibia 
#spatial_window = c(-69,-61,-3,0) #96 grid cells, works brilliant, Brazil 
spatial_window = c(-1,6,46,49) #84 grid cells, works brilliant, France

#define region
region <- "France"

#create output path for regional output 
regionpath <- paste0(getwd(), sep = "/", "Output/", region, sep="/")

#check if region extent is chosen correctly 
plot_spatialwindow(spatial_window)

```
![](png/Study_Region_France.pdf)

envDIR = temp_path <- gettemppath()

#create MadingleyR model inputs
#here, we just use standard spatial input of madingleyR to provide an example 
sptl_inp = madingley_inputs("spatial inputs")
chrt_def = madingley_inputs("cohort definition")
stck_def = madingley_inputs("stock definition")
mdl_prms = madingley_inputs("model parameters")

#check if hanpp layer is correctly loaded to madingley
plot(sptl_inp$hanpp, main="HANPP in gC/m^2/year")

#crop the raster to spatial window extent (safe resources while running the model)
sptl_inp = crop_spatial_rasters_to_window(inputs = sptl_inp, spatial_window = spatial_window)

#plot hanpp raster to see if cropping worked 
plot(sptl_inp$hanpp)

```


# References
Voldoire, A. 2019a. CNRM-CERFACS CNRM-CM6-1-HR model output prepared for CMIP6 CMIP historical (No. 20191021). doi: https://doi.org/http://doi.org/10.22033/ESGF/CMIP6.4067.

Voldoire, A. 2019b. CNRM-CERFACS CNRM-CM6-1-HR model output prepared for CMIP6 ScenarioMIP ssp126 (No. 20200127). doi: https://doi.org/http://doi.org/10.22033/ESGF/CMIP6.4067.

Voldoire, A. 2019c. CNRM-CERFACS CNRM-CM6-1-HR model output prepared for CMIP6 ScenarioMIP ssp585 (No. 20191202). doi: https://doi.org/http://doi.org/10.22033/ESGF/CMIP6.4067.

