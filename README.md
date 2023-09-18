---
The current project is used to identify unknown ligands of human CAR from gut micriobiome extracts
by Hui Peng (UofT), September 15, 2023 

## Example
This example uses the raw data to show how to identify unknown ligands. Only positive ion mode is analyzed
Multiple packages including the 'xcms' package need to be installed
---
Step 1: loading packages
```{r, message=FALSE, warning=FALSE}
  rm(list=ls())
  library(xcms)
  library(devtools)
  library(xlsx)
  library(isopat)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(ggrepel)
  library(openxlsx)
  path<-getwd()
  polarity<-1
```
Step 2: extracting peak features, please call the functions 'Peakextract' and 'LigandFeatureID'
```{r, message=FALSE, warning=FALSE}
  setwd(path)

  #extracting peaks, intensity cutoff is 10^6, ppm = 2.5
  peaks.raw<-Peakextract(10^5,2.5)
  
  #'identifying differentiated peaks,
  Control<-'NC'#the negative controls
  peaks<-LigandFeatureID(peaks.raw,'positive',Control)#This is the function extracting differential peaks compared to control
  
  ##Note that three documents including 'CAR_positive', 'CAR_MB_positive.csv' and 'CAR_MB_positive.tiff' are generated.
  #The 'CAR_positive.csv' document contains all raw data and differenated ratios and pvalues
  #The 'CAR_MB_positive.csv' document contains the peaks of putative compounds (fold>100, p<0.01)
  #The TIFF document showed the volcano plot
  
```