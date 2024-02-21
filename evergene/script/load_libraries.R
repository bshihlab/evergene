#### Load libraries
# Install libraries----------
if (!require("shinythemes", quietly = TRUE)) { install.packages("shinythemes") } 
if (!require("shinyvalidate", quietly = TRUE)) { install.packages("shinyvalidate") } 
if (!require("shinyWidgets", quietly = TRUE)) { install.packages("shinyWidgets") } 
if (!require("shinyhelper", quietly = TRUE)) { install.packages("shinyhelper") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") } 
if (!require("zip", quietly = TRUE)) { install.packages("zip") } 
if (!require("fst", quietly = TRUE)) { install.packages("fst") } 
if (!require("shinyBS", quietly = TRUE)) { install.packages("shinyBS") } 
if (!require("survival", quietly = TRUE)) { install.packages("survival") } 
if (!require("survminer", quietly = TRUE)) { install.packages("survminer") } 
if (!require("scales", quietly = TRUE)) { install.packages("scales") } 
if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") } 
if (!require("grid", quietly = TRUE)) { install.packages("grid") } 
if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") } 
if (!require("edgeR", quietly = TRUE)) { BiocManager::install("edgeR", quiet = TRUE) } 
#if (!require("svglite", quietly = TRUE)) { install.packages("svglite", quiet = TRUE) } 

# General
library(shiny)
library(stringr) # Makes working with strings simpler.
library(viridis) # Colour scheme.
library(zip) # Allows creation of zip files
library(fst) # Allows fast loading for saved r files
library(scales) # Allows fast loading for saved r files

# UI
library(shinythemes)
library(shinyvalidate)
library(shinyWidgets)
library(plyr)
library(shinyBS)
library(shinyhelper)
#

# Survival
library(survival) # For statistical analysis on survival data
library(survminer) # For plotting survival data

# PCA
library(factoextra) # For principal component analysis

# Plot
library(grid) # Allows fast loading for saved r files
library(ggplot2)
library(plotly) # Creates interactive graphs.
#library(svglight) # Creates save plots as svg files.

# Process count data
library(edgeR)
