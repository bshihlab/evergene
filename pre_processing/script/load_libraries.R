#### Load libraries

# Install libraries----------
if (!require("bslib", quietly = TRUE)) { install.packages("bslib") } 
if (!require("shinythemes", quietly = TRUE)) { install.packages("shinythemes") } 
if (!require("shinyWidgets", quietly = TRUE)) { install.packages("shinyWidgets") } 
if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("ggfortify", quietly = TRUE)) { install.packages("ggfortify") } 
if (!require("umap", quietly = TRUE)) { install.packages("umap") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") } 
if (!require("zip", quietly = TRUE)) { install.packages("zip") } 
if (!require("fst", quietly = TRUE)) { install.packages("fst") } 
if (!require("shinyBS", quietly = TRUE)) { install.packages("shinyBS") } 
if (!require("shinydashboard", quietly = TRUE)) { install.packages("shinydashboard") } 
if (!require("survival", quietly = TRUE)) { install.packages("survival") } 
if (!require("survminer", quietly = TRUE)) { install.packages("survminer") } 

# General
library(shiny)
library(stringr) # Makes working with strings simpler.
library(viridis) # Colour scheme.
library(zip) # Allows creation of zip files
library(fst) # Allows fast loading for saved r files

# UI
library(shinythemes)
library(forcats)
library(plyr)
library(shinydashboard)
library(shinyBS)

# Biological data
library(survival) # For statistical analysis on survival data
library(survminer) # For plotting survival data
library(ggplot2)
library(ggcorrplot) # Allows visulaisation of a correlation matrix using ggplot2.
library(plotly) # Creates interactive, publication-suitable graphs.

# PCA
library(factoextra) # For principal component analysis
library(ggfortify) # Allows plotting of PCA and survival analysis.
library(umap) # Algorithm for dimensional reduction
