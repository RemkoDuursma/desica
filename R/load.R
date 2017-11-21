
# Load custom functions
source("R/functions.R")

# Missing packages
if(!require("pacman"))install.packages("pacman")
pacman::p_load(plantecophys,dplyr,crayon,readxl)


if(!dir.exists("cache"))dir.create("cache")
