


library(plantecophys)
library(dplyr)
library(crayon)
library(readxl)
source("R/functions.R")


martin <- read_excel("data/DataBase.xlsx", sheet=4)
names(martin)[6] <- "Pgs50"


if(!dir.exists("cache"))dir.create("cache")
