
# Load functions, packages
source("R/load.R")


# Make met dataframe
# Timestep is 15min (default)
# Diurnals are simulated, but you can of course replace this
# dataframe with your own measurements.
met <- make_simdfr(Tmin=10, RH=30, ndays=200)



checkwatbal <- function(run, i=nrow(run)){
  
  # liters (1000 = nr liters for 1m3 soil volume, setting used)
  water_lost <- (run$sw[1] - run$sw[i]) * 1000
  
  # transpiration;
  water_used <- sum(run$Jrs[1:i] * 15 * 60 * 1E-03 * 18 * 1E-03, na.rm=TRUE)
  
  return(c(soilwater_decrease=water_lost, rootwateruptake=water_used))
  
}


# Run desica
# See desica() in R/functions.R for full list of arguments
run1 <- desica(met,
               psist0=0,
               AL=6,      # leaf area (m2)
               p50=-4,    # MPa
               psiv=-3,   # MPa (for Tuzet model)
               gmin=10,   # mmol m-2 s-1
               Cl=10000,  # Leaf capacitance (mmol MPa-1) (total plant)
               Cs=120000) # Stem capacitance (mmol MPa-1)


# Standard plot
plot_desica(run1)
checkwatbal(run1)


# Note some numerical noise in the above simulation.
# Redo by running each timestep twice.
run1_2 <- desica(met,
               AL=6,      # leaf area (m2)
               p50=-4,    # MPa
               psiv=-3,   # MPa (for Tuzet model)
               gmin=10,   # mmol m-2 s-1
               Cl=10000,  # Leaf capacitance (mmol MPa-1) (total plant)
               Cs=120000, # Stem capacitance (mmol MPa-1)
               runtwice=TRUE)

plot_desica(run1_2)
checkwatbal(run1_2)



