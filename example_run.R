
# Load functions, packages
source("R/load.R")


# Make met dataframe
# Timestep is 15min (default)
# Diurnals are simulated, but you can of course replace this
# dataframe with your own measurements.
met <- make_simdfr(Tmin=10, RH=30, ndays=200)



# Run desica
# See desica() in R/functions.R for full list of arguments
run1 <- desica(met,
               psist0=0,
               AL=6,      # leaf area (m2)
               p50=-4,    # MPa
               psiv=-3,   # MPa (for Tuzet model)
               gmin=10,   # mmol m-2 s-1
               Cl=10000,  # Leaf capacitance (mmol MPa-1) (total plant)
               Cs=120000, # Stem capacitance (mmol MPa-1)
               runtwice=TRUE,  # each timestep run twice, for numerical stability (default, can omit)
               stopsimdead=TRUE) # stop simulation when PLC>88 (default, can omit)

# Standard plot
plot_desica(run1)

# Note: run1 is a dataframe, with:
# - met variables (PPFD, Tair, vpd, precip)
# - Eleaf (mmol m-2 s-1)
# - psil MPa (leaf water potential)
# - psist (stem water potential)
# - psis (soil water potential)
# - sw (soil volumetric water content)
# - ks (soil conductance, mmol m-2 s-1 MPa-1) !!NOTE: ks is too sensitive to psis in the current setup!!
# - kp (plant hydr. conductance, mmol m-2 s-1 MPa-1)
# - Jsl (flux from stem to leaf, mmol s-1)
# - Jrs (root water uptake, mmol s-1)
# - krst (conductance from soil to stem water store)
# - kstl (conductance from stem to leaf)
# - plc (percent loss conductivity)
# - Eplant (mmol s-1)

# simple summary:
# The transition from Phase1 to Phase2 is (arbitrarily) when
# soil conductance has dropped to < 5% of plant conductance.
summarize_desica(run1)


# Parameter sensitivity
# A simple setup to test effects of parameters.

runs <- list()
simdfr <- data.frame(gmin=seq(10,40,by=10))
n <- nrow(simdfr)

for(i in 1:n){
  runs[[i]] <- desica(met,  #!! setup the met dataframe first
               psist0=0,
               AL=8,      # leaf area (m2)
               p50=-3,    # MPa
               psiv=-2,   # MPa (for Tuzet model)
               gmin=simdfr$gmin[i],   # mmol m-2 s-1
               Cl=10000,  # Leaf capacitance (mmol MPa-1) (total plant)
               Cs=120000)
}

# Plot one of the simulations
plot_desica(runs[[2]])

# Get dry-down lengths
dry_len <- summarize_desica(runs)

# Plot time to mortality (Phase 2 only) against parameters
plot(simdfr$gmin, dry_len$phase2, type='o',
     xlab="gmin", ylab="Time to mortality (days)")



