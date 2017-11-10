  
library(plantecophys)
library(fitplc) # not needed?

source("functions.R")



simdfr <- 
simdfr$precip <- 0

d <- desica(met=make_simdfr(Tmin=10, RH=30, ndays=50),
            psil0=-1, psist0=-0.5, psiv=-2.5, sf=1.5, 
            sw0=0.2,
              kpsat=1.5,
              Cl=500,
              b=6,
              Cs=8000,
              gmin=10,
              soildepth=1,
              AL=1,
              p50=-3,
              timestep=1*60)


with(d, {
  plot(t, psis, type='l', ylim=c(-4,0))
  lines(t, psil, col="forestgreen")
  lines(t, psist, col="blue2")
  abline(h=-3)
})


plot(d$sw)

with(d, {
  plot(t, Jsl, type='l', xlim=c(0,500))
  lines(t, Jrs, col="blue2")
})







