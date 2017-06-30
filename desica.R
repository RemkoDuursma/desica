
library(plantecophys)
library(fitplc)

source("functions.R")






simdfr <- make_simdfr(Tmin=10, RH=30, ndays=60)
simdfr$precip <- 0

d <- desicawb(psil0=0, psist0=0, psiv=-2.5, sf=1.5, 
            kpsat=1.5,
            Cl=100,
            b=6,
            Cs=8000,
            gmin=10,
            AL=30,
            met=simdfr,
            fracrootresist=0.01,
            timestep=15*60)



with(d, {
  plot(t, psis, type='l', ylim=c(-4,0))
  lines(t, psil, col="forestgreen")
  lines(t, psist, col="blue2")
})



with(d, {
  plot(t, Jsl, type='l', xlim=c(0,500))
  lines(t, Jrs, col="blue2")
})





