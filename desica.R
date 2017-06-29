
library(plantecophys)
library(fitplc)

source("functions.R")



dm <- diurnal_sim(PPFDmax=2000, RH=40, Tmin=10)
ndays <- 40
simdfr <- vector("list", length=ndays)
for(i in 1:length(simdfr))simdfr[[i]] <- dm
simdfr <- do.call(rbind, simdfr)
simdfr$psis <- seq(0, -3, length=nrow(simdfr))
simdfr$day <- rep(1:ndays, each=96)





d <- desicawb(psil0=0, psist0=0, psiv=-2.5, sf=1.5, 
            n=nrow(simdfr),
            kpsat=0.8,
            Cl=400,
            b=6,
            Cs=8000,
            gmin=0.005,
            PPFD=simdfr$PPFD,  timestep=15*60)
d <- cbind(d, simdfr)


with(d, {
  plot(t, psis, type='l', ylim=c(-4,0))
  lines(t, psil, col="forestgreen")
  lines(t, psist, col="blue2")
})


