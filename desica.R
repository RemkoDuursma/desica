  
library(plantecophys)
library(fitplc) # not needed?

source("functions.R")



d <- desica(met=make_simdfr(Tmin=10, RH=30, ndays=100),
            psil0=-1, 
            psist0=-0.5, 
            p50=-5,
            psiv=-2.5, 
            sf=5, 
            sw0=0.2,
            kpsat=1.5,
            Cl=500,
            b=6,
            Cs=8000,
            gmin=10,
            soildepth=1,
            AL=10,
            
            plcdead=88,
            timestep=1*60,
            stopsimdead=F)


# actual water uptake by roots has to be equal to this? (at steady state at least...)
# d$Jrs_a <- with(d, ks * (psis - psist))


plot_desica <- function(d){
  
  par(mfrow=c(2,2), mar=c(4,4,1,1))

  with(d, {
    plot(t, psis, type='l', ylim=c(-8,0))
    lines(t, psil, col="forestgreen")
    lines(t, psist, col="blue2")
    abline(h=-3)
  })
  
  plot(d$sw, type='l')
  
  with(d, {
    plot(t, Jsl, type='l')
    lines(t, Jrs, col="blue2")
  })
  
}
        
plot_desica(d)



# Random figures for ESA
  
  p50 <- -5
  s50 <- 30
  psiv <- -2
  sf <- 5
  
  plot(1, type='n', xlim=c(0,abs(1.8*p50)), ylim=c(0,1),
       ylab="Relative value",
       xlab="Water potential (-MPa)")
  curve(fsig_tuzet(-x, psiv=psiv, sf=sf), add=TRUE, col="blue2", lwd=2)
  curve(1-fsig_hydr(x, SX=s50, PX=p50), add=TRUE, col="red2", lwd=2)
  
  martin <- readxl::read_excel("data/DataBase.xlsx", sheet=4)
  names(martin)[6] <- "Pgs50"
  with(martin, plot(P50, Pgs50, pch=19, col=as.factor(Group),
                    xlim=c(-15,0), ylim=c(-5,0)))  
  abline(0,1)
  
  
  ge <- read.csv("data/ge_R.csv")
  plc <- read.csv("data/plc_R.csv")
  

  eme_plc <- subset(plc, species == "Eme" & !is.na(PLC..))
  f <- fitplc(eme_plc, varnames=c(PLC="PLC..", WP="Mpa"))
    
with(subset(ge, species == "Eme"), 
     plot(MPa, gssat/max(gssat), 
          pch=19, col="blue",
          xlim=c(0,8)))
curve(fsig_tuzet(-x, psiv=-1.5, sf=8), add=TRUE, col="blue2", lwd=2)

with(eme_plc, points(Mpa, PLC.. / 100, pch=15, col="red"))
plot(f, add=T, what="PLC", multiplier=0.01,linecol="red2",linelwd=2,
     plotPx=F, plotci=F,plotdata=F,px_ci="none")



baad_file <- "data/baad.rds"

if(!file.exists(baad_file)){
  url <- "https://github.com/dfalster/baad/releases/download/v1.0.1/baad.rds"
  download.file(url, baad_file, mode="wb")
}
baad <- readRDS(baad_file)$data

mss <- subset(baad, !is.na(m.ss))

with(mss, plot(log10(m.lf), log10(m.ss), pch=19, col=as.factor(studyName)))
abline(2,1)


library(smatr)
sfit <- sma(m.ss ~ m.lf * studyName, data=mss, log="xy")
plot(sfit, axes=FALSE, pch=19, 
     xlab="Leaf mass (kg)", ylab="Sapwood mass (kg)")
magicaxis::magaxis(side=1:2, unlog=1:2)

sfit0 <- sma(m.ss ~ m.lf, data=mss, log="xy")
coef(sfit0)





