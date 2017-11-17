  
library(plantecophys)
#library(fitplc) # not needed?

source("functions.R")


# stem volume

#1000*11.45*mf^1.33 * (1/500) * 10^3 * 0.5 * 10^3 / 18

d <- desica(met=make_simdfr(Tmin=10, RH=30, ndays=200),
            psil0=-1, 
            psist0=-0.5, 
            p50=-5,
            psiv=-2.5, 
            sf=5, 
            sw0=0.4,
            kpsat=1.5,
            
            b=6,
            
            gmin=10,
            soildepth=1,
            mf = 1,   # kg
            LMA = 100, # g m-2
            AL= mf / (LMA / 1000),  # m2
            
            Cl=4000,   # mmol MPa-1
            Cs=40000,  
            
            plcdead=88,
            timestep=1*60,
            stopsimdead=T)



dead_sim <- function()



plot_desica(d)


plot(cumsum(d$Jrs), type='l')
lines(cumsum(d$Eplant), col="red")

# somewhat arbitrary point where phase 2 begins?
# ks has dropped to 5% of plant conductance
# soil is now super limiting
x <- subset(d, ks < 0.05*max(kp))
with(d, plot(t, Eleaf, type='l'))
with(x, lines(t, Eleaf, col="red"))



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


library(smatr)
sfit <- sma(m.ss ~ m.lf * studyName, data=mss, log="xy")
plot(sfit, axes=FALSE, pch=19, 
     xlab="Leaf mass (kg)", ylab="Sapwood mass (kg)")
magicaxis::magaxis(side=1:2, unlog=1:2)

sfit0 <- sma(m.ss ~ m.lf, data=mss, log="xy")
coef(sfit0)


# Sapwood capacitance
m.lf <- 1   # kg
wood_dens <- 0.8  # g cm-3 = t m-3 = kg dm-3
pore_frac <- 0.6  
cs <- 0.15  # RWC MPa-1  - specific capacitance

m.ss <- 11.2 * m.lf^1.33  # kg of sapwood
v.ss <- m.ss / wood_dens  # liters of sapwood volume (dm3)

Vs <- pore_frac * v.ss    # liters of water in saturated tree

Cs_l <- Vs * cs   # liters = kg MPa-1
Cs <- Cs_l * 1000 * 1/18 * 1000  # mmol MPa-1

# Leaf capacitance
# ???
LMA <- 100

LWA <- 200 # see roderick script, leaf water content per unit area (g m-2)

cl <- 0.1  # specific capacitance, RWC MPa-1. I don't know the value, so let's assume
           # same as sapwood, based on Ximeng's comparison

Cl_g <- cl * LWA  # g m-2 MPa-1

Cl_m2 <- Cl_g * 1000 / 18  # mmol m-2 MPa-1

# using m.lf = 1 from above, and LMA=100, we get:
AL <- 1000 * m.lf / LMA

Cl <- Cl_m2 * AL













