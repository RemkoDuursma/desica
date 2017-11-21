

metsub <- met[1:(2*96),]
r1 <- desica(metsub, p50=-4, psiv=-3, 
                  gmin=10, Cl=0.1*150000,
                     Cs=0.9*150000)

     
#1
windows(6,5)
par(mar=c(4,4,1,1), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2)
with(r1, plot(24*t/96, Eleaf*2.5, 
              xlab="Time (hours)",
              ylab=expression(Water~flux~~(mmol~s^-1)),
              type='l', ylim=c(0,40)))
with(r1, lines(24*t/96, Jrs, type='l', col="red", lty=5))
legend("topleft", c("Transpiration","Root water uptake"),
                    lty=c(1,5), col=c("black","red"))
dev.copy2pdf(file="output/1.pdf")

#2
windows(6,5)
par(mar=c(4,4,1,1), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2)
with(r1, plot(24*t/96, psil, 
              xlab="Time (hours)",
              ylab=expression(Water~potential~~(MPa)),
              type='l', ylim=c(-3,0)))
with(r1, lines(24*t/96, psist, type='l', col="red", lty=5))
legend("bottomleft", c("Leaf","Stem"),lty=c(1,5), col=c("black","red"),
       inset=0.03)
dev.copy2pdf(file="output/2.pdf")



#3
library(fitplc)
ge <- read.csv("data/ge_R.csv")
plc <- read.csv("data/plc_R.csv")


eme_plc <- subset(plc, species == "Eme" & !is.na(PLC..))
f <- fitplc(eme_plc, varnames=c(PLC="PLC..", WP="Mpa"))

windows(7,5)
par(mar=c(4,4,1,4), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2)
with(subset(ge, species == "Eme"), 
     plot(MPa, gssat/max(gssat)*100, 
          xlab="Water potential (MPa)",
          ylab="",
          pch=19, col="blue",
          xlim=c(0,8)))
curve(100*fsig_tuzet(-x, psiv=-1.5, sf=8), add=TRUE, col="blue2", lwd=2)
mtext(side=4, line=2.2, cex=1.2, text="% Loss hydraulic conductance (PLC)", col="red2")
mtext(side=2, line=2.2, cex=1.2, text="Stomatal conductance (%)", col="blue2")
axis(4)
with(eme_plc, points(Mpa, 100*PLC.. / 100, pch=15, col="red"))
plot(f, add=T, what="PLC", linecol="red2",linelwd=2,
     plotPx=F, plotci=F,plotdata=F,px_ci="none")
dev.copy2pdf(file="output/3.pdf")




#4
library(nlshelper)


windows(6,5)
par(mar=c(4,4,1,1), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2)
with(martin, plot(P50, Pgs50, pch=19, cex=1.1, col="olivedrab",
                  xlab=expression(P["50"]~~(MPa)),
                  ylab=expression(Pgs["50"]~~(MPa)),
                  xlim=c(-15,0), ylim=c(-8,0)))  
abline(0,1)

fit0 <- loess(Pgs50 ~ P50, data=martin, span=0.8)
martin$Pgs50_fit <- predict(fit0, newdata=martin)
martin <- martin[order(martin$P50),]
with(martin, lines(P50, Pgs50_fit, lty=5))

dev.copy2pdf(file="output/4.pdf")



#5 drydown
r2 <- desica(met, p50=-3, psiv=-2, 
             gmin=15, Cl=0.1*180000,
             Cs=0.9*180000)


with(r2, plot(t/96, psil, type='l', xlim=c(5,20), ylim=c(-2.5,0)))
with(r2, lines(t/96, psist, col="red"))
with(r2, lines(t/96, psis, lwd=2, col="cornflowerblue"))


# iterate once
r2_2 <- desica(met, p50=-3, psiv=-2, sw0=0.4,
               gmin=18, Cl=0.1*150000, Lv=30000,
               Cs=0.9*150000, runtwice=TRUE)

windows(7,4)
par(mar=c(4,4,1,4), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2)
with(r2_2, plot(t/96, psil, type='l',
                xlab="Time (days)",
                ylab="Water potential (MPa)"))
with(r2_2, lines(t/96, psist, col="red"))
with(r2_2, lines(t/96, psis*2, lwd=2, col="cornflowerblue"))
legend(0,-2.5,c("Leaf","Stem","PLC"), lty=1, col=c("black","red","darkgrey"),
       inset=0.02)

par(new=TRUE)
with(r2_2, plot(t/96, plc, type='l', col="darkgrey",
                ann=FALSE,
                axes=FALSE))
axis(4)
mtext(side=4, line=2.2, cex=1.2, text="PLC (%)")

dev.copy2pdf(file="output/5.pdf")


# 8
baad_file <- "data/baad.rds"

if(!file.exists(baad_file)){
  url <- "https://github.com/dfalster/baad/releases/download/v1.0.1/baad.rds"
  download.file(url, baad_file, mode="wb")
}
baad <- readRDS(baad_file)$data

mss <- subset(baad, !is.na(m.ss))

library(smatr)
sfit <- sma(m.ss ~ m.lf * studyName, data=mss, log="xy")

windows(6,5)
par(mar=c(4,4,1,1), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2,
    family="Gotham Narrow Book")
plot(sfit, axes=FALSE, pch=19, 
     xlab="Leaf mass (kg)", ylab="Sapwood mass (kg)")
magicaxis::magaxis(side=1:2, unlog=1:2)

sfit0 <- sma(m.ss ~ m.lf, data=mss, log="xy")
coef(sfit0)

legend("topleft", expression(M[SS] == 11.5 %*% M[F] ^ 1.33), bty='n', cex=1.2)
legend("bottomright", "BAAD\nFalster et al. 2015", bty='n', cex=1.2, inset=0.03)
dev.copy2pdf(file="output/8.pdf")