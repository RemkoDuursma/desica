
library(plantecophys)
#library(fitplc
library(dplyr)
library(crayon)

source("functions.R")

dead_sim <- function(met, p50, psiv, gmin, capac, return=c("summary","simulation")){

  return <- match.arg(return)
  d <- desica2(met=met,
              psil0=-1, 
              psist0=-0.5, 
              p50=p50,
              psiv=psiv, 
              sf=5, 
              sw0=0.4,
              kpsat=1.5,
              
              b=6,
              
              gmin=gmin,
              soildepth=0.8,
              AL= 3, 
              
              Cl=0.1*capac,   # mmol MPa-1
              Cs=0.9*capac,  
              
              plcdead=88,
              
              stopsimdead=T)
  
  phase2 <- subset(d, ks < 0.05*max(kp))
  simtot <- nrow(d) / (4*24)
  sim2 <- nrow(phase2) / (4*24)
  sim1 <- simtot - sim2
  
  if(return == "summary"){
    return(c(phase1=sim1, phase2=sim2, p50=p50, psiv=psiv, gmin=gmin, 
      capac=capac, plcfinal=d$plc[nrow(d)]))
  } 
  if(return == "simulation"){
    return(d)
  }
}
  


# Take 1.
  met <- make_simdfr(Tmin=10, RH=30, ndays=200)
  dead_sim(met, p50=-5, psiv=-4, gmin=30, capac=300000)


  d <- dead_sim(met, p50=-5, psiv=-4, gmin=30, capac=300000, return="simulation")
  plot_desica(d)

  
  
  
# Take 2.

  met <- make_simdfr(Tmin=10, RH=30, ndays=200)
  p  <- expand.grid(p50 = -2:-8, sm = seq(0,2,length=5), capac=seq(25000,50000, by=5000), 
                         gmin=seq(5, 25, by=5))
  
  l <- list()
  for(i in seq_len(nrow(simdfr))){
    l[[i]] <- try(dead_sim(met,
                           p50=p$p50[i], psiv=p$p50[i] + p$sm[i], 
                           gmin=p$gmin[i], capac=p$capac[i]))
    
    cat(cyan("Simulation ") %+% chr(i) %+% cyan(" completed\n"))
  }

  
res <- do.call(rbind, l) %>% as.data.frame %>%
  filter(plcfinal > 88) # actual death!
  


# Take 3. 

# Sample psiv and p50 from Martin-StPaul
met <- make_simdfr(Tmin=10, RH=30, ndays=200)
mart <- martin[,c("P50","Pgs50")] %>% filter(!is.na(P50), !is.na(Pgs50))

mrt <- paste(mart[[1]], mart[[2]], sep="_")

p  <- expand.grid(comb=mrt, capac=seq(25000,50000, by=5000), 
                  gmin=seq(2.5, 20, length=8), stringsAsFactors = FALSE)

m <- strsplit(p$comb, "_")
p$p50 <- as.numeric(sapply(m, "[", 1))
p$psiv <- as.numeric(sapply(m, "[", 2))

g <- list()
for(i in seq_len(nrow(p))){
  g[[i]] <- try(dead_sim(met,
                         p50=p$p50[i], psiv=p$psiv[i], 
                         gmin=p$gmin[i], capac=p$capac[i],
                         return="simulation"))
  
  cat(cyan("Simulation ") %+% chr(i) %+% cyan(" completed\n"))
}

saveRDS(g, "g.rds")


res <- do.call(rbind,g) %>% as.data.frame

res <- subset(res, plcfinal > 88)
with(res, plot(-capac*p50/gmin, phase2, pch=19, col=rainbow(4)[cut(gmin,4)]))



x <- subset(res, gmin == 10 & capac==30000)
with(x, plot(p50, psiv, pch=19, col=rainbow(5)[cut(phase2,5)]))

with(x, plot(p50, phase2))


x <- subset(res, capac==30000)
with(x, plot(p50, phase2, pch=19, col=as.factor(gmin)))

with(x, plot(-p50/gmin, phase2))



system.time(r1  <- dead_sim(met, p50=-3.75, psiv=-3.3, gmin=15, capac=150000, return="simulation"))

plot_desica(r1, -3.75)




met <- make_simdfr(Tmin=10, RH=30, ndays=200)

d1 <- desica(met)
d2 <- desica2(met)

plot_desica(d1)
plot_desica(d2)
