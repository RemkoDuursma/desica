
source("load.R")

# Make simulation dataframe (combinations of inputs)
# Sample psiv and p50 from Martin-StPaul
met <- make_simdfr(Tmin=10, RH=30, ndays=200)
mart <- martin[,c("P50","Pgs50")] %>% filter(!is.na(P50), !is.na(Pgs50))

mrt <- paste(mart[[1]], mart[[2]], sep="_")

p  <- expand.grid(comb=mrt, capac=seq(100000,150000, by=10000), 
                  gmin=seq(5, 20, by=2.5), stringsAsFactors = FALSE)
m <- strsplit(p$comb, "_")
p$p50 <- as.numeric(sapply(m, "[", 1))
p$psiv <- as.numeric(sapply(m, "[", 2))


# run simulation for each row
g <- list()
for(i in seq_len(nrow(p))){
  t1 <- proc.time()
  g[[i]] <- try(desica(met,
                         p50=p$p50[i], psiv=p$psiv[i], 
                         gmin=p$gmin[i], 
                        Cl=0.1*p$capac[i],
                       Cs=0.9*p$capac[i]))
  t2 <- proc.time()
  elaps <- round((t2-t1)[3], 1)
  cat(cyan("Simulation ") %+% chr(i) %+% cyan(" completed in ") %+% silver(chr(elaps)) %+% silver(" sec\n"))
}
saveRDS(g, "cache/sim.rds")


tims <- do.call(rbind, lapply(g, summarize_desica) ) %>% as.data.frame
res <- cbind(p, tims)

res <- subset(res, plcfinal > 88)


# Figure 
# Test of Blackman et al. 2016
windows(8,4)
par(mar=c(5,4,1,1), tcl=0.2, mgp=c(4.5,0.4,0), cex.lab=1.2,
    family="Gotham Narrow Book")
cols <- RColorBrewer::brewer.pal(7, "YlOrRd")
par(mar=c(4,4,1,1), tcl=0.2, mgp=c(2.2,0.4,0), cex.lab=1.2, mfrow=c(1,2))
with(res, plot(-10^-3*capac*p50/gmin, phase2, pch=21, 
               ylab="Time to death (Phase 2) (days)",
               ylim=c(0,200),
               xlim=c(0,500),
               xlab=expression(C %.% P["50"] / g[min]),
               bg=cols[cut(gmin,7)]))
#dev.copy2pdf(file="output/6a.pdf")
legend("topleft", c("5","20"), fill=cols[c(1,7)], title=expression(g[min]))
# Get rid of non-linearity
with(res, plot((-10^-3*capac*p50/(gmin^1.5))^0.5, phase2, 
               xlim=c(0,15),
               ylab="Time to death (Phase 2) (days)",
               xlab=expression( (C %.% P["50"] / g[min]^1.5)^0.5),
               ylim=c(0,200),
               pch=21, bg=cols[cut(gmin,7)]))
dev.copy2pdf(file="output/6b.pdf")


# Figure

library(ppcor)
library(extrafont)

part_1 <- pcor(res[,c("phase1","capac","p50","psiv","gmin")])
sens_1 <- sort(abs(part_1$estimate[1,-1]),T)

part_2 <- pcor(res[,c("phase2","capac","p50","psiv","gmin")])
sens_2 <- sort(abs(part_2$estimate[1,-1]),T)

labfun <- function(x){
  
  if(length(x) > 1){
    return(sapply(x, labfun))
  }
  switch(x,
         gmin = expression(g[min]),
         capac = "Capac",
         p50 = expression(P[50]),
         psiv = expression(Pgs[50]))
}

windows(8,4)
par(mfrow=c(1,2), family="Gotham Narrow Book", mar=c(3,4,3,1),
    cex.axis=1.2)
barplot(sens_1, xlab="", ylab="Partial correlation", main="Phase 1", col="lightgrey",ylim=c(0,1),
        names.arg=labfun(names(sens_1)))
barplot(sens_2, xlab="", ylab="Partial correlation", main="Phase 2", col="dimgrey",ylim=c(0,1),
        names.arg=labfun(names(sens_2)))
dev.copy2pdf(file="output/7.pdf")





