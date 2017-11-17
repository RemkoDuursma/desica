
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
  cat(cyan("Simulation ") %+% chr(i) %+% cyan(" completed in ") %+% silver(chr(elaps)) %+% silver(" sec"))
}
saveRDS(g, "cache/sim.rds")

# redo using summarize_desica
# normally don't save all simulations but we need to
# res <- do.call(rbind,g) %>% as.data.frame
# 
# res <- subset(res, plcfinal > 88)
# with(res, plot(-capac*p50/gmin, phase2, pch=19, col=rainbow(4)[cut(gmin,4)]))



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
