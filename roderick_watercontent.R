
# rod <- read.csv("roderick1999.csv", stringsAsFactors = FALSE)
# 
# fixcol <- function(x)as.numeric(gsub("·",".",x))
# rod <- as.data.frame(lapply(rod[,c("ml","md","Vl","Va","As","An","Cd","Nd")], fixcol))

# write.csv(rod, "roderick1999_clean.csv", row.names=FALSE)

pacman::p_load(smatr, dplyr, plover)

rod <- read.csv("roderick1999_clean.csv") %>%
  mutate(LMA = 10^6 * md / As,
         LWA = 10^6 * ml /As)


fit <- sma(LWA ~ LMA, data=rod)

# scales with LWA ~ LMA^0.6
#fit2 <- sma(LWA ~ LMA, data=rod, log="xy")


par(mar=c(5,5,1,1), las=1, xaxs="i", yaxs="i",
    cex.axis=0.9, cex.lab=1.1)
with(rod, plot(LMA, LWA, 
               pch=19, col="slategrey",
               xlim=c(0,200),
               ylim=c(0,350),
               xlab=expression(LMA~~(gDM~m^-2)),
               ylab=expression(LWA~~(g~H[2]*O~m^-2))))
plot(fit, add=TRUE, type='l', col="black")
abline(0,1, lty=3)












