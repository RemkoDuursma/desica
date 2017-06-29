
library(plantecophys)


PhotosynWP <- function(..., 
                       VPD=2,
                       Patm=101,
                       kp = 2,
                       psis = 0,
                       psimin = -2){
  
  Emax <- kp * (psis - psimin)
  
  p <- Photosyn(VPD=VPD, Patm=Patm, ...)
  
  if(p$ELEAF > Emax){
    GS <- 1E-03 * Emax / (VPD / Patm)
    p <- Photosyn(GS=GS, VPD=VPD, Patm=Patm, ...)
  }
  
  p$Emax <- Emax
  
return(p)
  
}


  



  
  