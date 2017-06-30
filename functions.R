
diurnal_sim <- function(PPFDmax=2000, RH=30, 
                        Tmax=30, Tmin=15,
                        daylength=12, 
                        sunrise=8, 
                        timestep=15, 
                        lag=0.5){
  
  # PPFD
  crv <- function(relt, ymax)(ymax/2)*(1+sin((2*pi*relt)+1.5*pi))
  
  r <- seq(0, 24*60 - timestep, by=timestep)
  p <- rep(0, length(r))
  ii <- which(r >= (sunrise*60) & r < (sunrise*60 + daylength*60))
  p[ii] <- crv( (r[ii] - sunrise*60)/(daylength*60), PPFDmax)
  
  
  # during daylight
  partondiurnal <- function(timeh, Tmax,Tmin,daylen, lag){
    (Tmax-Tmin)*sin((pi*timeh)/(daylen+2*lag)) + Tmin
  }
  
  ta <- rep(NA, length(r))
  td <- partondiurnal(seq(0, daylength - timestep/60, by=timestep/60), Tmax, Tmin, daylength,lag)
  ta[ii] <- td
  
  nnight <- length(r)-length(td)
  tdec <- ta[max(ii)] + (Tmin - ta[max(ii)])* 1:nnight / nnight
  after <- (max(ii)+1):length(r)
  ta[after] <- tdec[1:length(after)]
  before <- 1:(min(ii)-1)
  ta[before] <- tdec[(length(after)+1):length(tdec)]
  
  vpd <- RHtoVPD(RH, ta)
  
  return(data.frame(PPFD=p, Tair=ta, VPD=vpd))
}


make_simdfr <- function(..., ndays=40){
  dm <- diurnal_sim(...)
  simdfr <- vector("list", length=ndays)
  for(i in 1:length(simdfr))simdfr[[i]] <- dm
  simdfr <- do.call(rbind, simdfr)
  simdfr$day <- rep(1:ndays, each=96)
  
return(simdfr)
}




ksoil_fun <- function(psis, 
                      Ksat, 
                      psie, 
                      b,
                      LAI,
                      Lv=50000,
                      rroot=1E-06,
                      soildepth=1
                      ){
  
  Ks <- Ksat*(psie/psis)^(2 + 3/b)
  Ks[psis == 0] <- Ksat
  
  rcyl <- 1/sqrt(pi*Lv)
  Rl <- Lv * soildepth
  
(Rl/LAI)*2*pi*Ks/log(rcyl/rroot)
}




desica_dt <- function(psist = -1,
                      psil=-2,
                      psis = 0,
                      psie= -0.8*1E-03,
                      psimin=-2.2,  # not used
                      psiv=-2,
                      sf=8,
                      g1=5,
                      Ksat=200,
                      kpsat=3,
                      b=6,
                      Ca=400,
                      Tair=25,
                      LAI=2,
                      soildepth=1,
                      Cs = 100 * 1000 / 18,  # mol (100 liters)
                      AL = 20,
                      VPD = 2,
                      Cl = 50 * 1000 / 18,  # mol (50 liters)
                      p50 = -4,
                      s50 = 30,
                      gmin = 10, # mmol m-2 s-1
                      fracrootresist=0.1,
                      ...
){

  ks <- ksoil_fun(psis, Ksat, psie, b, LAI, soildepth=soildepth)  # mmol m-2 s-1 MPa-1  # soil
  kp <- kpsat * fweibull(abs(psist), s50, abs(p50))           # plant
  ktot <- 1 / (1/ks + 1/kp)                                   # total pathway (not used)
  
  krst <- 1 / (1/ks + 1/(kp/fracrootresist))                               # from soil to stem pool
  kstl <- kp/(1 - fracrootresist)                                           # from stem pool to leaf
  
  
  fsig_tuzet <- function(psil, psiv, sf){
    (1 + exp(sf*psiv)) / (1 + exp(sf*(psiv - psil)))
  }
  
  # packageVersion("plantecophys") >= "1.2-6"
  p <- Photosyn(VPD=VPD, gsmodel="BBdefine", 
                BBmult=(g1/Ca)*fsig_tuzet(psil, psiv, sf), 
                Tleaf=Tair,
                ...)
  
  #p <- PhotosynWP(VPD=VPD, kp=kstl, psis=psist, psimin=psimin)
  Eleaf <- p$ELEAF + (VPD/101)*gmin  # mmol m-2 s-1
  
  Et <- 1E-03 * AL * Eleaf  # mol s-1
  
  Jrs <- 1E-03 * AL * krst * (psis - psist) # mol s-1
  Jsl <- 1E-03 * AL * kstl * (psist - psil) # mol s-1
  
  dpsist_dt <- (Jrs - Jsl)/Cs
  dpsil_dt <- (Jsl - Et)/Cl
  
  return(c(dpsist_dt=dpsist_dt, 
           dpsil_dt=dpsil_dt,
           Eleaf=Eleaf,
           kp=kp,
           ks=ks,
           krst=krst,
           kstl=kstl,
           Jsl=Jsl,
           Jrs=Jrs
           
  ))
}


desicawb <- function(psil0=-2, 
                     psist0=-1, 
                     timestep=900,   # seconds
                     thetasat=0.5, 
                     sw0=0.5,
                     AL=20,
                     soildepth=1,
                     groundarea=4,
                     b=6,
                     psie= -0.8*1E-03,
                     met=NULL,   # 
                     keepwet=FALSE,
                     ...){
  
  if(is.null(met)){
    stop("Must provide met dataframe with VPD, Tair, PPFD, precip (optional)")
  }
  n <- nrow(met)
  
  LAI <- AL / groundarea
  soilvolume <- groundarea * soildepth
  
  psil <- psist <- psis <- sw <- rep(NA, n)

  psil[1] <- psil0
  psist[1] <- psist0
  sw[1] <- sw0
  psis[1] <- psie*(sw0/thetasat)^-b
  
  d <- list()
  
  for(i in 2:n){
    
    d[[i]] <- desica_dt(psil=psil[i-1], psist=psist[i-1], 
                   b=b, psie=psie, LAI=LAI, soildepth=soildepth,
                   psis=psis[i-1], 
                   VPD=met$VPD[i], PPFD=met$PPFD[i], Tair=met$Tair[i], 
                   ...)
    
    # Simple difference equations.
    psil[i] <- psil[i-1] + timestep*d[[i]]["dpsil_dt"]
    psist[i] <- psist[i-1] + timestep*d[[i]]["dpsist_dt"]
    
    # Soil water increase: precip - transpiration (units kg total timestep-1)
    water_in <- groundarea * met$precip[i] - timestep*1E-03*18*d[[i]]["Jrs"]
    sw[i] <- pmin(1, sw[i-1] + water_in / (soilvolume * thetasat * 1E03))  # sw in units m3 m-3
    
    # for debugging
    if(keepwet)sw[i] <- sw[i-1]
    
    # update soil water potential
    psis[i] <- psie*(sw[i]/thetasat)^-b
  }
  
  
  d <- as.data.frame(do.call(rbind,d))
  d <- cbind(d, met[-1,], data.frame(psis=psis, psil=psil, psist=psist)[-1,])
  
  d$t <- 1:nrow(d)
  
return(d)
}



