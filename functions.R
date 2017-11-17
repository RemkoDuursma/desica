
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
  simdfr$precip <- 0
  
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


fsig_tuzet <- function(psil, psiv, sf){
  (1 + exp(sf*psiv)) / (1 + exp(sf*(psiv - psil)))
}

fsig_hydr <- function(P, SX, PX, X=50){
  
  P <- abs(P)
  PX <- abs(PX)
  X <- X[1] # when fitting; vector may be passed but X cannot actually vary.
  V <- (X-100)*log(1-X/100)
  p <- (P/PX)^((PX*SX)/V)
  relk <- (1-X/100)^p
  
  return(relk)
}

desica <- function(met=NULL,
                   psiv=-2,
                   sf=8,
                   g1=5,
                   
                   Ca=400,
                   Cs = 100000,
                   Cl = 10000,
                   
                   kpsat=3,
                   p50 = -4,
                   s50 = 30,
                   gmin = 10, # mmol m-2 s-1
                   
                   psil0=-1, 
                   psist0=-0.5, 
                   timestepmin=1,   # seconds
                   thetasat=0.5, 
                   sw0=0.4,
                   AL=2.5,
                   soildepth=1,
                   groundarea=1,
                   b=6,
                   psie= -0.8*1E-03,
                   Ksat=20,
                   Lv=10000,
                   
                   keepwet=FALSE,
                   stopsimdead=TRUE,
                   plcdead=88,
                   mf=NA,
                   LMA=NA
                   ){
  

  if(is.null(met)){
    stop("Must provide met dataframe with VPD, Tair, PPFD, precip (optional)")
  }
  
  # Increase length of met to conform with timestep of solution.
  # met will have 15min timestep, so if solution timestep (timestepmin) is 1,
  # repeat each row 15 times
  met$row <- 1:nrow(met)
  met$day <- rep(1:(nrow(met) / (24*4)), each=24*4)
  ntimes <- 15 / timestepmin
  ind <- rep(1:nrow(met), each=ntimes)
  met <- met[ind,]
  
  n <- nrow(met)

  timestep_sec <- 60*timestepmin
  
  LAI <- AL / groundarea
  soilvolume <- groundarea * soildepth
  
  Eleaf <- psil <- psist <- psis <- sw <- ks <- kp <- psirs <- Jsl <- Jrs <- rep(NA, n)
  krst <- kstl <- rep(NA,n)
  psil[1] <- psil0
  psist[1] <- psist0
  sw[1] <- sw0
  psis[1] <- psie*(sw0/thetasat)^-b
  Eleaf[1] <- 0
  ks[1] <- ksoil_fun(psis[1], Ksat, psie, b, LAI, soildepth=soildepth, Lv=Lv)
  
  psirs[1] <- psis[1]
  
  d <- list()
  
  for(i in 2:n){
    
    # Plant hydraulic conductance
    # Note how it depends on previous timestep stem water potential.
    kp[i] <- kpsat * fsig_hydr(psist[i-1], s50, p50)    
    
    # from soil to stem pool
    krst[i] <- 1 / (1/ks[i-1] + 1/(2*kp[i]))
    
    # from stem pool to leaf
    kstl[i] <- 2*kp[i]
    
    # packageVersion("plantecophys") >= "1.2-6"
    p <- Photosyn(VPD=met$VPD[i], 
                  gsmodel="BBdefine", 
                  BBmult=(g1/Ca)*fsig_tuzet(psil[i-1], psiv, sf), 
                  Tleaf=met$Tair[i],
                  PPFD=met$PPFD[i],
                  Ca=Ca)
    
    # Leaf transpiration (mmol m-2 s-1)
    Eleaf[i] <- p$ELEAF + (met$VPD[i]/101)*gmin
    
    # Xu method.
    # Can write the dynamic equation as: dPsil_dt = b + a*psil
    # Then it follows (Xu et al. 2016, Appendix, and Code).
    bp <- (AL * 2 * kstl[i]*psist[i-1] - AL*Eleaf[i])/Cl
    ap <- -(AL * 2 * kstl[i] / Cl)
    psil[i] <- ((ap*psil[i-1] + bp)*exp(ap*timestep_sec) - bp)/ap
    
    # Flux from stem to leaf= change in leaf storage, plus transpiration
    Jsl[i] <- (psil[i] - psil[i-1]) * Cl / timestep_sec + AL*Eleaf[i]
    
    # Update stem water potential
    # Also from Xu et al. 2016.
    bp <- (AL * 2 * krst[i] * psis[i-1] - Jsl[i])/Cs
    ap <- -(AL * 2 * krst[i] / Cs)
    psist[i] <- ((ap*psist[i-1] + bp)*exp(ap*timestep_sec) - bp)/ap
    
    # flux from soil to stem = change in stem storage, plus Jrl
    Jrs[i] <- (psist[i] - psist[i-1])*Cs/timestep_sec + Jsl[i]
    
    # Soil water increase: precip - transpiration (units kg total timestep-1)
    # (Note: transpiration is part of Jrs).
    water_in <- groundarea * met$precip[i] - timestep_sec*1E-06*18*Jrs[i]
    sw[i] <- pmin(1, sw[i-1] + water_in / (soilvolume * thetasat * 1E03))  # sw in units m3 m-3
    
    # for debugging / comparison.
    if(keepwet)sw[i] <- sw[i-1]
    
    # Update soil water potential and soil-to-root hydraulic conductance
    psis[i] <- psie*(sw[i]/thetasat)^-b
    ks[i] <- ksoil_fun(psis[i], Ksat, psie, b, LAI, soildepth=soildepth)
    
    # root surface water potential (not actually used?)
    psirs[i] <- psis[i] - Jrs[i]*1E03 / ks[i]
    
    if(stopsimdead){
      plc <- 100*(1 - kp[i]/kpsat)
      if(plc > plcdead)break
    }
  }
  
  d <- cbind(met[-1,], 
             data.frame(psis=psis, sw=sw, Eleaf=Eleaf, 
                        psirs=psirs, ks=ks, psil=psil, kp=kp,
                        psist=psist, Jsl=Jsl, Jrs=Jrs)[-1,])
  
  d$plc <- 100 * (1 - d$kp / kpsat)
  d$Eplant <- AL * Eleaf[-1]
  d$t <- 1:nrow(d)
  
  # stem water storage
  d$Wstem <- d$psist * Cs
  
  # when stopsimdead, last many rows are NA
  d <- d[!is.na(d$psis),]
  
  return(d)
}


plot_desica <- function(d, p50=-3){
  
  par(mfrow=c(2,2), mar=c(4,4,1,1))
  
  with(d, {
    plot(t, psis, type='l', ylim=c(-8,0))
    lines(t, psil, col="forestgreen")
    lines(t, psist, col="blue2")
    abline(h=p50)
  })
  
  plot(d$sw, type='l')
  
  with(d, {
    plot(t, Jsl, type='l')
    lines(t, Jrs, col="blue2")
  })
  
  plot(d$plc, type='l')
}



