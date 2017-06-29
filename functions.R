
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
                      Cs = 100 * 1000 / 18,  # mol (100 liters)
                      AL = 20,
                      VPD = 2,
                      Cl = 50 * 1000 / 18,  # mol (50 liters)
                      p50 = -4,
                      s50 = 30,
                      gmin = 10, # mmol m-2 s-1
                      ...
){
  
  # constant 8000 is made up :)  
  ksoil_fun <- function(psis, Ksat, psie, b){
    if(psis==0){
      Ksat * 8000
    } else {
      8000 *  Ksat*(psie/psis)^(2 + 3/b)
    }
  }
  #curve(ksoil_fun(200, -3*10^-3, x, 8), from=0, to=-3, ylim=c(0,2))
  
  ks <- ksoil_fun(psis, Ksat, psie, b)  # mmol m-2 s-1 MPa-1  # soil
  kp <- kpsat * fweibull(abs(psist), s50, abs(p50))           # plant
  ktot <- 1 / (1/ks + 1/kp)                                   # total pathway
  krst <- 1 / (1/ks + 1/(2*kp))                               # from soil to stem pool
  kstl <- 2*kp                                                # from stem pool to leaf
  
  
  fsig_tuzet <- function(psil, psiv, sf){
    (1 + exp(sf*psiv)) / (1 + exp(sf*(psiv - psil)))
  }
  
  # packageVersion("plantecophys") >= "1.2-6"
  p <- Photosyn(VPD=VPD, gsmodel="BBdefine", BBmult=(g1/Ca)*fsig_tuzet(psil, psiv, sf), ...)
  
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



# desica <- function(psil0=-2, psist0=-1, n=1000, timestep=900, 
#                    psis=-1,
#                    PPFD=1000,
#                    ...){
#   
#   psil <- psist <- ks <- kp <- Eleaf <- rep(NA, n)
#   psil[1] <- psil0
#   psist[1] <- psist0
#   
#   if(length(psis) < n)psis <- rep(psis[1],n)
#   if(length(PPFD) < n)PPFD <- rep(PPFD[1],n)
#   
#   for(i in 2:n){
#     d <- desica_dt(psil=psil[i-1], psist=psist[i-1], psis=psis[i], PPFD=PPFD[i],...)
#     psil[i] <- psil[i-1] + timestep*d["dpsil_dt"]
#     psist[i] <- psist[i-1] + timestep*d["dpsist_dt"]
#     kp[i] <- d["kp"]
#     ks[i] <- d["ks"]
#     Eleaf[i] <- d["Eleaf"]
#   }
#   
#   return(data.frame(t=1:n, psil=psil, psist=psist, ks=ks, kp=kp, Eleaf=Eleaf))
# }


desicawb <- function(psil0=-2, 
                     psist0=-1, 
                     timestep=900, 
                     n=960,
                     thetasat=0.5, 
                     sw0=0.4,
                     AL=20,
                     soilvolume=1,   # m3
                     b=6,
                     psie= -0.8*1E-03,
                     PPFD=1000,
                     ...){
  
  psil <- psist <- psis <- ks <- kp <- Eleaf <- sw <- rep(NA, n)
  

  psil[1] <- psil0
  psist[1] <- psist0
  sw[1] <- sw0
  psis[1] <- psie*(sw0/thetasat)^-b
  
  if(length(psis) < n)psis <- rep(psis[1],n)
  if(length(PPFD) < n)PPFD <- rep(PPFD[1],n)
  
  for(i in 2:n){
    
    d <- desica_dt(psil=psil[i-1], psist=psist[i-1], 
                   b=b, psie=psie,
                   psis=psis[i-1], PPFD=PPFD[i],...)
    psil[i] <- psil[i-1] + timestep*d["dpsil_dt"]
    psist[i] <- psist[i-1] + timestep*d["dpsist_dt"]
    kp[i] <- d["kstl"]
    ks[i] <- d["krst"]
    Eleaf[i] <- d["Eleaf"]
    
    # Soil water decreases by amount Jrs (mol s-1)
    sw[i] <- sw[i-1] - timestep*1E-03*18*d["Jrs"] / (soilvolume * thetasat * 1E03)
    psis[i] <- psie*(sw[i]/thetasat)^-b
  }
  
  return(data.frame(t=1:n, psil=psil, psist=psist, 
                    ks=ks, kp=kp, Eleaf=Eleaf,
                    sw=sw, psis=psis))
}



