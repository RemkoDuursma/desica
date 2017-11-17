
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
                    met_timestep = 15,
                    runs_per_timestep = 1,
                    
                    Ca = 400,  
                    sf=8,
                    g1=5,
                    
                    Cs = 100000,
                    Cl = 10000,
                    
                    kpsat=3,
                    p50 = -4,
                    psiv=-2,
                    s50 = 30,
                    gmin = 10, # mmol m-2 s-1
                    
                    psil0=-1, 
                    psist0=-0.5, 
                    
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
                    LMA=NA){
  
  if(is.null(met)){
    stop("Must provide met dataframe with VPD, Tair, PPFD, precip (optional)")
  }
  n <- nrow(met)
  
  timestep_sec <- 60*met_timestep / runs_per_timestep
  
  LAI <- AL / groundarea
  soilvolume <- groundarea * soildepth
  
  # Initial conditions
  na <- rep(NA, nrow(met))
  out <- data.frame(Eleaf=na,psil=na,psist=na,psis=na,sw=na,
                    ks=na,kp=na,Jsl=na,Jrs=na,krst=na,kstl=na)
  
  out$psil[1] <- psil0
  out$psist[1] <- psist0
  out$sw[1] <- sw0
  out$psis[1] <- psie*(sw0/thetasat)^-b
  out$Eleaf[1] <- 0
  
  # soil-to-root conductance
  out$ks[1] <- ksoil_fun(out$psis[1], Ksat, psie, b, 
                         LAI, soildepth=soildepth, Lv=Lv)
  
  pars <- list(timestep_sec=timestep_sec,psiv=psiv,sf=sf,g1=g1,Ca=Ca,Cs=Cs,
               Cl=Cl,kpsat=kpsat,p50=p50,s50=s50,gmin=gmin,
               thetasat =thetasat,AL=AL,soildepth=soildepth,
               soilvolume=soilvolume,groundarea=groundarea,
               LAI=LAI,b=b,psie=psie,Ksat=Ksat,
               Lv=Lv,keepwet=keepwet,stopsimdead=stopsimdead,
               plcdead=plcdead)
  
  for(i in 2:n){
    
    out <- calc_timestep(met, i, out, pars)
    
    if(stopsimdead){
      plc <- 100*(1 - out$kp[i]/pars$kpsat)
      if(plc > plcdead)break
    }
    
  }
  
  d <- cbind(met, out)
  d$plc <- 100 * (1 - d$kp / kpsat)
  d$Eplant <- AL * out$Eleaf
  d$t <- 1:nrow(d)
  
  # when stopsimdead, last many rows are NA
  d <- d[!is.na(d$psis),]
  
  return(d)
}




# Using previous timestep psist, psil and ks, calculate fluxes/pools for next timestep
# can only calculate timestep i when i-1 has been calculated!
calc_timestep <- function(met, i, out, pars){
  
  # Plant hydraulic conductance
  # Note how it depends on previous timestep stem water potential.
  out$kp[i] <- pars$kpsat * fsig_hydr(out$psist[i-1], pars$s50, pars$p50)    
  
  # from soil to stem pool
  out$krst[i] <- 1 / (1/out$ks[i-1] + 1/(2*out$kp[i]))
  
  # from stem pool to leaf
  out$kstl[i] <- 2*out$kp[i]
  
  # packageVersion("plantecophys") >= "1.2-6"
  p <- Photosyn(VPD=met$VPD[i], 
                gsmodel="BBdefine", 
                BBmult=(pars$g1/pars$Ca)*fsig_tuzet(out$psil[i-1], pars$psiv, pars$sf), 
                Tleaf=met$Tair[i],
                PPFD=met$PPFD[i],
                Ca=pars$Ca)
  
  # Leaf transpiration (mmol m-2 s-1)
  out$Eleaf[i] <- p$ELEAF + (met$VPD[i]/101)*pars$gmin
  
  # Xu method.
  # Can write the dynamic equation as: dPsil_dt = b + a*psil
  # Then it follows (Xu et al. 2016, Appendix, and Code).
  bp <- (pars$AL * 2 * out$kstl[i] * out$psist[i-1] - pars$AL * out$Eleaf[i])/pars$Cl
  ap <- -(pars$AL * 2 * out$kstl[i] / pars$Cl)
  out$psil[i] <- ((ap * out$psil[i-1] + bp) * exp(ap * pars$timestep_sec) - bp)/ap
  
  # Flux from stem to leaf= change in leaf storage, plus transpiration
  out$Jsl[i] <- (out$psil[i] - out$psil[i-1]) * pars$Cl / pars$timestep_sec + pars$AL * out$Eleaf[i]
  
  # Update stem water potential
  # Also from Xu et al. 2016.
  bp <- (pars$AL * 2 * out$krst[i] * out$psis[i-1] - out$Jsl[i]) / pars$Cs
  ap <- -(pars$AL * 2 * out$krst[i] / pars$Cs)
  out$psist[i] <- ((ap*out$psist[i-1] + bp) * exp(ap*pars$timestep_sec) - bp)/ap
  
  # flux from soil to stem = change in stem storage, plus Jrl
  out$Jrs[i] <- (out$psist[i] - out$psist[i-1]) * pars$Cs / pars$timestep_sec + out$Jsl[i]
  
  # Soil water increase: precip - transpiration (units kg total timestep-1)
  # (Note: transpiration is part of Jrs).
  water_in <- pars$groundarea * met$precip[i] - pars$timestep_sec * 1E-06*18 * out$Jrs[i]
  
  # soil water content (sw) in units m3 m-3
  out$sw[i] <- pmin(1, out$sw[i-1] + water_in / (pars$soilvolume * pars$thetasat * 1E03))  
  
  # for debugging / comparison.
  if(pars$keepwet){
    out$sw[i] <- out$sw[i-1]
  }
  
  # Update soil water potential and soil-to-root hydraulic conductance
  out$psis[i] <- pars$psie * (out$sw[i]/pars$thetasat)^-pars$b
  out$ks[i] <- ksoil_fun(out$psis[i], Ksat=pars$Ksat, psie=pars$psie, 
                         b=pars$b, LAI=pars$LAI, soildepth=pars$soildepth)
  
  
  return(out)
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



summarize_desica <- function(d){
  
  phase2 <- subset(d, ks < 0.05*max(kp))
  simtot <- nrow(d) / (4*24)
  sim2 <- nrow(phase2) / (4*24)
  sim1 <- simtot - sim2
  
  return(c(phase1=sim1, phase2=sim2, p50=p50, psiv=psiv, gmin=gmin, 
           capac=capac, plcfinal=d$plc[nrow(d)]))
  
}


