desica2 <- function(met=NULL,
                    met_timestep = 15,
                    runs_per_timestep = 1,
                    
                    Ca = 400,  
                   sf=8,
                   g1=5,
                   
                   Cs = 30000,
                   Cl = 3000,
                   
                   kpsat=3,
                   p50 = -4,
                   psiv=-3,
                   s50 = 30,
                   gmin = 10, # mmol m-2 s-1
                   
                   psil0=-2, 
                   psist0=-1, 
                   
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

