#!/usr/bin/env python

"""
Desica model

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.11.2017)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from generate_met_data import generate_met_data
from photosynthesis import FarquharC3

class Desica(object):

    def __init__(self, plc_dead=88.,soil_depth=1.0, ground_area=1.0,
                 met_timestep=15., Ca=400., sf=8., g1=10., Cs=100000.,
                 Cl=10000., kpsat=3., p50=-4., psiv=-2., s50=30., gmin=10,
                 psil0=-1., psist0=-0.5, theta_sat=0.5,sw0=0.5, AL=2.5, b=6.,
                 psie=-0.8*1E-03, Ksat=20., Lv=10000., F=None, keep_wet=False,
                 stop_dead=True, run_twice=True):

        self.keep_wet = keep_wet
        self.stop_dead = stop_dead
        self.plc_dead = plc_dead
        self.run_twice = run_twice
        self.soil_depth = soil_depth
        self.ground_area = ground_area
        self.soil_volume = self.ground_area * self.soil_depth
        self.met_timestep = met_timestep
        self.Ca = Ca
        self.sf = sf
        self.g1 = g1
        self.Cs = Cs
        self.Cl = Cl
        self.kpsat = kpsat
        self.p50 = p50
        self.psiv = psiv
        self.s50 = s50
        self.gmin = gmin
        self.psil0 = psil0
        self.psist0 =psist0
        self.theta_sat =theta_sat
        self.sw0 = sw0
        self.AL = AL
        self.lai = AL / self.ground_area
        self.b = b
        self.psie = psie
        self.Ksat = Ksat
        self.Lv = Lv
        self.F = F
        self.deg2kelvin = 273.15
        self.timestep_sec = 60. * self.met_timestep
        if self.run_twice:
            self.timestep_sec /= 2.


    def main(self, met=None):

        (n, out) = self.initial_model()

        for i in range(1, n):

            out = self.run_timestep(i, met, out)

            # save solutions, use as input for another run,
            # keeping everything else the same
            if self.run_twice:
                out2 = out
                out2.psil[i-1] = out.psil[i]
                out2.psist[i-1] = out.psist[i]
                out = self.run_timestep(i, met, out2)

            if self.stop_dead:
                plc = self.calc_plc(out.kp[i])
                if plc > self.plc_dead:
                    break

        #out["plc"] = self.calc_plc(out.kp)
        #out["Eplant"] = self.AL * out.Eleaf
        #out["t"] = np.arange(1, n+1)

    def initial_model(self):
        n = len(met)

        out = self.setup_out_df()
        out.psil[0] = self.psil0
        out.psist[0] = self.psist0
        out.sw[0] = self.sw0
        out.psis[0] = self.calc_swp(self.sw0)
        out.Eleaf[0] = 0.0

        # soil-to-root conductance
        out.ks[0] = self.calc_ksoil(out.psis[0])

        return n, out

    def setup_out_df(self):
        dummy = np.ones(len(met)) * np.nan
        out = pd.DataFrame({'Eleaf':dummy, 'psil':dummy, 'psist':dummy,
                            'psis':dummy, 'sw':dummy, 'ks':dummy, 'kp':dummy,
                            'Jsl':dummy, 'Jrs':dummy, 'krst':dummy,
                            'kstl':dummy})

        return out

    def run_timestep(self, i, met, out):

        # Plant hydraulic conductance
        # Note how it depends on previous timestep stem water potential.
        out.kp[i] = self.kpsat * self.fsig_hydr(out.psist[i-1])

        # from soil to stem pool
        out.krst[i] = 1.0 / (1.0 / out.ks[i-1] + 1.0 / (2.0 * out.kp[i]))

        # from stem pool to leaf
        out.kstl[i] = 2.0 * out.kp[i]

        Tleaf_K = met.tair[i] + self.deg2kelvin

        (An, Acn,
         Ajn, gsc) = self.F.calc_photosynthesis(Cs=Cs, Tleaf=Tleaf_K,
                                                Par=met.par[i], vpd=met.vpd[i],
                                                Rd25=0.92, Q10=1.92, Vcmax25=50,
                                                Jmax25=100., Eav=82620.87,
                                                deltaSv=645.1013, Eaj=39676.89,
                                                deltaSj=641.3615)

        print(An)
        sys.exit()

        # call photosynthesis ... add coupled code.
        # Use this for now ...
        Ci = 400.
        ALEAF = -0.4407474
        GS = 0.001
        ELEAF = 0.01102896
        Ac = 7.586432
        Aj = 0.0
        Ap = 3000.
        Rd = 0.4407474
        VPD = 1.102896
        Tleaf = 13.71879
        Ca = 400.
        Cc = 400.
        PPFD = 0.0
        Patm = 100.

        # Don't add gmin, instead use it as bottom value.
        gs = max(self.gmin, 1000. * GS)

        # Leaf transpiration (mmol m-2 s-1)
        out.Eleaf[i] = (met.vpd[i] / 101.0) * gs

        out.psil[i] = self.calc_xylem_water_potential(out.kstl[i],
                                                      out.psist[i-1],
                                                      out.psil[i-1],
                                                      out.Eleaf[i])

        # Flux from stem to leaf= change in leaf storage, plus transpiration
        out.Jsl[i] = self.calc_flux_to_leaf(out.psil[i], out.psil[i-1],
                                            out.Eleaf[i])

        # Update stem water potential
        out.psist[i] = self.update_stem_wp(out.krst[i], out.psis[i-1],
                                           out.Jsl[i], out.psist[i-1])

        # flux from soil to stem = change in stem storage, plus Jrl
        out.Jrs[i] = self.calc_flux_soil_to_stem(out.psist[i],
                                                 out.psist[i-1], out.Jsl[i])

        out.sw[i] = self.update_sw_balance(met.precip[i], out.Jrs[i],
                                           out.sw[i-1])

        # Update soil water potential
        out.psis[i] = self.calc_swp(out.sw[i])

        # Update soil-to-root hydraulic conductance
        out.ks[i] = self.calc_ksoil(out.psis[i])

        return out

    def calc_swp(self, sw):
        return self.psie * (sw / self.theta_sat)**-self.b

    def calc_ksoil(self, psis):
        rroot = 1E-06
        Ks = self.Ksat * (self.psie / psis)**(2. + 3. / self.b)
        if psis == 0.0:
            Ks = self.Ksat

        rcyl = 1.0 / np.sqrt(np.pi * self.Lv)
        Rl = self.Lv * self.soil_depth
        Ksoil = (Rl / self.lai) * 2. * np.pi * Ks / np.log(rcyl / rroot)

        return Ksoil

    def fsig_hydr(self, P):
        X = 50.
        SX = self.s50
        PX = self.p50

        P = np.abs(P)
        PX = np.abs(PX)
        V = (X - 100.) * np.log(1.0 - X / 100.)
        p = (P / PX)**((PX * SX) / V)
        relk = (1. - X / 100.)**p

        return (relk)

    def calc_xylem_water_potential(self, kstl, psist_prev, psil_prev, Eleaf):
        # Xu method.
        # Can write the dynamic equation as: dPsil_dt = b + a*psil
        # Then it follows (Xu et al. 2016, Appendix, and Code).
        bp = (self.AL * 2.0 * kstl * psist_prev - self.AL * Eleaf) / self.Cl
        ap = -(self.AL * 2.0 * kstl / self.Cl)
        psil = ((ap * psil_prev + bp) * np.exp(ap * self.timestep_sec) - bp)/ap

        return psil

    def calc_flux_to_leaf(self, psil, psil_prev, Eleaf):
        # Flux from stem to leaf = change in leaf storage, plus transpiration
        Jsl = (psil - psil_prev) * self.Cl / self.timestep_sec + self.AL * Eleaf
        return Jsl

    def update_stem_wp(self, krst, psis_prev, Jsl, psist_prev):
        # from Xu et al. 2016.
        bp = (self.AL * 2.0 * krst * psis_prev - Jsl) / self.Cs
        ap = -(self.AL * 2.0 * krst / self.Cs)
        psist = ((ap * psist_prev + bp) * np.exp(ap * self.timestep_sec)-bp)/ap

        return psist

    def calc_flux_soil_to_stem(self, psist, psist_prev, Jsl):
        return (psist - psist_prev) * self.Cs / self.timestep_sec + Jsl

    def update_sw_balance(self, precip, Jrs, sw_prev):

        # Soil water increase: precip - transpiration (units kg total tstep-1)
        # (Note: transpiration is part of Jrs).
        conv = 1E-06 * 18.
        water_in = self.ground_area * precip - self.timestep_sec * conv * Jrs

        # soil water content (sw) in units m3 m-3
        sw = min(1.0, sw_prev + water_in / (self.soil_volume * 1E03))

        return sw

    def calc_plc(self, kp):
        return 100.0 * (1.0 - kp / self.kpsat)


if __name__ == "__main__":

    met = generate_met_data(Tmin=10, RH=30, ndays=200)

    psist0 = 0.
    AL = 6.      # leaf area (m2)
    p50 = -4.    # MPa
    psiv = -3.   # MPa (for Tuzet model)
    gmin = 10.   # mmol m-2 s-1
    Cl = 10000.  # Leaf capacitance (mmol MPa-1) (total plant)
    Cs = 120000. # Stem capacitance (mmol MPa-1)


    Vcmax25 = 30.0
    Jmax25 = Vcmax25 * 2.0
    Rd25 = 2.0
    Eaj = 30000.0
    Eav = 60000.0
    deltaSj = 650.0
    deltaSv = 650.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 2.0
    gamma = 0.0
    g0 = 0.001
    g1 = 10.0
    theta_J = 0.85
    F = FarquharC3(peaked_Jmax=True, peaked_Vcmax=False, model_Q10=True,
                   gs_model="leuning", gamma=gamma, g0=g0,
                   g1=g1, theta_J=theta_J)

    D = Desica(psist0=psist0, AL=AL, p50=p50, psiv=psiv, gmin=gmin, Cl=Cl,
               Cs=Cs, F=F, run_twice=True, stop_dead=True)
    D.main(met)
