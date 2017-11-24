#!/usr/bin/env python

"""
Desica model: simple plant hydraulics model with mortality.

Stem water potential follows Xu et al.

Reference:
==========
* Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K.
  (2016), Diversity in plant hydraulic traits explains seasonal and
  inter-annual variations of vegetation dynamics in seasonally dry
  tropical forests. New Phytol, 212: 80–95. doi:10.1111/nph.14009.

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
                 met_timestep=30., sf=8., g1=10., Cs=100000., b=6.,
                 Cl=10000., kp_sat=3., p50=-4., psiv=-2., s50=30., gmin=10,
                 psi_leaf0=-1., psi_stem0=-0.5, theta_sat=0.5,sw0=0.5, AL=2.5,
                 psie=-0.8*1E-03, Ksat=20., Lv=10000., F=None, keep_wet=False,
                 stop_dead=True, nruns=1):

        self.keep_wet = keep_wet
        self.stop_dead = stop_dead
        self.plc_dead = plc_dead
        self.nruns = nruns
        self.soil_depth = soil_depth
        self.ground_area = ground_area
        self.soil_volume = self.ground_area * self.soil_depth
        self.met_timestep = met_timestep
        self.sf = sf
        self.g1 = g1
        self.Cs = Cs
        self.Cl = Cl
        self.kp_sat = kp_sat
        self.p50 = p50
        self.psiv = psiv
        self.s50 = s50
        self.gmin = gmin
        self.psi_leaf0 = psi_leaf0 # initial leaf water potential (MPa)
        self.psi_stem0 = psi_stem0 # initial stem water potential (MPa)
        self.theta_sat = theta_sat
        self.sw0 = sw0 # initial soil volumetric water content
        self.AL = AL
        self.lai = AL / self.ground_area
        self.b = b
        self.psie = psie
        self.Ksat = Ksat
        self.Lv = Lv
        self.F = F
        self.timestep_sec = 60. * self.met_timestep / self.nruns

    def main(self, met=None):

        (n, out) = self.initial_model()

        for i in range(1, n):

            out = self.run_timestep(i, met, out)

            # save solutions, use as input for another run,
            # keeping everything else the same
            #if self.run_twice:
            #    out2 = out
            #    out2.psi_leaf[i-1] = out.psi_leaf[i]
            #    out2.psi_stem[i-1] = out.psi_stem[i]
            #    out = self.run_timestep(i, met, out2)

            # save solutions, use as input for another run,
            # keeping everything else the same
            for j in range(1, self.nruns):
                out_temp = out
                out_temp.psi_leaf[i-1] = out.psi_leaf[i]
                out_temp.psi_stem[i-1] = out.psi_stem[i]
                out = self.run_timestep(i, met, out_temp)

            if self.stop_dead:
                # percent loss conductivity (%)
                plc = self.calc_plc(out.kp[i])
                if plc > self.plc_dead:
                    break

        out["plc"] = self.calc_plc(out.kp)
        # mmol s-1
        out["Eplant"] = self.AL * out.Eleaf
        out["t"] = np.arange(1, n+1)

        return (out)

    def initial_model(self):
        n = len(met)

        out = self.setup_out_df()
        out.psi_leaf[0] = self.psi_leaf0
        out.psi_stem[0] = self.psi_stem0
        out.sw[0] = self.sw0
        out.psi_soil[0] = self.calc_swp(self.sw0)
        out.Eleaf[0] = 0.0

        # soil-to-root hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ks[0] = self.calc_ksoil(out.psi_soil[0])

        return n, out

    def setup_out_df(self):
        dummy = np.ones(len(met)) * np.nan
        out = pd.DataFrame({'Eleaf':dummy, 'psi_leaf':dummy, 'psi_stem':dummy,
                            'psi_soil':dummy, 'sw':dummy, 'ks':dummy,
                            'kp':dummy, 'Jsl':dummy, 'Jrs':dummy, 'krst':dummy,
                            'kstl':dummy})

        return out

    def run_timestep(self, i, met, out):

        # Plant hydraulic conductance (mmol m-2 s-1 MPa-1). NB. depends on stem
        # water potential from the previous timestep.
        out.kp[i] = self.kp_sat * self.fsig_hydr(out.psi_stem[i-1])

        # Conductance from soil to stem water store (mmol m-2 s-1 MPa-1)
        out.krst[i] = 1.0 / (1.0 / out.ks[i-1] + 1.0 / (2.0 * out.kp[i]))

        # Conductance from stem water store to leaf (mmol m-2 s-1 MPa-1)
        out.kstl[i] = 2.0 * out.kp[i]

        mult = (self.g1 / met.Ca[i]) * self.fsig_tuzet(out.psi_leaf[i-1],
                                                       self.psiv, self.sf)

        # Calculate photosynthesis and stomatal conductance
        gsw = self.F.canopy(met.Ca[i], met.tair[i], met.par[i],
                            met.vpd[i], mult)

        # Don't add gmin, instead use it as bottom value.
        gsw = max(self.gmin, 1000. * gsw)

        # Leaf transpiration assuming perfect coupling (mmol m-2 s-1)
        out.Eleaf[i] = gsw * (met.vpd[i] / met.press[i])

        out.psi_leaf[i] = self.calc_xylem_water_potential(out.kstl[i],
                                                          out.psi_stem[i-1],
                                                          out.psi_leaf[i-1],
                                                          out.Eleaf[i])

        # Flux from stem to leaf (mmol s-1) = change in leaf storage,
        # plus transpiration
        out.Jsl[i] = self.calc_flux_to_leaf(out.psi_leaf[i], out.psi_leaf[i-1],
                                            out.Eleaf[i])

        # Update stem water potential
        out.psi_stem[i] = self.update_stem_wp(out.krst[i], out.psi_soil[i-1],
                                              out.Jsl[i], out.psi_stem[i-1])

        # flux from soil to stem, i.e. root water uptake (mmol s-1) = change in
        # stem storage, plus Jsl
        out.Jrs[i] = self.calc_flux_soil_to_stem(out.psi_stem[i],
                                                 out.psi_stem[i-1], out.Jsl[i])

        out.sw[i] = self.update_sw_balance(met.precip[i], out.Jrs[i],
                                           out.sw[i-1])

        # Update soil water potential
        out.psi_soil[i] = self.calc_swp(out.sw[i])

        # Update soil-to-root hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ks[i] = self.calc_ksoil(out.psi_soil[i])

        return out

    def calc_swp(self, sw):
        return self.psie * (sw / self.theta_sat)**-self.b

    def calc_ksoil(self, psi_soil):
        rroot = 1E-06
        Ks = self.Ksat * (self.psie / psi_soil)**(2. + 3. / self.b)
        if psi_soil == 0.0:
            Ks = self.Ksat

        rcyl = 1.0 / np.sqrt(np.pi * self.Lv)
        Rl = self.Lv * self.soil_depth
        Ksoil = (Rl / self.lai) * 2. * np.pi * Ks / np.log(rcyl / rroot)

        return Ksoil

    def fsig_hydr(self, P):
        X = 50.
        P = np.abs(P)
        PX = np.abs(self.p50)
        V = (X - 100.) * np.log(1.0 - X / 100.)
        p = (P / PX)**((PX * self.s50) / V)
        relk = (1. - X / 100.)**p

        return (relk)

    def calc_xylem_water_potential(self, kstl, psi_stem_prev, psi_leaf_prev,
                                   Eleaf):
        # Following Xu et al, see Appendix + code
        #
        # Reference:
        # ==========
        # * Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K.
        #   (2016), Diversity in plant hydraulic traits explains seasonal and
        #    inter-annual variations of vegetation dynamics in seasonally dry
        #    tropical forests. New Phytol, 212: 80–95. doi:10.1111/nph.14009.
        #
        # Can write the dynamic equation as: dpsi_leaf_dt = b + a*psi_leaf
        # Then it follows (Xu et al. 2016, Appendix, and Code).
        bp = (self.AL * 2.0 * kstl * psi_stem_prev - self.AL * Eleaf) / self.Cl
        ap = -(self.AL * 2.0 * kstl / self.Cl)
        psi_leaf = ((ap * psi_leaf_prev + bp) * \
                    np.exp(ap * self.timestep_sec) - bp) / ap

        return psi_leaf

    def calc_flux_to_leaf(self, psi_leaf, psi_leaf_prev, Eleaf):
        # Flux from stem to leaf = change in leaf storage, plus transpiration
        Jsl = (psi_leaf - psi_leaf_prev) * \
                self.Cl / self.timestep_sec + self.AL * Eleaf
        return Jsl

    def update_stem_wp(self, krst, psi_soil_prev, Jsl, psi_stem_prev):
        # Following Xu et al, see Appendix + code
        #
        # Reference:
        # ==========
        # * Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K.
        #   (2016), Diversity in plant hydraulic traits explains seasonal and
        #    inter-annual variations of vegetation dynamics in seasonally dry
        #    tropical forests. New Phytol, 212: 80–95. doi:10.1111/nph.14009.
        #
        bp = (self.AL * 2.0 * krst * psi_soil_prev - Jsl) / self.Cs
        ap = -(self.AL * 2.0 * krst / self.Cs)
        psi_stem = ((ap * psi_stem_prev + bp) * \
                    np.exp(ap * self.timestep_sec)-bp) / ap

        return psi_stem

    def calc_flux_soil_to_stem(self, psi_stem, psi_stem_prev, Jsl):
        return (psi_stem - psi_stem_prev) * self.Cs / self.timestep_sec + Jsl

    def update_sw_balance(self, precip, Jrs, sw_prev):

        # Soil water increase: precip - transpiration (units kg total tstep-1)
        # (Note: transpiration is part of Jrs).
        conv = 1E-06 * 18.
        water_in = self.ground_area * precip - self.timestep_sec * conv * Jrs

        # soil water content (sw) in units m3 m-3
        sw = min(1.0, sw_prev + water_in / (self.soil_volume * 1E03))
        return sw

    def calc_plc(self, kp):
        return 100.0 * (1.0 - kp / self.kp_sat)

    def fsig_tuzet(self, psi_leaf, psiv, sf):
        return (1.0 + np.exp(sf * psiv)) / (1.0 + np.exp(sf * (psiv - psi_leaf)))


class CanopySpace(object):

    def __init__(self, g0=0.001, gamma=0.0, g1=4.0, theta_J=0.85, Rd25=0.92,
                 Q10=1.92, Vcmax25=50, Jmax25=100., Eav=82620.87,
                 deltaSv=645.1013, Eaj=39676.89, deltaSj=641.3615):

        self.deg2kelvin = 273.15
        self.F = FarquharC3(peaked_Jmax=True, peaked_Vcmax=False,
                            model_Q10=True, gs_model="user_defined",
                            gamma=gamma, g0=g0, g1=g1, theta_J=theta_J)

        self.Rd25 = Rd25
        self.Q10 = Q10
        self.Vcmax25 = Vcmax25
        self.Jmax25 = Jmax25
        self.Eav = Eav
        self.Eaj = Eaj
        self.deltaSv = deltaSv
        self.deltaSj = deltaSj

    def canopy(self, Cs, tair, par, vpd, mult):

        tleaf_K = tair + self.deg2kelvin

        (An, gsc, gsw) = self.F.calc_photosynthesis(Cs=Cs, Tleaf=tleaf_K,
                                                    Par=par, vpd=vpd,
                                                    Rd25=self.Rd25,
                                                    Q10=self.Q10,
                                                    Vcmax25=self.Vcmax25,
                                                    Jmax25=self.Jmax25,
                                                    Eav=self.Eav,
                                                    deltaSv=self.deltaSv,
                                                    Eaj=self.Eaj,
                                                    deltaSj=self.deltaSj,
                                                    mult=mult)

        return (gsw)

def make_plot(out, timestep=15):

    if timestep == 15:
        ndays = out.t / 96
    elif timestep == 30:
        ndays = out.t / 96 * 2
    elif timestep == 60:
        ndays = out.t / 96 * 4

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ln1 = ax1.plot(ndays, out.psi_leaf, "k-", label="Leaf")
    ln2 = ax1.plot(ndays, out.psi_stem, "r-", label="Stem")
    ln3 = ax1.plot(ndays, out.psi_soil, "b-", label="Soil")

    ln4 = ax2.plot(ndays, out.plc, ls='-', color="darkgrey",
                   label="PLC")

    # added these three lines
    lns = ln1 + ln2 + ln3 + ln4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=(0.5,0.05), ncol=2)

    ax2.set_ylabel(r'PLC (%)')

    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Water potential (MPa)")
    #ax1.legend(numpoints=1, loc="best")
    fig.savefig("test_plot.pdf", bbox_inches='tight', pad_inches=0.1)



if __name__ == "__main__":

    time_step = 30

    met = generate_met_data(Tmin=10, RH=30, ndays=200, time_step=time_step)

    psi_stem0 = 0.
    AL = 6.      # leaf area (m2)
    p50 = -4.    # MPa
    psiv = -3.   # MPa (for Tuzet model)
    gmin = 10.   # mmol m-2 s-1
    Cl = 10000.  # Leaf capacitance (mmol MPa-1) (total plant)
    Cs = 120000. # Stem capacitance (mmol MPa-1)

    F = CanopySpace()

    D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psiv=psiv, gmin=gmin, Cl=Cl,
               Cs=Cs, F=F, nruns=3, stop_dead=True)
    out = D.main(met)

    make_plot(out, time_step)
