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
from canopy import Canopy, FarquharC3

class Desica(object):

    def __init__(self, plc_dead=88.,soil_depth=1.0, ground_area=1.0,
                 met_timestep=30., sf=8., g1=4., Cs=100000., b=6.,
                 Cl=10000., kp_sat=3., p50=-4., psi_f=-2., s50=30., gmin=10,
                 psi_leaf0=-1., psi_stem0=-0.5, theta_sat=0.5, sw0=0.5, AL=2.5,
                 psi_e=-0.8*1E-03, Ksat=20., Lv=10000., F=None, keep_wet=False,
                 stop_dead=True, nruns=1):

        self.keep_wet = keep_wet
        self.stop_dead = stop_dead
        self.plc_dead = plc_dead
        self.nruns = nruns
        self.soil_depth = soil_depth # depth of soil bucket, m
        self.ground_area = ground_area # m
        self.soil_volume = self.ground_area * self.soil_depth # m3
        self.met_timestep = met_timestep
        self.sf = sf # sensitivity parameter, MPa-1
        self.g1 = g1 # sensitivity of stomatal conductance to the assimilation rate, kPa
        self.Cs = Cs # stem capacitance, mmol MPa-1
        self.Cl = Cl # leaf capacitance, mmol MPa-1 (total plant)
        self.kp_sat = kp_sat # plant saturated hydraulic conductance (mmol m-2 s-1 MPa-1)
        self.p50 = p50 # xylem pressure inducing 50% loss of hydraulic conductivity due to embolism, MPa
        self.psi_f = psi_f # reference potential for Tuzet model, MPa
        self.s50 = s50
        self.gmin = gmin # minimum stomatal conductance, mmol m-2 s-1
        self.psi_leaf0 = psi_leaf0 # initial leaf water potential, MPa
        self.psi_stem0 = psi_stem0 # initial stem water potential, MPa
        self.theta_sat = theta_sat # soil water capacity at saturation (m3 m-3)
        self.sw0 = sw0 # initial soil volumetric water content (m3 m-3)
        self.AL = AL # plant leaf area, m2
        self.lai = AL / self.ground_area # leaf area index, m2 m-2
        self.b = b # empirical coefficient related to the clay content of the soil (Cosby et al. 1984).
        self.psi_e = psi_e # air entry point water potential (MPa)
        self.Ksat = Ksat
        self.Lv = Lv
        self.F = F
        self.timestep_sec = 60. * self.met_timestep / self.nruns


    def run_simulation(self, met=None):
        """
        Main wrapper to control everything

        Parameters:
        -----------
        met : object
            met forcing variables: day; Ca; par; precip; press; tair; vpd

        Returns:
        -------
        out : object
            output dataframe containing calculations for each timestep

        """
        (n, out) = self.initialise_model(met)

        for i in range(1, n):

            out = self.run_timestep(i, met, out)

            # save solutions, use as input for another run,
            # keeping everything else the same
            for j in range(1, self.nruns):
                out_temp = out
                out_temp.psi_leaf[i-1] = out.psi_leaf[i]
                out_temp.psi_stem[i-1] = out.psi_stem[i]
                out = self.run_timestep(i, met, out_temp)

            # Stop the simulation if we've died, i.e. reached P88
            if self.stop_dead:
                plc = self.calc_plc(out.kp[i])
                if plc > self.plc_dead:
                    break

        out["plc"] = self.calc_plc(out.kp)
        # mmol s-1
        out["Eplant"] = self.AL * out.Eleaf
        out["t"] = np.arange(1, n+1)

        return (out)

    def initialise_model(self, met):
        """
        Set everything up: set initial values, build an output dataframe to save
        things

        Parameters:
        -----------
        met : object
            met forcing variables: day; Ca; par; precip; press; tair; vpd

        Returns:
        -------
        n : int
            number of timesteps in the met file
        out : object
            output dataframe to store things as we go along

        """
        n = len(met)

        out = self.setup_out_df(met)
        out.psi_leaf[0] = self.psi_leaf0
        out.psi_stem[0] = self.psi_stem0
        out.sw[0] = self.sw0
        out.psi_soil[0] = self.calc_swp(self.sw0)
        out.Eleaf[0] = 0.0

        # soil-to-root hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ks[0] = self.calc_ksoil(out.psi_soil[0])

        return n, out

    def setup_out_df(self, met):
        """
        Create and output dataframe to save things

        Parameters:
        -----------
        met : object
            met forcing variables: day; Ca; par; precip; press; tair; vpd

        Returns:
        -------
        out : object
            output dataframe to store things as we go along.
        """
        dummy = np.ones(len(met)) * np.nan
        out = pd.DataFrame({'Eleaf':dummy, 'psi_leaf':dummy, 'psi_stem':dummy,
                            'psi_soil':dummy, 'sw':dummy, 'ks':dummy,
                            'kp':dummy, 'Jsl':dummy, 'Jrs':dummy, 'krst':dummy,
                            'kstl':dummy})

        return out

    def run_timestep(self, i, met, out):

        self.calc_conductances(out, i)

        mult = (self.g1 / met.Ca[i]) * self.fsig_tuzet(out.psi_leaf[i-1])

        # Calculate photosynthesis and stomatal conductance
        gsw = self.F.canopy(met.Ca[i], met.tair[i], met.par[i],
                            met.vpd[i], mult)

        # Don't add gmin, instead use it as bottom value.
        gsw = max(self.gmin, 1000. * gsw)

        # Leaf transpiration assuming perfect coupling (mmol m-2 s-1)
        out.Eleaf[i] = gsw * (met.vpd[i] / met.press[i])

        out.psi_leaf[i] = self.calc_leaf_water_potential(out.kstl[i],
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

        out.sw[i] = self.update_sw_bucket(met.precip[i], out.Jrs[i],
                                          out.sw[i-1])

        # Update soil water potential
        out.psi_soil[i] = self.calc_swp(out.sw[i])

        # Update soil-to-root hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ks[i] = self.calc_ksoil(out.psi_soil[i])

        return out

    def calc_conductances(self, out, i):
        """
        Update the simple bucket soil water balance

        Parameters:
        -----------
        out : object
            output dataframe to access previous calculations and update new
            states
        i : int
            current index

        """

        # Plant hydraulic conductance (mmol m-2 s-1 MPa-1). NB. depends on stem
        # water potential from the previous timestep.
        out.kp[i] = self.kp_sat * self.fsig_hydr(out.psi_stem[i-1])

        # Conductance from soil to stem water store (mmol m-2 s-1 MPa-1)
        out.krst[i] = 1.0 / (1.0 / out.ks[i-1] + 1.0 / (2.0 * out.kp[i]))

        # Conductance from stem water store to leaf (mmol m-2 s-1 MPa-1)
        out.kstl[i] = 2.0 * out.kp[i]

    def calc_swp(self, sw):
        """
        Calculate the soil water potential (MPa). The params The parameters b
        and psi_e are estimated from a typical soil moisture release function.

        Parameters:
        -----------
        sw : object
            volumetric soil water content (m3 m-3)

        Returns:
        -----------
        psi_swp : float
            soil water potential, MPa
        """
        return self.psi_e * (sw / self.theta_sat)**-self.b

    def calc_ksoil(self, psi_soil):
        rroot = 1E-06
        Ks = self.Ksat * (self.psi_e / psi_soil)**(2. + 3. / self.b)
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

    def calc_leaf_water_potential(self, kstl, psi_stem_prev, psi_leaf_prev,
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

    def update_sw_bucket(self, precip, water_loss, sw_prev):
        """
        Update the simple bucket soil water balance

        Parameters:
        -----------
        precip : float
            precipitation (kg m-2 s-1)
        water_loss : float
            flux of water out of the soil (transpiration (kg m-2 timestep-1))
        sw_prev : float
            volumetric soil water from the previous timestep (m3 m-3)
        soil_volume : float
            volume soil water bucket (m3)

        Returns:
        -------
        sw : float
            new volumetric soil water (m3 m-3)

        """
        M_2_MM = 1E03
        delta_sw = precip - (water_loss * 1E-06 * 18.0 * self.timestep_sec)
        sw = min(self.theta_sat, \
                 sw_prev + delta_sw / (self.soil_volume * M_2_MM))

        return sw

    def calc_plc(self, kp):
        """
        Calculates the percent loss of conductivity, PLC (-)

        Parameters:
        -----------
        kp : float
            plant hydraulic conductance (mmol m-2 s-1 MPa-1)

        Returns:
        -------
        plc : float
            percent loss of conductivity (-)

        """
        return 100.0 * (1.0 - kp / self.kp_sat)

    def fsig_tuzet(self, psi_leaf):
        """
        An empirical logistic function to describe the sensitivity of stomata
        to leaf water potential.

        Function assumes that stomata are insensitive
        to LWP at values close to zero and that stomata rapidly close with
        decreasing LWP.

        Parameters:
        -----------
        psi_leaf : float
            leaf water potential (MPa)

        Returns:
        -------
        fw : float
            sensitivity of stomata to leaf water potential

        Reference:
        ----------
        * Tuzet et al. (2003) A coupled model of stomatal conductance,
          photosynthesis and transpiration. Plant, Cell and Environment 26,
          1097–1116

        """
        num = 1.0 + np.exp(self.sf * self.psi_f)
        den = 1.0 + np.exp(self.sf * (self.psi_f - psi_leaf))
        fw = num / den

        return fw



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

def plot_swp_sw(out):

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

    ax1.plot(out.sw, out.psi_soil, "b.", label="Soil")

    ax1.set_xlabel("SW (m3 m-3)")
    ax1.set_ylabel("Soil Water potential (MPa)")
    #ax1.legend(numpoints=1, loc="best")
    fig.savefig("sw_swp.pdf", bbox_inches='tight', pad_inches=0.1)



if __name__ == "__main__":

    time_step = 30

    met = generate_met_data(Tmin=10, RH=30, ndays=200, time_step=time_step)

    psi_stem0 = 0. # initial stem water potential, MPa
    AL = 6.        # plant leaf area, m2
    p50 = -4.      # xylem pressure inducing 50% loss of hydraulic conductivity due to embolism, MPa
    psi_f = -3.    # reference potential for Tuzet model, MPa
    gmin = 10.     # minimum stomatal conductance, mmol m-2 s-1
    Cl = 10000.    # leaf capacitance, mmol MPa-1 (total plant)
    Cs = 120000.   # stem capacitance, mmol MPa-1
    g1 = 4.0       # sensitivity of stomatal conductance to the assimilation rate, kPa

    F = Canopy(g1=g1)
    D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
               Cl=Cl, Cs=Cs, F=F, g1=g1, nruns=3, stop_dead=True)
    out = D.run_simulation(met)

    make_plot(out, time_step)
    plot_swp_sw(out)
