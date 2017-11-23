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

class Desica(object):

    def __init__(self, keep_wet=False, stop_dead=True, plc_dead=88.,
                 run_twice=True, soil_depth=1.0, ground_area=1.0,
                 met_timestep=15., Ca=400., sf=8., g1=10., Cs=100000.,
                 Cl=10000., kpsat=3., p50=-4., psiv=-2., s50=30., gmin=10,
                 psil0=-1., psist0=-0.5, theta_sat=0.5,sw0=0.5, AL=2.5, b=6.,
                 psie=-0.8*1E-03, Ksat=20., Lv=10000., f=None, LMA=None):

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
        self.timestep_sec = 60. * self.met_timestep
        if self.run_twice:
            self.timestep_sec /= 2.


    def run_me(self, met=None):

        (n, out) = self.initial_model()

    def initial_model(self):
        n = len(met)

        out = self.setup_out_df()
        out.psil[0] = self.psil0
        out.psist[0] = self.psist0
        out.sw[0] = self.sw0
        out.psis[0] = self.psie * (self.sw0 / self.theta_sat)**-self.b
        out.Eleaf[0] = 0.0

        return n, out

    def setup_out_df(self):
        dummy = np.zeros(len(met))
        out = pd.DataFrame({'Eleaf':dummy, 'psil':dummy, 'psist':dummy,
                            'psis':dummy, 'sw':dummy, 'ks':dummy, 'kp':dummy,
                            'Jsl':dummy, 'Jrs':dummy, 'krst':dummy,
                            'kstl':dummy})

        return out

if __name__ == "__main__":

    met = generate_met_data(Tmin=10, RH=30, ndays=200)

    psist0 = 0.
    AL = 6.      # leaf area (m2)
    p50 = -4.    # MPa
    psiv = -3.   # MPa (for Tuzet model)
    gmin = 10.   # mmol m-2 s-1
    Cl = 10000.  # Leaf capacitance (mmol MPa-1) (total plant)
    Cs = 120000. # Stem capacitance (mmol MPa-1)

    D = Desica(psist0=psist0, AL=AL, p50=p50, psiv=psiv, gmin=gmin, Cl=Cl,
               Cs=Cs, run_twice=True, stop_dead=True)
    D.run_me(met)
