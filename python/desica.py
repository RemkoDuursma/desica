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
                 run_twice=True, soil_depth=1.0, ground_area=1.0):

        self.keep_wet = keep_wet
        self.stop_dead = stop_dead
        self.plc_dead = plc_dead
        self.run_twice = run_twice
        self.soil_depth = soil_depth
        self.ground_area = ground_area

    def run_me(self, met=None, met_timestep=15., Ca=400.,
               sf=8., g1=10., Cs=100000., Cl=10000., kpsat=3., p50=-4.,
               psiv=-2., s50=30., gmin=10, psil0=-1., psist0=-0.5, thetasat=0.5,
               sw0=0.5, AL=2.5, b=6., psie=-0.8*1E-03, Ksat=20., Lv=10000.,
               f=None, LMA=None):

        n = len(met)
        timestep_sec = 60. * met_timestep
        if self.run_twice:
            timestep_sec /= 2.

        lai = AL / self.ground_area
        soil_volume = self.ground_area * self.soil_depth

        


if __name__ == "__main__":

    met = generate_met_data(Tmin=10, RH=30, ndays=200)

    psist0 = 0.
    AL = 6.      # leaf area (m2)
    p50 = -4.    # MPa
    psiv = -3.   # MPa (for Tuzet model)
    gmin = 10.   # mmol m-2 s-1
    Cl = 10000.  # Leaf capacitance (mmol MPa-1) (total plant)
    Cs = 120000. # Stem capacitance (mmol MPa-1)

    D = Desica(run_twice=True, stop_dead=True)
    D.run_me(met, psist0=psist0, AL=AL, p50=p50, psiv=psiv, gmin=gmin, Cl=Cl,
             Cs=Cs)
