#!/usr/bin/env python

"""
Generate some fake met data for testing. Replication of Remko's R module

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

def generate_met_data(PPFDmax=2000, RH=30, Tmax=30, Tmin=10, day_length=12,
                      sunrise=8, time_step=15, lag=0.5, ndays=40):

    r = np.arange(0, (24 * 60 - time_step) + time_step, time_step)
    p = np.zeros(len(r))
    ta = np.zeros(len(r))
    ii = np.argwhere((r > (sunrise*60)) & (r <= (sunrise*60 + day_length*60)))

    relt = (r[ii] - sunrise * 60) / float(day_length * 60)
    p[ii] = generate_par(relt, PPFDmax)
    timeh = np.arange(0, (day_length - time_step / 60.) + 0.25, time_step / 60.)
    td = partondiurnal(timeh, Tmax, Tmin, day_length, lag)
    for i,ival in enumerate(ii):
        ta[ival-1] = td[i]

    nnight = len(r) - len(td)
    tdec = ta[max(ii)-1] + (Tmin - ta[max(ii)-1]) * \
            np.arange(1, nnight+1) / nnight
    after = np.arange((max(ii)+1), len(r))

    ta[after] = tdec[np.arange(len(after))+1]
    before = np.arange(min(ii)-1)
    ta[before] = tdec[np.arange(len(after)+1, len(tdec))]

    vpd = rh_to_vpd(RH, ta)

    day = np.repeat(np.arange(1, ndays+1), 96)
    new_size = int(len(day)/96)
    vpd = np.tile(vpd, new_size)
    tair = np.tile(ta, new_size)
    par = np.tile(p, new_size)
    precip = np.zeros(len(par))

    met = pd.DataFrame({'day':day, 'par':par, 'tair':tair,
                        'vpd':vpd, 'precip':precip})

    return met

def rh_to_vpd(rh, tair, pressure=101.0):
    esatval = esat(tair, pressure)
    e = (rh / 100.0) * esatval
    vpd = (esatval - e) / 1000.0
    return vpd

def esat(tair, pressure=101.0):
    a = 611.21
    b = 17.502
    c = 240.97
    f = 1.0007 + 3.46 * 10**-8 * pressure * 1000
    esatval = f * a * (np.exp(b * tair / (c + tair)))
    return esatval

def generate_par(relt, ymax):
    return (ymax / 2.0) * (1.0 + np.sin((2.0 * np.pi * relt) + 1.5 * np.pi))

def partondiurnal(timeh, Tmax, Tmin, daylen, lag):
    return (Tmax - Tmin) * np.sin((np.pi * timeh) / (daylen + 2 * lag)) + Tmin

if __name__ == "__main__":

    met = generate_met_data(Tmin=10, RH=30, ndays=200)