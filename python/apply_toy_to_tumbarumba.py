#!/usr/bin/env python

"""
Apply the desica model to tumbarumba

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
from desica import Desica, CanopySpace
from photosynthesis import FarquharC3
import netCDF4 as nc


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
    fig.savefig("tumbarumba_plot.pdf", bbox_inches='tight', pad_inches=0.1)

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
    fig.savefig("tumbarumba_sw_swp.pdf", bbox_inches='tight', pad_inches=0.1)


def qair_to_vpd(qair, tair, press):

    PA_TO_KPA = 0.001

    # saturation vapor pressure
    es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

    # vapor pressure
    ea = (qair * press / PA_TO_KPA) / (0.622 + (1.0 - 0.622) * qair)

    vpd = (es - ea) * PA_TO_KPA
    vpd = np.where(vpd < 0.0, 0.05, vpd)

    return vpd

if __name__ == "__main__":

    time_step = 60
    deg2kelvin = 273.15
    SW_2_PAR = 2.3
    PA_TO_KPA = 0.001
    MM_S_TO_MM_HR = 3600.0

    site = "TumbaFluxnet"
    met_dir = "/Users/mdekauwe/research/CABLE_runs/met_data/fluxnet2015"
    met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (site))
    f = nc.Dataset(met_fname)
    times = f.variables['time']
    date_time = nc.num2date(times[:], times.units)

    #print(date_time[0])
    #print(date_time[1])

    met = pd.DataFrame(f.variables["Tair"][:,0,0,0] - deg2kelvin,
                       columns=['tair'])

    met['par'] = f.variables["SWdown"][:,0,0] * SW_2_PAR
    met['qair'] = f.variables['Qair'][:,0,0]
    met['press'] = f.variables['PSurf'][:,0,0] * PA_TO_KPA
    met['precip'] = f.variables['Rainf'][:,0,0] * MM_S_TO_MM_HR
    #met['precip'] = np.zeros(len(met))

    met['Ca'] = np.ones(len(met)) * 380.0 # umol mol-1
    met['vpd'] = qair_to_vpd(met.qair, met.tair, met.press)

    # adding correct datetime information
    met['dates'] = date_time
    met = met.set_index('dates')
    met['year'] = met.index.year
    met['doy'] = met.index.dayofyear

    met = met[met.index.year < 2006]

    psi_stem0 = 0.
    AL = 2.      # leaf area (m2)
    p50 = -4.    # MPa
    psi_f = -3.  # Reference potential (MPa) for Tuzet model
    gmin = 10.   # mmol m-2 s-1
    Cl = 10000.  # Leaf capacitance (mmol MPa-1) (total plant)
    Cs = 120000. # Stem capacitance (mmol MPa-1)

    F = CanopySpace()
    D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
               Cl=Cl, Cs=Cs, F=F, nruns=1, stop_dead=True, met_timestep=60.)
    out = D.main(met)

    make_plot(out, time_step)
    plot_swp_sw(out)
