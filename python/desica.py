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

    def __init__(self, met):
        self.met = met
        

    def main(self):

        pass

if __name__ == "__main__":

    met = generate_met_data(Tmin=10, RH=30, ndays=200)
    D = Desica(met)
    D.main()
