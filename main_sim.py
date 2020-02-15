# -*- coding: utf-8 -*-
"""
SVV 2020- Structural Analysis Assignment


@author: vladg
"""

import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt

aircraft = "A320"

class Aileron(name):
    """The class containing the geomtric parameters of the aileron and the material properties.
    The name is a string that contains the name of the aircraft: A320, B737"""
    def __init__(self,name):
        self.la = 2.71  #m
        self.ca = 0.5
        self.ha = 0.22
        self.ts = 0.0011

