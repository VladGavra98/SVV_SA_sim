# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:42:31 2020
This meant to be used for the tests of units & modules.
@author: vladg
"""

from main_sim import calcCircum, calcStArea, calcStPose, calcInertia,calcCentroid
import unittest

import numpy as np

aircraft = "A320" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
ca = 0.547  # m
la = 2.771   # m
x1 = 0.172  # m change
x2 = 1.211  # m change
x3 = 2.591  # m cahnge
xa = 0.35   # m cahnge
ha = 0.225 # m
tsk = 1.1/1000  # m
tsp = 2.9/1000  # m
tst = 1.2/1000  # m
hst = 15./1000   # m
wst = 20./1000   # m
nst = 17  # -
d1 = 0.01154  # m change
d3 = 0.01840  # m change
theta = np.radians(26)  # rad
P = 97.4*1000  # N

class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_calcCircum(self):
        assert calcCircum(0,0) == 0

    def test_calcCentroid(self):
        zc_ver = -0.21577972811362234
        assert calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst) == zc_ver

    def test_calcInertia(self):
        Iyy_ver = 6.86413733566373e-05   #from verification model
        Izz_ver = 1.28074562408502e-05

if __name__ == '__main__':
    unittest.main()