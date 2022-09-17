import pyQuickBeam
import unittest
import os
import pytest
import numpy as np


class TestDSD(unittest.TestCase):
    def setUp(self):
        # Constants
        maxhclass = 50  # Number of hm classes
        nd = 85  # Number of diameters
        self.dmin, self.dmax = 0.1, 10000.  # Min and max particle size
        D = np.exp(np.linspace(np.log(self.dmin), np.log(self.dmax), nd))
        Dwidth = np.diff(D)
        self.Dwidth = np.concatenate((Dwidth, np.array([Dwidth[-1]])))
        Dmean = 0.5 * (D[1:]+D[:-1])
        Dmean = np.concatenate((Dmean, np.array([Dmean[-1]])))
        self.Dmean = D
        self.fc = np.zeros(D.shape)
        
        self.apm = (6/np.pi) * 100
        self.bpm = 2.
        self.Q = 1.34

    def test_DSD_F05(self):
        dsdN = pyQuickBeam.dsd(self.Q,  # Mixing ratio
                               self.Dmean,  # Size dist
                               6,  # Dist type
                               1.1,  # Density of air
                               -30,  # Temp
                               self.dmin,  # dmin
                               self.dmax,  # dmax
                               -1, # Size param 1
                               -1, # Size param 2
                               -1, # Size param 3
                               self.fc,  # fc
                               0,  # scaled?
                               self.apm,  # apm
                               self.bpm)  # bpm
        number = (self.Dwidth*dsdN).sum()*1e6/1.071
        assert(np.abs(number-6746) < 5)
        mass = (self.Dwidth*dsdN*self.apm*(self.Dmean)**self.bpm).sum()
        assert(np.abs(mass-1580) < 2)

    def test_DSD_MGamma(self):
        dsdN = pyQuickBeam.dsd(self.Q,  # Mixing ratio
                               self.Dmean,  # Size dist
                               1,  # Dist type
                               1.1,  # Density of air
                               -30,  # Temp
                               self.dmin,  # dmin
                               self.dmax,  # dmax
                               250000, # Size param 1
                               -2, # Size param 2
                               -3, # Size param 3
                               self.fc,  # fc
                               0,  # scaled?
                               self.apm,  # apm
                               self.bpm)  # bpm
        number = (self.Dwidth*dsdN).sum()*1e6/1.071
        assert(np.abs(number-250180) < 5)
        mass = (self.Dwidth*dsdN*self.apm*(self.Dmean)**self.bpm).sum()
        assert(np.abs(mass-1580) < 2)

    def test_DSD_ThompsonGraupel(self):
        dsdN = pyQuickBeam.dsd(self.Q,  # Mixing ratio
                               self.Dmean,  # Size dist
                               7,  # Dist type
                               1.1,  # Density of air
                               -30,  # Temp
                               self.dmin,  # dmin
                               self.dmax,  # dmax
                               100*self.Q, # Size param 1
                               10000, # Size param 2
                               -1, # Size param 3
                               self.fc,  # fc
                               0,  # scaled?
                               self.apm,  # apm
                               self.bpm)  # bpm
        number = (self.Dwidth*dsdN).sum()*1e6/1.071
        assert(np.abs(number-7) < 5)
        
        mass = (self.Dwidth*dsdN*self.apm*(self.Dmean)**self.bpm).sum()
        assert(np.abs(mass-1580) < 2)

