import pytest
import numpy as np
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils as util
from muldoon.read_data import *

#TODO: Needs missing values for this test to work
def test_vortex_wind(self):
    """
    Tests the vortex wind function
    """
    time = np.linspace(-1, 1, 1000)
    v_act = None
    r = None
    D_act = None
    expected = None

    v_r = util.vortex_wind(v_act, D_act, r)

    if (v_r == expected):
        self.assertTrue(True)

#TODO: Needs missing values for this test to work
def test_wind_profile(self):
    """
    Tests the vortex wind profile function
    """
    time = np.linspace(-1, 1, 1000)
    max_wind_speed = None
    U_1 = None
    U_2 = None
    slope = 1.
    t0 = 0.
    Gamma = 0.01 
    expected = None

    v_t = util.wind_vortex_profile(max_wind_speed, U_1, U_2, slope, time, t0,Gamma)

    if(v_t == expected):
        self.assertTrue(True) 

if __name__ == '__main__':
    pytest.main()
