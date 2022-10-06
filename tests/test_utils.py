import pytest
import numpy as np
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils as util
from muldoon.read_data import *

#TODO: Needs missing values for this test to work
def test_vortex_wind():
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
        assert True
    else:
        assert False

#TODO: Needs missing values for this test to work
def test_wind_profile():
    """
    Tests the vortex wind profile function
    """
    time = 22.1
    max_wind_speed = 23
    u_1 = 3.5
    slope = 1.
    t0 = 0.
    Gamma = 0.01 
    expected = 0.0000910710

    v_t = int(util.wind_vortex_profile(max_wind_speed, u_1, slope, time, t0,Gamma))

    if(v_t == expected):
        assert True 

if __name__ == '__main__':
    pytest.main()
