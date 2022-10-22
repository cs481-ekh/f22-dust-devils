import pytest
from muldoon import met_timeseries as met
from muldoon import utils as util
from muldoon.read_data import *

def test_vortex_wind():
    """
    Tests the vortex wind function
    """
    v_act = 25
    r = .75
    D_act = 0.01
    expected = 0.3333185

    # rounded to 7 digits after the decimal
    v_r = round(util.vortex_wind(v_act, D_act, r), 7)

    if (v_r == expected):
        assert True
    if (v_r != expected):
        assert False

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

    # rounded to 9 digits after the decimal
    v_t = round(util.wind_vortex_profile(max_wind_speed, u_1, slope, time, t0, Gamma), 9)

    if(v_t == expected):
        assert True 
    if(v_t != expected):
        assert False

if __name__ == '__main__':
    pytest.main()
