import pytest
import numpy as np
import pandas as pd
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils 
from muldoon.read_data import *

def test_init():
    time = np.linspace(-1, 1, 1000)
    baseline = 0.
    slope = 1.
    t0 = 0.
    DeltaP = 1.
    Gamma = 0.01
    right_answer = np.array([baseline, slope, t0, DeltaP, Gamma])
    profile = utils.modified_lorentzian(time, baseline, slope, t0, DeltaP, Gamma) +\
        normal(scale=slope/20., size=len(time))
   
    pts = met.PressureTimeseries(time, profile)
    assert pts is not None
