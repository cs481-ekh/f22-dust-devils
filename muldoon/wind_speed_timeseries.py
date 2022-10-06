import numpy as np
from numpy.random import choice
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths, boxcar
from astropy.convolution import convolve as astropy_convolve
from scipy.stats import mode
from statsmodels.robust import mad

import muldoon.utils as utils
import muldoon.met_timeseries as met
from emcee import *

class WindSpeedTimeseries(met.MetTimeseries):
    """
    A meteorological time-series tailored to wind and used for 
    looking for vortices
    """