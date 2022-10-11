import numpy as np
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils 
from muldoon.read_data import *
from muldoon.read_pds4_data import *

def test_read_input_file():
    data =  readInputFile('https://pdssbn.astro.umd.edu/holdings/pds4-gbo-kpno:hyakutake_spectra-v1.0/data/offset_0_arcsec.xml')
