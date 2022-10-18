import numpy as np
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils 
from muldoon.read_data import *
from pds4_tools import pds4_read
import pandas as pd

from muldoon.read_data import file_verification

def test_file_verification_pds4_file_exist():
    filename = 'https://pdssbn.astro.umd.edu/holdings/pds4-gbo-kpno:hyakutake_spectra-v1.0/data/offset_0_arcsec.xml'
    file_status, _ =  file_verification(filename)
    assert(file_status == 1)

def test_file_verification_pds4_file_does_not_exist():
    filename = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/sol_0000_0089/sol_0001/WE__0001___________CAL_ATS_________________P1.xml'
    file_status, _ =  file_verification(filename)
    # wrong link
    assert(file_status == 0)

def test_file_verification_csv_file_exist():
    filename = './tests/WE__0001___________DER_PS__________________P02.csv'
    file_status, _ =  file_verification(filename)
    assert(file_status == 1)

def test_file_verification_csv_file_does_not_exist():
    filename = './tests/WE__0001___________DER_PS__________________P0.csv'
    file_status, _ =  file_verification(filename)
    # wrong csv file
    assert(file_status == 0)