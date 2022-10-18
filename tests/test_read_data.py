import numpy as np
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils 
from muldoon.read_data import *

def test_file_verification_pds4_file_exist():
    data =  file_verification('https://pdssbn.astro.umd.edu/holdings/pds4-gbo-kpno:hyakutake_spectra-v1.0/data/offset_0_arcsec.xml')
    expect_file_status = 1
    assert(data == expect_file_status)

def test_file_verification_pds4_file_does_not_exist():
    data =  file_verification('https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/sol_0000_0089/sol_0001/WE__0001___________CAL_ATS_________________P1.xml')
    # wrong link
    expect_file_status = 0
    assert(data == expect_file_status)

def test_file_verification_csv_file_exist():
    data =  file_verification('./tests/WE__0001___________DER_PS__________________P02.csv')
    expect_file_status = 1
    assert(data == expect_file_status)

def test_file_verification_csv_file_does_not_exist():
    data =  file_verification('./tests/WE__0001___________DER_PS__________________P0.csv')
    # wrong csv file
    expect_file_status = 0
    assert(data == expect_file_status)