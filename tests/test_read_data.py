from unittest import result
import numpy as np
from numpy.random import normal
from muldoon import met_timeseries as met
from muldoon import utils 
from muldoon.read_data import *
from pds4_tools import pds4_read
import pandas as pd

from muldoon.read_data import read_data

from muldoon.read_data import time_conversion_to_seconds

def test_file_verification_pds4_file_exist():
    filename = 'https://pdssbn.astro.umd.edu/holdings/pds4-gbo-kpno:hyakutake_spectra-v1.0/data/offset_0_arcsec.xml'
    _, file_status =  read_data(filename)
    assert(file_status == 1)

def test_file_verification_pds4_file_does_not_exist():
    filename = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/sol_0000_0089/sol_0001/WE__0001___________CAL_ATS_________________P1.xml'
    _, file_status =  read_data(filename)
    # wrong link
    assert(file_status == 0)

def test_file_verification_csv_file_exist():
    filename = './tests/WE__0001___________DER_PS__________________P02.csv'
    _, file_status =  read_data(filename)
    assert(file_status == 1)

def test_file_verification_csv_file_does_not_exist():
    filename = './tests/WE__0001___________DER_PS__________________P0.csv'
    _, file_status =  read_data(filename)
    # wrong csv file
    assert(file_status == 0)

def test_time_conversion_to_seconds_ltst():
    times_str = ['15:28:01']
    delta_sols = [0.0]
    time_field = 'LTST'
    expected = [15.46694444]
    result = time_conversion_to_seconds(times_str, delta_sols,time_field)
    assert(round(result[0], 8) == expected[0])

def test_time_conversion_to_seconds_lmst():
    times_str = ['16:05:31.315']
    delta_sols = [0.0]
    time_field = 'LMST'
    expected = [16.09203194]
    result = time_conversion_to_seconds(times_str, delta_sols,time_field)
    assert(round(result[0], 8) == expected[0])

def test_make_seconds_since_midnight_ltst():
    filename = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/sol_0000_0089/sol_0001/WE__0001___________DER_PS__________________P02.xml'
    result = make_seconds_since_midnight(filename,'ltst', None)
    expected = [15.46694444]
    assert(round(result[0], 8) == expected[0])


def test_make_seconds_since_midnight_lmst():
    filename = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/sol_0000_0089/sol_0001/WE__0001___________DER_PS__________________P02.xml'
    result = make_seconds_since_midnight(filename,'lmst', None)
    expected = [16.09203194]
    assert(round(result[0], 8) == expected[0])

def test_read_Perseverance_PS_data_csv():
    filename = './tests/WE__0001___________DER_PS__________________P02.csv'
    _,result = read_Perseverance_PS_data(filename, None, 'LTST')
    expected = [715.96]
    assert(round(result[0], 2) == expected[0])

def test_read_Perseverance_PS_data_pds4():
    filename = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/sol_0000_0089/sol_0001/WE__0001___________DER_PS__________________P02.xml'
    _,result = read_Perseverance_PS_data(filename, None, 'LTST')
    expected = [715.96]
    assert(round(result[0], 2) == expected[0])