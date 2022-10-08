from ctypes import sizeof
from itertools import count
from operator import countOf
import numpy as np
import matplotlib.pyplot as plt
import unittest
import pytest

from numpy.random import normal

from scipy.stats import median_abs_deviation as mad

from muldoon.met_timeseries import MetTimeseries, PressureTimeseries, TemperatureTimeseries
from muldoon.utils import modified_lorentzian, fit_vortex, write_out_plot_data
from muldoon.read_data import *

# Temperature data
filename=r"docs\examples\WE__0089___________CAL_ATS_________________P01.CSV"
data = pd.read_csv(filename)
list = ['SCLK', 'LMST', 'LTST', 'ATS_LOCAL_TEMP1', 'ATS_LOCAL_TEMP2',
       'ATS_LOCAL_TEMP3', 'ATS_LOCAL_TEMP4', 'ATS_LOCAL_TEMP5']


# Unit testing for Temperature Time Series
class TestTTS(unittest.TestCase):
    def setUp(self):
        self.data = data
    
    def testHasHeaders(self):
        for i in list:
            self.assertIn(i,data.columns), "Please check all header values are in the CSV"

    def testNullValues(self):
        for i in data.columns:
            self.assertEqual(data[i].isnull().sum(), 0), "Please check for null values in CSV"


        


if __name__ == '__main__':
    unittest.main()


