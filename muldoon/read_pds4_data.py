"""
A collections of routines to read in data from various missions (PDS4 format)
"""
from urllib.request import urlretrieve
from pds4_tools import pds4_read
import matplotlib.pyplot as plt
import ssl

def readInputFile(filename:str):
    ssl._create_default_https_context = ssl._create_unverified_contextx
    data = pds4_read(filename)
    return data.info()
