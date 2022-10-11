"""
A collections of routines to read in data from various missions (PDS4 format)
"""
from urllib.request import urlretrieve
from pds4_tools import pds4_read
import matplotlib.pyplot as plt
import platform
import ssl

def readInputFile(filename:str):
    """
    Read in PDS4 file
    Args:
        filename (str): path to PDS4 xml/link
    Returns:
        Data information contained in PDS4 label
    """
    if(platform.system() == 'Windows'):
        print(platform.system())
        ssl._create_default_https_context = ssl._create_unverified_contextx
    data = pds4_read(filename)
    return data.info()
