"""
A collections of routines to read in data from various missions
"""
from urllib.error import URLError
import numpy as np
import pandas as pd
from urllib.request import urlretrieve
from pds4_tools import pds4_read
import matplotlib.pyplot as plt

def read_data(filename:str):
    """
    Verifies the existence of provide file

    Args:
        filename (str): path to CSV file/path to PDS4 file or link
        
    Returns:
        file_status: 1 if file exists, else expception is thrown
    """
    file_status = 0
    data = " "
    error_message = "\n===>Filenotfounderror: There is no file " + filename + " in the working directory. Please check file name or path\n"
    if filename.endswith('.xml'):
        try:
            data = pds4_read(filename)
            data.info()
            file_status = 1
        except Exception as e:
            print(error_message)
            file_status = 0
    else:
        print("Processing file: " + filename)
        try:
            data = pd.read_csv(filename)
            file_status = 1
        except Exception as e:
            print(error_message)
            file_status = 0

    return data, file_status
    

def read_Perseverance_PS_data(filename, sol=None, time_field='LTST'):
    """
    Read in Perseverance MEDA PS data - https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/

    Args:
        filename (str): path to CSV file

    Returns:
        time, pressure (float array): times and pressures, times in seconds
        since midnight of sol associated with filename
    """
    time_field = time_field.upper()
    time = make_seconds_since_midnight(filename, time_field=time_field)
    data, dummy = read_data(filename)
    pressure = []

    if filename.endswith('.xml'):
        pressure = data['TABLE']['PRESSURE']
    elif filename.endswith('.csv'):
        pressure = data['PRESSURE'].values

    return time, pressure

def read_Perseverance_ATS_data(filename, which_ATS=1, time_field='LTST',
        sol=None):
    """
    Read in Perseverance MEDA ATS data - https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/

    Args:
        filename (str): path to CSV file
        which_ATS (int or str): which of the five ATS sensors to read in;
        if "all", all of the ATS time series are returned
        time_field (str, optional): which time base to use

    Returns:
        time, pressure (float array): times and pressures, times in seconds
        since midnight of sol associated with filename

    """
    time_field=time_field.upper()
    # Note: ATS measures at 2 Hz, so there will be some duplicate LTST-values!
    time = make_seconds_since_midnight(filename, time_field=time_field)

    # Which ATS time-series to read in?
    if(isinstance(which_ATS, int)):
        which_ATS_str = "ATS_LOCAL_TEMP%i" % which_ATS
        temperature = pd.read_csv(filename)[which_ATS_str].values

        return time, temperature

    elif(isinstance(which_ATS, str) and (which_ATS == "all")):
        ATS1 = pd.read_csv(filename)["ATS_LOCAL_TEMP1"].values
        ATS2 = pd.read_csv(filename)["ATS_LOCAL_TEMP2"].values
        ATS3 = pd.read_csv(filename)["ATS_LOCAL_TEMP3"].values
        ATS4 = pd.read_csv(filename)["ATS_LOCAL_TEMP4"].values
        ATS5 = pd.read_csv(filename)["ATS_LOCAL_TEMP5"].values

        return time, [ATS1, ATS2, ATS3, ATS4, ATS5]
    
def read_Perseverance_WS_data(filename, sol=None, time_field='LTST'):
    """
    Read in WS data - (link)

    Args:
        filename (str): path to CSV file

    Returns:
        time, ws1-8 (float array): times and dimensions of wind speed, times in seconds
        since midnight of sol associated with filename
    """

    time = make_seconds_since_midnight(filename, time_field=time_field)

    ws1 = pd.read_csv(filename)['HORIZONTAL_WIND_SPEED'].values
    ws2 = pd.read_csv(filename)['HORIZONTAL_WIND_SPEED_UNCERTAINTY'].values
    ws3 = pd.read_csv(filename)['VERTICAL_WIND_SPEED'].values
    ws4 = pd.read_csv(filename)['VERTICAL_WIND_SPEED_UNCERTAINTY'].values
    ws5 = pd.read_csv(filename)['WIND_DIRECTION'].values
    ws6 = pd.read_csv(filename)['WIND_DIRECTION_UNCERTAINTY'].values
    ws7 = pd.read_csv(filename)['BOTH_BOOMS_USED_FOR_RETRIEVAL'].values
    ws8 = pd.read_csv(filename)['ROVER_STILL'].values

    return time, ws1, ws2, ws3, ws4, ws5, ws6, ws7, ws8

def make_seconds_since_midnight(filename, time_field='LTST', sol=None):
    """
    The MEDA data provide times in the LTST field in the format "sol hour:minute:second".

    Args:
        filename (str): name of the file
        time_field (str, optional): name of time field to analyze
        sol (int, optional): which is the primary sol; if not given, will
        determine from filename

    Returns:
        number of seconds in each row since midnight of the primary sol for that file
        """

    if(sol is None):
        primary_sol = which_sol(filename)

    data, _ = read_data(filename)
    time_field = time_field.upper()

    sols_str = [] #Sol is a solar day on Mars
    times_str = []
    delta_sols = []

    # Grab the sols and times associated with each row
    if(time_field == "LTST"):
        if(filename.endswith('.xml')):
            time_field_col = data['TABLE'][time_field]
            ltst_values = np.array(time_field_col)
            np.array(map(str, ltst_values))
            ltst_split = np.char.split(ltst_values,sep=' ') #splits fromchar 0 to space " "

            for x in range(len(ltst_values)):
                sols_str.append(float(ltst_split[x][0]))
                times_str.append(ltst_split[x][1])
            
            for x in range(len(ltst_values)):
                delta_sols.append(sols_str[x] - float(primary_sol))
        elif(filename.endswith('.csv')):
            sols_str = data[time_field].str.split(expand=True)[0].values #splits from char 0 to space " "
            times_str = data[time_field].str.split(expand=True)[1].values
            delta_sols = sols_str.astype(float) - float(primary_sol)

    elif(time_field == "LMST"):
        if(filename.endswith('.xml')):
            time_field_col = data['TABLE'][time_field]
            lmst_values = np.array(time_field_col)
            lmst_split = np.char.split(lmst_values,sep='M') #splits fromchar 0 to space " "
            
            for x in range(len(lmst_values)):
                times_str.append(lmst_split[x][1])
        elif(filename.endswith('.csv')):
            times_str = data[time_field].str.split("M", expand=True)[1].values

    time = time_conversion_to_seconds(times_str, delta_sols, time_field)

    return time

def time_conversion_to_seconds(times_str, delta_sols,time_field):
    """
    Time coversion for easy usability

    Args:
        times_str: time string for file
        time_field (str, optional): name of time field to analyze
        delta_sol (int, optional): difference between primary sol and sol from file
    Returns:
        converted time
    """
    time_field = time_field.upper()
    if(time_field=='LTST'):
        time = np.array([float(times_str[i].split(":")[0]) +\
                delta_sols[i]*24. + float(times_str[i].split(":")[1])/60 +\
                float(times_str[i].split(":")[2])/3600. for i in
                range(len(times_str))])
    elif(time_field=='LMST'):
        time = np.array([float(times_str[i].split(":")[0]) +\
                float(times_str[i].split(":")[1])/60 +\
                float(times_str[i].split(":")[2])/3600. for i in
                range(len(times_str))])
    else:
        raise Exception(time_field + " time field is not a valid option")

    return time


def which_sol(filename):
    """
    Based on the filename, returns the sol corresponding to a data file

    Args:
        filename (str): name of the file

    Returns:
        sol (int) associated with that file
    """

    ind = filename.find("WE__")

    return int(filename[ind+4:ind+8])