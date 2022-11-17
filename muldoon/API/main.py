from argparse import FileType
from urllib import response
from fastapi import *
import requests
from Create_Graphs import *



#create our fast api app
app = FastAPI()

#URL checker 
def checkUrl(PSurl):
    get = requests.get(PSurl)
    return get.status_code


def getUrlValues(val):
    
    if(int(val) < 90):
        val='00' + val
        val2 = 'sol_0000_0089'
    elif(int(val) > 89 and int(val) <180):
        if(len(val) < 3):
            val='00' + val
        else:
            val='0' + val
        val2 = 'sol_0090_0179'
    elif(int(val) >179 and int(val) < 300):
        val = '0' + val
        val2 = 'sol_0180_0299'
    else:
        val = '0' + val
        val2 = 'sol_0300_0419'
    return val,val2




@app.get("/PS/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]

    val,val2 = getUrlValues(val)
    

    PSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P01.CSV'
    PSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P02.CSV'


    if checkUrl(PSurl2) == 200:
        return create_PS_Graph(PSurl2)
    elif checkUrl(PSurl) == 200:
        return create_PS_Graph(PSurl)
    else: 
        raise HTTPException(status_code=404, detail="Item not found")  


@app.get("/PSGap/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]


    val,val2 = getUrlValues(val)

    PSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P01.CSV'
    PSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P02.CSV'


    if checkUrl(PSurl2) == 200:
        return create_PS_Graph_with_gaps(PSurl2)
    elif checkUrl(PSurl) == 200:
        return create_PS_Graph_with_gaps(PSurl)
    else: 
        raise HTTPException(status_code=404, detail="Item not found") 


@app.get("/PSDetrend/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]

    val,val2 = getUrlValues(val)

    PSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P01.CSV'
    PSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P02.CSV'


    if checkUrl(PSurl2) == 200:
        return create_PS_Graph_detrend(PSurl2)
    elif checkUrl(PSurl) == 200:
        return create_PS_Graph_detrend(PSurl)
    else: 
        raise HTTPException(status_code=404, detail="Item not found")


@app.get("/PSVortex/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]

    val,val2 = getUrlValues(val)

    PSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P01.CSV'
    PSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P02.CSV'


    if checkUrl(PSurl2) == 200:
        return create_PS_Vortex_Signals(PSurl2)
    elif checkUrl(PSurl) == 200:
        return create_PS_Vortex_Signals(PSurl)
    else: 
        raise HTTPException(status_code=404, detail="Item not found")

@app.get("/ATS/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]

    val,val2 = getUrlValues(val)

    ATSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________CAL_ATS_________________P01.CSV'
    ATSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/'+ val2 + '/sol_' + val + '/WE__' + val + '___________CAL_ATS__________________P02.CSV'

    if checkUrl(ATSurl) == 200:
        return create_ATS_Graph(ATSurl)
    elif checkUrl(ATSurl2) == 200:
        return create_ATS_Graph(ATSurl2)
    else: 
        raise HTTPException(status_code=404, detail="Item not found")


@app.get("/ATSDetrend/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]

    val,val2 = getUrlValues(val)

    ATSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________CAL_ATS_________________P01.CSV'
    ATSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/'+ val2 + '/sol_' + val + '/WE__' + val + '___________CAL_ATS__________________P02.CSV'


    if checkUrl(ATSurl) == 200:
        return create_ATS_Graph_Detrended(ATSurl)
    elif checkUrl(ATSurl2) == 200:
        return create_ATS_Graph_Detrended(ATSurl2)
    else: 
        raise HTTPException(status_code=404, detail="Item not found")   

         


@app.get("/ATS_PS/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"]    
    
    
    val,val2 = getUrlValues(val)

    PSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P01.CSV'
    PSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_PS__________________P02.CSV'

    ATSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_calibrated_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________CAL_ATS_________________P01.CSV'

    if(checkUrl(PSurl2) == 200 and checkUrl(ATSurl)) == 200:
        return create_ATS_PS_Graph(PSurl2,ATSurl)
    elif(checkUrl(PSurl) == 200 and checkUrl(ATSurl) == 200):
        return create_ATS_PS_Graph(PSurl,ATSurl)
    else: 
        raise HTTPException(status_code=404, detail="File path not found for {} or {}".format(PSurl,ATSurl))  


@app.get("/WS/{sol}")
async def root(sol: str):
    sol = {"sol": sol}
    val = sol["sol"] 

    val,val2 = getUrlValues(val)
    WSurl = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_WS__________________P01.xml'
    WSurl2 = 'https://pds-atmospheres.nmsu.edu/PDS/data/PDS4/Mars2020/mars2020_meda/data_derived_env/' + val2 + '/sol_' + val + '/WE__' + val + '___________DER_WS__________________P02.xml'


    if checkUrl(WSurl) == 200:
        return create_WS_Graph(WSurl)
    elif checkUrl(WSurl2) == 200:
        return create_WS_Graph(WSurl2)
    else: 
        print(WSurl2)
        raise HTTPException(status_code=404, detail="Item not found") 
    