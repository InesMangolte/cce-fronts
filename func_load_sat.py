"""

NB : The satellite data must be downloaded separately

"""
import sys

sys.path.append('./')
import func_sat as f_sat


import numpy as np
from netCDF4 import Dataset

import datetime as dt



path_submeso=''
path_ssh=''
path_fsle=''


def load_submesocolor(dataname, dates):
    '''
    data, lat, lon = load_submesocolor(date)
    dSST and CHL data was downloaded and HI was computed using the submeso-color GitHub
                (https://gitlab.in2p3.fr/clement.haeck/submeso-color) :
        SST and Chl are downloaded from Copernicus,
            CHL : L3=with clouds; L4=cloud-free (space-time interpolation)

    timestep = daily, 01/01/2008 to 31/12/2019
    resolution = 4 km

    unit of Chl = mg/m3

    name of HI file :
        - 20 = scale (size of window, in km)
        - 1 = number (eg to separate spatial resolution, I won't need to use another one)
        - 102 = coef (101 = cce_total, 102 = cce_cruise)

    arguments :
        dataname = 'hi101', 'hi_compV', 'sst', 'chl', 'chl_l3'
        date = datetime object or list of DTs
    returns :
        data = [lat, lon] or [time, lat, lon]
        lat = list of latitude
        lon = list of latitude
    '''
    if not isinstance (dates, list) :
        dates=[dates]


    date_str=[x.strftime('%Y%m%d') for x in dates] # string 'YYYYMMDD'

    if 'hi_comp' in dataname:
        path_sub=path_submeso +'HI/HI_20.0_1/1days/{}/HI'
        varname=dataname[-1]
    elif 'hi' in dataname:
        path_sub=path_submeso +'HI/HI_20.0_1/1days/{}/HI_final'+dataname[2:]
        varname='hi_final'
    elif dataname=='sst':
        path_sub= path_submeso + 'OSTIA/1days/{}/SST_processed'
        varname='sst'
    elif dataname=='chl':
        path_sub=path_submeso + 'GlobColour/L4/1days/{}/CHL'
        varname='CHL'
    elif dataname=='chl_l3':
        path_sub=path_submeso + 'GlobColour/L3/1days/{}/CHL'
        varname='CHL'
    else:
        print('dataname = sst,  chl, chl_l3,  hi+coef or hi_comp+V_S_B')
        return


    # only load the grid once
    t=0
    year=dates[t].strftime('%Y')
    filename=path_sub.format(year) + '_{}.nc'.format(date_str[t])
    lat=np.array(Dataset(filename,'r').variables['lat'])
    lon=np.array(Dataset(filename,'r').variables['lon'])

    out=np.zeros((len(dates), len(lat), len(lon)))

    for t in range(len(date_str)):
        year=dates[t].strftime('%Y')
        filename=path_sub.format(year) + '_{}.nc'.format(date_str[t])

        try:
            out[t]=np.array(Dataset(filename,'r').variables[varname])
        except (KeyError, ValueError): # some files are corrupted
            print('Couldn''t load {} in file {}'.format(varname, filename))
            out[t].fill(np.nan)

    out=np.where(out==-999, np.nan, out) # turn land into nans
    out=np.where(out==-32767, np.nan, out) # HI missing values

    return np.squeeze(out), lat, lon


def load_submeso_clouds(dates):
    '''
    loads a mask with the clouds (and land), from the CHL L3 data
    NB : the land are also nans, but in the plots it won't be visible (land in black, clouds in black hatches)

    arguments :
        date = datetime object or list of DTs

    '''
    C, lat, lon = load_submesocolor('chl_l3', dates)

    out=np.where(np.isnan(C), 1, 0)
    return out, lat, lon





def load_FLSE(dates, var):
    '''
    data, lat, lon = load_FSLE(date)

    loads FLSE computed by Pierre Chabert in the CCE region (cf Chabert et al., 2021)
    uses geostrophic (altimetry) + ageostrophic (Ekman currents computed from wind)
    timestep = 4 days
    resolution = 0.05° ~ 4 km

    argument :
        date = datetime object or list of DTs.
        var='lambda' (FSLE), 'touched' (Water age)
    returns :
        data = [lat, lon] or [time, lat, lon]
        lat = np.array of size [301]
        lon = np.array of size [341]

    '''

    if not isinstance (dates, list) :
        dates=[dates]

    #coordinates
    delta0 = .05
    lonv0, lonv1 = -135, -118
    latv0, latv1 = 30, 45
    lon = np.arange(lonv0, lonv1+ delta0, delta0)  #len=341, longitude
    lat = np.arange(latv0, latv1 + delta0, delta0)  #len=301, latitude



    out=np.zeros((len(dates), len(lat), len(lon)))

    for ii, d in enumerate(dates):

        # timestep = 4 days; year always starts on jan 01
        day_no = d.toordinal() - dt.date(d.year, 1, 1).toordinal() # number of the day in the year

        dd4=4*(day_no//4)+1  # day with a file

        d4=dt.date(d.year, 1, 1) + dt.timedelta(days=dd4-1)

        # print('day number : {}, previous day with a file : {}'.format(d, d4))

        filename=path_fsle+'ekman_CA_'+d4.strftime('%Y_%m_%d')+'.nc'

        out[ii]=np.array(Dataset(filename, 'r').variables[var])


    if var=='touched':
        out=np.where(out==0, np.nan, out)
        out=np.where(out==-1, np.nan, out)

    return np.squeeze(out), lat, lon



def load_ssh(dates):
    '''
    data, lat, lon = load_ssh(date)
    data is from COPERNICUS website : Global Reanalysis, sea_surface_height_above_geoid (m)
    https://resources.marine.copernicus.eu/product-download/GLOBAL_REANALYSIS_PHY_001_030

    timestep = daily 01/01/2008 - 31/12/2019
    resolution = 1/12° = 8 km

    argument :
        date = datetime object or list of DTs
    returns :
        data = [lat, lon] or [time, lat, lon]
        lat
        lon
    '''
    if not isinstance (dates, list) :
        dates=[dates]


    file_ssh=path_ssh+'global-reanalysis-phy-ssh.nc'
    data=Dataset(file_ssh,'r')

    time=list(np.array(data.variables['time'])+12) # there's a 12 hour delay in the .nc times
    ssh=np.array(data.variables['zos'])
    lat=np.array(data.variables['latitude'])
    lon=np.array(data.variables['longitude'])

    nh=f_sat.hours_since_1950(dates)
    inds=[time.index(x) for x in nh]

    out=ssh[inds,:,:]

    # turn the land into nans
    out=np.where(out==-32767, np.nan, out)

    data.close()

    return out, lat, lon
