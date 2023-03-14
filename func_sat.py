import sys
sys.path.append("./")
import func_load as fl
import func_cce as f_cce


import numpy as np
import pickle
import datetime as dt


def draw_landmask(ax, color='k'):
    '''
    draws the land mask (see compute_save.py)
    '''
    if ax.set_ylim()[0]<32: # for large scale maps
        fn='/Users/INES/Desktop/SIO/scripts/landmask.p'

        print('drawing landmask ...')
    else:
        # print('landmask zoom')
        fn='/Users/INES/Desktop/SIO/scripts/landmask_zoom1.p'

    with open(fn, 'rb') as file :
        mask=pickle.load(file)
        lat=pickle.load(file)
        lon=pickle.load(file)

    ax.contourf(lon, lat, mask, levels=[0,1], colors=color)
    # ax.pcolormesh(lon, lat, mask, cmap='Greys', vmin=0,vmax=1)
    # print('landmask done !')

def week_num(month, day):
    '''
    returns the week number
    '''
    mdays  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    ndays = np.sum(mdays[:month-1]) + day

    out =np.ceil(ndays/7)
    return int(out)


def week_date(N):
    '''
    returns a list of 7 'MMDD' for the given week number
    '''

    mdays  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    out=[[] for i in range(7)]

    for ii in range(7): # for each day of the week}
        ndays=7*(N-1)+ii+1 # number of day

        month=1
        while ndays > np.sum(mdays[:month]):
            month=month+1

        day=ndays - np.sum(mdays[:month-1])

        out[ii]='{:02d}{:02d}'.format(month, int(day))

    return out


def lim_zoom(zoom, cruise_only=False):
    '''
    defines xlim, ylim for various degrees of zooming in on the fronts, on satellite maps
    0 : entire CCS (HI statistics 'cce_total', 'hi101')
    1 : CCE cruises
    2 : transect scale
    3 : pixel scale


    for zoom=2 and 3, returns lists of 8 lims
    ['AFront', 'Upfront 2', 'Upfront 3', 'EFront1', 'EFront2',  'Transect1', 'Transect2', 'Transect3']

        if cruise_only == True, returns the 4 values for each cruise, without the duplicates for the transects
    '''

    if zoom==1: # regional scale = rectangle around cruises (=cce_cruise, hi102 in SUBMESOCOLOR)
        # xlim=[-125,-118]
        xlim=[-124,-119]
        # ylim=[31,36]
        ylim=[32,36]
    elif zoom==0: # very large scale = entire CCS (=cce_total, hi101 in SUBMESOCOLOR)
        xlim=[-133,-110]
        ylim=[20,50]
    elif zoom==2: # transect scale
        if not cruise_only:
            # xlim=[[-122,-119.5],[-122.5,-120.25],[-122.5,-120.25],[-124,-121.5],[-124,-121.5],[-123.3,-120.8],[-123.3,-120.8],[-123.3,-120.8]]
            # ylim=[[32,34],[33,34.8],[33,34.8],[33.6,35.6],[33.6,35.6],[33.5,35.5],[33.5,35.5],[33.5,35.5]]

            # wider region for F1-2-3
            xlim=[[-122,-119.5],[-122.5,-120.25],[-122.5,-120.25],[-124,-121.5],[-124,-121.5],[-124.5,-120],[-124.5,-120],[-124.5,-120]]
            ylim=[[32,34],[33,34.8],[33,34.8],[33.6,35.6],[33.6,35.6],[33,36.5],[33,36.5],[33,36.5]]
        else:
            xlim=[[-122,-119.5],[-122.5,-120.25],[-124,-121.5],[-123.3,-120.8]]
            ylim=[[32,34],[33,34.8],[33.6,35.6],[33.5,35.5]]

    elif zoom==3: # mega zoom to see individual pixels and stations
        if not cruise_only:
            xlim=[[-121,-120.5],[-121.9,-121.5],[-121.3,-120.7],[-123.3,-122.4],[-123.3,-122.4],[-121.8,-121],[-123.2,-121.6],[-123,-122.2]]
            ylim=[[32.5,33]    ,[33.65,34.05]    ,[33.3,33.8]    ,[34.2,35]      ,[34.2,35]      ,[34.7,35.4],[34,35.4],[33.6,34.4]]
        else:
            xlim=[[-121,-120.5],[-122,-120.6],[-123.4,-122.4],[-123.2,-120.8]]
            ylim=[[32.5,33],[33.3,34.2],[34.2,35],[33.5,35.5]]

    return xlim, ylim


def logclim(clim, ncont=25):
    '''
    creates levels with log spacing to use in contourf

    arguments :
        clim=[lo, hi]
        ncont = integer

    '''
    if clim[0]==0:
        clim[0]=1E-6
    lev_exp = np.linspace(np.log10(clim[0]), np.log10(clim[1]), ncont)
    levs = np.power(10, lev_exp)
    return levs


def remove_white_lines(cf):
    '''
    to get rid of the white contour lines in the pdf
    '''
    try:
        for c in cf.collections: # for contourf
                c.set_edgecolor("face")
    except AttributeError:
        cf.set_edgecolor("face") # for pcolormesh


def project_transects(data, lat, lon):
    '''
    find the pixel corresponding to the stations of the in-situ transects


    arguments :
        data = [tr, lat, lon] (8 transects = 8 images)
        lat, lon = lists or 1D arrays

    returns :
        data_tr = [tr, st]


    '''

    stations, var, transect_id, DT=fl.coord()
    stations=np.delete(stations, 1, axis=0)  #delete upfront1
    n_tr=stations.shape[0]
    n_st=f_cce.n_st(tr=None, n_tr=8)  #number of stations in each transect

    data_tr=np.zeros((n_tr,np.max(n_st))) # [tr, st]
    data_tr.fill(np.nan)

    for tr in range(n_tr):
        for st in range(n_st[tr]):
            c=stations[tr,st,:]  #station coord (lat, lon)
            i1=np.argmin(abs(c[0]-lat))
            j1=np.argmin(abs(c[1]-lon))
            data_tr[tr,st]=data[tr,i1,j1]
            print(i1,j1)

    return data_tr


def hours_since_1950(date):
    '''
    converts a datetime object to the number of hours since 01/01/1950
    argument :
        date : datetime object or list of DTs
    '''
    if type(date) is not list:
        date=[date]
    if type(date[0])==dt.date: #  create a datetime object and set the time to midnight
        date=[dt.datetime.combine(d, dt.time.min) for d in date]

    out=[[] for ii in range(len(date))]

    for ii,x in enumerate(date):

        out[ii]=int((x - dt.datetime(1950, 1, 1)).total_seconds() / 3600)


    return out

def date_range(start_date,end_date ):
    '''
    returns a list of dates between the start and end dates given
    '''
    delta = dt.timedelta(days=1)

    dates = []
    while start_date <= end_date:
        dates.append(start_date)
        start_date += delta
    return dates
