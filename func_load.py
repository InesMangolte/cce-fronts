
import sys
sys.path.append('./')
import func_cce as f

import numpy as np
from netCDF4 import Dataset
import geopy.distance
import datetime as dt


path_nc='./netcdf data files/'


def coord(withC1=False, datestr=False):
    '''
    data, var, transect_id, DT=coord(withC1=False)

    load the coordinates of the stations (latitude, longitude) from CTD file
    arguments:
        withC1 : if True, returns all 9 transects. Otherwise ignore C1
        datestr : if True, returns the dates a str, if False, as dt.datetime objects

    return :
        data_coord : numpy array [transect, station, var]  #2 var = lat, lon
        var : list of the names of the variables
        transect_id : name of the transects
        DT : date and time of the stations as a list of lists [tr][st]
    '''
    data=Dataset(path_nc+'CTD.nc', 'r')
    print(path_nc+'CTD.nc', data.description)

    trs=[0,1,2,3,4,5,6,7,8] if withC1 else [0,2,3,4,5,6,7,8]


    transect_id=[data.variables['Transect'][i] for i in trs]

    n_st=f.n_st(tr=None,n_tr=9) # for converting the datetimes, must take n_tr=9 since trs has the transect indices including C1
    DT=[list(data.variables['Date'][tr,:n_st[tr]]) for tr in trs] # str, all sublists have the same length


    data_coord=np.zeros((len(trs), np.max(f.n_st(tr=None, n_tr=len(trs))), 2))
    var_list=['Latitude', 'Longitude']
    for i, var in enumerate(var_list):
        data_coord[:,:,i]=data.variables[var][trs,:]


    data.close()
    if not datestr:
        DT=[[dt.datetime.strptime(x,'%Y-%m-%d %H:%M:%S').date() for x in xtr] for xtr in DT]

    return data_coord, var_list, transect_id, DT



def along():
    '''
    out=along()
    returns :
        along transect coordinates (in km), for the 8 main transects (excluding C1)
             dimensions = [transect, station] ; non-existant stations are nans
    '''

    data, var, transect_id, datetime=coord()  #read from CTD file; data = [tr, st, coord]
    # coord2, var, transect_id, datetime=fl.coord(withC1=True)  #read from CTD file; data = [tr, st, coord]


    ntr=data.shape[0]  #number of transects
    nst=data.shape[1]   #max number of stations

    dist=np.zeros((ntr,nst))   #[tr,st]
    dist.fill(np.nan)

    for tr in range(ntr):
        num_station=f.n_st(tr, n_tr=8)

        coords1=(data[tr,0,0], data[tr,0,1])   #1st station (lat,lon)
        for st in range(num_station):
            coords2=(data[tr,st,0], data[tr,st,1])
            dist[tr,st]=geopy.distance.distance(coords1, coords2).km

    return dist


def load_CTD_depth():
    '''
    z=load_CTD_depth()
    loads the vertical grid of CTD data, C1 is ecluded
    returns :
        numpy array [tr, st, z]
    '''


    data=Dataset(path_nc+'CTD.nc', 'r')
    print(path_nc+'CTD.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1

    out=np.zeros((len(trs), np.max(f.n_st(tr=None, n_tr=len(trs)))))
    var='Depth'
    out=data.variables[var][trs,:]
    data.close()

    return out

def load_dens(vave=None):
    '''
    out=load_temp()
    load composite density from the 2 CTD sensors
    arguments :
        vave : if None, returns all the vertical levels. if given, computes the vertical average.
                vave= number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns numpy array [tr, st, z] or [tr, st]
    NB : C1 is excluded
    '''

    data=Dataset(path_nc+'CTD.nc', 'r') 	#lit les données dans une classe python netcdf
    print(path_nc+'CTD.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1

    out=np.zeros((len(trs), np.max(f.n_st(tr=None, n_tr=len(trs)))))
    var='Density'
    out=data.variables[var][trs,:]

    if vave is not None:
        depths=data.variables['Depth'][trs,:]
        out=np.average(out[:,:,:vave], axis=2, weights=np.diff(depths[:,:,:vave+1])) #--> [tr,st]

    data.close()
    return out

def load_temp(vave=None):
    '''
    out=load_temp()
    load composite temperature from the 2 CTD sensors
    arguments :
        vave : if None, returns all the vertical levels. if given, computes the vertical average.
                vave= number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns numpy array [tr, st, z] or [tr, st]
    NB : C1 is excluded
    '''

    data=Dataset(path_nc+'CTD.nc', 'r') 	#lit les données dans une classe python netcdf
    print(path_nc+'CTD.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1

    out=np.zeros((len(trs), np.max(f.n_st(tr=None, n_tr=len(trs)))))
    var='Temperature'
    out=data.variables[var][trs,:]

    if vave is not None:
        depths=data.variables['Depth'][trs,:]
        out=np.average(out[:,:,:vave], axis=2, weights=np.diff(depths[:,:,:vave+1])) #--> [tr,st]

    data.close()
    return out



def load_sal(vave=None):
    '''
    out=load_salinity()
    load composite salinity from the 2 CTD sensors
    arguments :
        vave : if None, returns all the vertical levels. if given, computes the vertical average.
                vave= number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns numpy array [tr, st, z] or [tr, st]
    NB : C1 is excluded
    '''

    data=Dataset(path_nc+'CTD.nc', 'r') 	#lit les données dans une classe python netcdf
    print(path_nc+'CTD.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1

    out=np.zeros((len(trs), np.max(f.n_st(tr=None, n_tr=len(trs)))))
    var='Salinity'
    out=data.variables[var][trs,:]

    if vave is not None:
        depths=data.variables['Depth'][trs,:]
        out=np.average(out[:,:,:vave], axis=2, weights=np.diff(depths[:,:,:vave+1])) #--> [tr,st]

    data.close()
    return out



def load_fluo(vave=None):
    '''
    out=load_fluo()
    load fluorescence
    arguments :
        vave : if None, returns all the vertical levels. if given, computes the vertical average.
                vave= number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns numpy array [tr, st, z] or [tr, st]
    NB : C1 is excluded
    '''

    data=Dataset(path_nc+'CTD.nc', 'r') 	#lit les données dans une classe python netcdf
    print(path_nc+'CTD.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1

    out=np.zeros((len(trs), np.max(f.n_st(tr=None, n_tr=len(trs)))))
    var='Fluorescence'
    out=data.variables[var][trs,:]

    if vave is not None:
        depths=data.variables['Depth'][trs,:]
        out=np.average(out[:,:,:vave], axis=2, weights=np.diff(depths[:,:,:vave+1])) #--> [tr,st]

    data.close()
    return out



def load_flowcyt(vave=True, biomass=True, rel=False):
    '''
    out, var_list, data_name= fl.load_flowcyt(vave=True, biomass=True, rel=False)
    load flow cytometry data downloaded from DataZoo (Mike Landry; conversion to biomass was already done)
    4 groups : Pro, Syn, Pico-euk, heterotrophic microbes

    arguments :
        rel=True for relative (%), False for absolute
        vave = True for vertical average, False for all vertical levels
        biomass = True to load biomass, False to load abundance

    returns :
        out : [transect, station, depth, var] or [transect, station, var]
        var_list : name of the variables (depth + 3 plankton groups or only the 3 plankton groups)
        data_name : type of data and unit

    NB : unit of data in netcdf file : fc_unit= 'µgC/L' if biomass else 'nb cells/mL'

    '''
    if biomass:
        data=Dataset(path_nc+'flow_cytometry_biomass.nc', 'r')
        print(path_nc+'flow_cytometry_biomass.nc', data.description)
    else:
        data=Dataset(path_nc+'flow_cytometry_abundance.nc', 'r')
        print(path_nc+'flow_cytometry_abundance.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1
    nz=len(data.dimensions['z'])
    ntr=len(trs)
    nst=len(data.dimensions['station'])


    var_list=list(data.variables.keys())[4:] if biomass else list(data.variables.keys())[4:-1] # the first 4 are transect, lat, lon, and date; also don't take the pico-euk


    out_v=np.zeros((ntr, nst, nz,  len(var_list))) #[tr, st, z, var]
    for i, var in enumerate(var_list):
        out_v[:,:,:,i]=data.variables[var][trs,:,:]
    if biomass :
        out_v[:,:,:,1:]=out_v[:,:,:,1:]*1E3 # μg/L --> μg/m3

    if vave:
        zz=out_v[:,:,:,0]
        del var_list[0]
        out=np.zeros((ntr, nst, len(var_list))) #[tr, st, var]
        for i, var in enumerate(var_list):
            out[:,:,i]=vave_niskin(X=out_v[:,:,:,i+1], zz=zz)
    else:
        out=out_v


    # once all the taxa are loaded, compute the relative abundance or biomass
    if rel:
        if vave:
            out[:,:,:]=f.rel_abundance(out[:,:,:])
        else:
            out[:,:,:,1:]=f.rel_abundance(out[:,:,:,1:]) # first variable are the depths


    if rel:
        data_name='relative '
        if biomass:
            data_name = data_name + 'biomass (% of total bact. biomass)'
        else:
            data_name = data_name + 'abundance (% of total bact. abundance)'
    else:
        if biomass:
            data_name = 'biomass (µgC/m3)'
        else:
            data_name = 'abundance (nb cells/mL)'
    if vave:
        data_name = 'vertical average of ' + data_name

    data.close()

    return out, var_list, data_name


# def zooscan_names():
#     '''
#     returns the list of zooscan taxa from the .nc file, unnested

#     pretty names (for plot tick labels for instance)
#     var = ['Appendicularia', 'Chaetognaths', 'Cnidaria', 'Calanoids', 'Oithona', 'Other cop', 'other crust', 'Doliolids',
#            'Euphausiids', 'Ostracods', 'Polychaete', 'Pteropods', 'Pyrosomes', 'Rhizaria', 'Salps', 'Eggs',
#            'Nauplii', 'All Copepods']

#     '''
#     data=Dataset(path_nc+'ZooScan_biomass.nc', 'r')
#     var=data.variables.keys() # first 4 = tr, lat, lon, date
#     return list(var)[4:]

def load_zooscan(biomass=True, rel=False):
    '''
    out, var_list, data_name = fl.load_zooscan(biomass=True, rel=False)

    load all zooscan taxa (including eggs and nauplii)

    arguments :
        biomass = True to load biomass, False to load abundance. conversion was done by Mark Ohman
        rel= True for relative (% of total zoo. biomass), False for absolute

    returns :
        out = absolute or relative abundance/biomass [transect, station, groups]. units : μgC/m^3 or #/m^2
        var_list = list of taxa names
        data_name = data type + unit

    NB : unit of data in netcdf file : zs_unit='mgC/m^2' if biomass else '(#/m^2)'
    NB : there were only 12 stations sampled for EFRONT1 (station 13 is missing)

    '''
    if biomass:
        data=Dataset(path_nc+'ZooScan_biomass.nc', 'r')
        print(path_nc+'ZooScan_biomass.nc', data.description)
    else:
        data=Dataset(path_nc+'ZooScan_abundance.nc', 'r')
        print(path_nc+'ZooScan_abundance.nc', data.description)


    trs=[0,2,3,4,5,6,7,8] # all transects except C1
    ntr=len(trs)
    nst=len(data.dimensions['station'])
    # transect_id=[data.variables['Transect'][i] for i in trs]


    var_list=list(data.variables.keys())[4:] # the first 4 are transect, lat, lon, and date

    out=np.zeros((ntr, nst, len(var_list))) #[tr, st, var]
    for i, var in enumerate(var_list):
        out[:,:,i]=data.variables[var][trs,:] # [tr, st]
    if biomass:
        out= out* 1E3/100  # conversion mgC/m2 --> μgC/m3 (Bongo tow = 100 m depth)

    # group the copepods into 3 categories : Calanoids, Oithona (no modification necessary) and others
    ical = var_list.index('copepoda_calanoida')
    ieuc = var_list.index('copepoda_eucalanids')
    out[:,:,ical] = out[:,:,ical] + out[:,:,ieuc]
    var_list[ical] = 'Calanoids'

    ipoe = var_list.index('copepoda_poecilostomatoids')
    ihar = var_list.index('copepoda_harpacticoida')
    ioth = var_list.index('copepoda_others')
    out[:,:,ipoe] = out[:,:,ipoe] + out[:,:,ihar]+ out[:,:,ioth]
    var_list[ipoe] = 'Other copepods (Poec. + Harp.)'

    for ii in np.sort([ieuc, ihar, ioth])[::-1]: # delete the now extra taxa in decreasing order so the indices don't get confused
        out=np.delete(out, ii, axis=2)
        del var_list[ii]

    if rel :
        data_name = 'relative (%)'
    elif biomass:
        data_name = 'biomass (μgC/m^3)'
    else:
        data_name = 'abundance (#/m^2)'



    data.close()
    return out, var_list, data_name




def load_hplc(rel=False, Afront_ave=False, biomass=True):
    '''
    out, var_list, data_name = fl.load_hplc(rel=False, Afront_ave=False, biomass=True)

    arguments :
        rel= True for relative (%), False for absolute
        Afront_ave = True to compute vertical average of the A-front, False uses the surface value only
        biomass = True to convert from chl to carbon

    returns :
        out : [transect, station, var]
        var_list : list with the names of the taxa
        data_name : data type and unit

    NB : units in the netcdf file : hplc_unit= 'μg chla/L'
    NB : Chl:C ratios estimated from Liet al., 2010 --> diatoms are 0.02 and all other groups are 0.01
    '''
    ratios=np.array([1,0.010,0.020,0.010,0.010,0.010,0.010,0.010,0.010, 0.010])  #Chl:C, 1st value is for total chl a (no conversion).

    data=Dataset(path_nc+'HPLC.nc', 'r')
    print(path_nc+'HPLC.nc', data.description)

    var_list=list(data.variables.keys())[5:] # the first 5 are transect, date, lat, lon and depth


    trs=[0,2,3,4,5,6,7,8] # all transects except C1
    ntr=len(trs)
    nst=len(data.dimensions['station'])

    out=np.zeros((ntr, nst, len(var_list)))# [tr, st, var]
    out.fill(np.nan)

    for i, var in enumerate(var_list):
        X=data.variables[var][trs,:,:] # [tr, st, z]

        if Afront_ave:
            AA = vave_niskin(X, zz=data.variables['Depth'][trs,:,:]) # Afront will be the vertical average; [tr,st]
            out[0,:,i] = AA[0,:] # vertical average for Afront
            out[1:,:,i] = X[1:,:,0] # surface values for the others


        else:
            out[:,:,i] = X[:,:,0] # surface value for all transects


        if biomass:
             out[:,:,i]=out[:,:,i]/ratios[i]*1E3 # μg chla/L --> μgC/m3

    # once all the taxa are loaded, compute the relative abundance or biomass
    if rel:
        out=f.rel_abundance(out)
        data_name= 'relative biomass (% of total phyto. biomass)' if biomass else '% of total chl-a'
    else:
        data_name= 'biomass (μgC/m3)' if biomass else 'equivalent chl-a concentration (μg chla/L)'


    data.close()

    return out, var_list, data_name



def load_HPLC_Afront(biomass=False):
    '''
    out, z, var_list, hplc_unit = load_HPLC_Afront(biomass=False)
    loads HPLC A-front data for all vertical levels

    arguments:
        biomass : if True, convert from μg chl/L to μgC/m3
    returns :
        data : [st, z, var]
        z : [st, z] NB : are identical to FC depths since they were used while writing the .nc
        var : list of variables
        unit : unit of the data 'μgC/L' if biomass else 'μg chl/L'
    '''
    hplc_unit = 'μgC/m3' if biomass else 'μg chl/L'

    ratios=np.array([1,0.010,0.020,0.010,0.010,0.010,0.010,0.010,0.010, 0.010])  #Chl:C, 1st value is for total chl a (no conversion).

    data=Dataset(path_nc+'HPLC.nc', 'r')
    print(path_nc+'HPLC.nc', data.description)

    var_list=list(data.variables.keys())[5:]

    nst=len(data.dimensions['station'])
    nz=len(data.dimensions['z'])

    out=np.zeros((nst, nz, len(var_list)))# [st, z, var]
    z=data.variables['Depth'][0,:,:]

    for i, var in enumerate(var_list):
        out[:,:,i]=data.variables[var][0,:,:] # [st, z]
        if biomass:
             out[:,:,i]=out[:,:,i]/ratios[i]*1E3 # μg chla/L --> μgC/m3

    data.close()
    return out, z, var_list, hplc_unit



def load_nut(vave=None):
    '''
    out, var_list =load_nut(vave=None)
    load fluorescence
    arguments :
        vave : if None, returns all the vertical levels. if given, computes the vertical integral.
                vave= number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns :
        out : numpy array [tr, st, z, var] or [tr, st, var]
        var_list : list with names of variables

    NB : unit =  'µmol/L'
    NB : C1 is excluded
    '''

    data=Dataset(path_nc+'Nutrients.nc', 'r')
    print(path_nc+'Nutrients.nc', data.description)

    trs=[0,2,3,4,5,6,7,8] # all transects except C1
    nz=len(data.dimensions['z'])
    ntr=len(trs)
    nst=len(data.dimensions['station'])

    if vave:
        var_list=list(data.variables.keys())[5:] # the first 5 are transect, date, lat, lon and depth
        out=np.zeros((ntr, nst, len(var_list))) #[tr, st, var]
        for i, var in enumerate(var_list):
            data_v=data.variables[var][trs,:,:] # [tr, st, z]

            out[:,:,i]=vave_niskin(X=data_v, zz=data.variables['Depth'][trs,:,:])

    else:
        var_list=list(data.variables.keys())[4:] # the first 4 are transect, date, lat, and lon
        out=np.zeros((ntr, nst, nz,  len(var_list))) #[tr, st, z, var]
        for i, var in enumerate(var_list):
            out[:,:,:,i]=data.variables[var][trs,:,:]

    data.close()
    return out, var_list


def load_epi(vave=True):
    '''
    out, var_list, data_name= fl.load_epi(vave=True)
    load flow EPI data downloaded from DataZoo (Mike Landry; conversion to biomass was already done)

    arguments :
        vave = True for vertical average, False for all vertical levels

    returns :
        out : [transect, station, depth, var] or [transect, station, var]
        var_list : name of the variables (depth + plankton groups or only the 3 plankton groups)
        data_name : type of data and unit

    NB : unit of data in netcdf file : epi_unit= 'µgC/L'
    NB : the EPI analysis was only done in transects A and E1

    '''
    data=Dataset(path_nc+'EPI_biomass.nc', 'r')
    print(path_nc+'EPI_biomass.nc', data.description)

    trs=[0,1] # A and E1
    nz=len(data.dimensions['z'])
    ntr=len(trs)
    nst=len(data.dimensions['station'])


    if vave:
        var_list=list(data.variables.keys())[5:]# the first 5 are transect, lat, lon, date, and depth
        out=np.zeros((ntr, nst, len(var_list))) #[tr, st, var]
        for i, var in enumerate(var_list):
            data_v=data.variables[var][trs,:,:] # [tr, st, z]
            out[:,:,i]=vave_niskin(X=data_v, zz=data.variables['Depth'][trs,:,:])
            out=out*1E3 # μg/L --> μg/m3


    else:
        var_list=list(data.variables.keys())[4:]  # the first 4 are transect, lat, lon, and date
        out=np.zeros((ntr, nst, nz,  len(var_list))) #[tr, st, z, var]
        for i, var in enumerate(var_list):
            out[:,:,:,i]=data.variables[var][trs,:,:]
        out[:,:,:,1:]=out[:,:,:,1:]*1E3 # μg/L --> μg/m3


    data_name = 'biomass (µgC/m3)'
    if vave:
        data_name = 'vertical average of ' + data_name

    data.close()

    return out, var_list, data_name




def vave_niskin(X, zz):
    '''
    computes the vertical average of data in Niskin bottles
    The weights represents the thickness of the layer where the Niskin data is "representative" (halfway between successive bottles)

    arguments:
        X : [tr, st, z]
        zz : [tr, st, z]
    returns :
        out : [tr, st]

    NB : there are 2 stations where no bottle samples were analyzed :
        FC in transect F3 station 10;
        nutrients in E1 station 6
    '''


    out=np.zeros((X.shape[0],X.shape[1])) # [tr, st, z]
    out.fill(np.nan)

    for tr in range(X.shape[0]):
        stats_flag=np.any(np.isfinite(zz[tr,:,:]), axis=1) #[st]
        stats=np.arange(X.shape[1])[stats_flag]
        for st in stats:
            # compute the weights of each layer (=thickness, in m)
            z0=zz[tr, st] # all depths for station st
            nz=np.isfinite(z0).sum() # number of levels
            wei=np.zeros((nz))

            wei[0]=(z0[0] + z0[1])/2 # first layer
            wei[nz-1]=(z0[nz-1] - z0[nz-2])/2 # last layer
            for ii in np.arange(1,nz-1): #middle layers
                wei[ii]=(z0[ii+1] - z0[ii-1])/2

            # compute the vertical average (must do it this way because of the nans at depth, which cannnot be dealt with using np.average)
            out_int = np.nansum(np.multiply(X[tr, st,:nz], wei)) # integral --> [tr, st]
            out[tr, st] = np.divide(out_int, np.nansum(wei)) # average

    return out
