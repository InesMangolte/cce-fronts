"""
load biological data from original files (excel and csv) and convert them to netcdf.
The data in the .nc files can then be loaded using the functions in func_load.py
"""

import sys
sys.path.append('./')
import func_cce as f
import func_load as fl

import numpy as np
import pandas
import datetime as dt
from netCDF4 import Dataset


path_nc='./netcdf data files/'

path_original = './original data files/'

# add DOI of the article once submitted and published
description='{}of the overnight transects of CCE-LTER Process Cruises P0810 (transect A), P1107 (transects C1-C2-C3), P1208 (transects  E1-E2) and P1706 (transects F1-F2-F3). Each variable is organized by transect and station, with dimensions [transect, station, (depth)] {}. See Mangolte et al., 2023 for details.'

CTD_descr='Temperature (°C), Salinity (PSU), Density (kg/m3) and Fluorescence (V) measured during the CTD casts ', ''
nut_descr='Concentration of nutrients (Phosphate, Silicic Acid, Nitrite and Nitrate, {}) measured in the Niskin bottles ',''

hplc_descr='Concentration ({}) of phytoplankton taxa (equivalent Chl-a concentration) and total Chl-a determined with HPLC from the Niskin bottle samples ', 'Note that the A-front contains all vertical levels for the A-front, but the other transects only have the surface value.'

fc_descr='{} of pico-plankton taxa ({}) determined with flow cytometry from the Niskin bottle samples ', ''

zs_descr='{} of meso-zooplankton taxa ({}) determined with Zooscan from the Bongo tows ', ''

epi_descr='Carbon biomass of plankton taxa ({}) determined with Epifluorescence Microscopy from the Niskin bottle samples ', 'Note that the EPI analysis was only performed in transects A and E1'

#%% CTD  T, S, dens, fluo

loc = (path_original+'CTD Downcast Data (Process Cruise).csv')
df = pandas.read_csv(loc, index_col=1)  #the index is the transect name

# transect and variable names as they appear in the csv file
transect_id=['AFront', 'Upfront 1','Upfront 2', 'Upfront 3','EFront1','EFront2', 'Transect1', 'Transect2', 'Transect3']
var=['Depth (m)','Temperature 1 (C)','Temperature 2 (C)','Density 1 (kg/m³)','Density 2 (kg/m³)',
     'Potential Temp. 1 (C)','Potential Temp. 2 (C)', 'Surface PAR (µE/m²/s)','PAR (µE/m²/s)',
     'Corrected PAR (%)', 'Salinity 1 (PSU)','Salinity 2 (PSU)', 'Pressure (dbar)', 'Latitude (º)', 'Longitude (º)', 'Fluorescence (V)',]

data=np.zeros((len(transect_id), 13, 350, len(var)))  #13 stations max per transect, 307 depths
data.fill(np.nan)
for tr in range(len(transect_id)):
    df_t=df.loc[transect_id[tr]]  #get all data for the transect
    df_t=df_t.set_index('Cast Number')  #now the cast number is the index
    cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list shuffles the numbers, so have to sort them manually
    for st in range(len(cast)):  #for each station
        df_st=df_t.loc[cast[st], var]
        data[tr, st,0:len(df_st),:]=df_st  #take all the depths for all the variables

## compute composite T, S and D from the 2 sensors
t1=1 #temperature1
t2=2 #temperature2
d1=3 #density1
d2=4 #density2
s1=10 #salinity 1
s2=11 #salinity 2

temp_ctd=np.nanmean(data[:,:,:,[t1,t2]], axis=3) # [tr,st, z]
dens_ctd=np.nanmean(data[:,:,:,[d1,d2]], axis=3)
sal_ctd=np.nanmean(data[:,:,:,[s1,s2]], axis=3)



fluo_ctd=data[:,:,:,15]# fluorescence
z_ctd=data[:,:,:,0] #[tr, st, z]


## coordinates : latitude, longitude, date and time
var=['Latitude (º)','Longitude (º)']

#date and time : str, so put them in a list
DT=[[[] for ii in range(13)] for jj in range(len(transect_id))]

coord=np.zeros((len(transect_id), 13,  len(var)))  #13 stations max per transect, 15 depths
coord.fill(np.nan)
for tr in range(len(transect_id)):
    df_t=df.loc[transect_id[tr]]  #get all data for the transect
    df_t=df_t.set_index('Cast Number')  #now the cast number is the index
    cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list shuffles the numbers, so have to sort them manually
    for st in range(len(cast)):  #for each station
        df_st=df_t.loc[[cast[st]], var]    #take all the depths for all the variables (values are the same at all depths)
        coord[tr, st,:]=df_st.iloc[0]  #just take the first line of the station

        #get the date and time
        df_dt=df_t.loc[[cast[st]], 'Datetime GMT']
        DT[tr][st]=df_dt.iloc[0]


# manually correct the coordinates using the values in the ZooScan file
#correct the A-front coordinates : move all the stations up by 2, and add the last 2
coord[0,0:7,:]=coord[0,2:9]
coord[0,7]=[32.852, -120.710]
coord[0,8]=[32.897, -120.709]
# add missing values
coord[6,3,:]=[34.946,-121.457] # T1, station 4
coord[8,0,:]=[33.823,-122.715] # T3, station 1

# WRITE TO NETCDF FILE

tr_name=np.array(['A-front','C-front 1','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3'], dtype='object')


N_z=temp_ctd.shape[2] # max number of vertical levels
N=temp_ctd.shape[1] # max number of stations

fnc=Dataset(path_nc+'CTD.nc', 'w', format='NETCDF4')
fnc.description =description.format(CTD_descr[0], CTD_descr[1])

fnc.createDimension('transect', 9)
fnc.createDimension('station', N)
fnc.createDimension('z', N_z)

transect=fnc.createVariable('Transect', str, ('transect'))
transect[:]=tr_name
depth=fnc.createVariable('Depth', 'f4', ('transect','station','z'))
depth[:]=z_ctd
depth.units='meters'
density=fnc.createVariable('Density', 'f4', ('transect','station','z'))
density[:]=dens_ctd
density.units='kg/m³'
temperature=fnc.createVariable('Temperature', 'f4',  ('transect','station','z'))
temperature[:]=temp_ctd
temperature.units='Celsius'
salinity=fnc.createVariable('Salinity', 'f4',  ('transect','station','z'))
salinity[:]=sal_ctd
salinity.units='PSU'
fluorescence=fnc.createVariable('Fluorescence', 'f4',  ('transect','station','z'))
fluorescence[:]=fluo_ctd
fluorescence.units='V'

lat=fnc.createVariable('Latitude', 'f4',  ('transect','station'))
lat[:]=coord[:,:,0]
lat.units='°'
lon=fnc.createVariable('Longitude', 'f4',  ('transect','station'))
lon[:]=coord[:,:,1]
lon.units='°'
date=fnc.createVariable('Date', str,  ('transect','station'))
for tr in range(9): # only write the existing stations, otherwise error with the empty lists
    date[tr,:f.n_st(tr=tr, n_tr=9)]=np.array(DT[tr], dtype='object')[:f.n_st(tr=tr, n_tr=9)]

fnc.close() # save to disk



#%% Nutrients


loc = (path_original+'Dissolved Inorganic Nutrients (Process Cruise).csv')
df = pandas.read_csv(loc, index_col=4)  #the index is the transect name


transect_id=['5.5', 'Upfront 1','Upfront 2', 'Upfront 3','EFront1','EFront2', 'Transect1', 'Transect2', 'Transect3']   #as they appear in the csv file

var_nut=['Actual Depth (m)', 'Phosphate (µmol/L)', 'Silicic Acid (µmol/L)', 'Nitrite (µmol/L)', 'Nitrate (µmol/L)']

data_nut=np.zeros((len(transect_id), 13, 15, len(var_nut)))  #13 stations max per transect, 15 depths [tr, st, z, var]
data_nut.fill(np.nan)
for tr in range(len(transect_id)):
    df_t=df.loc[transect_id[tr]]  #get all data for the transect
    df_t=df_t.set_index('Cast')  #now the cast number is the index
    cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list shuffles the numbers, so have to sort them manually

    for st in range(len(cast)):  #for each station
        df_st=df_t.loc[[cast[st]], var_nut]    #take all the depths for all the variables

        #sort the get the depths in ascending order
        df_st=df_st.sort_values(by='Actual Depth (m)')
        data_nut[tr, st,0:len(df_st),:]=df_st
transect_id[0]='AFront'  #rename

z_nut=data_nut[:,:,:,0]
data_nut=data_nut[:,:,:,1:]
var_nut=var_nut[1:]
var_nut=[x[:-9] for x in var_nut] # remove the unit from the variable name

coord, var, transect_id, DT=fl.coord(withC1=True, datestr=True) # load from CTD.nc

# WRITE TO NETCDF FILE

nut_unit= 'µmol/L'

tr_name=np.array(['A-front','C-front 1','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3'], dtype='object')


N_z=data_nut.shape[2] # max number of vertical levels
N=data_nut.shape[1] # max number of stations

fnc=Dataset(path_nc+'Nutrients.nc', 'w', format='NETCDF4')
fnc.description = description.format(nut_descr[0].format(nut_unit), nut_descr[1])

fnc.createDimension('transect', 9)
fnc.createDimension('station', N)
fnc.createDimension('z', N_z)

# coordinates
transect=fnc.createVariable('Transect', str, ('transect'))
transect[:]=tr_name
lat=fnc.createVariable('Latitude', 'f4',  ('transect','station'))
lat[:]=coord[:,:,0]
lat.units='°'
lon=fnc.createVariable('Longitude', 'f4',  ('transect','station'))
lon[:]=coord[:,:,1]
lon.units='°'
date=fnc.createVariable('Date', str,  ('transect','station'))
for tr in range(9): # only write the existing stations, otherwise error with the empty lists
    date[tr,:f.n_st(tr=tr, n_tr=9)]=np.array(DT[tr], dtype='object')[:f.n_st(tr=tr, n_tr=9)]


# depth
depth=fnc.createVariable('Depth', 'f4', ('transect','station','z'))
depth[:]=z_nut
depth.units='meters'


for i, var in enumerate(var_nut):
    x=fnc.createVariable(var, 'f4',  ('transect','station','z'))
    x[:,:,:]=data_nut[:,:,:,i]
    x.units=nut_unit


fnc.close() # save to disk




#%% Flow cytometry
for biomass in [True, False]:

    if biomass:
        loc = (path_original+"flow_cytometry_biomass.csv")
        var_fc=['Depth (m)','Heterotrophic Bacteria (µg/L)','Prochlorococcus (µg/L)','Synechococcus (µg/L)']
        ii=5 # column with the transect name ("Cycle")
    else:
        loc = (path_original+"flow_cytometry.csv")
        var_fc=['Depth (m)','Heterotrophic Bacteria (number/ml)','Prochlorococcus (number/ml)','Synechococcus (number/ml)','Picoeukaryotes (number/ml)']
        ii=2 # column with the transect name ("Cycle (Event Log)")
    df = pandas.read_csv(loc, index_col=ii)  #the index is the transect name

    transect_id=['Afront', 'Upfront 1', 'Upfront 2', 'Upfront 3','EFRONT1','EFRONT2','Transect1','Transect2','Transect3'] #as they appear in the csv file


    data=np.zeros((len(transect_id), 13, 15, len(var_fc)))  #13 stations max per transect, 15 depths max : [tr, st, z, var]
    data.fill(np.nan)
    for tr in range(len(transect_id)):
        df_t=df.loc[transect_id[tr]]  #get all data for the transect
        df_t=df_t.set_index('Cast Number')  #now the cast number is the index
        cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list shuffles the numbers, so have to sort them manually
        for st in range(len(cast)):  #for each station
            df_st=np.array(df_t.loc[cast[st], var_fc])#take all the depths for all the variables [z, var]
            nz=df_st.shape[0]
            data[tr, st,0:nz,:]=df_st

    z_fc=data[:,:,:,0]

    data_fc=data[:,:,:,1:] # excludes depth
    var_fc=var_fc[1:] # excludes depth
    if biomass:
        var_fc=[x[:-7] for x in var_fc] # remove unit from variable name
    else:
        var_fc=[x[:-12] for x in var_fc] # remove unit from variable name

    coord, var, transect_id, DT=fl.coord(withC1=True, datestr=True) # load from CTD.nc

    # WRITE TO NETCDF FILE

    fc_unit= 'µgC/L' if biomass else 'nb cells/mL'
    filename='biomass' if biomass else 'abundance'


    tr_name=np.array(['A-front','C-front 1','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3'], dtype='object')


    N_z=data_fc.shape[2] # max number of vertical levels
    N=data_fc.shape[1] # max number of stations

    fnc=Dataset(path_nc+'flow_cytometry_'+filename+'.nc', 'w', format='NETCDF4')
    fnc.description = description.format(fc_descr[0].format(filename, fc_unit), fc_descr[1])

    fnc.createDimension('transect', 9)
    fnc.createDimension('station', N)
    fnc.createDimension('z', N_z)

    # coordinates
    transect=fnc.createVariable('Transect', str, ('transect'))
    transect[:]=tr_name
    lat=fnc.createVariable('Latitude', 'f4',  ('transect','station'))
    lat[:]=coord[:,:,0]
    lat.units='°'
    lon=fnc.createVariable('Longitude', 'f4',  ('transect','station'))
    lon[:]=coord[:,:,1]
    lon.units='°'
    date=fnc.createVariable('Date', str,  ('transect','station'))
    for tr in range(9): # only write the existing stations, otherwise error with the empty lists
        date[tr,:f.n_st(tr=tr, n_tr=9)]=np.array(DT[tr], dtype='object')[:f.n_st(tr=tr, n_tr=9)]


    # depth
    depth=fnc.createVariable('Depth', 'f4', ('transect','station','z'))
    depth[:]=z_fc
    depth.units='meters'


    for i, var in enumerate(var_fc):
        x=fnc.createVariable(var, 'f4',  ('transect','station','z'))
        x[:,:,:]=data_fc[:,:,:,i]
        x.units=fc_unit


    fnc.close() # save to disk




#%% HPLC

################# load surface data for C, E and F cruises

file=path_original+'/CCE-Process Cruise HPLC Transect Data - Nov 2019 original.xlsx'
df=pandas.read_excel(file, index_col=2, sheet_name='Data', header=3) #the index is the transect name (column "cycle"); column names are on line 4
df=df.iloc[2:] # skip the line with the units and an empty line

#as they appear in the csv file
var_hplc1=['Tchl a','Dinoflagellates','Diatoms','Prymnesiophytes','Pelagophytes','SYN','Chlorophytes','Cryptophytes','PRO']

transect_id=['Upfront 1', 'Upfront 2', 'Upfront 3','EFront1','EFront2','Transect1','Transect2','Transect3']


data_hplc1=np.zeros((len(transect_id), 13, len(var_hplc1)))  #13 stations max per transect : [tr, st, var]
data_hplc1.fill(np.nan)

for tr in range(len(transect_id)):
    df_t=df.loc[transect_id[tr]]  #get all data for the transect
    df_t=df_t.set_index('Cast')  #now the cast number is the index
    cast=list(df_t.index)
    for st in range(len(cast)):  #for each station
        df_st=np.array(df_t.loc[cast[st], var_hplc1]) #take all the variables
        data_hplc1[tr, st,:]=df_st # put them in the array

# replace negative values with 0
data_hplc1=np.where(data_hplc1<0, 0, data_hplc1)


################ load A-front vertical levels. 1 additional group = Prasino
'''
the surface bottle of station 2 is missing; there is a nan in the np.array
'''
file=path_original+'/A-Front HPLC PerContr.xlsx'
df_t=pandas.read_excel(file, index_col=2, sheet_name='Data').iloc[1:] #the index is the cast number; the line with the units must be skipped when loading the data

var_hplc_Afront=['Pressure', 'Tchl a','Dinoflagellates','Diatoms','Prymnesiophytes','Pelagophytes','SYN','Chlorophytes','Cryptophytes','PRO', 'Prasino']

data_load=np.zeros((9,8,len(var_hplc_Afront)))  #9 stations, 8 depths, var
data_load.fill(np.nan)
cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list shuffles the numbers, so have to sort them manually

for st in range(len(cast)):  #for each station
    df_st=np.array(df_t.loc[cast[st], var_hplc_Afront]) #take all the depths for all the variables [z,var]
    nz=df_st.shape[0]
    data_load[st,:nz,:]=df_st

data_load=data_load[:,::-1,:] # put the deepest bottle last

# replace negative values with 0
data_hplc_Afront=np.where(data_load<0, 0, data_load)

zz_hplc_A=data_hplc_Afront[:,:,0] # depths (Pressure, db)
data_hplc_Afront=data_hplc_Afront[:,:,1:] # take out the depth
del var_hplc_Afront[0]

coord, var, transect_id, DT=fl.coord(withC1=True, datestr=True) # load from CTD.nc

nzA=zz_hplc_A.shape[1]
z_hplc=np.zeros((9,13,15)) #[tr, st, z] same shape as the other files
z_hplc.fill(np.nan)
z_hplc[0,:9,:nzA]=zz_hplc_A




# WRITE TO NETCDF FILE

nst=f.n_st(tr=None, n_tr=9) # number of stations in each transect

hplc_unit= 'μg chla/L'

tr_name=np.array(['A-front','C-front 1','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3'], dtype='object')

N_z=z_hplc.shape[2] # max number of vertical levels
N=data_hplc1.shape[1] # max number of stations

fnc=Dataset(path_nc+'HPLC.nc', 'w', format='NETCDF4')
fnc.description = description.format(hplc_descr[0].format(hplc_unit), hplc_descr[1])

fnc.createDimension('transect', 9)
fnc.createDimension('station', N)
fnc.createDimension('z', N_z)

# coordinates
transect=fnc.createVariable('Transect', str, ('transect'))
transect[:]=tr_name
lat=fnc.createVariable('Latitude', 'f4',  ('transect','station'))
lat[:]=coord[:,:,0]
lat.units='°'
lon=fnc.createVariable('Longitude', 'f4',  ('transect','station'))
lon[:]=coord[:,:,1]
lon.units='°'
date=fnc.createVariable('Date', str,  ('transect','station'))
for tr in range(9): # only write the existing stations, otherwise error with the empty lists
    date[tr,:nst[tr]]=np.array(DT[tr], dtype='object')[:nst[tr]]


# depth
depth=fnc.createVariable('Depth', 'f4', ('transect','station','z'))
depth[:]=z_hplc
depth.units='meters'


for i, var in enumerate(var_hplc_Afront):
    x=fnc.createVariable(var, 'f4',  ('transect','station','z'))

    # A-front = 9 stations, 8 vertical levels
    x[0,:9,:8]=data_hplc_Afront[:,:,i]

    # other transects = surface
    for tr in np.arange(1,9):
        if var!='Prasino': # Prasino = only for Afront
            x[tr,:nst[tr],0]=data_hplc1[tr-1,:nst[tr],i] # no Afront in data_hplc1

    x.units=hplc_unit


fnc.close() # save to disk





#%% Zooscan


for biomass in [True, False]:

    filename='biomass' if biomass else 'abundance'
    zs_unit='mgC/m^2' if biomass else '(#/m^2)'




    file=path_original+'/zooscanVerticalBongoData_30Mar2022'
    df=pandas.read_csv(file + '.csv', index_col=0) #the index is the cruise name

    transect_id=['A-front','C-front 1','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3']
    cruise=['P0810', 'P1106','P1106', 'P1106', 'P1208', 'P1208', 'P1706', 'P1706', 'P1706'] # cruise name for each transect
    trnum=[1,1,2,3,1,2,1,2,3] # transect number in the cruise for each transect


    # create a new column with the transect name and make it the index
    df_tr=np.array(df[['Transect']]).squeeze()
    df_cr=np.array(df.index)
    trname=np.zeros(len(df_tr))

    for ii in range(len(transect_id)):
        boo=np.logical_and(df_cr==cruise[ii], df_tr==trnum[ii])
        trname=np.where(boo, transect_id[ii], trname)

    df['Transect name']=trname
    df=df.set_index('Transect name')  #now the transect name is the index

    var_zs=['appendicularia','chaetognatha','cnidaria_ctenophores','copepoda_calanoida','copepoda_eucalanids',
         'copepoda_oithona_like','copepoda_poecilostomatoids','copepoda_harpacticoida','copepoda_others',
         'crustacea_others','doliolids','euphausiids','ostracods','polychaete','pteropoda',
         'pyrosomes', 'rhizaria', 'salps', 'eggs', 'nauplii', 'All Copepods']


    data_name='Carbon_mg_sum' if biomass else'Abundance (#/m^2)'

    data_zs=np.zeros((len(transect_id), 13, len(var_zs)))  #13 stations max per transect, 15 depths max
    data_zs.fill(np.nan)

    for tr in range(len(transect_id)):
        df_t=df.loc[transect_id[tr]]  #get all data for the transect


        df_t=df_t.set_index('Station')  #now the station number is the index
        cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list sometimes shuffles the numbers, so have to sort them manually
        for st in range(len(cast)):  #for each station
            df_st=df_t.loc[cast[st]]  #take the abundance of all groups for this station

            df_st=df_st.set_index('Class')  #now the taxa name is the index

            if 'salps' not in df_st.index: # there is no salp category when salps where not detected; so it needs to be added (with the value 0)
                x = pandas.DataFrame([0], columns=[data_name], index=['salps'])
                df_st = pandas.concat([df_st,x])

            data_zs[tr, st,:]=df_st.loc[var_zs, data_name]  #take the values for all the categories


    coord, var, transect_id, DT=fl.coord(withC1=True, datestr=True)

    # WRITE TO NETCDF FILE

    tr_name=np.array(['A-front','C-front 1','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3'], dtype='object')


    N=data_zs.shape[1] # max number of stations

    fnc=Dataset(path_nc+'ZooScan_'+filename+'.nc', 'w', format='NETCDF4')
    fnc.description = description.format(zs_descr[0].format(filename, zs_unit), fc_descr[1])

    fnc.createDimension('transect', 9)
    fnc.createDimension('station', N)

    # coordinates
    transect=fnc.createVariable('Transect', str, ('transect'))
    transect[:]=tr_name
    lat=fnc.createVariable('Latitude', 'f4',  ('transect','station'))
    lat[:]=coord[:,:,0]
    lat.units='°'
    lon=fnc.createVariable('Longitude', 'f4',  ('transect','station'))
    lon[:]=coord[:,:,1]
    lon.units='°'
    date=fnc.createVariable('Date', str,  ('transect','station'))
    for tr in range(9): # only write the existing stations, otherwise error with the empty lists
        date[tr,:f.n_st(tr=tr, n_tr=9)]=np.array(DT[tr], dtype='object')[:f.n_st(tr=tr, n_tr=9)]

    for i, var in enumerate(var_zs):
        x=fnc.createVariable(var, 'f4',  ('transect','station'))
        x[:,:]=data_zs[:,:,i]
        x.units=zs_unit


    fnc.close() # save to disk


#%% EPI
'''
only for A and E1
already converted to carbon biomass, unit = (µgC/L)
'''


loc = (path_original+"EPI_biomass.csv") # there is no column with the transect name
var_fc=['Depth (m)','Heterotrophic Bacteria (µg/L)','Prochlorococcus (µg/L)','Synechococcus (µg/L)']

df = pandas.read_csv(loc, index_col=0)  #the index is the cruise name

cruise_names = ['P0810MV', 'P1208MV'] # in column 'StudyName'
station_names = [['FS{:02d}'.format(i) for i in np.arange(1,10)],['FS{:02d}'.format(i) for i in np.arange(1,14)] ] # in column 'Cycle'

var_names = list(df.keys()[8:]) # plankton taxa
var_names.insert(0, 'Depth (m)')

data=np.zeros((2, 13, 15, len(var_names)))  #13 stations max per transect, 15 depths max : [tr, st, z, var]
data.fill(np.nan)

for cr in range(2):
    df_t=df.loc[cruise_names[cr]]  #get all data for the cruise
    df_t=df_t.set_index('Cycle')  #now the cast number is the index
    for st in range(len(station_names[cr])):  #for each station
        df_st=np.array(df_t.loc[station_names[cr][st], var_names]) #take all the depths for all the variables
        data[cr, st,0:len(df_st),:]=df_st


data_epi=data[:,:,:,1:] # excludes depth [tr, st, z, var]
var_epi= [x[:-7] for x in var_names[1:]] # excludes depth and take out unit
z_epi=data[:,:,:,0] # [tr, st, z]

coord, var, transect_id, DT=fl.coord(withC1=True, datestr=True) # load from CTD.nc
# only take A and E1
coord=coord[[0,4]]
DT= [DT[ii] for ii in [0,4]]

# WRITE TO NETCDF FILE

epi_unit= 'µgC/L'


tr_name=np.array(['A-front', 'E-front 1'], dtype='object')
nst=[9, 13]

N_z=data_fc.shape[2] # max number of vertical levels
N=data_fc.shape[1] # max number of stations

fnc=Dataset(path_nc+'EPI_biomass.nc', 'w', format='NETCDF4')
fnc.description = description.format(epi_descr[0].format(epi_unit), epi_descr[1])

fnc.createDimension('transect', 2)
fnc.createDimension('station', N)
fnc.createDimension('z', N_z)

# coordinates
transect=fnc.createVariable('Transect', str, ('transect'))
transect[:]=tr_name
lat=fnc.createVariable('Latitude', 'f4',  ('transect','station'))
lat[:]=coord[:,:,0]
lat.units='°'
lon=fnc.createVariable('Longitude', 'f4',  ('transect','station'))
lon[:]=coord[:,:,1]
lon.units='°'
date=fnc.createVariable('Date', str,  ('transect','station'))
for tr in range(2): # only write the existing stations, otherwise error with the empty lists
    date[tr,:nst[tr]]=np.array(DT[tr], dtype='object')[:nst[tr]]


# depth
depth=fnc.createVariable('Depth', 'f4', ('transect','station','z'))
depth[:]=z_epi
depth.units='meters'


for i, var in enumerate(var_epi):
    x=fnc.createVariable(var, 'f4',  ('transect','station','z'))
    x[:,:,:]=data_epi[:,:,:,i]
    x.units=epi_unit


fnc.close() # save to disk


