'''
functions to load the data
'''
import geopy.distance


import numpy as np


import func_cce as f



# must be entirely rewritten with the new netcdf files


def coord():
    '''
    data, var, transect_id, DT=coord()

    load the coordinates of the stations (latitude, longitude)


    return :
        data : numpy array [transect, station, var]  #2 var = lat, lon
        var : name of the variables
        transect_id : name of the transects
        DT : date and time of the stations as a list of lists [tr][st]
    '''

    loc = (path+'data/CTD data/CTD Downcast Data (Process Cruise).csv')
    df = pandas.read_csv(loc, index_col=1)  #the index is the transect name

    transect_id=['AFront', 'Upfront 1','Upfront 2', 'Upfront 3','EFront1','EFront2', 'Transect1', 'Transect2', 'Transect3']   #as they appear in the csv file

    var=['Latitude (º)','Longitude (º)']

    #date and time
    DT=[[[] for ii in range(13)] for jj in range(len(transect_id))]

    data=np.zeros((len(transect_id), 13,  len(var)))  #13 stations max per transect, 15 depths
    data.fill(np.nan)
    for tr in range(len(transect_id)):
        df_t=df.loc[transect_id[tr]]  #get all data for the transect
        df_t=df_t.set_index('Cast Number')  #now the cast number is the index
        cast=sorted(list(set(df_t.index)))  #set() gives the unique values of the cast number. for some reason turning the set into a list shuffles the numbers, so have to sort them manually
        print(transect_id[tr], 'len(cast)=', len(cast))
        for st in range(len(cast)):  #for each station
            df_st=df_t.loc[[cast[st]], var]    #take all the depths for all the variables (values are the same at all depths)
            data[tr, st,:]=df_st.iloc[0]  #just take the first line of the station

            #get the date and time
            df_dt=df_t.loc[[cast[st]], 'Datetime GMT']
            DT[tr][st]=df_dt.iloc[0]


    # manually add missing values from UVP file
    data[6][3]=[34.9463,-121.4565]
    data[8][0]=[33.8235,-122.715167]
    #correct the A-front coordinates : move all the stations up by 2, and add the last 2 from UVP file
    data[0][0:7]=data[0][2:9]
    data[0][7]=[32.8515, -120.71066]
    data[0][8]=[32.89633, -120.70867]

    transect_id=['A-front', 'C-front 1','C-front 2', 'C-front 3','E-front 1','E-front 2', 'Transect 1', 'Transect 2', 'Transect 3']   #as they appear in the csv file

    return data, var, transect_id, DT



def along():
    '''
    out=along()
    returns :
        along transect coordinates (in km), for the 8 main transects (excluding C1)
             dimensions = [transect, station] ; non-existant stations are nans
    '''



    data, var, transect_id, datetime=coord()  #read from CTD file; data = [tr, st, coord]


    ntr=data.shape[0]  #number of transects
    nst=data.shape[1]   #max number of stations

    dist=np.zeros((ntr,nst))   #[tr,st]
    dist.fill(np.nan)

    for tr in range(ntr):
        num_station=f.n_st(tr, n_tr=9)

        coords1=(data[tr,0,0], data[tr,0,1])   #1st station (lat,lon)
        for st in range(num_station):
            coords2=(data[tr,st,0], data[tr,st,1])
            dist[tr,st]=geopy.distance.distance(coords1, coords2).km

    out=  np.delete(np.array(dist), 1, axis=0) # remove transect C1
    return out




def trophic(carb=True, cyan='fc', Afront_ave=False):
    '''
    out, names, transect_id = f.trophic()

    loads plankton data and computes the sum for each trophic levels :
        h-bact, cyanobact, other euk phy, diatoms, pico-grazers, meso-grazers, carnivores
    using ZooScan, HPLC and Flowcyt
        Zooscan = vertical integral, absolute abundance
        HPLC = surface values (or vertical average for A-front), absolute chlorophyll concentration
        flowcyt = vertical integral, absolute abundance


    also removes negative values and replaces them by 0

    arguments:
        carb : if True, converts all data to
        cyan : use flow cytometry (vertical integral) or HPLC (surface only) data
        Afront_ave : if True, uses vertical average for A-front; if False, uses surface values

    returns :
        out= np.array [tr, st, troph_lev]
        names = name of each trophic level + instrument + unit

    '''


    #zooscan
    if carb:
        data_zs, var, transect_id , data_name= fl.load_zooscan(rel=False, data_type='cbiomass', longnames=True)  #[tr, st, ig]
        data_zs=data_zs*1E3/100 # conversion mgC/m2 --> μgC/m3
    else:
        data_zs, var, transect_id , data_name= fl.load_zooscan(rel=False, data_type='abundance', longnames=True)  #[tr, st, ig]
        # unit='(#/m^2)'  #absolute


    #trophic groups
    ig4 = [var.index('rhizaria'),var.index('doliolids') ,var.index('appendicularia'),var.index('pyrosomes')]
    ig5 = [var.index('Calanoids'),var.index('copepoda_oithona_like') ,var.index('Other copepods (Poec. + Harp.)'),var.index('euphausiids'),var.index('pteropoda')]
    ig6 = [var.index('chaetognatha'),var.index('cnidaria_ctenophores') ,var.index('ostracods') ,var.index('polychaete') ]

    g4=np.nansum(data_zs[:,:,ig4],axis=2)
    g5=data_zs[:,:,ig5].sum(axis=2)
    g6=data_zs[:,:,ig6].sum(axis=2)


    #flow cyt
    data_fc, var, transect_id=fl.load_flowcyt(rel=False, vint=True, biomass=carb)  #[transect, station, species]
    # unit = μg/m3

    ig0=[0]  #heterotrophic bacteria
    ig1_fc=[1,2]#PRO+SYN
    g0=np.squeeze(data_fc[:,:,ig0])

    #HPLC
    data_hplc, var, transect_id=fl.load_hplc(rel=False, Afront_ave=Afront_ave, carb=carb)  #[transect, station, species]
    # unit = μg/m3


    #trophic groups
    ig1_hplc=[5,8]#PRO+SYN
    ig2=[1,3,4,6,7]#other euk-phyto
    ig3=[2] #diatoms

    g2=data_hplc[:,:,ig2].sum(axis=2)
    g3=data_hplc[:,:,ig3].squeeze()

    if cyan=='hplc':
        g1=data_hplc[:,:,ig1_hplc].sum(axis=2)
    elif cyan=='fc':
        g1=data_fc[:,:,ig1_fc].sum(axis=2)

    if carb:
        if cyan=='hplc':
            names = ['h-bact (flowcyt, μgC/m3 )','cyanobact (HPLC, μgC/m3 ','other euk-phyto (HPLC, μgC/m3','diatoms (HPLC, μgC/m3','pico-grazers (ZooScan, μgC/m3)','meso-grazers (ZooScan, μgC/m3)','carnivores (ZooScan, μgC/m3)']
        elif cyan=='fc':
            names = ['h-bact (flowcyt, μgC/m3 )','cyanobact (flowcyt, μgC/m3 ','other euk-phyto (HPLC, μgC/m3','diatoms (HPLC, μgC/m3','pico-grazers (ZooScan, μgC/m3)','meso-grazers (ZooScan, μgC/m3)','carnivores (ZooScan, μgC/m3)']

    else:
        if cyan=='hplc':
            names = ['h-bact (flowcyt, #/m^2 )','cyanobact (HPLC, μg chla/m3 ','other euk-phyto (HPLC, μg chla/m3','diatoms (HPLC, μg chla/m3','pico-grazers (ZooScan, #/m^2)','meso-grazers (ZooScan, #/m^2)','carnivores (ZooScan, #/m^2)']
        elif cyan=='fc':
            names = ['h-bact (flowcyt, #/m^2 )','cyanobact (flowcyt, #/m^2 ','other euk-phyto (HPLC, μgC/m3','diatoms (HPLC, μgC/m3','pico-grazers (ZooScan, μgC/m3)','meso-grazers (ZooScan, μgC/m3)','carnivores (ZooScan, μgC/m3)']

    out=np.stack([g0, g1, g2, g3, g4, g5, g6], axis=-1)
    print('end of f.trophic() : ', names)

    out=np.where(out<0, 0, out)

    return out, names, transect_id


def trophic_detail(carb=True, cyan='hplc', Afront_ave=False):
    '''
    out, names, var_names, transect_id = f.trophic_detail(carb=True, cyan='fc', Afront_ave=False)

    loads plankton data and returns it organized by trophic levels :
        h-bact, cyanobact, other euk phy, diatoms, pico-grazers, meso-grazers, carnivores
    using ZooScan, HPLC and Flowcyt
        Zooscan = vertical integral, absolute abundance
        HPLC = surface values (or vertical average for A-front), absolute chlorophyll concentration
        flowcyt = vertical integral, absolute abundance


    also removes negative values and replaces them by 0

    arguments:
        carb : if True, converts all data to
        cyan : use flow cytometry (vertical integral) or HPLC (surface only) data
        Afront_ave : if True, uses vertical average for A-front; if False, uses surface values

    returns :
        out= list of 7 np.array [tr, st, var]
        var_names = list of lists with names of taxa in each trophic level
        names = name of each trophic level + instrument + unit
        transect_id = name of transects

    old zooscan file (before april 2022) :
    ig4=[0,7,13,14]     #  'pico-grazers '  #rhizaria, doliolids, appendicularians, salps
    ig5=[3,4,5,6,8,11,12]   # 'meso-grazers '  #copepods, euphausiids, pteropoda, pyrosoma
    ig6=[1,2,9,10]  #  'carnivorous '  #chaetognaths, cnidarians, ostracods, polychaetes


    '''


    #zooscan
    if carb:
        data_zs, var, transect_id , data_name= fl.load_zooscan(rel=False, data_type='cbiomass', longnames=True)  #[tr, st, ig]
        data_zs=data_zs*1E3/100 # conversion mgC/m2 --> μgC/m3
    else:
        data_zs, var, transect_id , data_name= fl.load_zooscan(rel=False, data_type='abundance', longnames=True)  #[tr, st, ig]
        # unit='(#/m^2)'  #absolute

    # get the short names
    var_zs= fl.load_zooscan(rel=False, data_type='abundance', longnames=False, names_only=True)  #[tr, st, ig]


    #trophic groups
    ig4 = [var.index('rhizaria'),var.index('doliolids') ,var.index('appendicularia'),var.index('pyrosomes'), var.index('salps')]
    ig5 = [var.index('Calanoids'),var.index('copepoda_oithona_like') ,var.index('Other copepods (Poec. + Harp.)'),var.index('euphausiids'),var.index('pteropoda'),var.index('crustacea_others')]
    ig6 = [var.index('chaetognatha'),var.index('cnidaria_ctenophores') ,var.index('ostracods') ,var.index('polychaete') ]


    g4=data_zs[:,:,ig4]
    g5=data_zs[:,:,ig5]
    g6=data_zs[:,:,ig6]
    var4=np.array(var_zs)[ig4]
    var5=np.array(var_zs)[ig5]
    var6=np.array(var_zs)[ig6]


    #flow cyt
    data_fc, var_fc, transect_id=fl.load_flowcyt(rel=False, vint=True, biomass=carb)  #[transect, station, species]
    # unit=μgC/m3

    ig0=[0]  #heterotrophic bacteria
    ig1_fc=[1,2]#PRO+SYN
    g0=np.squeeze(data_fc[:,:,ig0])
    var0=np.array(var_fc)[ig0]

    #HPLC
    data_hplc, var_hplc, transect_id=fl.load_hplc(rel=False, Afront_ave=Afront_ave, carb=carb)  #[transect, station, species]
    # unit= μgC/m3

    #trophic groups
    ig1_hplc=[5,8]#PRO+SYN
    ig2=[1,3,4,6,7]#other euk-phyto
    ig3=[2] #diatoms

    g2=data_hplc[:,:,ig2]
    var2=np.array(var_hplc)[ig2]
    g3=data_hplc[:,:,ig3].squeeze()
    var3=np.array(var_hplc)[ig3]

    if cyan=='hplc':
        g1=data_hplc[:,:,ig1_hplc]
        var1=np.array(var_hplc)[ig1_hplc]
    elif cyan=='fc':
        g1=data_fc[:,:,ig1_fc]
        var1=np.array(var_fc)[ig1_fc]

    if carb:
        if cyan=='hplc':
            names = ['h-bact (flowcyt, μgC/m3 )','cyanobact (HPLC, μgC/m3 ','other euk-phyto (HPLC, μgC/m3','diatoms (HPLC, μgC/m3','pico-grazers (ZooScan, μgC/m3)','meso-grazers (ZooScan, μgC/m3)','carnivores (ZooScan, μgC/m3)']
        elif cyan=='fc':
            names = ['h-bact (flowcyt, μgC/m3 )','cyanobact (flowcyt, μgC/m3 ','other euk-phyto (HPLC, μgC/m3','diatoms (HPLC, μgC/m3','pico-grazers (ZooScan, μgC/m3)','meso-grazers (ZooScan, μgC/m3)','carnivores (ZooScan, μgC/m3)']

    else:
        if cyan=='hplc':
            names = ['h-bact (flowcyt, #/m^2 )','cyanobact (HPLC, μg chla/m3 ','other euk-phyto (HPLC, μg chla/m3','diatoms (HPLC, μg chla/m3','pico-grazers (ZooScan, #/m^2)','meso-grazers (ZooScan, #/m^2)','carnivores (ZooScan, #/m^2)']
        elif cyan=='fc':
            names = ['h-bact (flowcyt, #/m^2 )','cyanobact (flowcyt, #/m^2 ','other euk-phyto (HPLC, μgC/m3','diatoms (HPLC, μgC/m3','pico-grazers (ZooScan, μgC/m3)','meso-grazers (ZooScan, μgC/m3)','carnivores (ZooScan, μgC/m3)']

    out=[g0, g1, g2, g3, g4, g5, g6]
    var_names=[var0, var1, var2, var3, var4, var5, var6]
    print('end of f.trophic() : ', names)

    for ii, x in enumerate(out):
        out[ii]=np.where(x<0, 0, x)

    return out, names, var_names, transect_id












