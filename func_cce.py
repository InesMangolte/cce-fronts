
import sys
sys.path.append('./')
import func_load as fl

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

cm=1/2.54  # centimeters to inches


def transect_id(short=False, add_date=False):
    '''
    names of the 8 transects analysed
    arguments :
        short : True = short or False = long version of transect name
        add_date : if True, adds the date of the transect
    '''
    if short:
        out= ['A','C2', 'C3', 'E1', 'E2', 'F1', 'F2', 'F3']
    else:
        out=['A-front','C-front 2', 'C-front 3', 'E-front 1', 'E-front 2', 'F-front 1', 'F-front 2', 'F-front 3']
    if add_date:
        tr_date=['oct 24-25 2008', 'july 2-3 2011', 'july 15-16 2011', 'august 4-5 2012', 'august 20-21 2012', 'june 7-8 2017', 'june 17-18 2017', 'june 22-23 2017']
        for i,x in enumerate(tr_date):
            out[i]=out[i] + '  ('+x + ')'
    return out



def n_st(tr, n_tr):
    '''
    number of stations in each transect
    arguments :
        tr = transect index. if None, returns a list with all the numbers
        n_tr = total number of transects (9 = all ; 8 = C1 is removed)
    '''
    if n_tr==9:
        st=[9, 5, 10, 10, 13, 10, 11, 11, 11]
    elif n_tr==8:
        st=[9, 10, 10, 13, 10, 11, 11, 11]

    if tr is None :
        return np.array(st)
    else:
        return np.array(st)[tr]


######## definitions of fronts and biomass peaks

def def_fronts():
    '''
    Definition of the front/warm/cold stations

    warm, front, cold, tr_number, name = f.def_fronts()
    returns :
        warm, front, cold : list of station numbers (not indices !) on each side of the front
        tr_number : index of the transect of each front
        name : name of the fronts

    '''
    # vérifier quelle def est la plus récente dans le texte de l'article. voir pourquoi les diatomées apparaissent en pic dans E1
    warm =  [[2]      , [1]              , [10]          , [1,2,3,4]      , [10]     , [2]    , [1]    , [11]    , [11]    , [1,2,3,4] ]
    front = [[3,4,5]  , [2,3,4,5,6,7,8,9], [7,8,9]       , [5,6,7,8,9,10] , [6,7,8,9], [3,4]  , [2,3]  , [9,10]  , [9,10]  , [5,6,7, 8, 9,10]]
    cold =  [[6,7,8]  , [10]             , [1,2,3,4,5,6] , [11,12,13]     , [5]      , [5]    , [4,5,6], [6,7,8] , [7,8]   , [11] ]
    name = ['A'       , 'C2'             ,'C3'           , 'E1'           , 'E2a'    , 'E2b'  , 'F1a'  ,'F1b'    , 'F2'    , 'F3' ]
    tr_number = [0,1,2,3,4,4,5,5,6,7 ]


    return warm, front, cold, tr_number, name


def biomass_peaks():
    '''
    Definition of the biomass peaks
    returns :
        width=dict with width of peaks (km)
        pos=dict with center of the peak (km)
        p = dict with list of peak station numbers
        var : list of variables (density grad + fluo + all taxa names, unnested)


    in the dictionaries, front ids are the keys of the dict, each dict value is a list with the variables

    In the initial creation of the dict, the variables are in the same order as in f.trophic_detail but not nested, with fluo first:
    fluo
    [array(['Het Bact'], dtype='<U8'),
     array(['PRO', 'SYN'], dtype='<U8'),
     array(['Dino', 'Prym', 'Pelago', 'Chloro', 'Crypto'], dtype='<U7'),
     array(['Diatoms'], dtype='<U7'),
     array(['Rhizaria', 'Doliolids', 'Appendicularia', 'Pyrosomes', 'Salps'],
           dtype='<U14'),
     array(['Calanoids', 'Oithona', 'Other cop', 'Euphausiids', 'Pteropods',
            'other crust'], dtype='<U14'),
     array(['Chaetognaths', 'Cnidaria', 'Ostracods', 'Polychaete'],
           dtype='<U14')]

    '''
    # just to get the names of taxa
    data_detail, taxa_names, group_names= trophic_detail(Afront_ave=True)

    var=[val for sublist in taxa_names for val in sublist]

    var.insert(0,'fluo')

    p=dict() # station number in the peak. If empty, it's a transition


    ## define the stations with fluo and biomass peaks
    #A
    p['A']=[[4,5,6,7],
            [],
            [],[4,5],
            [4,5], [4,5], [4,5],[4,5], [],
            [4,5],
            [2,3],[], [3,4,5], [], [],
            [3,4], [5], [], [3,4], [], [],
            [3,4,5], [], [3], [5]
            ]


    p['C2']=[[3,4,5,6,7],
            [3,4],
             [], [5,6,7],
             [2,3,4], [5], [5], [], [],
             [3,4],
             [4,5],[3,4,5],[4,5,6], [],[],
             [4,5,6,7],[5,6,7],[4,5,6,7], [3,4,5,6],[],[5],
             [4,5,6,7,8,9], [3,4,5,6,7], [], [],
             ]



    p['C3']=[[],
            [],
             [], [],
             [], [7,8], [], [7], [],
             [],
             [], [6,7], [5,6,7], [],[],
             [], [], [], [], [], [],
             [], [], [],[7]
             ]


    p['E1']=[[3,4,5,6,7,8,9,10,11],
            [],
             [], [],
             [11], [11], [10,11], [11], [11],
             [10,11],
             [10], [11], [6,7,8], [],[],
             [5,6,7,8,9,10,11], [6,7,8,9,10], [6,7,8,9,10], [8,9], [10], [6,7],
             [6,7,8,9], [6,7,8], [6],[7,8,9]
             ]

    p['E2b']=[[],
            [],
              [], [2,3],
              [2,3,4],[3],[3],[3],[3],
              [],
              [3], [], [2,3,4],[],[],
              [3], [2,3],[3], [2,3,4], [3],[4],
              [2,3],[3],[], [3]
              ]

    p['E2a']=[[],
            [],
             [], [],
             [], [], [], [], [],
             [],
             [7], [], [7], [],[],
             [], [5,6,7], [7], [], [], [],
             [7], [], [],[7]
             ]

    p['F1a']=[[4,5],
            [3,4,5],
             [], [],
             [3], [3], [2], [3], [4],
             [3,4],
             [2], [], [], [],[],
             [3], [], [], [], [3], [],
             [2], [], [],[]
             ]


    p['F1b']=[[9],
            [],
             [], [],
             [], [6,7,8], [], [], [6,7],
             [7],
             [7], [7], [], [],[7],
             [], [], [], [], [], [],
             [], [], [],[]
             ]

    p['F2']=[[9],
            [],
             [], [],
             [], [10], [], [9], [9],
             [9],
             [9], [9,10], [], [],[10],
             [10], [9,10], [9], [10], [10], [],
             [], [10], [9],[]
             ]

    p['F3']=[[8,9,10],
            [7,8],
             [], [],
             [], [7,8], [7,8], [7], [7],
             [7,8],
             [8], [7], [7,8,9], [9],[],
             [8,9], [8], [8], [9], [9], [9],
             [7,8,9,10], [9], [8,9],[8]
             ]





    # add frontal stations (based on density gradient)
    var.insert(0,'density grad')
    warm, front, cold, tr_number, f_id = def_fronts()
    f_st=dict(zip(f_id, front))
    for l in p.keys():
        p[l].insert(0,f_st[l])

    N=len(p['A']) #26 = 24 taxa + fluo + density


    ############ compute position and widths ############
    #station coordinates (km) [tr, st]
    a=fl.along()

    # transect resolution of each front
    fr_res = [round(np.nanmax(a[tr_number[i]])/n_st(tr_number[i], n_tr=8),1) for i in range(10)]
    fr_res=dict(zip(f_id, fr_res))


    pos=dict()
    nst=np.zeros((N)) #number of peak stations
    w=dict()


    for iif,flab in enumerate(f_id):
        print(iif, flab)
        tr=tr_number[iif]
        pos_vars=np.zeros((N))
        for ivar in range(N):
            nst[ivar]=len(p[flab][ivar]) # number of of peak stations
            if nst[ivar]>0: # there is a peak
                st_inds=np.array(p[flab][ivar])-1
                pos_vars[ivar]=np.nanmean(a[tr][st_inds])
            else: # transition
                pos_vars[ivar]=np.nan

        w[flab]=nst*fr_res[flab]
        pos[flab]=pos_vars

    return pos, w, p, var


###### COMPUTE THE PHYSICAL GRADIENTS

def compute_grads(do_abs, z):
    '''
    gradt, gradd, grads = fl.compute_grads(do_abs, z)
    computes the gradients of the vertically averaged density, temperature, salinity
    arguments:
        do_abs : if True, returns the absolute value (positive); is only used
        z : number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns :
        gradt, gradd, grads : [tr, st] if do_norm

    NB : np.gradient() documentation :
        The gradient is computed using second order accurate central differences in the interior points
        and either first or second order accurate one-sides (forward or backwards) differences at the boundaries.
        The returned gradient hence has the same shape as the input array.

    '''

    along=fl.along()
    #load vertically average variables from netcdf files
    tm=fl.load_temp(vave=z)         # [tr, st]
    dm=fl.load_dens(vave=z)
    sm=fl.load_sal(vave=z)

    gradt=np.full(along.shape, np.nan)
    gradd=np.full(along.shape, np.nan)
    grads=np.full(along.shape, np.nan)


    for tr in range(len(along)):
        tmtr=tm[tr]
        dmtr=dm[tr]
        smtr=sm[tr]
        atr=along[tr]
        ind_st=np.where(np.isfinite(tmtr))
        gradt[tr, ind_st]=np.gradient(tmtr[ind_st], atr[ind_st] )
        gradd[tr, ind_st]=np.gradient(dmtr[ind_st], atr[ind_st])
        grads[tr, ind_st]=np.gradient(smtr[ind_st], atr[ind_st])

    if do_abs:
        return abs(gradd), abs(gradt), abs(grads)
    else:
        return gradd, gradt, grads


def compute_norm_grads(z):

    '''
    computes the gradients of the vertically averaged density, temperature, salinity normalized by the transect average value

    arguments:
        z : number of vertical levels, equivalent to depth in m (in CTD casts, 1 level  1m)

    returns :
        rel_grad : normalized gradients, [tr, st, var] ; var = dens, temp, sal
        ave_grad : average gradients for each transect, [tr, st, var] ; var = dens, temp, sal

    NB : ave_grad is not actually the average of all the gradients but the average of the min and the max, to get an estimate of the gradient variation amplitude

    '''
    grad_list = compute_grads(do_abs=True, z=z)

    rel_grad=np.zeros((grad_list[0].shape[0], grad_list[0].shape[1], len(grad_list))) #[tr, st, var]
    ave_grad=np.zeros((grad_list[0].shape[0], len(grad_list))) # [tr, var]
    for ivar, data_var in enumerate(grad_list):
        for tr in range(8):
            ave_grad[tr, ivar]=(np.nanmin(abs(data_var[tr])+np.nanmax(abs(data_var[tr]))))/2 # average of min and max
            rel_grad[tr,:,ivar]=abs(data_var[tr])/ave_grad[tr, ivar]


    return rel_grad, ave_grad






###### computations based on the front definitions


def def_front_in_transect(numbers=False):
    '''
    returns the indices of warm/front/cold stations in each transect (to use as markers in lineplots)
        NB : organized by transect, contrary to def_fronts() where they are organized by front

    arguments :
        numbers : if True, returns the numbers ; if False, returns the indices
    returns:
        fr_tr,w_tr,c_tr : list of size [tr]
        trid : list of the names of the transects [tr]
    '''

    warm, front, cold, tr_number, name = def_fronts()
    ii=0 if numbers else 1

        ########### indices (not numbers !) of front stations (for lineplots) in transect order (excluding Up1)
    # station numbers of each category for each TRANSECT
    trid=transect_id()
    fr_no=[[] for i in range(len(trid))]
    w_no=[[] for i in range(len(trid))]
    c_no=[[] for i in range(len(trid))]
    for iif in range(len(name)):
        tr=tr_number[iif]
        fr_no[tr].extend(np.array(front[iif])-ii)
        w_no[tr].extend(np.array(warm[iif])-ii)
        c_no[tr].extend(np.array(cold[iif])-ii)

    return w_no, fr_no,c_no, trid



def core_front(var):
    '''
    core_front=core_front(var)
    returns a dictionary of the core position for each front or the along coordinates centered on the core ascending order (negative values = first stations on the left)
    for instance core_front(var='core')['A'] is a single value
    arguments :
        var = 'core' (returns a single value for each front) or 'coord' (returns a list of coords for each front)

    '''


    along=fl.along() #[tr, st]

    warm, front, cold, tr_number, f_id = def_fronts() # station numbers

    core=[]
    core_front=[[] for i in range(len(f_id))]

    for iif, flab in enumerate(f_id):
        tr=tr_number[iif]
        x=np.nanmean(along[tr][np.array(front[iif])-1])
        core.append(x)

        # compute the coordinates centered on the front
        core_front[iif]=along[tr] - x
        # if flab in ['E2a','F1b', 'F2', 'C3']: # reverse the sign to put the cold side on the left
        #     core_front[iif] = - core_front[iif]

    if var=='core':
        return dict(zip(f_id,core))
    elif var=='coord':
        return dict(zip(f_id,core_front))




def peak_intensity(data, maxF):
    '''
    computes the intensity of the biomass peaks (as percentage variation compared to the background biomass). also returns a mask indicating transitions

    arguments :
        data = [tr, st, var]
        maxF : if True, computes the max in the front (aka the biomass peak intensity). if False, computes the average

    returns :
        X = [nf, var]
        mask_trans = [nf, var] nans in peaks, 1 in transitions


    '''
    nvar=data.shape[2]


    warm, front, cold, tr_number, f_id = def_fronts()
    nf=len(f_id)


    X=np.zeros((nf,nvar))  # [fronts, var]
    mask_trans=np.zeros((nf,nvar))  # [fronts, var]
    mask_trans.fill(np.nan)


    for iif in range(nf):
        tr=tr_number[iif]

        wi=np.array(warm[iif])-1
        fi=np.array(front[iif])-1
        ci=np.array(cold[iif])-1


        for ivar in range(nvar): # for each variable

            if maxF: # detect a min or max in each zone
                fwc=np.concatenate([fi,wi,ci])
                x=data[tr, fwc, ivar]
                ind_max=fwc[list(x).index(np.nanmax(x))]
                ind_min=fwc[list(x).index(np.nanmin(x))]
                if ind_max in fi : # positive value : compute amplitude of the peak (compared to minn or average background)
                    # b=np.nanmin(data_i[tr, np.append(wi,ci), ivar]) # --> (max(front) - min(background))/min(background)
                    b=np.nanmean(data[tr, np.append(wi,ci), ivar]) # --> (max(front) - mean(background))/mean(background)
                    X[iif, ivar]= (np.nanmax(data[tr, fi, ivar])- b )/ b


                elif ind_min in fi:# negative value : compute amplitude of peak (compared to min or average background)
                    # b= np.nanmax(data_i[tr, np.append(wi,ci), ivar])  # --> (min(front) - max(background))/max(background)
                    b= np.nanmean(data[tr, np.append(wi,ci), ivar]) # --> (min(front) - mean(background))/mean(background)
                    X[iif, ivar]= (np.nanmin(data[tr, fi, ivar]) -b) / b

                else:
                    mask_trans[iif, ivar]=1

            else: # compare averages in each zone
                dw=np.nan_to_num(np.nanmean(data[tr, fi, ivar])) - np.nan_to_num(np.nanmean(data[tr, wi, ivar]))  # front - warm
                dc=np.nan_to_num(np.nanmean(data[tr, fi, ivar])) - np.nan_to_num(np.nanmean(data[tr, ci, ivar]))  # front - cold


                if dw*dc < 0: # dw and dc have opposite signs = front is a transition
                    mask_trans[iif, ivar]=1
                else: # # dw and dc have the same sign = front is a max or min
                    X[iif, ivar]= (np.nanmean(data[tr, fi, ivar]) - np.nanmean(data[tr, np.append(wi,ci), ivar]))/ np.nanmean(data[tr, np.append(wi,ci), ivar])


    return X*100, mask_trans







########################### LOAD DATA

def trophic(Afront_ave=False):
    '''
    out, group_names= f.trophic()

    loads plankton data and computes the total carbon biomass (µgC/m3) for each of the 7 functional groups :
        h-bact + cyanobact                       : flowcyt; vertical integral
        other euk phy + diatoms                  : HPLC ; surface values (or vertical average for A-front)
        pico-grazers + meso-grazers + carnivores : Zooscan; vertical integral


    arguments:
        Afront_ave : if True, uses vertical average for A-front; if False, uses surface values

    returns :
        out= np.array [tr, st, var]
        group_names = name of each trophic level

    NB : to check if there is data : np.any(np.isnan(data_zs), axis=2)
    '''

    group_names = ['h-bact','cyanobact','other euk-phyto','diatoms','pico-grazers','meso-grazers','carnivores']

    #flow cyt
    data_fc, var_fc, data_name= fl.load_flowcyt(vave=True, biomass=True, rel=False)#[transect, station, var] ; (µgC/m3)

    ig0=[0]  #heterotrophic bacteria
    ig1=[1,2]#PRO+SYN

    g0=np.squeeze(data_fc[:,:,ig0])
    g1=np.sum(data_fc[:,:,ig1],axis=2)


    #HPLC
    data_hplc, var_hplc, data_name = fl.load_hplc(rel=False, Afront_ave=Afront_ave, biomass=True)#[transect, station, var] ; (µgC/m3)

    ig2=[1,3,4,6,7]#other euk-phyto
    ig3=[2] #diatoms

    g2=np.nansum(data_hplc[:,:,ig2],axis=2)
    g3=data_hplc[:,:,ig3].squeeze()

    #zooscan
    data_zs, var_zs, data_name = fl.load_zooscan(biomass=True, rel=False) # [tr, st, var] ;μgC/m^3

    ig4 = [var_zs.index('rhizaria'),var_zs.index('doliolids') ,var_zs.index('appendicularia'),var_zs.index('pyrosomes'),var_zs.index('salps')]
    ig5 = [var_zs.index('Calanoids'),var_zs.index('copepoda_oithona_like') ,var_zs.index('Other copepods (Poec. + Harp.)'),var_zs.index('euphausiids'),var_zs.index('pteropoda'),var_zs.index('crustacea_others')]
    ig6 = [var_zs.index('chaetognatha'),var_zs.index('cnidaria_ctenophores') ,var_zs.index('ostracods') ,var_zs.index('polychaete') ]


    g4=np.sum(data_zs[:,:,ig4],axis=2)
    g5=np.sum(data_zs[:,:,ig5],axis=2)
    g6=np.sum(data_zs[:,:,ig6],axis=2)

    out=np.stack([g0, g1, g2, g3, g4, g5, g6], axis=-1)

    return out, group_names



def trophic_detail(Afront_ave=False):
    '''
    out, taxa_names, group_names= f.trophic_detail()

    loads all the plankton data and organizes them by functional groups :
        h-bact + cyanobact                       : flowcyt; vertical integral
        other euk phy + diatoms                  : HPLC ; surface values (or vertical average for A-front)
        pico-grazers + meso-grazers + carnivores : Zooscan; vertical integral


    arguments:
        Afront_ave : if True, uses vertical average for A-front; if False, uses surface values

    returns :
        out= list of np.arrays [tr, st, var]
        taxa_names= list of lists with the taxa names
        group_names= name of each functional group

    '''

    group_names = ['h-bact','cyanobact','other euk-phyto','diatoms','pico-grazers','meso-grazers','carnivores']

    #flow cyt
    data_fc, var_fc, data_name= fl.load_flowcyt(vave=True, biomass=True, rel=False)#[transect, station, var] ; (µgC/m3)

    ig0=[0]  #heterotrophic bacteria
    ig1=[1,2]#PRO+SYN

    g0=np.squeeze(data_fc[:,:,ig0])
    g1=data_fc[:,:,ig1]

    #HPLC
    data_hplc, var_hplc, data_name = fl.load_hplc(rel=False, Afront_ave=Afront_ave, biomass=True)#[transect, station, var] ; (µgC/m3)

    ig2=[1,3,4,6,7]#other euk-phyto
    ig3=[2] #diatoms

    g2=data_hplc[:,:,ig2]
    g3=data_hplc[:,:,ig3].squeeze()

    #zooscan
    data_zs, var_zs, data_name = fl.load_zooscan(biomass=True, rel=False) # [tr, st, var] ;μgC/m^3

    ig4 = [var_zs.index('rhizaria'),var_zs.index('doliolids') ,var_zs.index('appendicularia'),var_zs.index('pyrosomes'),var_zs.index('salps')]
    ig5 = [var_zs.index('Calanoids'),var_zs.index('copepoda_oithona_like') ,var_zs.index('Other copepods (Poec. + Harp.)'),var_zs.index('euphausiids'),var_zs.index('pteropoda'),var_zs.index('crustacea_others')]
    ig6 = [var_zs.index('chaetognatha'),var_zs.index('cnidaria_ctenophores') ,var_zs.index('ostracods') ,var_zs.index('polychaete') ]

    g4=data_zs[:,:,ig4]
    g5=data_zs[:,:,ig5]
    g6=data_zs[:,:,ig6]


    out=[g0, g1, g2, g3, g4, g5, g6]

    # because no fancy indexing in lists
    var_fc=np.array(var_fc)
    var_hplc=np.array(var_hplc)
    var_zs=np.array(var_zs)

    taxa_names=[list(var_fc[ig0]), list(var_fc[ig1]), list(var_hplc[ig2]), list(var_hplc[ig3]), list(var_zs[ig4]), list(var_zs[ig5]), list(var_zs[ig6])]

    return out, taxa_names, group_names


def load_naup():
    '''
    loads the egg and nauplii abundances and computes the ratio to adults
    eggs = all copepods + euphausiids
    nauplii = all copepods

    '''

    data_zs, var_zs, data_name = fl.load_zooscan(biomass=False, rel=False) # [tr, st, var] ;#/m^2

    data_cop=data_zs[:,:,var_zs.index('All Copepods')]
    data_eggs = data_zs[:,:,var_zs.index('eggs')]
    data_nauplii = data_zs[:,:,var_zs.index('nauplii')]
    data_euph = data_zs[:,:,var_zs.index('euphausiids')]
    r_eggs=data_eggs / (data_cop + data_euph)
    r_naup=data_nauplii / data_cop

    out=np.stack([r_eggs, r_naup], axis=-1) # [tr, st, var]

    return out







################### PLOTS


def layout_stack(nrows, lA, hA, org, fs=13):
    '''
    fig_list, axes_list, tr_list, along_list, f_id = layout_stack(nrows, lA, hA, org, fs=13)

    prepares the layout of the figures for stack plots
    --> 4 figures (1 for each cruise); in each figure 1 column for each transect (or each front)

        all the x-axis scales (the distance along the transects) will be the same for all transects, even in separate figures
        some x-axis are reversed so that all the warm sides are on the left, and E1 and E2 have the same orientation


    arguments :
        nrows = number of rows
        lA= length of A-front figure. The length of the other figures are determined automatically to keep the same x-axis scale
        hA= height of figures
        org= which layout. 1 = full transects; 2 = centered on each front
        fs= fontsize of axes labels

    returns :
        fig_list = list of 4 fig handles [cr]
        axes_list = list of 4 arrays with the axes handles [cr][nrows, ncol]
        tr_list = list of 4 lists with the index of the transect in each figure [cr][ncol]
        along_list = list of 4 arrays with the x-axis for each column [cr][ncol]


    NB : the only things that should be changed are the paddings (hpad, padv1 and padv2) and the xlims for org=2. EVERYTHING ELSE is automatic
        pad_xlab and pad_stat can also be adjusted if needed
    '''


    padv1=0.6*cm/hA # vertical padding at top and bottom of figure (portion of fig size). NB : the numerical value (0.4) is in cm
    padv2=0.4*cm/hA # vertical padding between subplots (portion of fig size). NB : the numerical value (0.1) is in cm


    pad_xlab=padv1*hA*0.8 # vertical pad for the x label; in inches; 0.2 = fraction of bottom vertical pad


    fig_list=[[] for i in range(4)]
    axes_list=[[] for i in range(4)]
    along_list=[[] for i in range(4)]

    xlab='Station number'

    if org==1: # full transects
        tr_list=[[0],[1,2],[3,4],[5,6,7]] # indices of the transects in each figure
        f_id=[['A'],['C2', 'C3'], ['E1', 'E2'], ['F1', 'F2', 'F3']]

        # lengths of each plot in axis coordinates (km)

        along= fl.along()  # [tr, st]
        a=np.nanmax(along[0])-np.nanmin(along[0])
        c2,c3=np.nanmax(along[1])-np.nanmin(along[1]), np.nanmax(along[2])-np.nanmin(along[2])
        e1,e2=np.nanmax(along[3])-np.nanmin(along[3]), np.nanmax(along[4])-np.nanmin(along[4])
        f1,f2,f3=np.nanmax(along[5])-np.nanmin(along[5]), np.nanmax(along[6])-np.nanmin(along[6]),  np.nanmax(along[7])-np.nanmin(along[7])


        lkm=[[a], [c2,c3], [e1,e2], [f1,f2,f3]] # length of each subplot in data coord (km)


        along_list=[[along[tr] for tr in tri] for tri in tr_list ]

        flipx=['C3','E2', 'F2']


        hpad=0.5*cm*lA # horizontal pad between subplots in inches. NB : the numerical value (0.5) is in cm


    elif org==2:
        tr_list=[[0],[1,2],[3,4,4],[5,5,6,7]] # indices of the transects in each subplot
        f_id=[['A'],['C2','C3'], ['E1', 'E2a', 'E2b'],['F1a', 'F1b', 'F2', 'F3']]

        xl=[[9.3], [10, 10], [15, 13, 13], [10, 10,15,20]] # xlims, plots are symmetrical around 0

        a, c2, c3=xl[0][0]*2, xl[1][0]*2, xl[1][1]*2
        e1, e2a, e2b=xl[2][0]*2, xl[2][1]*2, xl[2][2]*2
        f1a, f1b, f2, f3 = xl[3][0]*2, xl[3][1]*2, xl[3][2]*2 , xl[3][3]*2

        lkm=[[a], [c2,c3], [e1,e2a, e2b], [f1a, f1b,f2,f3]] # length of each subplot in data coord (km)

        along=core_front(var='coord') # dict
        along_list=[[along[flab] for flab in tri] for tri in f_id ]


        flipx=['C3','E2a','F1b', 'F2']

        hpad=1.2*cm*lA # horizontal pad between subplots in inches. NB : the numerical value (0.5) is in cm


    # common to both org=1 and org=2
    X=0.9*lA/a # factor of conversion between data coordinates (km) and axis length (inches)

    lfig=[(len(x)+1)*hpad + X*np.sum(x) for x in lkm]# length of each figure (inches)
    h=(1-2*padv1 - (nrows-1)*padv2)/nrows # height of each subplot (portion of fig height)

    bottoms=[padv1+(nrows-irow)*padv2 + (nrows-irow+1)*h for irow in range(nrows)] # bottoms of each subplot (proportion of fig height)


    for ifig in range(4): # for each figure

        fig=plt.figure(figsize=(lfig[ifig], hA))

        fig_list[ifig]=fig


        ncol=len(tr_list[ifig])
        axes=[[[] for i in range(ncol)] for j in range(nrows)]


        w=[X*L/lfig[ifig] for L in lkm[ifig]]# width of the subplots (proportion of fig width)


        for icol in range(ncol):
            for irow in range(nrows):


                left=(icol+1)*hpad/lfig[ifig] + np.sum(w[:icol]) # left position of the subplot (proportion of fig width)

                ax=plt.axes([left, bottoms[irow], w[icol], h]) # [left, bottom, width, height]

                axes[irow][icol] = ax

                tr=tr_list[ifig][icol]

                # don't display the xticks (the station ticks will be enough)
                ax.set_xticks([])
                ax.set_xticklabels([])


                ###############
                if org==1:
                    ax.set_xlim(np.nanmin(along[tr]), np.nanmax(along[tr]))

                elif org==2:
                    ax.set_xlim(-xl[ifig][icol], xl[ifig][icol])

                    ax.axvline(0, linestyle='-', color='k', linewidth=0.5) # vertical line at core of front

                ###############
                pad_stat = pad_xlab*0.99/h/hA  # distance from axis to station labels, fraction of subplot height;

                if irow == nrows-1: # bottom row
                    ax.set_xlabel(xlab, fontsize=fs, labelpad=pad_xlab*72) # pad must be given in points; pad_xlab is in inches;
                    # ax.tick_params(axis='x', pad=pad_xlab*72) #  useless since the tick display is off
                    draw_stations(f_id=f_id[ifig][icol], ax=ax, vpad=pad_stat, fs=fs-1, ticks_only=False) #
                else:  #don't display the station numbers (only the station ticks)
                    draw_stations(f_id=f_id[ifig][icol], ax=ax, vpad=pad_stat, fs=fs-1, ticks_only=True)

                if f_id[ifig][icol] in flipx:
                    ax.invert_xaxis()


                axes_list[ifig]=axes
                fig.subplots_adjust(top=bottoms[0]+h, bottom=padv1)


    return fig_list, axes_list, tr_list, along_list, f_id





def draw_stations(f_id, ax, ticks_only=False, vpad=0.15, fs=13, tcolor='k'):
    '''
    draws vertical ticks at the bottom of the plot to show the location of the stations
    can also write labels with the station numbers in red, blue and black

    arguments :
        f_id= str with the id of the transect or the front
        ax : ax handle (must already have a x-axis in km)
        ticks_only : if True, only draws the ticks at each station but not the numbers
        vpad : vertical padding between the station numbers and the bottom of the axis (fraction of subplot height)
        fs : fontsize
        tcolor : color of the vertical line
    '''



    core_flag= ax.get_xlim()[0]<0


    # pos = where to put the ticks (in km, coordinates along the transect or centered on the front)
    # d= move the numbers horizontally so they are centered on the ticks rather than start on their edge (in transect or front chronological order)
            # negative values when they were flipped

    if core_flag:
        warm, front, cold, tr_number, ff_id=def_fronts() # numbers of front stations, organized by front
        ii=ff_id.index(f_id) # front index, to use in front warm cold
        tr=tr_number[ii] # for the station numbers
        pos=core_front(var='coord')[f_id]
        d=[0.55,0.8,-0.8,0.6,-1.2,1.2,1,-1,-0.5,1]



    else:
        tr=transect_id(short=True).index(f_id)
        warm, front, cold, tr_number = def_front_in_transect(numbers=True) # numbers of front stations, organized by transect
        pos=fl.along()[tr]
        d=[0.55,0.8,-0.8,0.6,-1.2,1,-0.5,1]

        ii=tr # transect index, to use in front warm cold


    num_s=np.arange(n_st(tr=tr,n_tr=8))+1   #list of all station numbers


    ytext=- vpad # in axis coord

    xtext=pos-d[ii] # x-position of the labels, in data coord

    for s,  sn in enumerate(num_s): # s = station index; sn=station number
        if np.min(ax.set_xlim())<=pos[s]<=np.max(ax.set_xlim()): # only draw stations between the xlims
            ax.axvline(pos[s], color=tcolor, linestyle='--', linewidth=0.5, alpha=0.7, ymin=-0.2, ymax=0.07) # dashed line indicating the position of the station, ymin and ymax are in axis coord

            if not ticks_only: # write the station numbers
                if sn in front[ii]: # front = bold
                    ax.text(xtext[s], ytext, str(sn), fontsize=fs+1, alpha=1, fontweight='bold', transform=ax.get_xaxis_transform()) # fronts in bold
                elif sn in warm[ii]: # warm = red
                    ax.text(xtext[s], ytext, str(sn), fontsize=fs, alpha=1, color='red', transform=ax.get_xaxis_transform()) # fronts in bold
                elif sn in cold[ii]: # cold = blue
                    ax.text(xtext[s], ytext, str(sn), fontsize=fs, alpha=1, color='blue', transform=ax.get_xaxis_transform()) # fronts in bold

                else:
                    ax.text(xtext[s], ytext, str(num_s[s]), fontsize=fs, alpha=1, transform=ax.get_xaxis_transform()) # other stations


def front_colors():
    '''
    colors associated with each front (for lineplots)
    '''
    colors=['red', 'orange', 'chocolate', 'gold', 'palegreen', 'green', 'lime', 'turquoise', 'cyan', 'deepskyblue', 'royalblue', 'darkviolet', 'violet']
    name = ['A'       , 'E1'                 , 'E2a'   , 'E2b' , 'C2'       ,'C3'   , 'F1a' ,'F1b'  , 'F2'   , 'F3' ]

    return dict(zip(name, colors))

############# other usefull functions ###################

def interp_z(data, zz, step_z=20, method='linear'):
    '''
    interpolates data on a regular vertical grid (useful for averaging vertical profiles for instance)

    arguments :
        data = [tr, st, z] or [tr, st, z, var]
        zz = [tr, st, z]
        step_z = interpolation step (in m)
        method : interpolation method, see values in "griddata" function
    returns :
        out = interpolated data, same dimensions as input
        grid_z = interpolated grid [z]

    '''

    grid_z=-np.arange(0, -np.nanmin(zz)+step_z, step_z)
    sh=list(data.shape)
    sh[2]=len(grid_z)
    out=np.zeros(sh)
    out_near=np.zeros(sh)

    for tr in range(data.shape[0]):
        for st in range(data.shape[1]):
            if len(data.shape)==4: # multiple variables
                for ivar in range(data.shape[3]):
                    out[tr, st, :,ivar]=griddata(zz[tr,st], data[tr,st,:, ivar], grid_z, method='linear')
                    out_near[tr, st, :,ivar]=griddata(zz[tr,st], data[tr,st,:, ivar], grid_z, method='nearest')
            else: # 1 variable
                out[tr, st]=griddata(zz[tr,st], data[tr,st,:], grid_z, method='linear')
                out_near[tr, st]=griddata(zz[tr,st], data[tr,st,:], grid_z, method='nearest')

    out=np.where(np.isnan(out), out_near, out) # fill the edges with nearest neighbor
    return out, grid_z


def flat(x):
    '''
    flattens a list of lists or an array of list into a single list
    '''

    return [val for sublist in x for val in sublist]

