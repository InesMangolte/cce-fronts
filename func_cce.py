'''
functions for various computations and definitions
'''

import numpy as np


import func_load_cce as fl




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
        n_tr = total number of transects (9 = all ; 8 = upfront1 is removed)
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
    warm =  [[1,2]    , [1]              , [10]          , [1,2,3,4]      , [9,10]   , [2]    , [1]        , [11]        , [11]    , [1,2,3,4] ]
    front = [[3,4,5]  , [2,3,4,5,6,7,8,9], [7,8,9]       , [5,6,7,8,9,10] , [6,7,8]  , [3]    , [2,3]      , [9,10]      , [9,10]  , [5,6,7, 8, 9,10]]
    cold =  [[6,7,8]  , [10]             , [1,2,3,4,5,6] , [11,12,13]     , [4,5]    , [4,5]  , [4,5,6,7,8], [4,5,6,7,8] , [7,8]   , [11] ]
    name = ['A'       , 'C2'             ,'C3'           , 'E1'           , 'E2a'    , 'E2b'  , 'F1a'      ,'F1b'        , 'F2'    , 'F3' ]
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
    data_detail, names, var_names, transect_id = trophic_detail(carb=True, cyan='fc', Afront_ave=True)

    var=[val for sublist in var_names for val in sublist]
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
    warm, front, cold, tr_number, f_id = def_fronts(order=False)
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


###### computations based on the front definitions


def def_front_in_transect():
    '''
    returns the indices of front/warm/cold stations in each transect (to use as markers in lineplots)

    returns:
        fr_tr,w_tr,c_tr : list of station indices in each transect
        trid : list of the names of the transects
    '''

    warm, front, cold, tr_number, name = def_fronts()

        ########### indices (not numbers !) of front stations (for lineplots) in transect order (excluding Up1)
    # station numbers of each category for each TRANSECT
    trid=transect_id()
    fr_no=[[] for i in range(len(trid))]
    w_no=[[] for i in range(len(trid))]
    c_no=[[] for i in range(len(trid))]
    for iif in range(len(name)):
        tr=tr_number[iif]
        fr_no[tr].extend(np.array(front[iif])-1)
        w_no[tr].extend(np.array(warm[iif])-1)
        c_no[tr].extend(np.array(cold[iif])-1)

    return fr_no,w_no,c_no, trid



def PEF(data):
    '''
    computes the peak enhancement factor

    arguments :
        data = [tr, st, var]

    returns :
        X = [nf, var] or [nfp, var] (all fronts or peaks only)
        mask_trans= [nf, var] nans in peaks, 1 in transitions



    '''

    warm, front, cold, tr_number, f_id = def_fronts()

    nf=len(f_id)
    nvar=data.shape[2]


    X=np.zeros((nf,nvar))  # fronts, var
    mask_trans=np.zeros((nf,nvar))  #fronts, var
    mask_trans.fill(np.nan)


    for iif in range(nf):
        tr=tr_number[iif]

        wi=np.array(warm[iif])-1
        fi=np.array(front[iif])-1
        ci=np.array(cold[iif])-1

        fwc=np.concatenate([fi,wi,ci])

        for ivar in range(nvar): # for each variable

            x=data[tr, fwc, ivar]

            ind_max=fwc[list(x).index(np.nanmax(x))]
            ind_min=fwc[list(x).index(np.nanmin(x))]

            if ind_max in fi : #(max(front) - mean(background))/mean(background)
                b=np.nanmean(data[tr, np.append(wi,ci), ivar])
                X[iif, ivar]= (np.nanmax(data[tr, fi, ivar])- b )/ b


            elif ind_min in fi: #(min(front) - mean(background))/mean(background)
                b= np.nanmean(data[tr, np.append(wi,ci), ivar])
                X[iif, ivar]= (np.nanmin(data[tr, fi, ivar]) -b) / b

            else: # transition
                mask_trans[iif, ivar]=1


    return X, mask_trans, f_id




def core_front():
    '''
    returns a dictionary with the position of the center of each front

    '''

    along=fl.along() #station coordinates (km) [tr, st]

    warm, front, cold, tr_number, f_id = def_fronts() # station numbers

    core=[]

    for iif, flab in enumerate(f_id):
        tr=tr_number[iif]
        x=np.nanmean(along[tr][np.array(front[iif])-1])
        core.append(x)

    return dict(zip(f_id,core))

def core_coord_front(var):
    '''
    core_front=core_front(var)
    returns a dictionary with the along coordinates centered on the core

    '''

    along=fl.along() #station coordinates (km) [tr, st]

    warm, front, cold, tr_number, f_id = def_fronts() # station numbers

    core=core_front()
    core_coord=[[] for i in range(len(f_id))]

    for iif, flab in enumerate(f_id):
        tr=tr_number[iif]
        core_coord[iif]=along[tr] - core[flab]
        if flab in ['E2a','F1b', 'F2', 'C3']: # reverse the sign to put the cold side on the left
            core_coord[iif] = - core_coord[iif]

    return dict(zip(f_id,core_front))


############# other usefull functions ###################

def flat(x):
    '''
    flattens a list of lists or an array of list into a single list
    '''

    return [val for sublist in x for val in sublist]


def spiciness1(theta, sal):
    '''
    computes the spiciness from Jackett, D. R., and T. J. McDougall, 1985: An oceanographic variable for the characterization of intrusions and water masses. Deep-Sea Res., 32, 1195–1207.
    for 1 point
    arguments :
        theta = potential temperature (°C)
        sal = salinity (PSU)
    '''
    A = np.array([[1.609705E-1 , 6.542397E-1, 5.222258E-4, -2.586742E-5, 7.565157E-7],
                  [-8.007345E-2, 5.309506E-3, -9.612388E-5, 3.211527E-6, -4.610513E-8],
                  [1.081912E-2, -1.561608E-4,3.774240E-6, -1.150394E-7, 1.146084E-9],
                  [-1.451748E-4, 3.485063E-6,-1.387056E-7, 3.737360E-9, -2.967108E-11],
                  [1.219904E-6, -3.591075E-8,1.953475E-9,-5.279546E-11, 4.227375E-13]])

    if not np.isfinite(theta) or not np.isfinite(sal):
        return np.nan

    theta_power=np.array([[np.power(theta, i) for j in range(5)] for i in range(5)])
    sal_power=np.array([[np.power(sal, j) for j in range(5)] for i in range(5)])

    prod=np.multiply(np.multiply(A, theta_power), sal_power)

    return np.nansum(prod)

def spiciness(theta, sal):
    '''
    computes the spiciness from Jackett, D. R., and T. J. McDougall, 1985: An oceanographic variable for the characterization of intrusions and water masses. Deep-Sea Res., 32, 1195–1207.
    for arrays of up to 3 dimensions
    arguments :
        theta = potential temperature (°C)
        sal = salinity (PSU)
    '''
    theta=np.squeeze(theta)
    sal=np.squeeze(sal)
    out=np.zeros(theta.shape)
    dims=len(theta.shape)
    if dims==0:
        out=spiciness1(theta, sal)
    elif dims==1:
        for ii in range(theta.shape[0]):
            out[ii]=spiciness1(theta[ii], sal[ii])

    elif dims==2:
        for ii in range(theta.shape[0]):
            for jj in range(theta.shape[1]):
                out[ii,jj]=spiciness1(theta[ii,jj], sal[ii,jj])

    elif dims==3:
        for ii in range(theta.shape[0]):
            for jj in range(theta.shape[1]):
                for kk in range(theta.shape[2]):
                    out[ii,jj,kk]=spiciness1(theta[ii,jj,kk], sal[ii,jj,kk])
    else:
        print('too many dimensions')
        return

    return out

