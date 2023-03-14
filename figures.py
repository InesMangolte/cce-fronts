"""

the intermediate figures are saved in path_fig.
path_fig_bg contains the final figures (assembled with Illustrator of directly from python in a couple of cases)
the individual images for the videos are saved in path_vid (in separate folders for each cruise) and must then be assembled

"""

import sys

sys.path.append("./")
import func_cce as f
import func_load as fl
import func_load_sat as fl_sat
import func_sat as f_sat

import matplotlib.pyplot as plt
import numpy as np


import datetime as dt
from matplotlib.patches import Rectangle
from matplotlib import colors as mpl_colors



path_fig_bg = './figures/'
path_fig    = path_fig_bg+ '/from_python/' # raw figures, then processed with Illustrator to get the final images
path_vid    = path_fig_bg+ '/videos/' # all the individual images, which must then be assembled into videos externally

cm=1/2.54  # centimeters to inches # pour la conversion de la taille des figures

#%% fig 1 : SST and CHL maps climatology (COPERNICUS 4km daily)


dates=f_sat.date_range(start_date= dt.date(2008, 1, 1), end_date = dt.date(2017, 1, 31))

# sort the filenames by month
month_list=[[] for i in range(12)]
for it, date in enumerate(dates):
    m=int(date.month)-1 # index of the month
    month_list[m].append(date)


# load each month, do average
data_sst=[[] for i in range(12)]
data_chl=[[] for i in range(12)]

# for m in range(12): # full climatology
for m in [6]: # just do the computation for a single month = July
    X, lat, lon=fl_sat.load_submesocolor(dataname='sst', dates=month_list[m]) # COPERNICUS
    data_sst[m] = np.nanmean(X, axis=0)

    X, lat, lon=fl_sat.load_submesocolor(dataname='chl_l3', dates=month_list[m]) # COPERNICUS
    data_chl[m] = np.nanmean(X, axis=0)




#%%

# PLOT
fs=9 # size of axes labels, colorbar labels, subplot titles and fig title

ncont=20
xlim, ylim=f_sat.lim_zoom(zoom=1) # regional CCE

stations, var, tr_id, DT=fl.coord()
warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations


data=[data_sst, data_chl]
var_list=['SST (°C)', 'Chlorophyll-a (mg/m3)']
cmaps=['plasma','viridis']

clim_sst=[12,20]
cti_sst=[12, 14, 16, 18, 20]
levs_sst=np.linspace(clim_sst[0],clim_sst[1], ncont)

clim_chl=[0.03,4]
levs_chl=f_sat.logclim(clim_chl, ncont)
cti_chl=[0.03,0.05,0.1,0.2,0.5,1, 2, 4]


clim=[clim_sst, clim_chl]
cti=[cti_sst, cti_chl]
levs=[levs_sst, levs_chl]


m=6 # JULY

fig, axes=plt.subplots(1,2, figsize=(12*cm, 6*cm))
# fig.suptitle('2007-2018 climatology, july', fontsize=fs)

for ivar, var in enumerate(var_list):
    ax=axes[ivar]
    ax.set_title(var, fontsize=fs)
    cf=ax.contourf(lon, lat,data[ivar][m],levs[ivar] , extend='both', cmap=cmaps[ivar])
    cf.set_clim(clim[ivar])
    f_sat.remove_white_lines(cf)


    f_sat.draw_landmask(ax=ax, color='k')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect(aspect=1)
    ax.tick_params(axis='both', labelsize=fs-1)


    cb=fig.colorbar(cf, ax=ax, shrink=0.5)
    cb.set_ticks(cti[ivar])
    cb.set_ticklabels(cti[ivar], fontsize=fs)


    # draw the stations
    for tr in range(8):
        for st in range(stations.shape[1]): # for each station
            if st in fr[tr]:
                ax.plot(stations[tr,st,1],stations[tr,st,0],'x', color='red',markersize=3) # front = red cross
            elif st in warm[tr] or st in cold[tr]:
                ax.plot(stations[tr,st,1],stations[tr,st,0],'o', color='black',markersize=2) # back = black circle
            else: # other
                ax.plot(stations[tr,st,1],stations[tr,st,0], 'o', color='black', markersize=1) # other = small black circle

fig.savefig(path_fig+'clim_CCEpy.pdf')




#%% fig 1 inlays : North Pacific climatology with rectangle around CCE

start_date = dt.date(2008, 1, 1)
end_date = dt.date(2017, 1, 31)
delta = dt.timedelta(days=1)

dates = []
while start_date <= end_date:
    dates.append(start_date)
    start_date += delta

X, lat, lon=fl_sat.load_submesocolor(dataname='sst', dates=dates) # COPERNICUS
data_sst = np.nanmean(X, axis=0)


X, lat, lon=fl_sat.load_submesocolor(dataname='chl_l3', dates=dates) # COPERNICUS
data_chl = np.nanmean(X, axis=0)



#% PLOT
fs=9 # size of axes labels, colorbar labels, subplot titles and fig title

ncont=20
xlim, ylim= [-133,-110],[20,50] # North Pacific


data=[data_sst, data_chl]
var_list=['SST', 'CHL']
cmaps=['plasma','viridis']


clim_sst=[5,30]
cti_sst=[5,10,15,20,25,30]
levs_sst=np.linspace(clim_sst[0],clim_sst[1], ncont)

clim_chl=[1E-2,1.5]
levs_chl=f_sat.logclim(clim_chl, ncont)
cti_chl=[0.01,0.05,0.1,0.5,1, 1.5]



clim=[clim_sst, clim_chl]
cti=[cti_sst, cti_chl]
levs=[levs_sst, levs_chl]


fig, axes=plt.subplots(1,2, figsize=(12*cm, 6*cm))
# fig.suptitle('2007-2018 Annual average')

for ivar, var in enumerate(var_list):
    ax=axes[ivar]
    # ax.set_title(var)
    cf=ax.contourf(lon, lat,data[ivar],levs[ivar] , extend='both', cmap=cmaps[ivar])
    cf.set_clim(clim[ivar])
    f_sat.remove_white_lines(cf)


    f_sat.draw_landmask(ax=ax, color='k')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect(aspect=1)
    # don't display the ticks
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])



    cb=fig.colorbar(cf, ax=ax,shrink=0.9)
    cb.set_ticks(cti[ivar])
    cb.set_ticklabels(cti[ivar], fontsize=fs)



    # rectangle for CCE
    lon_r, lat_r=f_sat.lim_zoom(zoom=1)
    rect_CCE = Rectangle((lon_r[0], lat_r[0]), np.diff(lon_r), np.diff(lat_r), color='k',linewidth=1, fill=False)
    ax.add_patch(rect_CCE)


    fig.savefig(path_fig+'clim_North_Pacificpy.pdf')


#%% fig 2 : schéma stations

along=True
fs=9

warm, front, cold, tr_number, f_id = f.def_fronts() # stations numbers
nf=len(f_id)
n_st=f.n_st(tr=None,n_tr=8)

along_coord=fl.along()  #distance along transect (km) [tr, st].
if along:
    x=along_coord
else:
    x=np.tile(np.arange(0,along_coord.shape[1]),along_coord.shape[0]).reshape(along_coord.shape)


# yy=np.arange(nf,0,-1) # even vertical spacing
yy=[10,  9, 8, 7, 5.8, 5.2, 3.8,  3.2,  2,  1] # As and Bs closer
y=np.repeat(yy,along_coord.shape[1]).reshape([nf,along_coord.shape[1]]) # 10 fronts, 13 stations

# cruise chronological order. Blank stations from cutting the filament transects in 2
blank=[[],[], [], [], [1,2,3,4], [6,7,8,9,10], [7,8,9,10,11],[1,2,3,4,5], [],[]]


fig, ax=plt.subplots(1,1, figsize=(12*cm,10*cm))
ax.tick_params(axis='both', labelsize=fs)
ax.set_xlabel('Distance along transect (km)', fontsize=fs)

# ax.axis('off')
ax.set_yticks(y[:,0])
ax.set_yticklabels(f_id)
ax.set_ylim(0.2,10.8)
# ax.set_xlabel('distance along transect (km)')

for fr in range(len(f_id)):
    tr=tr_number[fr]

    ax.plot(x[tr, np.arange(n_st[tr])],y[fr, np.arange(n_st[tr])], '.', color='k', markersize=2)

    l1,=ax.plot(x[tr, np.array(front[fr])-1],y[fr, np.array(front[fr])-1], 'x', color='k', markersize=3, label='frontal station')
    l2,=ax.plot(x[tr, np.array(warm[fr])-1],y[fr, np.array(warm[fr])-1], 'o', color='red', markersize=3, label='warm background station')
    l3,=ax.plot(x[tr, np.array(cold[fr])-1],y[fr, np.array(cold[fr])-1], 'o', color='blue', markersize=3, label='cold background station')
    if len(blank[fr])>0:
        ax.plot(x[tr, np.array(blank[fr])-1],y[fr, np.array(blank[fr])-1], 'o', color='white', markeredgecolor='white', markersize=3)

    if fr in [0,1,2,3,5,7,8,9]:
        for st in range(1,n_st[tr]+1):
            ax.text(x[tr, st-1]-0.8,y[fr, st-1]-0.5,str(st), fontsize=fs-2)

ax.legend([l1,l2,l3],[l1.get_label(),l2.get_label(),l3.get_label()] , fontsize=fs)

fig.savefig(path_fig+'stations_classpy.pdf') # no numbers for as and bs



#%% fig 3 : schéma instruments (Illustrator)
#%% fig 4 : stack phy, grad, troph


troph, group_names = f.trophic(Afront_ave=True) #[transects, stations, troph_lev]; ZS and flowcyt = vertical integral; HPLC = surface


# physics
v=50 # in meters (CTD : 1 level1m)
tm=fl.load_temp(vave=v)         # [tr, st, z]
sm=fl.load_sal(vave=v)
dm=fl.load_dens(vave=v)

rel_grads, ave_grads=f.compute_norm_grads(z=v) # [tr, st, var]

fluo=fl.load_fluo(vave=v) # [tr, st]


X=[[dm, tm, sm],[rel_grads[:,:,0], rel_grads[:,:,1], rel_grads[:,:,2]],[troph[:,:,0],troph[:,:,1]], [troph[:,:,2],troph[:,:,3] , fluo], [troph[:,:,4],troph[:,:,5], troph[:,:,6] ]]


var=[[ 'Dens (0-50m)','Temp (0-50m)', 'Sal (0-50m)'],['Relative dens gradient','Relative temp gradient','Relative sal gradient'],['h-bact', 'cyano-bact'], ['other phyto', 'diatoms', 'fluo (0-100m)'],['pico-grazers','meso-grazers','carnivores'] ]

colors=[['k','blue', 'orange'],['k','blue', 'orange'],['dodgerblue', 'turquoise'], ['palegreen', 'darkgreen','limegreen'],['hotpink', 'darkviolet', 'red'] ]

delta=[1.1, 3.6, 0.7] # range of y-axis for density and temp and sal

warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations

fil_center=[[], [], [], [], [5], [6], [], []] # station number of the center of the filament, for drawing a vertical line

#%
fs=9
lA=4*cm # length of A-front figure
hA=16*cm # height of figures
org=1



fig_list, axes_list, tr_list, along_list, f_id=f.layout_stack(nrows=len(X), lA=lA, hA=hA, org=org, fs=fs)

for cr in range(4): # 1 figure for each cruise
    ntr=len(tr_list[cr])
    for ip_col in range(ntr): # for each column = each transect
        tr=tr_list[cr][ip_col]
        along=along_list[cr][ip_col]

        for ip in range(len(X)): # for each subplot = each row
            ax=axes_list[cr][ip][ip_col]

            # if ip_col==0: # first column
                # ax.set_ylabel(titl[ip] + '\n μg C /m3')

            nl=len(X[ip]) # number of lines on the subplot (=number of variables)
            lns=[]

            axes=[ax]
            for il in range(nl-1):
                if ip in [0,2,3,4]:
                    axes.append(ax.twinx())
                else:# no twinx for relative gradients
                    axes.append(ax)



            ls = ['-','-','--'] if ip==3 else ['-' for i in range(nl)] # fluo in dashed lines
            lw = [2,1,1] if ip in [0,1] else [1 for i in range(nl)] # density in thicker line


            for ivar in range(nl):
                axi=axes[ivar]
                axi.tick_params(axis='both', labelsize=fs)

                axi.ticklabel_format(axis='y', style='sci', scilimits=(0,2))
                axi.yaxis.offsetText.set_fontsize(fs-4) # fontsize of the exponent


                if ip in [0,2,3,4]: # no twinx for the gradients
                    axi.tick_params(axis='y', labelcolor=colors[ip][ivar])

                l,=axi.plot(along,X[ip][ivar][tr,:], label=var[ip][ivar], color=colors[ip][ivar], linestyle=ls[ivar], linewidth=lw[ivar])

                axi.plot(along[fr[tr]],X[ip][ivar][tr, fr[tr]], 'x', color=colors[ip][ivar], markersize=8, markeredgewidth=lw[ivar]*2) # front

                lns.append(l)


                if ip==0: # dens and temp ans sal
                    lo=np.nanmin(X[ip][ivar][tr,:]) - 0.08*abs(np.nanmin(X[ip][ivar][tr,:]) - np.nanmax(X[ip][ivar][tr,:]))
                    axi.set_ylim(lo,lo+delta[ivar] )
                elif ip==1:
                    axi.set_ylim(0,2.2 )

            if tr in [4,5]: # vertical lines at the center of the filaments E2 and F1
                ax.axvline(along[np.array(fil_center[tr])-1], linestyle='-',  linewidth=2, color='k')
            if ip in [1]: # red horizontal line at 1 for relative gradients
                ax.axhline(1, linestyle='-',  linewidth=2, color='red')

            labs = [l.get_label() for l in lns]

    fig=fig_list[cr]
    fig.savefig(path_fig + 'stack_phy_troph{}.pdf'.format(cr), bbox_inches = 'tight')


#%% fig 5 : recap detail
fs=9 # fontsize of tick labels
hfigs=10*cm

data_detail, taxa_names, group_names = f.trophic_detail(Afront_ave=True) #list of [transects, stations, var]; μgC/m^3
data_tot, group_names = f.trophic(Afront_ave=True) #[transects, stations, troph_lev]; μgC/m^3

data_naup=f.load_naup() # [tr, st, var], ratios


ba_titl='Biomass'

warm, front, cold, tr_number, f_id = f.def_fronts()
nf=len(f_id)

fig_titles=['bacteria', 'euk-phyto', 'pico-grazers', 'meso-grazers', 'carnivores', 'nauplii']

xticks=[['hbact', 'total \n cyano', 'pro', 'syn'],
        ['total phyto','diatoms', 'euk-phyto', 'dino', 'prym', 'pelago', 'chloro', 'crypto'],
        ['total','rhizaria', 'doliolids', 'append.','salps', 'pyro.' ],
        ['total', 'Calan.', 'Oithona','Other \n cop.', 'euph.', 'pterop.','other \n crust.'],
        ['total','chaeto.', 'cnid. + \n cten.', 'ostra.', 'poly.'],
        ['egg ratio', 'nauplii ratio']]


ncols=[len(ii) for ii in xticks] # number of columns in each figure
nst=data_tot.shape[1]
ntr=data_tot.shape[0]

for ifig in range(6): # 6 figures
    nvar=ncols[ifig]
    data_ifig=np.zeros((ntr, nst, nvar))  #[tr, st, ncol];


    if ifig==0: # bact
        data_ifig[:,:,0]=data_tot[:,:,0] # hbact
        data_ifig[:,:,1]=data_tot[:,:,1] # total cyano-bact
        data_ifig[:,:,2:]=data_detail[1][:,:,:] # pro and syn


    elif ifig==1: # phytoplankton
        data_ifig[:,:,0]=data_tot[:,:,3] + data_tot[:,:,2]   # total phyto
        data_ifig[:,:,1]=data_tot[:,:,3] # diatoms
        data_ifig[:,:,2]=data_tot[:,:,2] # total other euk-phyto
        data_ifig[:,:,3:]=data_detail[2][:,:,:] # all other euk-phyto



    elif ifig in [2,3,4]: # zoo
        data_ifig[:,:,0]=data_tot[:,:,ifig+2] # total
        data_ifig[:,:,1:]=data_detail[ifig+2][:,:,:] # all other taxa

    elif ifig==5: # eggs and nauplii
        data_ifig[:,:,:]=data_naup



    ## COMPUTE PEAK INTENSITY
    X, mask_trans =f.peak_intensity(data_ifig, maxF=True)


    fig, ax=plt.subplots(1,1, figsize=(X.shape[1]*3*cm,hfigs))
    # fig.suptitle(titls[ifig])

    cf=ax.pcolormesh(X[::-1], cmap='seismic') # A-front at the top
    ax.pcolormesh(mask_trans[::-1], cmap='Greys', vmin=0,vmax=2) # adjust vmax to determine the shade of grey plotted


    plt.xticks(np.arange(X.shape[1])+0.5, xticks[ifig], rotation=20, fontsize=fs)
    plt.yticks(np.arange(X.shape[0])+0.5, f_id[::-1] , fontsize=fs)

    cb=fig.colorbar(cf)
    cf.set_clim(-300,300)
    cti=np.arange(-300,400,100)
    clab=[str(x)+' %' for x in cti]
    clab[3]='0 %'
    cb.set_ticks(cti)
    cb.set_ticklabels(clab, fontsize=fs)

    fig.savefig(path_fig+'recap_detail{}.pdf'.format(ifig))



#############################################################################################
# figure with the group totals, fluo and the total phyto and zoo

data_tot, group_names = f.trophic(Afront_ave=True) #[transects, stations, troph_lev]; μgC/m^3


tot_phy=data_tot[:,:,group_names.index('cyanobact')] + data_tot[:,:,group_names.index('other euk-phyto')] + data_tot[:,:,group_names.index('diatoms')]
data_tot=np.insert(data_tot, 0, tot_phy, axis=2)
group_names.insert(0, 'total \n phyto')

tot_zoo=data_tot[:,:,group_names.index('pico-grazers')] + data_tot[:,:,group_names.index('meso-grazers')] + data_tot[:,:,group_names.index('carnivores')]
data_tot=np.insert(data_tot, 1, tot_zoo, axis=2)
group_names.insert(1, 'total zoo')


fluo_vint=fl.load_fluo(vave=100) # [tr,st]
data_tot=np.insert(data_tot, 0, fluo_vint, axis=2)
group_names.insert(0, 'fluo')


X, mask_trans=f.peak_intensity(data_tot, maxF=True)



fig, ax=plt.subplots(1,1, figsize=(X.shape[1]*3*cm,hfigs))
cf=ax.pcolormesh(X[::-1], cmap='seismic') # peak fronts at the top
ax.pcolormesh(mask_trans[::-1], cmap='Greys', vmin=0,vmax=2) # adjust vmax to determine the shade of grey plotted

x_ticks= ['fluo', 'total \n phyto', 'total \n zoo','h-bact','cyano.','euk-phyto','diatoms','pico-graz','meso-graz','carn']

plt.xticks(np.arange(X.shape[1])+0.5, x_ticks, rotation=0, fontsize=fs )

plt.yticks(np.arange(X.shape[0])+0.5, f_id[::-1] , fontsize=fs)

cb=fig.colorbar(cf)
cf.set_clim(-300,300)
cti=np.arange(-300,400,100)
clab=[str(x)+' %' for x in cti]
clab[3]='0 %'
cb.set_ticks(cti)
cb.set_ticklabels(clab, fontsize=fs)

fig.savefig(path_fig+'recappy.pdf')


#%% fig 6 : peak examples
'''
examples of peaks to show intra-front heterogeneities

autres subplots possibles :
['total cyan', 'PRO', 'SYN']  [data_detail[1].sum(axis=2),data_detail[1][:,:,0],data_detail[1][:,:,1]]    ['turquoise','green', 'blue']
[tot_euk_phyto,data_detail[3],data_detail[2][:,:,1]]        ['total euk phyto', 'diatoms', 'Prym']         ['green','palegreen', 'blue']
[data_detail[5][:,:,0],data_detail[5][:,:,1]]           ['calanoids', 'oithona']               ['darkorange', 'magenta']




 array(['Het Bact'], dtype='<U8'),
 array(['PRO', 'SYN'], dtype='<U8'),
 array(['Dino', 'Prym', 'Pelago', 'Chloro', 'Crypto'], dtype='<U7'),
 array(['Diatoms'], dtype='<U7'),
 array(['Rhizaria', 'Doliolids', 'Appendicularia', 'Pyrosomes', 'Salps'],
       dtype='<U14'),
 array(['Calanoids', 'Oithona', 'Other cop', 'Euphausiids', 'Pteropods',
        'other crust'], dtype='<U14'),
 array(['Chaetognaths', 'Cnidaria', 'Ostracods', 'Polychaete'],

#colors
cols_phys= ['k', 'blue']                                                    # dens          , temp
cols_bact= [ 'green', 'blue', 'red']                                        # hbact         , pro        , syn
cols_phyto=[ 'green', 'blue', 'red', 'darkorange', 'magenta', 'turquoise']  # 'Diatoms'     , 'Dino'     , 'Prym'          , 'Pelago'     , 'Chloro', 'Crypto'
cols_pico= [ 'green', 'blue', 'red', 'darkorange', 'magenta']               # 'Rhizaria'    , 'Doliolids', 'Appendicularia', 'Pyrosomes'  , 'Salps'
cols_meso= [ 'green', 'blue', 'red', 'darkorange', 'magenta', 'turquoise']  # 'Calanoids'   , 'Oithona'  , 'Other cop'     , 'Euphausiids', 'Pteropods'
cols_carn= [ 'green', 'blue', 'red', 'darkorange']                          # 'Chaetognaths', 'Cnidaria' , 'Ostracods'     , 'Polychaete'

'''
fs=13 # fontsize
v=50 # temp and dens averages, in meters (CTD : 1 level1m)

tm=fl.load_temp(vave=v)        # [tr, st, z]
dm=fl.load_dens(vave=v)

data_detail,taxa_names, group_names= f.trophic_detail(Afront_ave=True) # list of [tr, st, var]

# fluorescence
fluo=fl.load_fluo(vave=v) # [tr, st]


for cruise in [[0], [2], [3]]  :
    if cruise==[0]:# A
        X=[[dm, tm]             , [data_detail[4][:,:,2],data_detail[6][:,:,0]] , [data_detail[1][:,:,0], data_detail[1][:,:,1]]]
        var=[['density', 'temp'], ['appendicularians', 'chaetognaths']          , ['PRO', 'SYN']]
        colors= [['k', 'blue']  , ['darkviolet', 'red']                         , [ 'dodgerblue', 'turquoise']  ]

    if cruise==[2]:# E
        X=[[dm, tm]             , [fluo,data_detail[5][:,:,0],data_detail[5][:,:,3]], [data_detail[1][:,:,0], data_detail[1][:,:,1]]]
        var=[['density', 'temp'], ['fluo','calanoids', 'euphausiids']                    , ['PRO', 'SYN']                  ]
        colors= [['k', 'blue']  , ['palegreen','darkorange', 'hotpink']                  , [ 'dodgerblue', 'turquoise']                    ]


    if cruise==[3]:# F
        X=[[dm, tm]             ,[fluo,data_detail[3],data_detail[5][:,:,0]] , [data_detail[1][:,:,0], data_detail[1][:,:,1]]]
        var=[['density', 'temp'], ['fluo', 'diatoms', 'calanoids']                , ['PRO', 'SYN']        ]
        colors= [['k', 'blue']  , ['palegreen','green', 'darkorange']             , [ 'dodgerblue', 'turquoise']  ]



    warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations


    lA=6*cm # length of A-front figure
    hA=12*cm # height of figures
    org=2
    fig_list, axes_list, tr_list, along_list, f_id=f.layout_stack(nrows=len(X), lA=lA, hA=hA, org=org, fs=fs)

    for cr in cruise:
        ntr=len(tr_list[cr])
        fig_lns=[]
        fig=fig_list[cr]
        for ip_col in range(ntr): # for each column = each transect
            tr=tr_list[cr][ip_col]
            along=along_list[cr][ip_col]

            for ip in range(len(X)): # for each subplot
                ax=axes_list[cr][ip][ip_col]
                # if ip_col==0: # first column
                    # ax.set_ylabel(titl[ip])


                nl=len(X[ip]) # number of lines on the subplot (=number of variables)
                lns=[]
                ax_list=[ax]
                for il in range(nl-1):
                    ax_list.append(ax.twinx())

                lw = [2,1] if ip in [0] else [1 for i in range(nl)] # density in thicker line

                for ivar in range(nl):
                    Xi=X[ip][ivar][tr,:] # [st]
                    axi=ax_list[ivar]
                    axi.tick_params(axis='y', labelcolor=colors[ip][ivar])
                    axi.tick_params(axis='both', labelsize=fs)
                    axi.ticklabel_format(axis='y', style='sci', scilimits=(0,2))
                    axi.yaxis.offsetText.set_fontsize(fs-4) # fontsize of the exponent

                    l,=axi.plot(along,Xi, label=var[ip][ivar], color=colors[ip][ivar], linewidth=lw[ivar])

                    axi.plot(along[fr[tr]],Xi[fr[tr]], 'x', color=colors[ip][ivar], markersize=8, markeredgewidth=3) # front
                    lns.append(l)


                labs = [l.get_label() for l in lns]
                if ip==1: #middle row
                    fig.legend(lns, labs, fontsize=fs, loc=(0.2,0.01))
                else: # first and third row : only add legend for last panel
                    if ip_col==2:
                        ax.legend(lns, labs, fontsize=fs, loc='right')
                fig_lns.append(lns)

        fig.savefig(path_fig+'stack_expeak{}.pdf'.format(cr), bbox_inches = 'tight')


#%% fig 7 : peak width and position
fs=6
# med=True # True for median, False for mean
med=True


pos_km, w, peak_dict, var=f.biomass_peaks()

# taxa names in var are in the order defined in trophic and trophic_detail, but unnested and shorter, and with density grad and fluo at the beginning
var=['Phys', 'Fluo',
    'Het-bact',
    'PRO', 'SYN',
    'Dino', 'Prym', 'Pelago', 'Chloro', 'Crypto',
    'Diatoms',
    'Rhiz', 'Dolio', 'Append', 'Pyro', 'Salps',
    'Calanoids', 'Oithona', 'Other cop', 'Euphau', 'Ptero','Others',
    'Chaeto', 'Cnid', 'Ostra', 'Polychaete']



nf=len(pos_km)
nv=len(var)
f_id=list(pos_km.keys())

core_km=f.core_front(var='core') # core of the density peak


# compute the offset between the biomass peaks and the density gradient peaks
offset=np.zeros((nf, nv)) # [iif, ivar]
width=np.zeros((nf, nv))

# 1 if the km increase from warm to cold
warm_flag={'A':1, 'C2':1, 'C3':-1, 'E1':1, 'E2a':-1, 'E2b':1, 'F1a':1, 'F1b':-1, 'F2':-1, 'F3':1}

for iif, flab in enumerate(f_id):
    offset[iif, :]=(core_km[flab] - pos_km[flab][:])*warm_flag[flab]
    width[iif, :]=w[flab][:] # turn winto a  np.array to use in the plot


# add median by front and by trophic level. [front, var]
width=np.where(width==0,np.nan, width)# mask zeros in widths

if med :
    width=np.append(width, np.expand_dims(np.nanmedian(width, axis=0), axis=0), axis=0) # between the fronts for each variable, ignore density
    width=np.append(width, np.expand_dims(np.nanmedian(width[:,1:], axis=1), axis=1), axis=1) # between the variables for each front

    offset=np.append(offset, np.expand_dims(np.nanmedian(offset, axis=0), axis=0), axis=0) # between the fronts for each variable, ignore density
    offset=np.append(offset, np.expand_dims(np.nanmedian(offset[:,1:], axis=1), axis=1), axis=1) # between the variables for each front

    med_name='Median'
else:

    width=np.append(width, np.expand_dims(np.nanmean(width, axis=0), axis=0), axis=0) # between the fronts for each variable, ignore density
    width=np.append(width, np.expand_dims(np.nanmean(width[:,1:], axis=1), axis=1), axis=1) # between the variables for each front

    offset=np.append(offset, np.expand_dims(np.nanmean(offset, axis=0), axis=0), axis=0) # between the fronts for each variable, ignore density
    offset=np.append(offset, np.expand_dims(np.nanmean(offset[:,1:], axis=1), axis=1), axis=1) # between the variables for each front
    med_name='Average'

# between the variables for each front
width_std=np.nanstd(width[:,1:], axis=1)


nf=nf+1
nv=nv+1
f_id = np.insert(f_id,nf-1,med_name)
var = np.insert(var,nv-1,med_name)

for iif, fn in enumerate(f_id):
    print(fn, round(width[iif,-1],1), ' +/- ', round(width_std[iif], 1), '({} %)'.format(round(width_std[iif]/width[iif,-1]*100)), 'density : {}'.format(round(width[iif,0],1)))


fig, ax=plt.subplots(1,1, figsize=(12*cm,8*cm))
x=np.append([0.5],np.arange(nv-1)+2) # leave space between phy and bio
x[-1]=x[-1]+1 # leave space with the median
xx=np.tile(x, [nf,1])
y=np.arange(nf)
y[-1]=y[-1]+1 # leave space with the median
yy=np.transpose(np.tile(y, [nv,1]))

plt.scatter(xx, yy, c=offset , s=width*2, cmap='seismic', edgecolors='k', linewidths=0.5)



plt.ylabel('Front', fontsize=fs+2)
plt.xlabel('Plankton taxa', fontsize=fs+2)
plt.xticks(x, var, rotation=60 , fontsize=fs-1)
plt.yticks(y, f_id , fontsize=fs)
plt.xlim(-0.5,x[-1]+1)
plt.ylim(-0.5,y[-1]+1)
cti=np.arange(-10,15,5)
plt.clim(-10,10)
cb=plt.colorbar(shrink=0.6)
cb.set_ticks(cti)
cb.set_ticklabels(cti,fontsize=fs)
ax.invert_yaxis()

# vlines between phy and bio
ax.axvline(1.5, linestyle='-', color='k', linewidth=1)
# vlines between trophic levels
ytroph=[0.5,1.5,3.5,8.5,9.5,14.5,20.5]
for ii in ytroph:
    ax.axvline(ii+2, linestyle=':', color='k', linewidth=0.5)

# lines to separate the median
ax.axvline(x[-1]-1, linestyle='-', color='k', linewidth=1.5)
ax.axhline(y[-1]-1, linestyle='-', color='k', linewidth=1.5)

fig.tight_layout()
fig.savefig(path_fig + 'peak_widthspy.pdf', bbox_inches = 'tight')



#%% fig 8 : schéma processus (Illustrator)

#%% fig S1-S4 : satellite maps snapshots

## LOAD DATA
data, var, transect_id, dates=fl.coord(withC1=False)
dates=[x[0] for x in dates] # take the first station of each transect


data=[[] for i in range(7)] # empty lists to contain all data
lat=[[] for i in range(7)]
lon=[[] for i in range(7)]


data[0], lat[0], lon[0]=fl_sat.load_submesocolor(dataname='sst', dates=dates) # COPERNICUS SST; daily, 4 km; L4 ("cloud-free", with interpolation)
data[1], lat[1], lon[1]=fl_sat.load_submesocolor(dataname='chl', dates=dates) # COPERNICUS CHL (mg/m3); daily, 4 km; L4 ("cloud-free", with interpolation)

data[2], lat[2], lon[2]=fl_sat.load_FLSE( dates=dates, var='lambda') # [time, lat, lon] FLSE
data[3], lat[3], lon[3]=fl_sat.load_FLSE( dates=dates, var='touched') # [time, lat, lon] Water Age

data[4], lat[4], lon[4] = fl_sat.load_submeso_clouds( dates=dates) # cloud mask
data[5], lat[5], lon[5]=fl_sat.load_submesocolor(dataname='hi101', dates=dates) # [time, lat, lon] HI; daily, 4 km; window size = 20 km

data[6], lat[6], lon[6]=fl_sat.load_ssh( dates=dates) # [time, lat, lon] SSH ; m above geoid; daily, 8 km; global reanalysis

#% PLOT
# cruise=['2008-Afront','2011-Cfront','2012-Efront','2017-transect']

zoom=1
fs=10 # size of axes labels, colorbar labels, subplot titles and fig title
ncont=25
cmap_list=['plasma','plasma', 'plasma','plasma']

clim=[None,[0.05,5],[0,0.6], [0,90]] # - , Chl, FSLE, Water Age
clim_sst=[[14,18], [14,16], [14,16], [13,17], [13,17], [11,16], [13,16], [12,15]] # SST clim
cti=[np.arange(12,20+1,2),[0.01,0.1,0.5,1,5],np.round(np.arange(0,0.61,0.1), 1),np.round(np.arange(0,91,20),1)] # colorbar ticks
var_name = ['SST (°C)', 'Chl (mg/m3)', 'FSLE', 'Water Age (days)']

stations, var, transect_id, DT=fl.coord()

xlim, ylim=f_sat.lim_zoom(zoom=zoom, cruise_only=False)

warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations

transect_id=f.transect_id(short=True, add_date=True)
for ic in range(len(dates)): # for each transect
    fig, axes=plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12*cm,10*cm))
    fig.suptitle(transect_id[ic], fontsize=fs)
    axes=axes.flatten()


    for iv in range(4): # for each variable
        ax=axes[iv]
        ax.set_title(var_name[iv], fontsize=fs)
        if ic==0:
            ax.set_title(var_name[iv])
        clim_iv=clim_sst[ic] if iv==0 else clim[iv]

        if iv==1: # chl = logscale contourf
            levs=f_sat.logclim(clim_iv, ncont) if clim_iv is not None else ncont
            cf=ax.contourf(lon[iv], lat[iv], data[iv][ic], cmap=cmap_list[iv], levels=levs, extend='both', norm=mpl_colors.LogNorm())

        elif iv==0: # sst = contourf
            levs=np.linspace(clim_iv[0], clim_iv[1], ncont) if clim_iv is not None else ncont
            cf=ax.contourf(lon[iv], lat[iv], data[iv][ic], cmap=cmap_list[iv], levels=levs, extend='both')
        else : # fsle and water age = pcolor
            cf=ax.pcolormesh(lon[iv], lat[iv], data[iv][ic], cmap=cmap_list[iv])

        f_sat.remove_white_lines(cf)

        if iv==0: # contours HI sur SST
            ax.contour(lon[5], lat[5], data[5][ic], levels=[10], colors='k', linewidths=1)

        if iv==2: # contours SSH sur fsle
            ax.contour(lon[6], lat[6], data[6][ic], levels=[0,0.05,0.1,0.15,0.2,0.25,0.3], colors='white', linestyles=['--', '--', '--', '-', '-', '-', '-'], linewidths=1)

        # colorbar
        cf.set_clim(clim_iv)
        cb=fig.colorbar(cf, ax=ax, shrink=0.8)
        cti_iv=np.arange(clim_iv[0], clim_iv[1]+1, 1).astype(int) if iv==0 else cti[iv]
        cb.set_ticks(cti_iv)
        cb.set_ticklabels(cti_iv, fontsize=fs)


        if iv in [0,1]: # for submesocolor data : add hatches to show the location of clouds
            ax.contourf(lon[4], lat[4], data[4][ic], levels=[0.5,1.5], hatches=['.'], colors='none')


        # plot limits
        if len(xlim)>2: # different lims for each cruise
            ax.set_xlim(xlim[ic])
            ax.set_ylim(ylim[ic])
        else: # large scale = same lims for all cruises
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        ax.set_aspect(aspect=1)
        ax.tick_params(axis='both', labelsize=fs)

        for st in range(stations.shape[1]): # for each station
            if st in fr[ic]:
                ax.plot(stations[ic,st,1],stations[ic,st,0],'x', color='red',markersize=4) # front = red cross
            elif st in warm[ic] or st in cold[ic]:
                ax.plot(stations[ic,st,1],stations[ic,st,0],'o', color='black',markersize=3) # back = black circle
            else: # other
                ax.plot(stations[ic,st,1],stations[ic,st,0], 'o', color='black', markersize=2) # other = small black circle

        f_sat.draw_landmask(ax=ax, color='k')

    fig.tight_layout()

    fig.savefig(path_fig+'transects_maps{}.png'.format(ic))




#%% fig S5 : vertical sections fluo

v=100 # CTD : 1level = 1m

fluo=fl.load_fluo(vave=None) # [tr, st, z]
fluo_vave=fl.load_fluo(vave=v) # [tr, st]
depths=fl.load_CTD_depth()

data_hplc, var_hplc, transect_id= fl.load_hplc(rel=False, biomass=False, Afront_ave=False) #data=[transect, station, var]
chla_tot=data_hplc[:,:,0]

X=[fluo, [fluo_vave, chla_tot]]
titl=['fluo', ['fluorescence vertical integral', 'total chl a surface (HPLC)']]
colors=[None, ['limegreen', 'black'] ]

clim=[[4], [0.15]*2, [3.5]*2, [9,7,2]] # vertical sections, for each transect. same values for C2-3 and E1-2
ctis=[[[0,1,2,3,4]], [[0,0.05,0.1,0.15]]*2, [[0,1,2,3]]*2, [[0,2,4,6,8], [0,2,4,6], [0,0.5,1,1.5,2]]]

ncont=8
cmap='plasma'

warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations

lA=5*cm # length of A-front figure
hA=9*cm # height of figures
org=1 # full transects
fs=12

fig_list, axes_list, tr_list, along_list, f_id=f.layout_stack(nrows=len(X), lA=lA, hA=hA, org=org, fs=fs)



for cr in range(4): # 1 figure for each cruise
    fig=fig_list[cr]
    ntr=len(tr_list[cr])
    for ip_col in range(ntr): # for each column = each transect
        tr=tr_list[cr][ip_col]
        along=along_list[cr][ip_col]

        for ip in range(len(X)): # for each subplot
            ax=axes_list[cr][ip][ip_col]

            if ip==0: # vertical sections

                ax.set_ylim(-100,0)
                ax.tick_params(axis='both', labelsize=fs)
                levs=np.linspace(0,clim[cr][ip_col], ncont)
                if ip_col==0:
                    # ax.set_ylabel(titl[ip])
                    ax.set_ylabel('Depth (m)', fontsize=fs)


                ind_st=np.arange(f.n_st(tr=tr, n_tr=8))
                grid_y=-depths[tr,ind_st,:]
                grid_z=X[ip] [tr, ind_st,:]
                grid_x=np.reshape(np.repeat(along[ind_st],grid_z.shape[1]),grid_z.shape) # add a depth dimension to get the same shape

                cf=ax.contourf(grid_x, grid_y,grid_z,levels=levs, cmap=cmap,extend='both')

                if ip_col>0: # only display depth label for left subplot
                    ax.tick_params(labelleft=False)

                if not (cr in [1,2] and ip_col in [0]): # no cb for C2 and E1 since they are identical to C3 and E2, respectively
                    cb_width=0.25*cm/fig.get_size_inches()[0] # fraction of figsize; numerical value (0.3) is in cm
                    cb_left=ax.get_position().x1 + 0.1*cm/fig.get_size_inches()[0] # fraction of figsize; numerical value (0.3) is in cm

                    fig.subplots_adjust(right=cb_left)


                    cbar_ax = fig.add_axes([cb_left, 1.45, cb_width, 0.3]) # [left, bottom, width, height]
                    cb=fig.colorbar(cf,ax=ax, cax=cbar_ax, shrink=0.8)

                    cti= ctis[cr][ip_col]
                    cb.set_ticks(cti)
                    cb.set_ticklabels(cti, fontsize=fs-2)



            else: # lineplots

                axes = [ax, ax.twinx()]
                nvar=len(X[ip])
                lns=[]

                for ivar in range(nvar):
                    axi=axes[ivar]

                    l,=axi.plot(along,X[ip][ivar][tr,:], label=titl[ip][ivar], color=colors[ip][ivar])

                    axi.plot(along[fr[tr]],X[ip][ivar][tr, fr[tr]], 'x', color=colors[ip][ivar], markersize=8, markeredgewidth=2)
                    axi.plot(along[warm[tr]],X[ip][ivar][tr, warm[tr]], 'o', color=colors[ip][ivar], markersize=6, markerfacecolor='none') # warm
                    axi.plot(along[cold[tr]],X[ip][ivar][tr, cold[tr]], 'o', color=colors[ip][ivar], markersize=6, markerfacecolor='none') # cold

                    axi.tick_params(axis='y', labelcolor=colors[ip][ivar])
                    axi.tick_params(axis='both', labelsize=fs)

                    axi.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) # scientific notation with exponent at the top of the axis
                    axi.yaxis.offsetText.set_fontsize(fs-5)





                    lns.append(l)


                    labs = [l.get_label() for l in lns]
                    # ax.legend(lns, labs) # added in Illustrator


    fig.savefig(path_fig+'/fluo_sections{}py.pdf'.format(cr), bbox_inches = 'tight')



#%% fig S6 : vertical sections PRO and SYN A-front
'''
save directly, no Illustrator

'''


data_hplc, z_hplc, var_hplc, hplc_unit = fl.load_HPLC_Afront(biomass=True) #[st, z, var]

data_fc, var_fc, data_name_fc= fl.load_flowcyt(vave=False, biomass=True, rel=False) #[tr, st, z, var]
z_fc=data_fc[0,:,:,0]

nst=9 # 9 stations
nz=7 # 7 niskin bottles


# PLOT

data=[[data_hplc[:nst,:nz,8],data_fc[0,:nst,:nz,2] ], [ data_hplc[:nst,:nz,5], data_fc[0,:nst,:nz,3]]]
tit=[['PRO hplc'            , 'PRO fc']             , ['SYN hplc'            , 'SYN fc']]

ind_st=np.arange(f.n_st(tr=0, n_tr=8))
along=fl.along()[0,0:nst]

ncont=20
clim=[13E3, 23E3] # PRO, SYN (=first row, second row)
cb_ticks=[3E3, 5E3] # step between cb ticks

fs=10

Z=[-z_hplc[:nst,:nz].data, -z_fc[:nst,:nz]] #[st,z] # z_hplc is a masked array

X=np.reshape(np.repeat(along,Z[0].shape[1]),Z[0].shape) # add a depth dimension to get the same shape as the others [st,z]

fig, axes=plt.subplots(2,2, sharey=True, sharex=True, figsize=(12*cm,8*cm))
for ii in range(2): #row
    for jj in range(2): #column
        ax=axes[ii,jj]
        # ax.set_title(tit[ii][jj], fontsize=fs)

        levs=np.linspace(0,clim[ii], ncont)

        C=data[ii][jj] #[st, z]

        cf=ax.contourf(X, Z[jj], C,levels=levs, cmap='plasma',extend='both')
        f_sat.remove_white_lines(cf)
        ax.set_xlim(along[0], along[-1])

        cb=fig.colorbar(cf,ax=ax)
        cb.formatter.set_powerlimits((0, 0))
        tc_bi=np.arange(0,clim[ii],cb_ticks[ii])
        cb.set_ticks(tc_bi)
        cb.ax.tick_params(labelsize=fs)
        # don't display the xticks (the station ticks will be enough)
        ax.set_xticks([])
        ax.set_xticklabels([])

        f.draw_stations(f_id='A', ax=ax, fs=fs, ticks_only=False, vpad=0.15, tcolor='white')

axes[0,0].set_ylabel('PRO \n depth (m)')
axes[1,0].set_ylabel('SYN \n depth (m)')
axes[0,0].set_title('HPLC')
axes[0,1].set_title('FC')
for ii in [0,1]:
    axes[1,ii].set_xlabel('Station Number', labelpad=0.4*cm*72)

fig.tight_layout()

fig.savefig(path_fig_bg+'figA06.pdf')


#%% fig S7 : HPLC vs FC for PRO and SYN, all transects

data,  var_fc, data_name=fl.load_flowcyt(rel=False, vave=False, biomass=True)  #[tr, st, z, var] μg/m3


fc0=data[:,:,0,:] #[tr, st, var]  # surface, take out depth so fc0 and fc1 have the same variables

fc1, var_fc, transect_id=fl.load_flowcyt(rel=False, vave=True, biomass=True)  #[tr, st, var] # VERTICAL INTEGRAL, μg/m3


hplc0, var_hplc, transect_id=fl.load_hplc(rel=False, Afront_ave=False, biomass=True)  #[tr, st, var] μg/m3



X=  [[hplc0[:,:,8]      , fc0[:,:,1]      , fc1[:,:,1]]        , [hplc0[:,:,5]      , fc0[:,:,2]      , fc1[:,:,2]]]
var=[['PRO hplc surface', 'PRO fc surface', 'PRO fc vertical'] , ['SYN hplc surface', 'SYN fc surface', 'SYN fc vertical']]

colors=[['blue', 'turquoise', 'green'],['blue', 'turquoise', 'green'] ]

fs=9

warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations


lA=5*cm # length of A-front figure
hA=9*cm # height of figures
org=1
fig_list, axes_list, tr_list, along_list, f_id=f.layout_stack(nrows=len(X), lA=lA, hA=hA, org=org, fs=fs)

for cr in range(4): # 1 figure for each cruise
    ntr=len(tr_list[cr])
    for ip_col in range(ntr): # for each column = each transect
        tr=tr_list[cr][ip_col]
        along=along_list[cr][ip_col]

        for ip in range(len(X)): # for each subplot
            axi=axes_list[cr][ip][ip_col]

            # if ip_col==0: # first column
                # ax.set_ylabel(titl[ip] + '\n μg C /m3')

            nvar=len(X[ip])
            lns=[]


            for ivar in range(nvar):
                l,=axi.plot(along,X[ip][ivar][tr,:], label=var[ip][ivar], color=colors[ip][ivar])

                axi.plot(along[fr[tr]],X[ip][ivar][tr, fr[tr]], 'x', color=colors[ip][ivar], markersize=8, markeredgewidth=2)


                lns.append(l)


                labs = [l.get_label() for l in lns]
                # axi.legend(lns, labs)

                axi.tick_params(axis='both', labelsize=fs)
                axi.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) # scientific notation with exponent at the top of the axis
                axi.yaxis.offsetText.set_fontsize(fs-3)

    fig=fig_list[cr]
    fig.savefig(path_fig + 'stack_cyano{}.pdf'.format(cr), bbox_inches = 'tight')






#%% fig S8 : HPLC and EPI vertical sections A-front
'''
save directly, no Illustrator
'''

nst=9 # 9 stations
nz=7 # 7 niskin bottles


data_hplc, z, var_hplc, hplc_unit = fl.load_HPLC_Afront(biomass=True) #[st, z, var]

z_hplc=z[:nst,:nz]
ig_hplc=[2,1,3] # diatoms, dino, prym

X_hplc=[data_hplc[:nst,:nz,ii] for ii in ig_hplc]

data_epi, var_epi, data_name = fl.load_epi(vave=False) #[tr, st, z, var]
data_epi = data_epi[0] # Afront
z_epi=data_epi[:nst,:nz,0]
ig_epi=[40,41,42]# diatoms, dino, prym

X_epi=[data_epi[:nst,:nz,ii] for ii in ig_epi]

data=[[X_hplc[ii], X_epi[ii]] for ii in range(3) ]
tit =[['Diatoms HPLC', 'Diatoms EPI'],['Dino HPLC', 'Dino EPI'] ,['Prym HPLC', 'Prym EPI'] ]

#PLOT

ind_st=np.arange(f.n_st(tr=0, n_tr=8))
along=fl.along()[0,0:nst]

ncont=20
clim=[6E4, 2E4, 2E4] # fo each taxa (=each row)
cb_ticks=[1E4, 5E3, 5E3] # step between cb ticks

fs=10

Z=[-z_hplc,-z_epi] #[st,z]

X=np.reshape(np.repeat(along,nz),(nst,nz)) # add a depth dimension to get the same shape as the others

fig, axes=plt.subplots(len(data),2, sharey=True, sharex=True, figsize=(12*cm,12*cm))
for ii in range(len(data)): #row
    for jj in range(2): #column
        ax=axes[ii,jj]
        # ax.set_title(tit[ii][jj], fontsize=fs)

        levs=np.linspace(0,clim[ii], ncont)

        C=data[ii][jj] #[st, z]

        cf=ax.contourf(X, Z[jj], C,levels=levs, cmap='plasma',extend='both')
        f_sat.remove_white_lines(cf)
        ax.set_xlim(along[0], along[-1])
        ax.set_ylim(np.nanmin(Z[jj]),np.nanmax(Z[jj]))


        cb=fig.colorbar(cf,ax=ax)
        cb.formatter.set_powerlimits((0, 0))
        tc_bi=np.arange(0,clim[ii],cb_ticks[ii])
        cb.set_ticks(tc_bi)
        cb.ax.tick_params(labelsize=fs)
        # don't display the xticks (the station ticks will be enough)
        ax.set_xticks([])
        ax.set_xticklabels([])

        f.draw_stations(f_id='A', ax=ax, fs=fs, ticks_only=False, vpad=0.15, tcolor='white')

axes[0,0].set_ylabel('Diatoms \n depth (m)')
axes[1,0].set_ylabel('Dino \n depth (m)')
axes[2,0].set_ylabel('Prym \n depth (m)')
axes[0,0].set_title('HPLC')
axes[0,1].set_title('EPI')
for ii in [0,1]:
    axes[2,ii].set_xlabel('Station Number', labelpad=0.4*cm*72)

fig.tight_layout()


fig.savefig(path_fig_bg+'figA08.pdf')


#%% fig S9 : vertical sections T, S, dens

depths=fl.load_CTD_depth()

sal=fl.load_sal(vave=None)
temp=fl.load_temp(vave=None)
dens=fl.load_dens(vave=None)


X=[temp, sal, dens]
titl=['Temperature \n (°C)', 'Salinity \n (PSU)', 'Density \n (kg/m3)']

clim=[[5,17], [33,33.9], [24.25,26.25]]
r=[0,1,0] # how many decimals to round the colorbar ticklabels
ncont=15
cmap='plasma'

lA=5*cm # length of A-front figure
hA=10*cm # height of figures
fs=11
org=1
fig_list, axes_list, tr_list, along_list, f_id=f.layout_stack(nrows=len(X), lA=lA, hA=hA, org=org, fs=fs)

for cr in range(4): # 1 figure for each cruise
    fig=fig_list[cr]
    ntr=len(tr_list[cr])
    for ip_col in range(ntr): # for each column = each transect
        tr=tr_list[cr][ip_col]
        along=along_list[cr][ip_col]

        for ip in range(len(X)): # for each row = each variable
            ax=axes_list[cr][ip][ip_col]
            ax.set_ylim(-100,0)
            ax.tick_params(axis='both', labelsize=fs)
            # if ip_col==0 and cr in [0,3]:
            #     ax.set_ylabel(titl[ip], fontsize=fs+4)


            ind_st=np.arange(f.n_st(tr=tr, n_tr=8))
            grid_y=-depths[tr,ind_st,:]
            grid_z=X[ip] [tr, ind_st,:]
            grid_x=np.reshape(np.repeat(along[ind_st],grid_z.shape[1]),grid_z.shape) # add a depth dimension to get the same shape

            levs=np.linspace(clim[ip][0],clim[ip][1], ncont)
            cf=ax.contourf(grid_x, grid_y,grid_z,levels=levs, cmap=cmap,extend='both')
            cf.set_clim(clim[ip])
            ax.contour(grid_x, grid_y,grid_z,levels=levs, colors='k', linewidths=0.8)

            if (cr==2 and ip_col==1) or (cr==3 and ip_col==2): # cb for E2 and F3
                cb=fig.colorbar(cf,ax=ax, shrink=0.7)
                cti=[round(x, r[ip]) for x in levs[::4]]
                cb.set_ticks(cti)
                cb.set_ticklabels(cti, fontsize=fs-2)

            if (cr==1 and ip_col==1) or (cr==2 and ip_col==1)or (cr==3 and ip_col in [1,2]): # no yticks for C3, E2, F2 and F3
                # ax.set_yticks([])
                ax.set_yticklabels([])

    fig.savefig(path_fig+'v_sections{}py.pdf'.format(cr), bbox_inches = 'tight')





#%% fig S10 : vertical profiles nutrients
'''
no modifications in Illustrator, save directly in BG folder
'''


data, var_list=fl.load_nut() # [transect, station, depth, var]

transect_id=f.transect_id()

zz=-data[:, :,:,0]
data=data[:, :,:,1:] # take out depth
var=var_list[1:]


nutri=np.full((data.shape[0],data.shape[1],data.shape[3]),np.nan ) # tr, st, ivar
thresh=[-0.01, -0.1, None, -0.2, None] # gradient thresh to detect the nutricline

warm, front, cold, tr_number, f_id = f.def_fronts()

ave_nutri=np.zeros((len(f_id), 3, len(var))) # nutricline depth


# project on a regular z grid to be able to average the stations
data_grid, grid_z=f.interp_z(data, zz, step_z=5) #[tr, st, z, var]

ave_data=np.zeros((len(f_id), 3, len(grid_z), len(var)))
# can't use f.warm_front_cold because there are vertical levels
for fi in range(len(f_id)):
    tr=tr_number[fi]
    ii=[np.array(warm[fi])-1, np.array(front[fi])-1, np.array(cold[fi])-1]
    for iiz in range(3):
        for ivar in range(len(var)):
            if len(ii[iiz])>1: #more than 1 station : do average
                ave_data[fi, iiz, :,ivar] = np.nanmean(data_grid[tr, ii[iiz], :,ivar], axis=0)
            else: # 1 station : take value directly
                ave_data[fi, iiz,:, ivar] = data_grid[tr, ii[iiz], :,ivar]


labs=['warm', 'front', 'cold']
colors=['red', 'k', 'blue']
fs=10

fig, axes=plt.subplots(len(var),len(f_id), figsize=(16*cm,20*cm), sharey=True)
for ivar in range(len(var)): # for each row = each nutrient
    for fi in range(len(f_id)):
        tr=tr_number[fi]

        ax=axes[ivar,fi]
        if ivar==0:
            ax.set_title(f_id[fi], fontsize=fs)
        if fi==0:
            ax.set_ylabel(var[ivar], fontsize=fs)
        ax.set_ylim(-100,0)
        ax.tick_params(axis='both', labelsize=fs-3)

        for iiz in range(3):
            ax.plot(ave_data[fi, iiz, :,ivar], grid_z, color=colors[iiz], label=labs[iiz])


        iiz=1 # at the front : compute nutricline depth to draw horizontal line
        x=ave_data[fi, iiz,:, ivar]
        grad=np.gradient(x,grid_z )
        if thresh[ivar] is not None:
            ii=0
            while grad[ii]>thresh[ivar] and ii < len(grad)-1:
                ii=ii+1
            ave_nutri[fi, iiz, ivar]=grid_z[ii]
            ax.axhline(ave_nutri[fi, 1, ivar], linestyle='-',  linewidth=0.8, color='green')

fig.savefig(path_fig_bg+'figA10.pdf')



#%% fig S11-12-13 : full biological transect data, with biomass peaks
'''
no modifications in Illustrator, save directly in BG folder
'''

# vérifier valeurs de biomasse cyano avec les fichiers originaux


fs=10 # titles and legend
fsa=7 # axes labels

for flab, fnum in zip([['A', 'C2', 'E1'], ['E2b', 'F1a', 'F1b', 'F3'],['E2a', 'F2', 'C3'] ], [11,12,13]): # fnum = figure number A11-A12-A13


    nl = len(flab) # number of fronts = number of lines


    data_detail, taxa_names, group_names = f.trophic_detail(Afront_ave=True)
    data_tot, group_names = f.trophic(Afront_ave=True) #[transects, stations, troph_lev]; ZS and flowcyt = vertical integral; HPLC = surface
    troph_id= ['h-bact','cyanobact','other euk-phyto','diatoms','pico-grazers','meso-grazers','carnivores']

    ind_detail=[[], [0,1], [0,1,2,3,4], [], [0,1,2,3,4], [0,1,2,3,4,5], [0,1,2,3]]

    warm, front, cold, tr_number, f_id = f.def_fronts() # to get tr_number and f_id
    f_ind= [list(f_id).index(x) for x in flab]
    ind_tr= [tr_number[x] for x in f_ind]
    cols=[f.front_colors()[x] for x in flab]
    core_front=[f.core_front(var='coord')[x] for x in flab] # coordinates centered on the fronts


    pos_km, width, peak_dict , var= f.biomass_peaks()
    peak_list=[peak_dict[x][1:] for x in flab] # select the fronts, and ignore density

    list_ii_taxa=[1,np.nan, 2,3,np.nan, 4,5,6,7,8,9, np.nan, 10,11,12,13,14,np.nan, 15,16,17,18,19,20, np.nan, 21,22,23,24]
    # index of the taxa in peak_list to use for each subplot

    warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations in each transect

    fig, axes=plt.subplots(len(group_names), 7, figsize=(30*cm,20*cm), sharex=True)
    lines=[[] for i in range(len(flab))]
    ii_sub=-1 # index of the subplot (start at -1 because the increment is before the computations)
    for ii in range(len(group_names)): # for each row = trophic group
        ncol=len(ind_detail[ii]) +1 # total + individual taxa

        for jj in range(ncol): # for each column
            ax1=axes[ii, jj]
            ax1.set_xlim(-15,15)
            ax1.axvline(0, linestyle='-', color='k', linewidth=0.5) # core of the front
            ii_sub=ii_sub+1


            ax_list=[ax1]
            crosses=[[] for i in range(nl)] # lines in the subplot
            for il in range(nl-1):
                ax_list.append(ax1.twinx())

            for il in range(nl):
                # decide where to draw crosses

                if np.isfinite(list_ii_taxa[ii_sub]):
                    crosses[il]=np.array(peak_list[il][list_ii_taxa[ii_sub]])-1 # index of the stations with the peak

                    if len(crosses[il])==0:
                        crosses[il]=None
                elif jj==0: # no crosses in the totals because it's confusing
                    crosses[il]=None
                else:
                    crosses[il]=fr[tr]

                # print(flab[il],ii_sub, crosses)
            if jj==0: # first column = total trophic group

                ax1.set_title('total', fontsize=fs)
                ax1.set_ylabel(group_names[ii], fontsize=fs)


                for il in range(nl):
                    ax=ax_list[il]
                    ax.tick_params(axis='y', labelcolor=cols[il])
                    ax.tick_params(axis='both', labelsize=fsa)
                    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                    ax.yaxis.offsetText.set_fontsize(fsa)

                    tr=ind_tr[il]


                    l1,=ax.plot(core_front[il],data_tot[tr, :, ii] , color=cols[il], label=flab[il])

                    if crosses[il] is not None:
                        l2, = ax.plot(core_front[il][crosses[il]], data_tot[tr, :, ii][crosses[il]], 'x', label=flab[il], color=cols[il], markersize=8, markeredgewidth=2) # symbols for front stations for the legend
                    lines[il]=l1

            else:# trophic detail
                tit=taxa_names[ii][ind_detail[ii][jj-1]]
                ax1.set_title(tit, fontsize=fs)


                for il in range(nl):
                    ax=ax_list[il]
                    ax.tick_params(axis='y', labelcolor=cols[il])
                    ax.tick_params(axis='both', labelsize=fsa)
                    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                    ax.yaxis.offsetText.set_fontsize(fsa)

                    tr=ind_tr[il]
                    data_il=data_detail[ii][tr,:, ind_detail[ii][jj-1]]
                    xil=core_front[il]

                    ax.plot(xil,data_il , label=flab[il], color=cols[il])
                    if crosses[il] is not None:
                        ax.plot(xil[crosses[il]], data_il[crosses[il]], 'x', label=flab[il], color=cols[il], markersize=8, markeredgewidth=2)


        for j in range(len(ind_detail[ii]) +1, 7):
            axes[ii, j].axis('off')
        axes[0, 3].legend(lines, flab, fontsize=fs)

    fig.tight_layout(w_pad=0.1, h_pad=0.6)
    fig.align_ylabels()

    fig.savefig(path_fig_bg + 'figA{}.pdf'.format(fnum))

#%% videos

ndays=90 # number of days before and after the cruise (for submesocolor data, timestep = 1 day)
delta=dt.timedelta(days=ndays)

## LOAD DATA
data, var, transect_id, dates=fl.coord(withC1=False)
dates=[x[0] for x in dates] # take the first station of each transect

# find the range for each cruise
cruise_name=['Afront','Cfront','Efront','Ffront'] # titles of the figures, also used in the path to save the images
stat=[[0],[1,2],[3,4], [5,6,7]]  #transects in each cruise
date_cruise=[ [dates[np.min(x)]-delta,dates[np.max(x)]+delta ] for x in stat ]


# plot parameters (must be the same as the snapshots Figs S1-4)
zoom=1
fs=10 # size of axes labels, colorbar labels, subplot titles and fig title
ncont=25
cmap_list=['plasma','plasma', 'plasma','plasma']

clim=[None,[0.05,5],[0,0.6], [0,90]] # SST, Chl, FSLE, Water Age
clim_sst=[[14,18], [13,17], [13,17], [13,17]]
cti=[np.arange(10,20+1,2),[0.01,0.1,0.5,1,5],np.round(np.arange(0,0.61,0.1), 1),np.round(np.arange(0,91,20),1)] # colorbar ticks
var_name = ['SST (°C)', 'Chl (mg/m3)', 'FSLE', 'Water Age (days)']


stations, var, transect_id, DT=fl.coord()
transect_id=f.transect_id(short=False, add_date=False)
tr_id=list(transect_id)  # list with elements deleted as they are used. This is the simplest way to deal with the multiple transects in the cruises

xlim, ylim=f_sat.lim_zoom(zoom=zoom, cruise_only=True)
warm, fr , cold , transect_id =f.def_front_in_transect() # indices of front stations


data=[[] for i in range(7)] # empty list to contain all data
lat=[[] for i in range(7)]
lon=[[] for i in range(7)]


for ic in range(4): # one video for each cruise (must do them all, otherwise transect names are displayed wrong)
    # find all the list of dates to use in the video

    DD=f_sat.date_range(start_date=date_cruise[ic][0], end_date=date_cruise[ic][1])

    print('cruise {} : start = {}, end = {}'.format( ic, DD[0],DD[-1]))


    # LOAD DATA

    data[0], lat[0], lon[0]=fl_sat.load_submesocolor(dataname='sst', dates=DD) # COPERNICUS SST; daily, 4 km; L4 ("cloud-free", with interpolation)
    data[1], lat[1], lon[1]=fl_sat.load_submesocolor(dataname='chl', dates=DD) # COPERNICUS CHL (mg/m3); daily, 4 km; L4 ("cloud-free", with interpolation)

    data[2], lat[2], lon[2]=fl_sat.load_FLSE(dates=DD, var='lambda') # [time, lat, lon] FLSE
    data[3], lat[3], lon[3]=fl_sat.load_FLSE(dates=DD, var='touched') # [time, lat, lon] Water Age

    data[4], lat[4], lon[4] = fl_sat.load_submeso_clouds(dates=DD) # cloud mask
    data[5], lat[5], lon[5]=fl_sat.load_submesocolor(dataname='hi101', dates=DD) # [time, lat, lon] HI; daily, 4 km; window size = 20 km

    data[6], lat[6], lon[6]=fl_sat.load_ssh(dates=DD) # [time, lat, lon] SSH ; m above geoid; daily, 8 km; global reanalysis





    #% PLOT

    for it,d in enumerate(DD):

        fig, axes=plt.subplots(2, 2, sharex=True, sharey=True, figsize=(22*cm,16*cm))
        if d in dates:
            fig.suptitle(d.strftime('%d %b %Y') +' : '+ tr_id[0] , fontweight='bold')
            del tr_id[0]
        else:
            fig.suptitle(cruise_name[ic] +' : '+ d.strftime('%d %b %Y'))
        axes=axes.flatten()



        for iv in range(4): # for each variable
            ax=axes[iv]
            ax.set_title(var_name[iv], fontsize=fs)
            if ic==0:
                ax.set_title(var_name[iv])
            clim_iv=clim_sst[ic] if iv==0 else clim[iv]

            if iv==1: # chl = logscale contourf
                levs=f_sat.logclim(clim_iv, ncont) if clim_iv is not None else ncont
                cf=ax.contourf(lon[iv], lat[iv], data[iv][it], cmap=cmap_list[iv], levels=levs, extend='both', norm=mpl_colors.LogNorm())

            elif iv==0: # sst = contourf
                levs=np.linspace(clim_iv[0], clim_iv[1], ncont) if clim_iv is not None else ncont
                cf=ax.contourf(lon[iv], lat[iv], data[iv][it], cmap=cmap_list[iv], levels=levs, extend='both')
            else : # fsle and water age = pcolor
                cf=ax.pcolormesh(lon[iv], lat[iv], data[iv][it], cmap=cmap_list[iv])

            f_sat.remove_white_lines(cf)

            if iv==0: # contours HI sur SST
                ax.contour(lon[5], lat[5], data[5][it], levels=[10], colors='k', linewidths=1)

            if iv==2: # contours SSH sur fsle
                ax.contour(lon[6], lat[6], data[6][it], levels=[0,0.05,0.1,0.15,0.2,0.25,0.3], colors='white', linestyles=['--', '--', '--', '-', '-', '-', '-'], linewidths=1)

            # colorbar
            cf.set_clim(clim_iv)
            cb=fig.colorbar(cf, ax=ax, shrink=0.8)
            cti_iv=np.arange(clim_iv[0], clim_iv[1]+1, 1).astype(int) if iv==0 else cti[iv]
            cb.set_ticks(cti_iv)
            cb.set_ticklabels(cti_iv, fontsize=fs)


            if iv in [0,1]: # for submesocolor data : add hatches to show the location of clouds
                ax.contourf(lon[4], lat[4], data[4][it], levels=[0.5,1.5], hatches=['.'], colors='none')


            # plot limits
            if len(xlim)>2: # different lims for each cruise
                ax.set_xlim(xlim[ic])
                ax.set_ylim(ylim[ic])
            else: # large scale = same lims for all cruises
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
            ax.set_aspect(aspect=1)
            ax.tick_params(axis='both', labelsize=fs)

            trs=stat[ic] # indices of all the transects in the cruise
            for tr in trs:
                for st in range(stations.shape[1]): # for each station
                    if st in fr[tr]:
                        ax.plot(stations[tr,st,1],stations[tr,st,0],'x', color='red',markersize=4) # front = red cross
                    elif st in warm[tr] or st in cold[tr]:
                        ax.plot(stations[tr,st,1],stations[tr,st,0],'o', color='black',markersize=3) # back = black circle
                    else: # other
                        ax.plot(stations[tr,st,1],stations[tr,st,0], 'o', color='black', markersize=2) # other = small black circle

            f_sat.draw_landmask(ax=ax, color='k')

        # fig.tight_layout()
        fig.savefig(path_vid+cruise_name[ic]+'/_{}'.format(it)+'.jpg')
        # fig.show()

