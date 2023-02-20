'''
functions for plotting
'''



import matplotlib.pyplot as plt
import numpy as np


import func_cce as f
import func_load_cce as fl



def front_colors():
    '''
    returns :
        dictionary with the color of each front (to use in line plots
    '''
    colors =['red', 'orange', 'chocolate', 'gold', 'palegreen', 'green', 'lime', 'turquoise', 'cyan', 'deepskyblue']
    name = ['A'   , 'E1'    , 'E2a'      , 'E2b' , 'C2'       ,'C3'    , 'F1a' ,'F1b'       , 'F2'  , 'F3' ]
    return dict(zip(name, colors))




def draw_stations(pos, ind_st=None, ax=None, tr=None, fil_flag=False, ticks_only=False):
    '''
    draws vertical lines on an existing (sub)plot at the positions specified, and adds labels with station number
    arguments :
        pos = list of positions along the transect (in km)
        ind_st = station number, must be passed for omics, but for the others they are sequential and can be guessed
        ax : draw on ax if specified, or on the current ax
        tr = transect index. if given, writes the fronts in bold
            If not given, the plot is probably a line plot and has the fr/fil as markers
        fil_flag : if True, draws rectangles around filament station numbers
        ticks_only : if True, only draws the ticks at each station but not the numbers
    '''
    if ax is None:
        ax=plt.gca()
    # print('z : ', z)
    if ind_st is None:
        num_s=np.arange(np.isfinite(pos).sum())+1   #list of station numbers
    else:  #for omics, when not all stations are sampled
        num_s=ind_st
        # print('draw at num_s : ', num_s)

    warm, front, cold, tr_number, f_id = f.def_fronts()
    fr_ind=[np.where(tr_number==x)[0] for x in range(8)] # index of the front in each transect

    ytext=np.min(ax.set_ylim()) - 0.10*abs(np.diff(ax.set_ylim()))# data coord ; under the plot

    for s in range(len(num_s)):
        if np.min(ax.set_xlim())<=pos[s]<=np.max(ax.set_xlim()): # only draw stations between the xlims
            ax.axvline(pos[s], color='k', linestyle='--', linewidth=0.5, alpha=0.7, ymax=0.07) # dashed line indicating the position of the station

            if not ticks_only: # write the station numbers
                if tr is not None:
                    if s+1 in f.flat(front[fr_ind[tr]]): # front = bold
                        ax.text(pos[s], ytext, str(num_s[s]), fontsize=8, alpha=1, fontweight='bold') # fronts in bold
                    elif s+1 in f.flat(warm[fr_ind[tr]]): # warm = red
                        ax.text(pos[s], ytext, str(num_s[s]), fontsize=7, alpha=1, color='red') # fronts in bold
                    elif s+1 in f.flat(cold[fr_ind[tr]]): # cold = blue
                        ax.text(pos[s], ytext, str(num_s[s]), fontsize=7, alpha=1, color='blue') # fronts in bold

                    else:
                        ax.text(pos[s], ytext, str(num_s[s]), fontsize=7, alpha=1) # other stations





def xlims(tr, n_tr):
    '''
    defines the xlimits of transects subplots
    comment here to return values or None (autoscale)
    arguments :
        tr=index of transect
        n_tr = number of transects
    '''
    # print('in xlims, n_tr=', n_tr)
    if n_tr==3:  #P1706 omics
        left=[0,42,0]
        right=[52,94,52]
    elif n_tr==8:
#        left=[0,0,0,0,0,0,42,0]  #truncate P1706 transect 2
        left=[0,0,0,0,0,0,0,0]
#        right=[26,26,26,52,52,52,94,52] #3 different ranges
        right=[52,52,52,52,52,52,94,52] #all at 50km (except transect2)
    return [left[tr],right[tr] ]
#    return [None,None]  #--> autoscale





def lim_zoom(zoom, cruise_only=False):
    '''
    defines xlim, ylim for various degrees of zooming in on the fronts, on satellite maps
    arguments :
        zoom=0,1,2,3
                0 : entire CCS (HI statistics 'cce_total', 'hi101')
                1 : CCE cruises
                2 : transect scale
                3 : pixel scale
        cruise_only :   if False, zoom=2 and 3 returns lists of 8 lims
                            ['AFront', 'Upfront 2', 'Upfront 3', 'EFront1', 'EFront2',  'Transect1', 'Transect2', 'Transect3']
                        if True, returns only the 4 values for each cruise, without the duplicates for the transects

    '''

    if zoom==0: # very large scale = entire CCS
        xlim=[-133,-110]
        ylim=[20,50]
    elif zoom==1: # regional scale = rectangle around cruises
        xlim=[-124,-119]
        # ylim=[31,36]
        ylim=[32,36]
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

    elif zoom==3: # closer zoom to see individual pixels and stations
        if not cruise_only:
            xlim=[[-121,-120.5],[-121.9,-121.5],[-121.3,-120.7],[-123.3,-122.4],[-123.3,-122.4],[-121.8,-121],[-123.2,-121.6],[-123,-122.2]]
            ylim=[[32.5,33]    ,[33.65,34.05]    ,[33.3,33.8]    ,[34.2,35]      ,[34.2,35]      ,[34.7,35.4],[34,35.4],[33.6,34.4]]
        else:
            xlim=[[-121,-120.5],[-122,-120.6],[-123.4,-122.4],[-123.2,-120.8]]
            ylim=[[32.5,33],[33.3,34.2],[34.2,35],[33.5,35.5]]

    return xlim, ylim



def logclim(clim, ncont=25):
    '''
    creates levels with log spacing

    arguments :
        clim=[lo, hi]
        ncont = integer

    '''
    if clim[0]==0:
        print('clim[0] must be >0')
        return
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



########### prepare the layout of plots


def layout_stack(nrows, lA, hA, org):
    '''
    prepares the layout of the figures for stack plots
    4 figures (1 for each transects); in each figure 1 column for each transect
    the figsize and the position of the axes are set to account for the different lengths of the transects (so they are displayed on the same scale)


    arguments :
        nrows = number of rows in the subplots
        lA= length of A-front figure
        hA= height of figures
        org :  1 = full transects; 2 = centered on each front
    returns :
        fig_list = list of 4 fig handles
        axes_list = list of 4 arrays [nrows, ntr]
        tr_list = index of the transect in each column subplot [cr][ncol]
        along_list = xaxis for each subplot column [cr][ncol]
    '''
    # cr_id=['A-front 2008', 'C-front 2011', 'E-front 2012', '2017 cruise']

    # padhs=0.08 # horizontal padding around subplots (portion of fig size)
    padhs=[0.08,0.22,0.14,0.1] # different padding for each figure

    padv1=0.05 # vertical padding at top and bottom of figure (portion of fig size)
    padv2=0.035 # vertical padding between subplots (portion of fig size)


    fig_list=[[] for i in range(4)]
    axes_list=[[] for i in range(4)]
    along=fl.along()

    core=f.core_front()

    if org==1:
        tr_list=[[0],[1,2],[3,4],[5,6,7]] # indices of the transects in each figure
        transect_id=f.transect_id(short=False)
        f_id=[['A'],['C2', 'C3'], ['E1', 'E2'], ['F1', 'F2', 'F3']]

        a=np.nanmax(along[0])-np.nanmin(along[0])
        c2,c3=np.nanmax(along[1])-np.nanmin(along[1]), np.nanmax(along[2])-np.nanmin(along[2])
        e1,e2=np.nanmax(along[3])-np.nanmin(along[3]), np.nanmax(along[4])-np.nanmin(along[4])
        t1,t2,t3=np.nanmax(along[5])-np.nanmin(along[5]), np.nanmax(along[6])-np.nanmin(along[6]),  np.nanmax(along[7])-np.nanmin(along[7])

        along_list= [[[]],[[],[]],[[],[]],[[],[],[]]] #[cr][col]



        wfact=[[1],[1,c3/c2], [1,e2/e1], [1,t2/t1,t3/t1]] # adjust the width of subplots within each figure

        h=(1-(nrows-1)*(padv1+padv2))/nrows # height of the subplots
        bottoms=[padv1 +(nrows-irow)*(h+padv2) for irow in range(nrows)] # bottoms of the subplots



        for ifig in range(4): # for each figure

            padh=padhs[ifig] if len(padhs)==4 else padhs # same padding for all figures or different values
            padkm=a*padh/(1-2*padh)
            widths_left=[a/(a+2*padkm), c2/(c2+c3+3*padkm), e1/(e1+ e2+3*padkm), t1/(t1+ t2+ t3+4*padkm)]

            fact=[(a+2*padkm)/a,(c2+c3+3*padkm)/a, (e1+e2+3*padkm)/a, (t1+t2+t3+4*padkm)/a] # adjust figsize to account for the differing lengths of transects


            fig=plt.figure(figsize=(lA*fact[ifig], hA))
            fig_list[ifig]=fig
            # fig.suptitle(cr_id[ifig])

            ntr=len(tr_list[ifig])
            w=np.zeros((ntr))
            axes=[[[] for i in range(ntr)] for j in range(nrows)]

            for irow in range(nrows):
                for icol in range(ntr):

                    w[icol]=wfact[ifig][icol]*widths_left[ifig] # width of the subplot

                    left=(icol+1)*padh + w[:icol].sum() # left of the subplot

                    tr=tr_list[ifig][icol]
                    along_list[ifig][icol]=along[tr] # store the xaxis (for org==1 it is trivial)


                    ax=plt.axes([left, bottoms[irow], w[icol], h])
                    # if irow==0:
                    #     ax.set_title(transect_id[tr])

                    ax.set_xlim(np.nanmin(along[tr]), np.nanmax(along[tr]))
                    axes[irow][icol] = ax

                    if irow == nrows-1: # bottom row
                        ax.set_xlabel('Distance along transect (km)')
                        ax.tick_params(axis='x', pad=12)
                    else: # don't display the xticks
                        ax.tick_params(labelbottom=False)

                    #flip xaxis
                    if tr in [2,4,6]:  #Up3,E2, T2
                        # print('invert xaxis : ', cr,tr, ip, ip_col)
                        ax.invert_xaxis()

            axes_list[ifig]=axes
            fig.subplots_adjust(top=bottoms[0]+h, bottom=padv1)

    elif org==2:
        tr_list=[[0],[1,2],[3,4,4],[5,5,6,7]] # indices of the transects in each figure (excluding Up1)
        transect_id=f.transect_id(short=False)
        f_id=[['A'],['C2','C3'], ['E1', 'E2a', 'E2b'],['F1a', 'F1b', 'F2', 'F3']]

        xl=[[9.3], [10, 10], [15, 13, 13], [10, 10,15,20]] # xlims, plots are symmetrical around 0

        along_list= [[[]],[[],[]],[[],[],[]],[[],[],[],[]]] #[cr][col]
        a, c2, c3a=xl[0][0]*2, xl[1][0]*2, xl[1][1]*2
        e1, e2a, e2b=xl[2][0]*2, xl[2][1]*2, xl[2][2]*2
        t1a, t1b, t2a, t3 = xl[3][0]*2, xl[3][1]*2, xl[3][2]*2 , xl[3][3]*2


        wfact=[[1],[1,c3a/c2], [1,e2a/e1, e2b/e1], [1,t1b/t1a,t2a/t1a,t3/t1a]] # adjust the width of subplots within each figure

        h=(1-(nrows-1)*(padv1+padv2))/nrows # height of the subplots
        bottoms=[padv1 +(nrows-irow)*(h+padv2) for irow in range(nrows)] # bottoms of the subplots

        for ifig in range(4): # for each figure
            padh=padhs[ifig] if len(padhs)==4 else padhs # same padding for all figures or different values
            padkm=a*padh/(1-2*padh)
            fact=[1,(c2+c3a+3*padkm)/(a+2*padkm), (e1+ e2a+ e2b+4*padkm)/(a+2*padkm), (t1a+ t1b+ t2a+ t3+5*padkm)/(a+2*padkm)] # adjust figsize to account for the differing lengths of transects

            widths_left=[a/(a+2*padkm), c2/(c2+c3a+3*padkm), e1/(e1+ e2a+ e2b+4*padkm), t1a/(t1a+ t1b+ t2a+ t3+5*padkm)]

            fig=plt.figure(figsize=(lA*fact[ifig], hA))
            fig_list[ifig]=fig
            # fig.suptitle(cr_id[ifig])

            ntr=len(tr_list[ifig])
            w=np.zeros((ntr))
            axes=[[[] for i in range(ntr)] for j in range(nrows)]

            for irow in range(nrows):
                for icol in range(ntr):

                    w[icol]=wfact[ifig][icol]*widths_left[ifig] # width of the subplot

                    left=(icol+1)*padh + w[:icol].sum() # left of the subplot

                    tr=tr_list[ifig][icol]
                    along0=along[tr] - core[ifig][icol]

                    # print('ifig, icol', ifig, icol)
                    along_list[ifig][icol]= along0 # store the xaxis (for org==1 it is trivial)


                    ax=plt.axes([left, bottoms[irow], w[icol], h])
                    if irow==0:
                        ax.set_title(transect_id[tr]+' - '+ f_id[ifig][icol])

                    ax.set_xlim(-xl[ifig][icol], xl[ifig][icol])
                    axes[irow][icol] = ax

                    if irow == nrows-1: # bottom row
                        ax.set_xlabel('Distance to core of the front (km)')
                        ax.tick_params(axis='x', which='major', pad=12)
                    else: # don't display the xticks
                        ax.tick_params(labelbottom=False)

                    ax.axvline(0, linestyle='-', color='k', linewidth=0.5) # vertical line at core of front

                    #flip xaxis
                    if f_id[ifig][icol] in ['E2a','F1b', 'F2', 'C3']:
                        ax.invert_xaxis()

            axes_list[ifig]=axes

    return fig_list, axes_list, tr_list, along_list, f_id


