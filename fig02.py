#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

import func_cce as f
import func_load_cce as fl
import func_plot_cce as fp



warm, front, cold, tr_number, f_id = f.def_fronts() # stations numbers
nf=len(f_id)
n_st=f.n_st(tr=None,n_tr=8)

x=fl.along()  #distance along transect (km) [tr, st].


yy=np.arange(nf,0,-1) # even vertical spacing
yy=[10,  9, 8, 7, 5.8, 5.2, 3.8,  3.2,  2,  1] # A and B closer
y=np.repeat(yy,x.shape[1]).reshape([nf,x.shape[1]]) # 10 fronts, 13 stations


fig, ax=plt.subplots(1,1, figsize=(10,7))

# ax.axis('off')
ax.set_yticks(y[:,0])
ax.set_yticklabels(f_id)
ax.set_ylim(0.2,10.8)
# ax.set_xlabel('distance along transect (km)')

for fr in range(len(f_id)):
    tr=tr_number[fr]

    ax.plot(x[tr, np.arange(n_st[tr])],y[fr, np.arange(n_st[tr])], '.', color='k', markersize=2)

    ax.plot(x[tr, np.array(front[fr])-1],y[fr, np.array(front[fr])-1], 'X', color='k')
    ax.plot(x[tr, np.array(warm[fr])-1],y[fr, np.array(warm[fr])-1], 'o', color='red')
    ax.plot(x[tr, np.array(cold[fr])-1],y[fr, np.array(cold[fr])-1], 'o', color='blue')
    for st in range(1,n_st[tr]+1):
            ax.text(x[tr, st-1]-0.6,y[fr, st-1]-0.4,str(st))

