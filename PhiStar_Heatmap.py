# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:26:48 2024

@author: tanne
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Default Parameters
x0 = 0.1
r = 2.5
kappa = 68040
beta = 0.81
T = 10
c = 1/1744
seed = 140000

def q_effective(a, phi):
    return np.log(phi + 1) / np.log(a)

def get_mstar(beta,s,kappa,phi,r):
    return (beta*s/kappa)*(1-(1/((1-phi)*r)))

def get_profit(phi, r, kappa, beta, seed, c, theta, m, d, b):
    m_0 = get_mstar(beta,seed,kappa,0,r)
    profit = b*(c*beta*seed)*(1-theta*(m/(m_0))) - d*phi*beta*seed
    return profit

def get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, d, c, theta):
    x = [x0]
    x_old = x0

    profit = 0

    for i in range(T):
        x_new = (1-phi)*r*x_old*(1-kappa*(x_old/(beta*seed)))
        x.append(x_new)
        profit += get_profit(phi, r, kappa, beta, seed, c, theta, x_new, d, b)
        x_old = x_new

    return profit

accuracy = 41

# sensetivity parameters
a_values = [2,2.718,10]
a_values_str = ['2','e','10']
theta_values = [0.3,0.7]
cost_values = np.linspace(0, 0.002, accuracy)
bushel_values = np.linspace(0, 20, accuracy)
phi_values = np.linspace(0.0, 1, 50)
cmap = cm.hot

f, ax = plt.subplots(3, 2, figsize=(5, 7.5), sharex='col', sharey='row')

all_res = []
#iterate through all parameter conditions to run simulations
for iiii,a in enumerate(a_values):
    all_res.append([])
    for jjjj,theta in enumerate(theta_values):
        res_all = []
        for i,b in enumerate(bushel_values):
            res_all.append( [] )
            for ii,d in enumerate(cost_values):
                res_all[-1].append( [] )
                for j,phi_1 in enumerate(phi_values):
                    phi = q_effective(a, phi_1)
                    res_all[-1][-1].append(get_timeseries_profit(x0, phi, r, kappa, beta, seed, T, d, c, theta))
        res_all = np.array(res_all)
        #find highest profit from all simulations
        res = res_all.argmax(2)/(len(phi_values)-1)
        all_res[-1].append(res)
        images=[]
        images.append(ax[iiii,jjjj].imshow(res,cmap=cmap,origin='lower',vmin=0,vmax=1))

        #Tilde_C_values_ticks = [0,0.02,0.04,0.06,0.08]
        yticks = [0,(accuracy-1)/2,accuracy-1]
        ylabels = ['0',str(np.round(bushel_values[int((accuracy-1)/2)])),bushel_values[-1]]
        ax[iiii,jjjj].set_yticks(yticks)
        if jjjj==0:
            ax[iiii,jjjj].set_yticklabels(ylabels)
            ax[iiii,jjjj].set_ylabel('b - Bushel Price ($)')
            
       
        xticks = [0,(accuracy-1)/2,accuracy-1]
        xlabels = ['0',str((cost_values[int((accuracy-1)/2)])),cost_values[-1]]
        ax[iiii,jjjj].set_xticks(xticks)
        if iiii==2:
            ax[iiii,jjjj].set_xticklabels(xlabels)
            ax[iiii,jjjj].set_xlabel('c - Fungicide Cost ($/plant)')
            
all_res = np.array(all_res)

# Get the position of the first and last subplots
pos1 = ax[0, 0].get_position()
pos2 = ax[2, 1].get_position()

# Set the colorbar position based on the combined height of the subplots
cbar_ax = f.add_axes([pos2.x1 + 0.02, pos2.y0, 0.02, pos1.y1 - pos2.y0])

# Add the colorbar
cbar = f.colorbar(images[0], cax=cbar_ax)
cbar.set_label(r'$\phi^*$ - Optimal Fungicide Application Rate')

# adding exta labels to the plots 
ax_left = f.add_axes([pos1.x0 - 0.18, pos2.y0, 0.02, pos1.y1 - pos2.y0])
ax_left.set_xticks([])
ax_left.set_yticks([])
for spine in ax_left.spines.values():
    spine.set_visible(False)
ax_left.plot([0,0],[0,1],'k-',lw=0.5)
ax_left.set_ylim([0,1])
ax_left.text(-0.4,0.5,'Maximal Fungicide Efficiency',ha='center',va='center',rotation=90)    
ax_left.text(-0.1,0.16666,'Low ('+r'$a=$ '+a_values_str[2]+')',ha='center',va='center',rotation=90)    
ax_left.text(-0.1,0.5,'Medium ('+r'$a =$ '+a_values_str[1]+')',ha='center',va='center',rotation=90)    
ax_left.text(-0.1,0.8333,'High ('+r'$a =$ '+a_values_str[0]+')',ha='center',va='center',rotation=90)    

ax_top = f.add_axes([pos1.x0, pos1.y1 + 0.01, pos2.x1 - pos1.x0, 0.02])

ax_top.set_xticks([])
ax_top.set_yticks([])
for spine in ax_top.spines.values():
    spine.set_visible(False)
ax_top.plot([0,1],[0,0],'k-',lw=0.5)
ax_top.set_xlim([0,1])
ax_top.text(0.5,0.52,'Low ($r = 2$)',ha='center',va='center') 
ax_top.text(0.5,0.72,'Mold Growth Rate',ha='center',va='center') 
ax_top.text(0.5,0.28,'Maximal Mold Damage to Crops',ha='center',va='center')    
ax_top.text(0.25,0.1,'Low ('+r'$\theta =$ '+str(theta_values[0])+')',ha='center',va='center',rotation=0)    
ax_top.text(0.75,0.1,'High ('+r'$\theta =$ '+str(theta_values[1])+')',ha='center',va='center',rotation=0)    


ax_top_2 = f.add_axes([pos1.x0, pos1.y1 + 0.085, pos2.x1 - pos1.x0, 0.02])

ax_top_2.set_xticks([])
ax_top_2.set_yticks([])
for spine in ax_top_2.spines.values():
    spine.set_visible(False)
ax_top_2.plot([0,1],[0,0],'k-',lw=0.5)

# Adjust the spacing between subplots
f.subplots_adjust(wspace=0.09, hspace=0.09)

#saving plots
f.savefig('fig3_02.png',dpi=400,bbox_inches = "tight")   

