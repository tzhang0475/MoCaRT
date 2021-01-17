#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : 1dMoCaRT.py
# Author            : tzhang
# Date              : 02.05.2020
# Last Modified Date: 06.01.2021
# Last Modified By  : tzhang

# libs
import numpy as np
import random 
import sys
from matplotlib import pyplot as plt

#################################################################################
# a demonstration 1d Monte Carlo solving neutron transport for education purpose#
#################################################################################

# Basic Parameters

# Benchmark problem UD2O-1-0-SL, 1-D slab of uranium d2o mixture
# dimension
x_r = 10.371065             # up boundary of the geometryi, in cm
x_l = -10.371065            # lower boundary of the geometry, in cm

# marco cross section set: total, scattering, fission
s_t = 0.54628    # total 
s_s = 0.464338   # scattering
s_c = 0.027314  # capture
s_f = 0.054628  # fission

S = [s_t, s_s, s_c, s_f]

# average fission neutron prodcution 
nu = 1.70

# neutron population data
N = 5000
C = 50
C_pre = 20

# flux detector mesh number
N_det = 11

# number of maximum loops per cycle
n_loop = 1e3

#####################
# start calculation #
#####################

# discrete fission neutron production
I_nu = int(nu)
R_nu = nu - I_nu

# generate meshes
dx = x_r - x_l # length of the slab

# generate detector mesh
dm = dx/N_det
mesh = [x_l]
for i in range(N_det):
    dm_new = mesh[-1] + dm
    mesh.append(dm_new)

mesh_plot = []
for i in range(N_det):
    mesh_plot.append((mesh[i]+mesh[i+1])/2)
print (mesh_plot)

# plot structure
dy = 0.1*dx 

# initialize neutron position
pos_ini = x_l + dx * np.random.rand(N)

# start cycling
pos_new = pos_ini

# array for plot
pos_fiss_axial = []

# colision counter
n_col = 0
# track length counter
n_length = 0

# initialize detector
det_col = np.zeros(N_det)

fiss_sum = []
N_sum = 0
k_set = []
for cycle in range(C+C_pre):
    pos_fiss = []
    loop = 0
    N_neu = N
    # initialize neutron flying angle in cosine
    ang = 2*np.random.rand(N) - 1   # sample cosine of scattering angle

    while (N_neu > 0 and loop<=n_loop):
        # generate flying path length
        xi = np.random.rand(N_neu)
        r= -np.log(xi)
        l_fly = r/S[0]
        
        # calculate flying distance to 1-D
        d_fly = l_fly*ang
        
        # calculate new position of neutrons
        pos_new = pos_new + d_fly
        
        # determine if there are neutrons out of boundary
        pos_new = pos_new[(x_l<=pos_new)&(pos_new<=x_r)]
        
        # account number of collisions in this loop
        if cycle > C_pre:
            n_col = n_col + len(pos_new)
            # collision in each detector mesh
            for dN in range(N_det):
                det_col[dN] = det_col[dN] + len(pos_new[(mesh[dN]<=pos_new)&(pos_new<mesh[dN+1])])

        # nuclear reaction
        # reaction mesh
        t_s = S[1]/S[0]
        t_c = (S[1]+S[2])/S[0]
        eta = np.random.rand(len(pos_new))
    
        # account position for fission
        pos_fiss_curr = (pos_new[(t_c<eta)&(eta<=1)])
        # emitting fission neutron
        lmbda = np.random.rand(len(pos_fiss_curr))
        # number of position I+1 neutrons are emitted
        pos_plus = pos_fiss_curr[(lmbda>=0)&(lmbda<=R_nu)]
        # number of position I neutrons are emitted
        pos_int = pos_fiss_curr[(lmbda>R_nu)&(lmbda<=1)]
        # account fission neutron position 
        for pos in pos_plus:
            for i in range(I_nu+1):
                pos_fiss.append(pos)
        for pos in pos_int:
            for i in range(I_nu):
                pos_fiss.append(pos)
    
        # account position for scattering
        pos_new=pos_new[(0<eta)&(eta<=t_s)]
        N_neu = len(pos_new)
        ang = 2*np.random.rand(N_neu) - 1   # sample cosine of scattering angle
    
        # check if continue looping
        if N_neu == 0 or loop>n_loop:
            if loop>n_loop:
                print ('!!! WARNING: looping not converged !!!')
            break
        loop = loop + 1

    fiss_sum.append(len(pos_fiss))
    N_sum = N_sum+N
    # fission neutron generated 
    N_gen = sum(fiss_sum)
    
    # calculate analog k_eff
    k_eff = N_gen/N_sum
    sig = (1/N_sum*(N_gen/N_sum)**2 - (1/N_sum*k_eff)**2)**(1./2.)
    sig_nd = sig/N_sum**(1./2.)
    k_set.append(k_eff)

    # collision esitimator
    phi_col = 1/(N*(cycle+1))*1/(dx/100*S[0])*n_col

    # collision estimator for 1-D phi
    phi_det = 1/(N*(cycle+1))*1/(dx/100*S[0])*det_col
    
    if cycle < C_pre:
        print ('( pre cycle ',cycle,'/',C_pre,')','pre_cycle k_eff: ',"%.6f"%k_eff)
    else:
        print ('(',cycle-C_pre,'/',C,')','analog simulation k_eff: ',"%.6f"%k_eff,'\u00B1',"%.6f"%sig_nd)
    
    # sampling new fission position
    if len(pos_fiss) > N:
        pos_new = random.sample(pos_fiss,N)
    elif len(pos_fiss) == N:
        pos_new = pos_fiss
    else:
        pos_new = pos_fiss
        pos_extend = x_l + dx * np.random.rand(N-len(pos_fiss))
        for pos in pos_extend:
            pos_new.append(pos)

# plot normalized phi
norm_phi_det = phi_det/(phi_col/N_det)
plt.figure(figsize = (14,8))
plt.plot(mesh_plot,norm_phi_det, color = 'firebrick',linewidth = '3',marker = 'o', markersize = '5')
plt.xlabel('Position',fontsize = '16')
plt.ylabel('Normalized Phi', fontsize = '16')
plt.grid(linestyle='--',linewidth = '1')

pltName = 'phi'+'_'+'det.png'
plt.savefig(pltName,dpi = 100)


# plot fission distribution
pos_fiss_axial = pos_fiss
pos_fiss_horizon = np.random.uniform(low=-1,high=1,size=len(pos_fiss_axial))
fig,ax = plt.subplots()
ax.set_aspect(aspect=0.2)
plt.axis("off")
ax.scatter(pos_fiss_horizon,pos_fiss_axial,color='coral',marker='o',s=1)
fig.show()
fname = 'fission.png'
plt.savefig(fname,dpi=300)

# plot k_eff in generations
Cycle_tot = C + C_pre
Cycle = np.arange(Cycle_tot)
plt.figure(figsize = (14,8))
plt.plot(Cycle,k_set,color='darkgreen',linewidth= '3')
plt.axvline(C_pre,color='black',linewidth='4',linestyle='--')
plt.grid(linestyle='--',linewidth = '1')
fig.show()
fname = 'keffCycle.png'
plt.savefig(fname,dpi=300)
