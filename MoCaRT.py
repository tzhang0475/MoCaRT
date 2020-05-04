#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : MoCaRT.py
# Author            : tzhang
# Date              : 02.05.2020
# Last Modified Date: 05.05.2020
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
s_t = 0.5628    # total 
s_s = 0.464338   # scattering
s_c = 0.027314  # capture
s_f = 0.054628  # fission

S = [s_t, s_s, s_c, s_f]

# average fission neutron prodcution 
nu = 1.70

# neutron population data
N = 5000
C = 500
C_pre = 50

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
#print (pos_ini)    

# initialize neutron flying angle in cosine
ang = np.random.uniform(low=-1,high=1,size=N)

# start cycling
pos_fiss = []
pos_new = pos_ini

# array for plot
pos_fiss_axial = []

# colision counter
n_col = 0
# track length counter
n_length = 0

# initialize detector
det_col = np.zeros(N_det)

for cycle in range(C+C_pre):
    pos_fiss = []
    loop = 0
    N_neu = N
    ang = np.random.uniform(low=-1,high=1,size=N)
    while (N_neu > 0 and loop<=n_loop):
        # generate flying path length
        xi = np.random.rand(N_neu)
        r= -np.log(xi)
        l_fly = r/S[0]
        
        # calculate flying distance to 1-D
        d_fly = l_fly*ang
     #   print ('flying distance', d_fly)
        
        # calculate new position of neutrons
        pos_new = pos_new + d_fly
     #   print ('new postion of neutron', pos_new)
        
        # determine if there are neutrons out of boundary
        pos_new = pos_new[(x_l<=pos_new)&(pos_new<=x_r)]
     #   print (pos_new,len(pos_new))
        
        # account number of collisions in this loop
        if cycle > C_pre:
            n_col = n_col + len(pos_new)
            # collision in each detector mesh
            for dN in range(N_det):
                det_col[dN] = det_col[dN] + len(pos_new[(mesh[dN]<=pos_new)&(pos_new<mesh[dN+1])])

#            print (pos_new)
#            print (mesh)
#            print (det_col)
#
        # nuclear reaction
        # reaction mesh
        t_s = S[1]/S[0]
        t_c = (S[1]+S[2])/S[0]
     #   print (t_s,t_c)
        eta = np.random.rand(len(pos_new))
     #   print (eta)
    
        # account position for fission
        pos_fiss_curr = (pos_new[(t_c<eta)&(eta<=1)])
       # for pos in pos_fiss_curr:
       #     pos_fiss.append(pos)
       #     pos_fiss_axial.append(pos)
        # emitting fission neutron
        lmbda = np.random.rand(len(pos_fiss_curr))
#        print (lmbda)
        # number of position I+1 neutrons are emitted
        pos_plus = pos_fiss_curr[(lmbda>=0)&(lmbda<=R_nu)]
        # number of position I neutrons are emitted
        pos_int = pos_fiss_curr[(lmbda>R_nu)&(lmbda<=1)]
#        print (pos_plus,pos_int)
        # account fission neutron position 
        for pos in pos_plus:
            for i in range(I_nu+1):
                pos_fiss.append(pos)
        for pos in pos_int:
            for i in range(I_nu):
                pos_fiss.append(pos)
#        print (pos_fiss) 


     #   print ('fission location',pos_fiss)
    
        # account position for scattering
        pos_new=pos_new[(0<eta)&(eta<=t_s)]
        N_neu = len(pos_new)
        ang = np.random.uniform(low=-1,high=1,size=N_neu)   # sample cosine of scattering angle
     #   print ('neutron remaining',len(pos_new))
    
        # check if continue looping
        if N_neu == 0 or loop>n_loop:
            if loop>n_loop:
                print ('!!! WARNING: looping not converged !!!')
            break
        loop = loop + 1

            
    
#    print ('final fission location',pos_fiss)
#    print ('number of loops', loop)
    
    # fission neutron generated 
    N_gen = len(pos_fiss)
    
    # calculate analog k_eff
    k_eff = N_gen/N

    # collision esitimator
    phi_col = 1/(N*(cycle+1))*1/(dx/100*S[0])*n_col

    # collision estimator for 1-D phi
    phi_det = 1/(N*(cycle+1))*1/(dx/100*S[0])*det_col
    
    if cycle < C_pre:
        print ('( pre cycle ',cycle,'/',C_pre,')','pre_cycle k_eff: ',"%.6f"%k_eff)
    else:
        print ('(',cycle-C_pre,'/',C,')','analog simulation k_eff: ',"%.6f"%k_eff)
    
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
plt.savefig(fname)

