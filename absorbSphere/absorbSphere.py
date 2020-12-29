#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : absorbSphere.py
# Author            : tzhang
# Date              : 28.12.2020
# Last Modified Date: 29.12.2020
# Last Modified By  : tzhang

'''

An absorb sphere as a Monte Carlo example

Miller and E.E. Lewis Computational Methods pf Neutron Transport p.321

'''

import numpy as np
import random
import sys
from matplotlib import pyplot as plt


# basic parameters 

R_sph = 0.5   # radius of the sphere 
S_t = 1 # cross section

# Number of particles 
N = 10000
PI = 3.1415926

# volume of the sphere
V = 4./3. * PI * R_sph**3

# sampling source position 
theta = np.random.rand(N)
r = R_sph*theta**(1./3.)

# sampling transport distrance 
dr = -np.log(np.random.rand(N))

# angle sampling
eta = np.random.rand(N)
mu = 2 * eta - 1

# distance to surface 
dy = (R_sph**2 - r**2*(1-mu**2))**(1./2.) - r*mu

# net distance 
net_d = dy - dr

# collision esitimator
# number of collisions in shpere
n_c = len(list(filter(lambda x: (x >= 0), net_d)))
c_cap = n_c/N
# error estimation
c_sd = (n_c/N - c_cap**2)**(1./2.)/N**(1./2.)
print ('collision estimator',"%.5f"%c_cap,'\u00B1',"%.6f"%c_sd)


# track length estimator
#escape = list(filter(lambda x: (x < 0), net_d))
#l_total = sum(dr)+sum(escape)
l = np.asarray([min(dr[i],dy[i]) for i in range(len(dr))],dtype=float)
l_cap = 1/N * sum(l)
#error estimation
l_sd = (1/N*sum(l**2) - l_cap**2)**(1./2.)/N**(1./2.)
print ('track length estimator',"%.5f"%l_cap,'\u00B1',"%.6f"%l_sd)
