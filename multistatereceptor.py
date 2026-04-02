#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 10:08:49 2025

@author: moor
"""

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import matplotlib

from General_Algebra import Reaction, System, calculateCovariance, MatrixEquation

# for plotting 
lwd=3
matplotlib.rc('font',size=30)
plt.rcParams['ps.useafm']=True
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['FreeSans']})
plt.rcParams['pdf.fonttype'] = 42

t_max = 1000 #maximal integration time
dt= 0.01 #step size of the integraton

#%% Systems' definition - Multi-State

def fill_x0(n, base_z0, Mstar): #species vector depending on number of receptor states
    z0 = [Mstar/(n+1)] + base_z0[1:4] + [Mstar/(n+1)]*n
    return z0

def fill_stoich(n, numrec, vec0, vec_add): #stoichiometry depending on number of receptor states
    '''
    n : number of states
    numcell : index of receptor state 
    vec0 : vector of entries for species that are not constantly added
    vec_add : the vector that we would to continuously add
    -------
    vec : final stoichiometry vector 

    '''
    vec = vec0
    for i in range(n):
        if i==numrec:
            vec = vec + vec_add
        if i!=numrec:   
            vec = vec + [0]
    return vec

def ligrec_downstream(const, x, obsindx): #function that defines the reaction network
    roff0=const[-1]
    r1 = Reaction("r1", fill_stoich(n, 0, [0, 0, 0], [0]), fill_stoich(n, 0, [0, 1, 0], [0]), const[0], const[0])
    r2 = Reaction("r2", fill_stoich(n, 0, [0, 1, 0], [0]), fill_stoich(n, 0, [0, 0, 0], [0]), const[1], const[1]*x[1])
    r3 = Reaction("r3", fill_stoich(n, 0, [0, 1, 0], [0]), fill_stoich(n, 0, [1, 1, 0], [0]), const[2]*roff0, const[2]*roff0*x[1])
    r4 = Reaction("r4", fill_stoich(n, 0, [1, 0, 0], [0]), fill_stoich(n, 0, [0, 0, 0], [0]), const[3], const[3]*x[0])
    r5 = Reaction("r5", fill_stoich(n, 0, [1, 0, 0], [0]), fill_stoich(n, 0, [1, 0, 1], [0]), const[4], const[4]*x[0])
    r6 = Reaction("r6", fill_stoich(n, 0, [0, 0, 1], [0]), fill_stoich(n, 0, [0, 0, 0], [0]), const[5], const[5]*x[2])
    r7 = Reaction("r7", fill_stoich(n, 0, [1, 0, 0], [0]), fill_stoich(n, 0, [0, 0, 0], [1]), const[6], const[6]*x[0]) 
    r = []
    for i in range(n-1):
        r.append(Reaction("r7_%d"%i, fill_stoich(n, i, [0, 0, 0], [1]), fill_stoich(n, i+1, [0, 0, 0], [1]), const[6], const[6]*x[3+i]))
        
    r.append(Reaction("r8", fill_stoich(n, n-1, [0, 0, 0], [1]), fill_stoich(n, 0, [0, 0, 0], [0]), const[6], const[6]*x[n-1]))
    for i in range(n):
        r.append(Reaction("r9_%d"%i, fill_stoich(n, i, [0, 0, 0], [1]), fill_stoich(n, i, [0, 0, 1], [1]), const[4], const[4]*x[3+i]))
        
    S = System(len(x), obsindx)
    S.addReactions([r1, r2, r3, r4, r5, r6, r7]+r)
    return S

#%% Mutual Information - Multi-State

def MutualInformationRate_downstream(const_vec, variances_x, variances_lx, x_eq, n): 
    mirate = const_vec[4]/2*(variances_x[-1,0,0] - variances_lx[-1,0,0])/(x_eq[0])
    for i in range(n-1):
        mirate = mirate + const_vec[4]/2*(variances_x[-1,3+i,3+i] - variances_lx[-1,2+i,2+i])/(x_eq[4+i])
    return mirate

#%% Pipeline rate of multi-state receptor depending on n

kswitch = 1
kmin = 1

const={'k1': 1, 'k2': 0.01, 'koff': 1, 'kon': kmin, 'k4': 500, 'k5': 50, 'kswitch': kswitch, 'roff': 1} #reaction constants
const=list(const.values())

M = 10 #number of receptors
roff0 = (const[1]*const[3]*M)/(const[0]*const[2]+const[1]*const[3]) #initial condition for unbound state
const[-1] = roff0
ron0 = M-roff0  #initial condition for bound state
lig0 = const[0]/const[1]  #initial condition for ligand
x0 = const[4]/const[5]*ron0  #initial condition for downstream signal 
n = 2 #initialising number of receptor states
const[3] = kmin*(n+1) #initialising constant 

iniconds = [ron0, roff0, lig0, x0] #initial conditions of the species

nlist = [1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40] #list of receptor states
rate_ms = np.zeros(len(nlist)) #information rate 

for i in range(len(nlist)):
    n = nlist[i]
    const[3] = kmin/(n)
    x0vec = fill_x0(n, iniconds, ron0)
    const[-1] = roff0

    obsindx = [2]
    dim = 2+n
    ligrec_downstream_sol_x = calculateCovariance(const, [x0vec[0]]+x0vec[2:], ligrec_downstream, obsindx, dim, t_max, dt)
    obsindx = [1,2]
    dim = 1+n
    ligrec_downstream_sol_lx = calculateCovariance(const, [x0vec[0]]+x0vec[2:], ligrec_downstream, obsindx, dim, t_max, dt)
    rate_ms[i] = MutualInformationRate_downstream(const, ligrec_downstream_sol_x, ligrec_downstream_sol_lx, x0vec, n)
    print(rate_ms[i])
    
np.savetxt(f'mi-ms-k1-{const[0]}-k2-{const[1]}-koff-{const[2]}-kon-{const[3]}-k4-{const[4]}-k5-{const[5]}-kswitch-{const[6]}-M-{M}-t_max-{t_max}.out', [rate_ms, nlist])
    
#%% Plotting 

kswitch = 1
kmin = 0.025
const={'k1': 1, 'k2': 0.01, 'koff': 1, 'kon': 1, 'k4': 500, 'k5': 50, 'kswitch': kswitch} 
const=list(const.values())

misb = -const[1]/2 + const[1]/2*np.sqrt(const[3]*const[2]*M/(const[1]*const[3]+const[0]*const[2])+1) #state-based mi rate
mirb = const[1]/2*(-1+(np.sqrt(2*const[3]*const[2]*M/(const[0]*const[2]+const[1]*const[3])+1))) #reaction-based mi rate 

const[3] = kmin

rate_ms, xaxis = np.loadtxt(f'mi-circle-k1-{const[0]}-k2-{const[1]}-koff-{const[2]}-kon-{const[3]}-k4-{const[4]}-k5-{const[5]}-kswitch-{const[6]}-M-{M}-t_max-{t_max}.out')

plt.figure(1, figsize= (10,10))
plt.plot(xaxis, rate_ms, 'o', color = 'blue', linewidth = lwd, markersize = 8)
plt.plot(xaxis, np.full(len(xaxis), mirb), color = 'green', linestyle = '--', linewidth = lwd)
plt.plot(xaxis, np.full(len(xaxis), misb), color = 'grey', linestyle = '--', linewidth = lwd)
plt.xlabel(r'Number of Receptor States $n$')
plt.ylabel(r'Mutual Information Rate $i^{lr}$')
plt.legend(('multi-state receptor', 'reaction-based', 'state-based'), loc = 'best')
plt.axis([1,40,0.01,0.019])
plt.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
# plt.savefig('mirate_multistate_dots.pdf',dpi=250)
plt.show()