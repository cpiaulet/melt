#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 12:59:20 2021

@author: caroline

Utility functions
Energy balance model to calculate tidal energy flux
and temperature in the convecting rock mantle
following Barr+ 2018
"""

# Import modules ---------- 
import numpy as np
import matplotlib.pyplot as plt
from auxbenneke.constants import Mearth, Rearth, G, day
import astropy.constants as const

#%% Constants
sec_in_a_day = 24. * 3600.

# rock properties
Qstar = 333*1e3 # J/mol
ktherm = 3.2 # W/m/K
alpha = 3e-5 # [K^-1]
Cp = 1200 # J/kg
rho_rock = 5000 # kg/m3

# rock temperatures
Ts = 1600. # K, solidus point
Tb = 1800. # K, breakdown point, volume fraction of solid rock crystals = volume fraction of melted rock
Tl = 2000. # K, liquidus point

eta0 = 2.13e6 # Pa s

Rg = 8.314 #J /mol /K gas constant

mu1 = 8.2e4 # K
mu2 = -40.6

#%% Functions

# functions for eta and mu as a function of T

# rock
def calc_f(T):
    """ 
    volume fraction of rock melt present
    following Barr+ 2018
    """
    f = (1./(Tl-Ts)) * (T-Ts)
    return f

def eta_visc_rock(T, B = 25.):
    """ 
    viscosity of rock as a function of T
    Eqns (13-16) Barr+ 2018
    B dimensionless between 10 and 40
    """
    if T <= Ts:
        eta = eta0 * np.exp(Qstar/(Rg * T))
    elif (T > Ts) and (T <= Tb):
        f = calc_f(T)
        eta = eta0 * np.exp(Qstar/(Rg * T)) * np.exp(-B * f)
        
    elif (T > Tb) and (T <= Tl):
        f = calc_f(T)
        eta = 1e-7 * np.exp(40000./T) * (1.35 * f - 0.35)**(-5./2.)
        
    elif T > Tl:
        f = 1.
        eta = 1e-7 * np.exp(40000./T) * (1.35 * f - 0.35)**(-5./2.)
    return eta

def mu_shearmod_rock(T):
    """ 
    shear modulus of rock as a function of T
    Section 3.4.1 Barr+ 2018
    """
    if T <= Ts:
        mu = 50e9
    elif (T > Ts) and (T <= Tb):
        mu = 10**(mu1/T + mu2)
        ## I add a constant for continuity at Ts
        A0 = 50e9/(10**(mu1/Ts + mu2))
        mu = A0 * mu
        
    elif T > Tb:
        mu = 1e-7

    return mu  

def get_Ft_Fc_intersect(Ft, Fc):
    '''
    get intersection of profiles
    '''
    
    idx = np.argwhere(np.diff(np.sign(Ft - Fc))).flatten()
    return idx, Fc[idx]

    
