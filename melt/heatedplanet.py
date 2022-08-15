#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 12:59:20 2021

@author: caroline

HeatedPlanet Class definition

Energy balance model to calculate tidal energy flux
and temperature in the convecting rock mantle
following Barr+ 2018
"""

# Import modules ---------- 
import numpy as np
import astropy.constants as const
import melt.heatutils as hu

#%% Constants
sec_in_a_day = 24. * 3600.

#%% Object definition
class HeatedPlanet(object):
    '''
    Create a model of a tidally-heated planet
    '''
    def __init__(self, pla_dict):
        '''
        pla_dict should contain: 
            planet or run name
            mass in Earth masses
            radius in Earth radii
            orbital period in days
            eccentricity
        '''
        self.name = pla_dict["name"]
        self.mp = pla_dict["mp_earth"] * const.M_earth.value
        self.rp = pla_dict["rp_earth"] * const.R_earth.value
        self.per = pla_dict["per_day"] * sec_in_a_day
        self.ecc = pla_dict["ecc"]
    
    def calc_im_k2(self, T, B):
        """
        Calculate the imaginary part of the k2 Love number for a Maxwell
        viscoelastic rheolgy
        inputs in SI units
        Eq. (6) in Barr et al. 2018
        """
        
        omega = 2.*np.pi/self.per 
        rhomoy = self.mp/((4./3.)* np.pi * self.rp**3.)
        g = const.G.value * self.mp /(self.rp**2.)
        
        rho_bar = rhomoy # alternatively, use rho_rock        
        
        eta_visc = hu.eta_visc_rock(T, B=B)
        mu_shearmod = hu.mu_shearmod_rock(T)
        
        num = 57. * eta_visc * omega 
        denom_parenth1 = 1. + (19. * mu_shearmod/(2. * rho_bar * g * self.rp))
        denom_parenth2 = 1. + (denom_parenth1**2. * ((eta_visc**2. * omega**2.)/(mu_shearmod**2.)))
        denom = 4. * rho_bar * g * self.rp * denom_parenth2
        
        im_k2 = num/denom
        
        return im_k2

    def calc_dEdt_tidal(self, T, B):
        """ 
        energy produced by tidal dissipation in a synchronously rotating body
        Eq. (5) in Barr et al. 2018    
        """
        im_k2 = self.calc_im_k2(T, B)
        omega = 2.*np.pi/self.per 
        
        term = (self.rp**5. * omega**5. * self.ecc**2.)/const.G.value
        dEdt = (-21./2.) * im_k2 * term 
        
        return dEdt
    
    def calc_F_tidal(self, T, B):
        """ 
        Globally averaged tidal heat flux
        Eq. (7) in Barr et al. 2018        
        """        
        dEdt_tidal = self.calc_dEdt_tidal(T,B)
        F_tidal = dEdt_tidal/(4.*np.pi * self.rp**2.)
        return -F_tidal

    
    def calc_F_conv(self, T, B):
        """ 
        Convective heat flux
        Eq. (11) in Barr et al. 2018  
        """
        g = const.G.value * self.mp /(self.rp**2.)
        kappatherm = hu.ktherm/(hu.rho_rock * hu.Cp)
        
        term1 = hu.Qstar/(hu.Rg * T**2.)
        term2 = (hu.rho_rock * g * hu.alpha * hu.ktherm**3./(kappatherm * hu.eta_visc_rock(T, B=B)))
        
        Fconv = 0.53 * term1**(-4./3.) * term2**(1./3.)
        return Fconv
    
    def calc_F_tidal_profile(self, T_arr, B):
        Ftidal = np.zeros_like(T_arr)
        for i, t in enumerate(T_arr):
            Ftidal[i] = self.calc_F_tidal(T=t, B=B)
        
        return Ftidal

    def calc_F_conv_profile(self, T_arr, B):
        Fconv = np.zeros_like(T_arr)
        for i, t in enumerate(T_arr):
            Fconv[i] = self.calc_F_conv(T=t, B=B)
        
        return Fconv
        
