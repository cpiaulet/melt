#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:12:41 2022

@author: caroline
Calculate and plot energy balance profiles using melt
"""

import numpy as np
import matplotlib.pyplot as plt
import melt.heatedplanet as hp
import melt.heatutils as hu

# Test planet

pla_dict = dict()

pla_dict["name"] = "TEST"
pla_dict["mp_earth"] = 0.9
pla_dict["rp_earth"] = 1.032
pla_dict["per_day"] =  2.7534360
pla_dict["ecc"] = 0.019

hp1 = hp.HeatedPlanet(pla_dict)

# Array of temperatures
temp = np.linspace(1400, 1800, 1000)
Ft_B20 = hp1.calc_F_tidal_profile(temp, B=20)
Ft_B30 = hp1.calc_F_tidal_profile(temp, B=30)
Fc_B20 = hp1.calc_F_conv_profile(temp, B=20)
Fc_B30 = hp1.calc_F_conv_profile(temp, B=30)


#%% plot results 


fig, ax = plt.subplots(1,1)  
# ax.axvline(1671, color="k", ls="--")      
ax.axvline(hu.Ts, color="r", ls="-")      
ax.axvline(hu.Tb, color="r", ls="-")      
 
# ax.axhline(np.log10(0.38), color="k", ls="--")      
ax.plot(temp, np.log10(Ft_B20), color="g", lw=2, label = "Tidal flux B=20")
ax.plot(temp, np.log10(Ft_B30), color="b", label = "Tidal flux B=30")
ax.plot(temp, np.log10(Fc_B20), color="g", lw=2, ls="--",label = "Convective flux B=20")
ax.plot(temp, np.log10(Fc_B30), color="b", ls="--", label = "Convective flux B=30")

ax.axhspan(np.log10(1.),np.log10(2.), color="gray", alpha=0.5, lw=0)
ax.text(1450,np.log10(1.)+0.05, "Io's tidal heat flux", fontsize=12)

# get intersection index
id20, _ = hu.get_Ft_Fc_intersect(Ft_B20, Fc_B20)
id30, _ = hu.get_Ft_Fc_intersect(Ft_B30, Fc_B30)
ax.scatter(temp[id20], np.log10(Fc_B20[id20]), color="g")
ax.scatter(temp[id30], np.log10(Fc_B30[id30]), color="b")

# ax.legend()
ax.set_xlabel("Mantle temperature [K]")
ax.set_ylabel(r"log$_{10}$ F [W/m$^2$]")
ax.set_xlim(1400, 1800)
# ax.set_ylim(9, 18)

ax.text(1420, -3.3, "B=30", color="b")
ax.text(1420, -3.0, "B=20", color="g")

ax.text(1500., -1.1, r"F$_\mathrm{conv}$", rotation=10)
ax.text(1500., 0.75, r"F$_\mathrm{tidal}$", rotation=20)

fig.savefig("../melt_results/energy_balance_example.pdf")