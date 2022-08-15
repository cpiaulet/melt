# melt
Energy balance calculation for rocky planet interiors under the influence of tides

*melt* (siMple Energy baLance with Tides) is a tool to calculate equilibrium between tidal and convective energy profiles in the interior of a rocky planet. This code is based on a simple Maxwell rheology prescription as described in Barr et al. (2018), which follows Moore (2003) and  Dobos and Turner (2015).

If you use this code, please cite Caroline Piaulet as well as the following papers: 
* https://ui.adsabs.harvard.edu/abs/2018A%26A...613A..37B/abstract
* https://ui.adsabs.harvard.edu/abs/2003JGRE..108.5096M/abstract
* https://ui.adsabs.harvard.edu/abs/2015ApJ...804...41D/abstract


## Installation
You can install *melt* from GitHub:

    git clone https://github.com/cpiaulet/melt.git
    cd melt
    python setup.py install

### Dependencies
The only dependencies of *melt* are *NumPy* and *astropy*!

### Inputs
All the required inputs are set in a dictonary which should contain:
* the planet name, 
* mass in Earth masses, 
* radius in Earth radii, 
* orbital period in days,
* eccentricity

See ```melt_analysis/python plot_energy_balance_example.py``` for an example of how to create a HeatedPlanet object from this dictionary, calculate convective and tidal flux profiles, and plot them! 

### Example
Copy the two folders in ```melt/example/``` wherever in your installation you want to run your code;

In ```melt_analysis/```, run ```python plot_energy_balance_example.py``` and inspect the results in ```melt_results/```.
