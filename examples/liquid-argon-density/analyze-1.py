#!/usr/bin/env python

"""
Analyze multiple replicates of a simulation of liquid argon at constant pressure.

"""

import numpy as np
import os.path
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
#import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import netCDF4

#from equilibration import statisticalInefficiency_geyer as statisticalInefficiency
from equilibration import statisticalInefficiency_multiscale as statisticalInefficiency

netcdf_filename = 'data.nc'
figure_directory = 'figures'

#
# Ensure figure directory exists
#

if not os.path.exists(figure_directory):
    print "Creating figures directory..."
    os.makedirs(figure_directory)
#
# READ DATA
#

# Open NetCDF file
ncfile = netCDF4.Dataset(netcdf_filename, 'r')
[nreplicates, niterations] = ncfile.variables['reduced_density'].shape
observation_interval = ncfile.variables['observation_interval'].getValue() # ps
ps_per_iteration = observation_interval

# Select data to analyze.
A_it = np.array(ncfile.variables['reduced_density'][:,:], np.float64)

# Compute true expectation using mean of second half of all simulations.
t0 = int(niterations/2)
true_expectation = A_it[:,t0:].mean(1).mean(0)

#
# PLOT RELAXATION OF REDUCED DENSITIES
#

# Write initial reduced densities.
print "Initial reduced densities:"
print A_it[0,0]

x = np.arange(niterations+1) * ps_per_iteration
A_t = A_it.mean(0)
dA_t = A_it.std(0)

# Save plot to PDF.
filename = os.path.join(figure_directory, 'argon-density.pdf') # PDF file to write

pp = PdfPages(filename)

#
# FIGURE 1
# Plot average density decay from initial conditions, along with cumulative average and average with some data discarded to equilibrium.
#

fontsize = 10
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : fontsize}

plt.rc('font', **font)

figure = pylab.figure(figsize=(7.5, 3.5), dpi=300)

nskip = 10 # number of samples to skip in uncertainty shading
nmax = 4000 # maximum samples to show

#
# TOP PANEL: AVERAGE DENSITY AND 95% CONFIDENCE INTERVAL
#

subplot = pylab.subplot(211)
pylab.hold(True)

pylab.fill_between(x[:nmax:nskip], A_t[:nmax:nskip]+2*dA_t[:nmax:nskip], A_t[:nmax:nskip]-2*dA_t[:nmax:nskip], facecolor='grey', edgecolor='grey', alpha=0.5, linewidth=0)
pylab.plot(x[:nmax], A_t[:nmax], 'k-')

pylab.ylabel(r'reduced density $\rho^*$', fontsize=fontsize)

# Adjust axes.
oldaxis = pylab.axis()
pylab.axis([0, x[:nmax].max(), oldaxis[2], oldaxis[3]])

subplot.set_xticks([]) # no x-tick labels

#
# BOTTOM PANEL: CUMULATIVE AVERAGE WITH INITIAL BURN-IN PERIOD DISCARDED
#

Nequil = 100
Acumavg_it = np.zeros([nreplicates, niterations+1], np.float64)
Aburnin_it = np.zeros([nreplicates, niterations+1], np.float64)
for i in range(nreplicates):
    for t in range(niterations+1):
        Acumavg_it[i,t] = A_it[i,0:(t+1)].mean()

        Ninit = Nequil
        if t < Nequil:
            Ninit = 0
        Aburnin_it[i,t] = A_it[i,Ninit:(t+1)].mean()

# Compute reverse cumulative average.
Arevcum_it = np.zeros([nreplicates, niterations+1], np.float64)
for i in range(nreplicates):
    for t0 in range(niterations+1):
        Arevcum_it[i,t0] = A_it[i,t0:].mean()

Acumavg_mean_t = Acumavg_it.mean(0)
Aburnin_mean_t = Aburnin_it.mean(0)
Arevcum_mean_t = Arevcum_it.mean(0)

Acumavg_std_t = Acumavg_it.std(0)
Aburnin_std_t = Aburnin_it.std(0)
Arevcum_std_t = Arevcum_it.std(0)

# Create a subplot.
subplot = pylab.subplot(212)
pylab.subplots_adjust(hspace=0.001)
pylab.hold(True)

# Draw shaded confidence intervals
pylab.fill_between(x[:nmax:nskip], Acumavg_mean_t[:nmax:nskip]+2*Acumavg_std_t[:nmax:nskip], Acumavg_mean_t[:nmax:nskip]-2*Acumavg_std_t[:nmax:nskip], facecolor='red', edgecolor='red', alpha=0.5, linewidth=0)
pylab.fill_between(x[Nequil:nmax:nskip], Aburnin_mean_t[Nequil:nmax:nskip]+2*Aburnin_std_t[Nequil:nmax:nskip], Aburnin_mean_t[Nequil:nmax:nskip]-2*Aburnin_std_t[Nequil:nmax:nskip], facecolor='blue', edgecolor='blue', alpha=0.5, linewidth=0)

pylab.plot([x[0], x[nmax]], [true_expectation, true_expectation], 'k--') # true expectation (from all data)
pylab.plot(x[:nmax], Acumavg_mean_t[:nmax], 'r:') # cumulative average from beginning
pylab.plot(x[Nequil:nmax], Aburnin_mean_t[Nequil:nmax], 'b-') # cumulative average after discarding some initial data to equilibrium

pylab.legend(['true expectation', 'cumulative average', 'discarding first %d samples to equilibration' % Nequil], fontsize=fontsize, frameon=False)

pylab.xlabel('simulation time / ps', fontsize=fontsize)
pylab.ylabel(r'reduced density $\rho^*$', fontsize=fontsize)

# Adjust axes.
oldaxis = pylab.axis()
pylab.axis([0, x[:nmax].max(), oldaxis[2], oldaxis[3]])

figure.tight_layout()

# Write figure to PDF
pp.savefig()
pp.close()
