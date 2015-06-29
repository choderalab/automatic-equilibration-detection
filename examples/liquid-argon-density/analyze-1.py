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

#from equilibration import statisticalInefficiency_geyer as statisticalInefficiency
from equilibration import statisticalInefficiency_multiscale as statisticalInefficiency

data_directory = 'data'
figure_directory = 'figures'

nreplicates = 100
niterations = 10000

# TODO: Read this in?
timestep = 0.002846993438 # ps
nsteps_per_iteration = 25
iterations_per_ns = timestep * nsteps_per_iteration * 1000 # ns per iteration

#
# Ensure figure directory exists
#

if not os.path.exists(figure_directory):
    print "Creating figures directory..."
    os.makedirs(figure_directory)
#
# READ DATA
#

reduced_density_it = np.zeros([nreplicates, niterations+1], np.float64)
reduced_potential_it = np.zeros([nreplicates, niterations+1], np.float64)

for replicate in range(nreplicates):
    filename = os.path.join(data_directory, '%05d.output' % replicate)
    lines = open(filename, 'r').readlines()
    for iteration in range(niterations+1):
        elements = lines[iteration].split()
        reduced_density_it[replicate, iteration] = float(elements[1])
        reduced_potential_it[replicate, iteration] = float(elements[2])

A_it = reduced_density_it

# Compute true expectation using mean of second half of all simulations.
t0 = int(niterations/2)
true_expectation = A_it[:,t0:].mean(1).mean(0)

# Compute true expectation using all data. # DEBUG
#true_expectation = A_it.mean().mean()
#print "true expectation (from all data) : %f" % true_expectation

#
# PLOT RELAXATION OF REDUCED DENSITIES
#

# Write initial reduced densities.
print "Initial reduced densities:"
print A_it[0,0]

x = np.arange(niterations+1) / iterations_per_ns # ns
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

pylab.xlabel('simulation time / ns', fontsize=fontsize)
pylab.ylabel(r'reduced density $\rho^*$', fontsize=fontsize)

# Adjust axes.
oldaxis = pylab.axis()
pylab.axis([0, x[:nmax].max(), oldaxis[2], oldaxis[3]])

figure.tight_layout()

# Write figure to PDF
pp.savefig()
pp.close()
