#!/usr/bin/env python

"""
Analyze multiple replicates of a simulation of liquid argon at constant pressure.

"""

import numpy as np
import os.path
import matplotlib.pyplot as plt

data_directory = 'data'

nreplicates = 35
niterations = 10000

reduced_density_it = np.zeros([nreplicates, niterations+1], np.float64)
reduced_potential_it = np.zeros([nreplicates, niterations+1], np.float64)


# Read data
for replicate in range(nreplicates):
    filename = os.path.join(data_directory, '%05d.output' % replicate)
    lines = open(filename, 'r').readlines()
    for iteration in range(niterations+1):
        elements = lines[iteration].split()
        reduced_density_it[replicate, iteration] = float(elements[1])
        reduced_potential_it[replicate, iteration] = float(elements[2])

print reduced_density_it

# Plot

A_it = reduced_density_it

# TODO: Read this in?
timestep = 0.002846993438 # ps
nsteps_per_iteration = 25
iterations_per_ns = timestep * nsteps_per_iteration * 1000 # ns per iteration

x = np.arange(niterations+1) / iterations_per_ns # ns
A_t = A_it.mean(0)
dA_t = A_it.std(0)

print "<A> = "
print A_t
print "d<A> = "
print dA_t

# Save plot to PDF.
filename = 'argon-density.pdf' # PDF file to write

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(filename)

import pylab
#import seaborn as sns

#
# FIGURE
#

fontsize=10

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
#pylab.errorbar(range(niterations+1), A_t, yerr=dA_t)
#for replicate in range(nreplicates):
#    pylab.plot(x, A_it[replicate,:])
#pylab.errorbar(x, A_t, yerr=2*dA_t, fmt='ko')
#pylab.fill_between(x, A_t+2*dA_t, A_t-2*dA_t, facecolor='grey', alpha=0.5)
pylab.fill_between(x[:nmax:nskip], A_t[:nmax:nskip]+2*dA_t[:nmax:nskip], A_t[:nmax:nskip]-2*dA_t[:nmax:nskip], facecolor='grey', edgecolor='grey', alpha=0.5, linewidth=0)
pylab.plot(x[:nmax], A_t[:nmax], 'k-')

#pylab.xlabel('simulation time / ns')
pylab.ylabel(r'density $\rho^*$', fontsize=fontsize)
# Adjust axes.
oldaxis = pylab.axis()
pylab.axis([0, x[:nmax].max(), oldaxis[2], oldaxis[3]])

subplot.set_xticklabels([]) # no x-tick labels
subplot.set_yticklabels([0.8, 0.9, 1.0, 1.1, 1.2]) # no x-tick labels

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

Acumavg_mean_t = Acumavg_it.mean(0)
Aburnin_mean_t = Aburnin_it.mean(0)

Acumavg_std_t = Acumavg_it.std(0)
Aburnin_std_t = Aburnin_it.std(0)



subplot = pylab.subplot(212)
pylab.subplots_adjust(hspace=0.001)
pylab.hold(True)

pylab.fill_between(x[:nmax:nskip], Acumavg_mean_t[:nmax:nskip]+2*Acumavg_std_t[:nmax:nskip], Acumavg_mean_t[:nmax:nskip]-2*Acumavg_std_t[:nmax:nskip], facecolor='red', edgecolor='red', alpha=0.5, linewidth=0)
pylab.fill_between(x[Nequil:nmax:nskip], Aburnin_mean_t[Nequil:nmax:nskip]+2*Aburnin_std_t[Nequil:nmax:nskip], Aburnin_mean_t[Nequil:nmax:nskip]-2*Acumavg_std_t[Nequil:nmax:nskip], facecolor='blue', edgecolor='blue', alpha=0.5, linewidth=0)

pylab.plot(x[:nmax], Acumavg_mean_t[:nmax], 'r-')
pylab.plot(x[Nequil:nmax], Aburnin_mean_t[Nequil:nmax], 'b-')
pylab.legend(['cumulative average', 'discarding first %d samples to equilibration' % Nequil], fontsize=fontsize)

pylab.xlabel('simulation time / ns', fontsize=fontsize)
pylab.ylabel(r'reduced density $\rho^*$', fontsize=fontsize)

subplot.set_yticklabels([0.8, 0.9, 1.0, 1.1, 1.2]) # no x-tick labels

# Adjust axes.
oldaxis = pylab.axis()
pylab.axis([0, x[:nmax].max(), oldaxis[2], oldaxis[3]])

#pylab.show()
figure.tight_layout()

# Write figure to PDF
pp.savefig()
pp.close()
