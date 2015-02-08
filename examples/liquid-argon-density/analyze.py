#!/usr/bin/env python

"""
Analyze multiple replicates of a simulation of liquid argon at constant pressure.

"""

import numpy as np
import os.path

data_directory = 'data'

nreplicates = 35
niterations = 10000

reduced_density_it = np.zeros([nreplicates, niterations+1], np.float64)
reduced_potential_it = np.zeros([nreplicates, niterations+1], np.float64)


# Read data
for replicate in range(nreplicates):
    filename = os.path.join(data_directory, 'output.%05d' % replicate)
    lines = open(filename, 'r').readlines()
    for iteration in range(niterations+1):
        elements = lines[iteration].split()
        reduced_density_it[replicate, iteration] = float(elements[1])
        reduced_potential_it[replicate, iteration] = float(elements[2])

print reduced_density_it

# Plot

A_it = reduced_density_it

x = range(niterations+1)
A_t = A_it.mean(0)
dA_t = A_it.std(0)

print "<A> = "
print A_t
print "d<A> = "
print dA_t

import pylab
pylab.figure()

pylab.subplot(211)
pylab.hold(True)
#pylab.errorbar(range(niterations+1), A_t, yerr=dA_t)
#for replicate in range(nreplicates):
#    pylab.plot(x, A_it[replicate,:])
#pylab.errorbar(x, A_t, yerr=2*dA_t, fmt='ko')
pylab.plot(x, A_t, 'ko')
pylab.fill_between(x, A_t+2*dA_t, A_t-2*dA_t, facecolor='grey', alpha=0.5)

pylab.xlabel('number of samples')
pylab.ylabel(r'reduced density $\rho^*$')

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

pylab.subplot(212)
pylab.hold(True)
pylab.plot(x, Acumavg_it.mean(0), 'k.')
pylab.plot(x, Aburnin_it.mean(0), 'r.')
pylab.legend(['cumulative average', 'discarding first %d samples to equilibration' % Nequil])

pylab.xlabel('number of samples')
pylab.ylabel(r'reduced density $\rho^*$')

pylab.show()
