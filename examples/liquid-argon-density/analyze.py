#!/usr/bin/env python

"""
Analyze multiple replicates of a simulation of liquid argon at constant pressure.

"""

import numpy as np
import os.path

data_directory = 'data'

nreplicates = 100
niterations = 1000

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
x = range(niterations+1)
A_t = reduced_density_it.mean(0)
dA_t = reduced_density_it.std(0)

print "<A> = "
print A_t
print "d<A> = "
print dA_t

import pylab
pylab.figure()
pylab.hold(True)
#pylab.errorbar(range(niterations+1), A_t, yerr=dA_t)
for replicate in range(nreplicates):
    pylab.plot(x, reduced_density_it[replicate,:])

pylab.errorbar(x, A_t, yerr=2*dA_t, fmt='ko')
pylab.show()
