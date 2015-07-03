#!/usr/bin/env python

"""
Run a simulation of liquid argon at constant pressure.

"""

import os, os.path, copy

from simtk import openmm
from simtk import unit
from simtk.openmm import app

from openmmtools import testsystems, integrators

import utils # simulation and analysis utilities

# Set conditions for simulation.
# Roughly corresponds to conditions from http://www.cstl.nist.gov/srs/LJ_PURE/mc.htm
nparticles = 500
mass = 39.9 * unit.amu
sigma = 3.4 * unit.angstrom
epsilon = 0.238 * unit.kilocalories_per_mole
#reduced_density = 0.860     # reduced_density = density * (sigma**3)
reduced_density = 0.960     # reduced_density = density * (sigma**3)
reduced_temperature = 0.850 # reduced_temperature = kB * temperature / epsilon
reduced_pressure = 1.2660   # reduced_pressure = pressure * (sigma**3) / epsilon

netcdf_filename = 'data.nc' # NetCDF file to store data in

r0 = 2.0**(1./6.) * sigma   # minimum potential distance for Lennard-Jones interaction
characteristic_timescale = unit.sqrt((mass * r0**2) / (72 * epsilon)) # characteristic timescale for bound Lennard-Jones interaction
                                                            # http://borisv.lk.net/matsc597c-1997/simulations/Lecture5/node3.html
timestep = 0.01 * characteristic_timescale # integrator timestep

# From http://www.cstl.nist.gov/srs/LJ_PURE/md.htm
#characteristic_timescale = unit.sqrt(mass * sigma**2 / epsilon)
#timestep = 0.05 * characteristic_timescale

print "characteristic timescale = %.3f ps" % (characteristic_timescale / unit.picoseconds)
print "timestep = %.12f ps" % (timestep / unit.picoseconds)

#collision_rate = 5.0 / unit.picoseconds # collision rate for Langevin thermostat
collision_rate = 1.5 / characteristic_timescale # collision rate for Langevin thermostat
barostat_frequency = 100 # number of steps between barostat updates

# Set parameters for number of simulation replicates, number of iterations per simulation, and number of steps per iteration.
nreplicates = 500
niterations = 10000
nsteps_per_iteration = barostat_frequency
observation_interval = timestep * nsteps_per_iteration

# Compute real units.
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
density = reduced_density / (sigma**3)
temperature = reduced_temperature * epsilon / kB
pressure = reduced_pressure * epsilon / (sigma**3)
kT = kB * temperature

# Create the Lennard-Jones fluid.
testsystem = testsystems.LennardJonesFluid(nparticles=nparticles, mass=mass, sigma=sigma, epsilon=epsilon, reduced_density=reduced_density)

# Construct initial positions by minimization.
print "Minimizing positions..."
initial_positions = utils.minimize(testsystem.system, testsystem.positions)

# Write initial positions.
print "Writing initial positions to initial.pdb"
utils.write_pdb('initial.pdb', initial_positions)

# Create NetCDF file to store data.
ncfile = utils.create_netcdf_datastore(netcdf_filename, testsystem.system, initial_positions, nreplicates, niterations, observation_interval)

# Run replicates of the simulation.
for replicate in range(nreplicates):
    # Reconstitute System object.
    print "Reconstituting System object..."
    system = openmm.System()
    system.__setstate__(str(ncfile.variables['system'][0]))

    # Add a barostat to the system.
    # NOTE: This is added to a new copy of the system to ensure barostat random number seeds are unique.
    print "Adding barostat..."
    barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency)
    system.addForce(barostat)

    # Create integrator
    print "Creating integrator..."
    integrator = integrators.VVVRIntegrator(temperature, collision_rate, timestep)

    # Create context.
    print "Creating Context..."
    context = openmm.Context(system, integrator)

    # Set initial conditions.
    print "Setting initial positions..."
    context.setPositions(initial_positions)
    print "Setting initial velocities appropriate for temperature..."
    context.setVelocitiesToTemperature(temperature)

    # Record initial data.
    state = context.getState(getEnergy=True)
    reduced_volume = state.getPeriodicBoxVolume() / (nparticles * sigma**3)
    reduced_density = 1.0 / reduced_volume
    reduced_potential = state.getPotentialEnergy() / kT
    print "replicate %5d / %5d : initial                 : density %8.3f | potential %8.3f" % (replicate, nreplicates, reduced_density, reduced_potential)
    ncfile.variables['reduced_density'][replicate,0] = reduced_density
    ncfile.variables['reduced_potential'][replicate,0] = reduced_potential

    # Run simulation.
    for iteration in range(1,niterations+1):
        # Integrate the simulation.
        integrator.step(nsteps_per_iteration)

        # Record data.
        state = context.getState(getEnergy=True)
        reduced_volume = state.getPeriodicBoxVolume() / (nparticles * sigma**3)
        reduced_density = 1.0 / reduced_volume
        reduced_potential = state.getPotentialEnergy() / kT

        if (iteration % 1000) == 0:
            print "replicate %5d / %5d : iteration %5d / %5d : density %8.3f | potential %8.3f" % (replicate, nreplicates, iteration, niterations, reduced_density, reduced_potential)
            ncfile.sync()

        # Store data.
        ncfile.variables['reduced_density'][replicate,iteration] = reduced_density
        ncfile.variables['reduced_potential'][replicate,iteration] = reduced_potential

    # Clean up.
    del context, integrator

    # Ensure all data is flushed to NetCDF file.
    ncfile.sync()

    print ""

# Clean up.
ncfile.close()

