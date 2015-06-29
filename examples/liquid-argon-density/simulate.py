#!/usr/bin/env python

"""
Run a simulation of liquid argon at constant pressure.

"""

import os, os.path, copy

from simtk import openmm
from simtk import unit
from simtk.openmm import app

from openmmtools import testsystems, integrators

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

data_directory = 'data'     # Directory in which data is to be stored

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
barostat_frequency = 25 # number of steps between barostat updates

# Set parameters for number of simulation replicates, number of iterations per simulation, and number of steps per iteration.
nreplicates = 100
niterations = 10000
nsteps_per_iteration = 25

# Compute real units.
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
density = reduced_density / (sigma**3)
temperature = reduced_temperature * epsilon / kB
pressure = reduced_pressure * epsilon / (sigma**3)
kT = kB * temperature

# Create the Lennard-Jones fluid.
testsystem = testsystems.LennardJonesFluid(nparticles=nparticles, mass=mass, sigma=sigma, epsilon=epsilon, reduced_density=reduced_density)

# Construct initial positions by minimization.
print "Minimizing to obtain initial positions..."
#integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
integrator = integrators.VVVRIntegrator(temperature, collision_rate, timestep)
context = openmm.Context(testsystem.system, integrator)
context.setPositions(testsystem.positions)
openmm.LocalEnergyMinimizer.minimize(context)
state = context.getState(getPositions=True)
initial_positions = state.getPositions(asNumpy=True)
del context, integrator, state

# DEBUG: Write initial positions.
outfile = open("initial.pdb", 'w')
for particle in range(nparticles):
    outfile.write("ATOM  %5d  AR   AR     1    %8.3f%8.3f%8.3f\n" % (particle, initial_positions[particle,0]/unit.angstrom, initial_positions[particle,1]/unit.angstrom, initial_positions[particle,2]/unit.angstrom))
outfile.close()

# Make directory to store data.
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

# Run replicates of the simulation.
for replicate in range(nreplicates):
    # Make a new copy of the system.
    print "Making a deep copy of the system..."
    system = copy.deepcopy(testsystem.system)

    # Add a barostat to the system.
    # NOTE: This is added to a new copy of the system to ensure barostat random number seeds are unique.
    print "Adding barostat..."
    barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency)
    system.addForce(barostat)

    # Open output file.
    print "Opening output file..."
    output_filename = os.path.join(data_directory, '%05d.output' % replicate)
    outfile = open(output_filename, 'w')

    # Create integrator
    print "Creating LangevinIntegrator..."
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

    # Create context.
    print "Creating Context..."
    context = openmm.Context(system, integrator)

    # Set initial conditions.
    print "Setting initial positions..."
    context.setPositions(initial_positions)
    print "Setting initial velocities appropriate for temperature..."
    context.setVelocitiesToTemperature(temperature)

    # DEBUG: Write out initial conditions.
    print "Serializing..."
    def write_file(filename, contents):
        file = open(filename, 'w')
        file.write(contents)
        file.close()
    from simtk.openmm import XmlSerializer
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    write_file(os.path.join(data_directory, '%05d.system.xml' % replicate), XmlSerializer.serialize(system))
    write_file(os.path.join(data_directory, '%05d.state.xml' % replicate), XmlSerializer.serialize(state))
    write_file(os.path.join(data_directory, '%05d.integrator.xml' % replicate), XmlSerializer.serialize(integrator))

    # Record initial data.
    state = context.getState(getEnergy=True)
    reduced_volume = state.getPeriodicBoxVolume() / (nparticles * sigma**3)
    reduced_density = 1.0 / reduced_volume
    reduced_potential = state.getPotentialEnergy() / kT
    print "replicate %5d / %5d : initial                 : density %8.3f | potential %8.3f" % (replicate, nreplicates, reduced_density, reduced_potential)
    outfile.write('%8d %12.6f %12.6f\n' % (0, reduced_density, reduced_potential))

    # Run simulation.
    for iteration in range(niterations):
        # Integrate the simulation.
        integrator.step(nsteps_per_iteration)

        # Record data.
        state = context.getState(getEnergy=True)
        reduced_volume = state.getPeriodicBoxVolume() / (nparticles * sigma**3)
        reduced_density = 1.0 / reduced_volume
        reduced_potential = state.getPotentialEnergy() / kT

        if ((iteration + 1) % 100) == 0:
            print "replicate %5d / %5d : iteration %5d / %5d : density %8.3f | potential %8.3f" % (replicate, nreplicates, iteration+1, niterations, reduced_density, reduced_potential)
        outfile.write('%8d %12.6f %12.6f\n' % (iteration+1, reduced_density, reduced_potential))

    # Clean up.
    del context, integrator
    outfile.close()

    print ""


