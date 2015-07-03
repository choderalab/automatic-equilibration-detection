#!/usr/bin/env python

"""
Run a simulation of liquid argon at constant pressure.

"""

import os, os.path, copy
import netCDF4

from simtk import openmm, unit
from simtk.openmm import app

from openmmtools import testsystems, integrators

def minimize(system, positions):
    """
    Minimize the specified testsystem.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to minimize
    positions : simtk.unit.Quantity of size (nparticles,3) with units compatible with angstroms
        The initial positions to be minimized.
    
    Returns
    -------
    minimized_positions : simtk.unit.Quantity of size (nparticles,3) with units compatible with angstroms
        Minimized positions.

    """
    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    openmm.LocalEnergyMinimizer.minimize(context)
    final_positions = context.getState(getPositions=True).getPositions(asNumpy=True)
    del context, integrator
    return final_positions

def write_pdb(filename, positions):
    """
    Write PDB file for argon particles.

    Parameters
    ----------
    filename : str
        Filename to write PDB file to.
    positions : simtk.unit.Quantity of size (nparticles,3) with units compatible with angstroms
        Positions to write.

    """
    nparticles = positions.shape[0]
    outfile = open(filename, 'w')
    for particle in range(nparticles):
        outfile.write("ATOM  %5d  AR   AR     1    %8.3f%8.3f%8.3f\n" % (particle, positions[particle,0]/unit.angstrom, positions[particle,1]/unit.angstrom, positions[particle,2]/unit.angstrom))
    outfile.close()
    return

def create_netcdf_datastore(filename, system, positions, nreplicates, observation_interval):
    """
    Create (or resume from) NetCDF data storage file.

    Parameters
    ----------
    filename : str
        Filename of NetCDF file.
    system : simtk.openmm.System
        The system to minimize
    positions : simtk.unit.Quantity of size (nparticles,3) with units compatible with angstroms
        The initial positions used for all simulations
    nreplicates : int
        The number of simulation replicates to be performed
    obervation_interval : simtk.unit.Quantity with units compatible with ps
        Observation interval between frames.
    
    Returns
    -------
    ncfile : netCDF4.Dataset

    """

    if os.path.exists(filename):
        raise Exception("Datafile '%s' already exists." % filename)

    # Create a new file.
    ncfile = netCDF4.Dataset(filename, 'w', version='NETCDF4')

    # Determine some extra dimensions
    nparticles = positions.shape[0]
    
    # Initialize NetCDF file.
    ncfile.createDimension('replicate', nreplicates)
    ncfile.createDimension('iteration', 0) # unlimited number of iterations
    ncfile.createDimension('atom', nparticles) # number of atoms in system
    ncfile.createDimension('spatial', 3) # number of spatial dimensions
    ncfile.createDimension('singleton', 1)

    # Set global attributes.
    import time
    setattr(ncfile, 'title', 'liquid argon simulation density data')
    setattr(ncfile, 'CreationDate', time.ctime())

    # Store global data.
    ncvar = ncfile.createVariable('observation_interval', 'f4')
    ncvar.assignValue(observation_interval / unit.picoseconds)
    setattr(ncvar, 'units', 'ps')

    # Store initial positions.
    ncvar_positions = ncfile.createVariable('initial_positions', 'f4', ('atom','spatial'), zlib=True, chunksizes=(nparticles,3))
    setattr(ncvar_positions, 'units', 'nm')
    setattr(ncvar_positions, "long_name", "initial_positions[atom][spatial] is initial position of coordinate 'spatial' of atom 'atom' used for all simulations.")
    x = positions / unit.nanometers
    ncfile.variables['initial_positions'][:,:] = x[:,:]

    # Store system.
    ncvar_system = ncfile.createVariable('system', str, ('singleton',), zlib=True)
    setattr(ncvar_system, 'long_name', "system is the serialized OpenMM System used for all simulations")
    ncvar_system[0] = system.__getstate__()

    # Create storage for simulation data.
    ncvar_densities = ncfile.createVariable('reduced_density', 'f4', ('replicate','iteration'), zlib=True, chunksizes=(nreplicates,1))
    setattr(ncvar_densities, "long_name", "reduced_density[replicate][iteration] is the density (in reduced, dimensionless units) of iteration 'iteration' of replicate 'replicate'")
    ncvar_potential = ncfile.createVariable('reduced_potential', 'f4', ('replicate','iteration'), zlib=True, chunksizes=(nreplicates,1))
    setattr(ncvar_potential, "long_name", "reduced_potential[replicate][iteration] is the density (in kT) of iteration 'iteration' of replicate 'replicate'")

    ncfile.sync()

    return ncfile
