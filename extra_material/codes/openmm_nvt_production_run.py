
#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os

# OpenMM Imports
import openmm as mm
import openmm.app as app

# ParmEd Imports
from parmed import load_file, unit as u
from parmed.openmm import StateDataReporter
from mdtraj import reporters

import argparse

parser = argparse.ArgumentParser(description='production run with openmm, NVT system')
parser.add_argument("-o","--out_dir", default=os.path.abspath(os.curdir), help="output directory, default: current directory")
parser.add_argument("-t","--topology_file", required=True, help="name of the amber topology file")
parser.add_argument("-i","--inpcrd_file", required=True, help="name of the amber initial coordinate file")
parser.add_argument("-p","--prefix", required=True, help="prefix for names of the output files")
parser.add_argument("-r","--run_time", default=100, type=int, help="time of the simulation in ns. Only accept ")

args = parser.parse_args()

prefix = args.prefix
topology_file = os.path.abspath(args.topology_file)
inpcrd_file = os.path.abspath(args.inpcrd_file)
output_dir = os.path.abspath(args.out_dir)

# Load the Amber files
print('Loading AMBER files...')
pep_solv = load_file(topology_file, inpcrd_file)

# Create the OpenMM system
print('Creating OpenMM System')
system = pep_solv.createSystem(nonbondedMethod=app.PME,
                                nonbondedCutoff=8.0*u.angstroms,
                                constraints=app.HBonds,
)

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        300*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
                        4.0*u.femtoseconds, # Time step
)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('CUDA')
prop = dict(CudaPrecision='mixed') # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(pep_solv.topology, system, integrator, platform, prop)

# Set the particle positions
sim.context.setPositions(pep_solv.positions)

# Minimize the energy
print('Minimizing energy')
sim.minimizeEnergy(maxIterations=500)

sim.saveState('%s/%s_state_000.xml' % (output_dir, prefix))
# Save the state, coordinates, and stats into run_num increments (10 ns each run)
run_num = int(args.run_time / 10)
for run_id in range(0,run_num + 1):
    if not os.path.exists('%s/%s_state_%.3d.xml' % (output_dir, prefix, run_id)):
        sim.loadState('%s/%s_state_%.3d.xml' % (output_dir, prefix, run_id - 1))
        dcf_reporter = reporters.NetCDFReporter('%s/%s_%.3d.nc' % (output_dir, prefix, run_id), 25000)
        sim.reporters.append( dcf_reporter )
        sim.reporters.append(
            StateDataReporter("%s/%s_%.3d_log.txt" % (output_dir, prefix, run_id), 25000, step=True,
                              potentialEnergy=True, temperature=True, volume=True)
        )
        # Run dynamics
        print('Running dynamics')
        sim.step(2500000)
        sim.saveState('%s/%s_state_%.3d.xml' % (output_dir, prefix, run_id))
        sim.reporters.clear() # remove the reporters from the current run_num to move to the next one
