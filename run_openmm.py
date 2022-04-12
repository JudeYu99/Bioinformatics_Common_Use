#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pdbfixer import PDBFixer
from openmm.app.pdbfile import PDBFile
import os

def fix_pdb(pdb_id):
    path = os.getcwd()
    if len(pdb_id) != 4:
        print("Creating PDBFixer...")
        fixer = PDBFixer(pdb_id)
        print("Finding missing residues...")
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                print("ok")
                del fixer.missingResidues[key]

        print("Finding nonstandard residues...")
        fixer.findNonstandardResidues()
        print("Replacing nonstandard residues...")
        fixer.replaceNonstandardResidues()
        print("Removing heterogens...")
        fixer.removeHeterogens(keepWater=True)

        print("Finding missing atoms...")
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(7)
        print("Writing PDB file...")

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)),
                 "w"),
            keepIds=True)
        return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)

fix_pdb('protein.pdb')


from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout

# Importing pdb files, also supports PDBx/mmCIF format
pdb = app.PDBFile('protein_fixed_pH_7.pdb')


# Determine the force field used, Amber14: amber14-all.xml; TIP3P-FB water model: amber14/tip3pfb.xml
# Use XML files to define standard force fields
# Available force fields to view from: http://docs.openmm.org/latest/userguide/application.html#force-fields
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')


# Combine the force field with the molecular topology file to create a complete mathematical description of this simulated system
# Some options：
# Use the particle lattice Ewald for long-range electrostatic interactions (nonbondedMethod = PME)
# Set a cutoff of 1 nm for direct spatial interactions (nonbondedCutoff = 1 * nm)
# Constrain the length of all bonds involving hydrogen atoms (constraint = HBonds).
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,nonbondedCutoff=1*nanometer, constraints=HBonds)


# Langevin Dynamics, simulated temperature (300 K), friction coefficient (1 ps-1) and step size (0.004 ps)
# Many other integration methods can also be used. For example, if the system is to be simulated at constant energy rather than constant temperature, the VerletIntegrator can be used.
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

# Combine topology, system, and integrator to create a Simulation object, and assign it to a variable named Simulation.
# Simulation objects manage all the processes involved in running a simulation
simulation = Simulation(pdb.topology, system, integrator)

# Specify the initial atomic positions to be used for the simulation.
simulation.context.setPositions(pdb.positions)

# Perform local energy minimization, recommended before the start of the simulation to reduce stress
simulation.minimizeEnergy()

# Create a reporter to output the results during the simulation, PDBReporter writes the structure to the PDB file.
# Specify the output file as output.pdb and should write a structure every 1000 steps.
simulation.reporters.append(PDBReporter('output.pdb', 1000))

# Get regular status reports while the simulation is running and you can monitor its progress
# Add another reporting program that prints some basic information at every 1000 time steps: the current step index, the potential energy of the system and the temperature.
# Specify stdout (without quotes) as the output file, which means to display the results in the terminal.
# You can also give the file name (caused by quotes) and write the information to the file.
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,potentialEnergy=True, temperature=True))

# Run simulation
simulation.step(10000)
  
