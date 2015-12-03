#!/usr/bin/env python
"""
Force Generation
Date: 2015/11/16
Author: Junhao Li
"""
from ase.lattice import bulk
from kimcalculator import KIMCalculator
from string import Template
from ase.optimize import FIRE
from ase.lattice.cubic import Diamond
from kimservice import KIM_API_get_data_double
from ase import Atoms
import numpy as np
import os, sys
import random

# Configurations
elem = 'Si'
lattice = 'diamond'
model = 'EDIP_BOP_Bazant_Kaxiras_Si__MO_958932894036_001'
latticeConstant = '5.43049777059'
randomSize = 0.1
size = 5

def main():
    # Make them global for ipython debugging
    global atoms, superAtoms
    global length, positions, nAtoms
    global latticeConstant
    
    # Create single conventional unit cell
    latticeConstant = float(latticeConstant)
    calc = KIMCalculator(model)
    atoms = bulk(elem, lattice, a = latticeConstant, cubic = True)

    # Create bulk supercell
    superAtoms = atoms.copy()
    superAtoms *= (size, size, size)

    superCell = superAtoms.get_cell()
    length = superCell[0][0]
    nAtoms = superAtoms.get_number_of_atoms()

    # Randomly purturb positions
    positions = superAtoms.get_positions().copy()
    positions += np.random.uniform(0, randomSize, (nAtoms, 3))
    superAtoms.set_positions(positions)

    # Copy 27 times to simulate periodic boundary condition
    superAtoms *= (3, 3, 3)

    # move real atoms to the beginning of the list
    positionsOriginal = superAtoms.get_positions().copy()
    tmp = np.hstack([positionsOriginal >= length, positionsOriginal < length * 2])
    isReal = np.all(tmp, axis = 1)
    realPositions = positionsOriginal[isReal]
    ghostPositions = positionsOriginal[np.logical_not(isReal)]
    positions = np.vstack((realPositions, ghostPositions))
    superAtoms.set_positions(positions)

    print superAtoms.cell
    print len(superAtoms)
    # Setup superAtoms object
    superAtoms.set_calculator(calc)
    superAtoms.set_pbc([0, 0, 0])
    ghosts = np.ones(nAtoms * 27, dtype = 'int8')
    ghosts[0: nAtoms] = 0

    print len(superAtoms) - ghosts.sum()
    print 'superAtoms:', superAtoms.get_positions().shape
    print 'ghosts:', ghosts
    np.set_printoptions(suppress = True)
    for i in range(27):
        tmp = realPositions - positionsOriginal[nAtoms * i: nAtoms * (i + 1)]
        print '#', i, 'avg:', tmp.mean(axis = 0)
        # print 'std:', tmp.std(axis = 0)
    # tmp = realPositions - ghostPositions[nAtoms: nAtoms * 2]
    # print tmp
    # print (np.abs(tmp - length) > 1.0e-10).sum()
    calc.set_ghosts(ghosts)
    energy = superAtoms.get_potential_energy()
    forces = superAtoms.get_forces()
    
    # Generate output
    output = []
    output.append(str(energy))
    output.append(str(length))
    for i in range(nAtoms * 0, nAtoms * 1):
        row = []
        for j in range(3):
            row.append(str(positions[i][j] - length))
        for j in range(3):
            row.append(str(forces[i][j]))
        output.append('\t'.join(row))
    output = '\r\n'.join(output)
    print 'cohesive-potential-energy:', energy

    with open(os.path.abspath('forces.txt'), 'w') as f:
        f.write(output)

main()
