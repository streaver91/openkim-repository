#!/usr/bin/env python
"""
Force Generation

Date: 2015/11/03
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

elem = 'Si'
lattice = 'diamond'
model = 'EDIP_BOP_Bazant_Kaxiras_Si__MO_958932894036_001'
latticeConstant = '5.43049777059'
randomSize = 0.2

latticeConstant = float(latticeConstant)
calc = KIMCalculator(model)
atoms = bulk(elem, lattice, a = latticeConstant, cubic = True)
atoms.set_calculator(calc)

size = 5;
superAtoms = atoms.copy()
superAtoms *= (size, size, size)

positions = superAtoms.get_positions()
superCell = superAtoms.get_cell()
length = superCell[0][0]
nAtoms = len(positions)

for i in range(nAtoms):
    cnt = 0
    for j in range(3):
        coord = positions[i][j]
        coord = coord + (random.random() - 0.5) * randomSize * 2
        positions[i][j] = coord

superAtoms.set_positions(positions)
# superAtoms *= (3, 3, 3)
superAtoms.set_calculator(calc)
superAtoms.set_pbc([1, 1, 1])
# ghosts = np.ones(nAtoms * 27, dtype = 'int8')
# for i in range(nAtoms * 13, nAtoms * 14):
    # ghosts[i] = 0
# calc.set_ghosts(ghosts)
energy = superAtoms.get_potential_energy()
forces = superAtoms.get_forces()
positions = superAtoms.get_positions()

output = []
output.append(str(energy))
output.append(str(length))
for i in range(nAtoms):
    row = []
    for j in range(3):
        row.append(str(positions[i][j]))
    for j in range(3):
        row.append(str(forces[i][j]))
    output.append('\t'.join(row))

output = '\r\n'.join(output)
print 'cohesive-potential-energy:', energy

with open(os.path.abspath('forces.txt'), 'w') as f:
	f.write(output)

