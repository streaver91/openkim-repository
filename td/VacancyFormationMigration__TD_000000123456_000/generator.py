"""
Generate test_generator.json for tests generation

Date: 2015/10/30
Author: Junhao Li <streaver91@gmail.com>
"""

import json
from ase.data import chemical_symbols
from ase.data import reference_states
from random import randint

lattices = ['fcc', 'hcp']
chemical_symbols = ['Ac', 'Ag', 'Al', 'Ar', 'Au', 'Ba', 'Be', 'Bi', 'Ca', 'Cd', 'Ce', 'Co', 'Cr', 'Cs', 'Cu', 'Dy', 'Er', 'Fe', 'Ga', 'Ge', 'He', 'Hf', 'Ho', 'In', 'Ir', 'K', 'Kr', 'La', 'Li', 'Mg', 'Mn', 'Mo', 'Na', 'Nb', 'Ne', 'Ni', 'Os', 'P', 'Pa', 'Pb', 'Pd', 'Pr', 'Pt', 'Rb', 'Re', 'Rh', 'Ru', 'Sc', 'Si', 'Sn', 'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'V', 'W', 'Xe', 'Y', 'Zn', 'Zr']

# For Debug
# chemical_symbols = ['Al']

with open('test_generator.json', 'w') as f:
    for pk, elem in enumerate(chemical_symbols):
        for lattice in lattices:
            kimnum = '{:012d}'.format(randint(0, 10**12 - 1))
            
            # hcp lattice will use a different lattice constant driver
            if lattice == 'hcp':
                latticeProperty = 'structure-hexagonal-crystal-npt'
            else:
                latticeProperty = 'structure-cubic-crystal-npt'
            
            f.write(json.dumps({
                'elem': elem,
                'lattice': lattice, 
                'kimnum': kimnum,
                'latticeProperty': latticeProperty,
            }) + '\n')
