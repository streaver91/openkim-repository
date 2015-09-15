"""
Fitting EPID Potential With Python Minimization

Date: 2015/09/12
Author: Junhao Li
"""

import numpy as np
import scipy as sp
import ctypes
import os
from copy import copy
import ase
from ase import Atoms
from kimcalculator import KIMCalculator
import kimservice as km
import sys

MODEL_DIR = '/home/openkim/openkim-repository/mo'
MODEL = 'EDIP_BOP_Bazant_Kaxiras_Si__MO_958932894036_001'
INPUT_FILE = 'forces_small.txt'
INITIAL_PARAM = {
    'PARAM_FREE_A': 0.99,
    # ...
}

paramNames = []

def getParamNames():
    desFile = os.path.join(MODEL_DIR, MODEL, 'descriptor.kim')
    with open(desFile, 'r') as f:
        for line in f:
            if 'PARAM_FREE' in line:
                paramName = line.split()[0]
                paramNames.append(paramName)

def setZeroInitial():
    for name in paramNames:
        INITIAL_PARAM[name] = 0
    print INITIAL_PARAM

def readData():
    with open(INPUT_FILE, 'r') as f:
        energy = float(f.readline())
        cellSize = float(f.readline())
        positions = []
        forces = []
        for line in f:
            line = line.split()
            if len(line) == 6:
                for i in range(6):
                    line[i] = float(line[i])
                # print line
                positions.append([line[0], line[1], line[2]])
                forces.append([line[3], line[4], line[5]])
    nAtoms = len(positions)
    atoms = Atoms(
        'Si' + str(nAtoms),
        positions = positions,
        cell = (cellSize, cellSize, cellSize),
        pbc = (1, 1, 1),
    )
    ghosts = np.zeros(nAtoms, dtype = 'int8')
    calc = KIMCalculator(MODEL)
    atoms.set_calculator(calc)
    calc.set_ghosts(ghosts)
    return atoms, forces, energy
    
def getResiduals(params):
    calc = KIMCalculator(MODEL)
    atoms.set_calculator(calc)
    for name in params:
        pkim = calc.pkim
        p = km.KIM_API_get_data_double(calc.pkim, name)
        print p
        raw_input("Press Enter to continue...")
        # calc.free_kim()
        # km.KIM_API_set_data_double(pkim, name, [])
        km.KIM_API_model_reinit(calc.pkim)
        # km.KIM_API_set_data_double(pkim, name, [hex(id(1.1))])
        p = km.KIM_API_get_data_double(pkim, name)
        print p

    return 0

getParamNames()
# setZeroInitial()
atoms, forces, energy = readData()
getResiduals(INITIAL_PARAM)
# res = sp.optimize.leastsq(getResiduals, INITIAL_PARAM, full_output = 1)
# print res