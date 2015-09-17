"""
Fitting EPID Potential With Python Minimization

Date: 2015/09/12
Author: Junhao Li
"""
import ase
from ase import Atoms
from ase.lattice import bulk

from kimcalculator import KIMCalculator
import kimservice as km

import numpy as np
from scipy import optimize
import ctypes
import os
from copy import copy
import sys
import math
import random

LATTICE = 'diamond'
ELEM = 'Si'
MODEL_DIR = '/home/openkim/openkim-repository/mo'
MODEL = 'EDIP_BOP_Bazant_Kaxiras_Si__MO_958932894036_001'
INPUT_FILE = 'forces_small.txt'
PARAMS_FILE = 'params.txt'
FORCE_WEIGHT = 0.9
ENERGY_WEIGHT = 0.1

def getParamNames():
    paramNames = []
    desFile = os.path.join(MODEL_DIR, MODEL, 'descriptor.kim')
    with open(desFile, 'r') as f:
        for line in f:
            if 'PARAM_FREE' in line:
                paramName = line.split()[0]
                paramNames.append(paramName)
    return paramNames

def readParams():
    # For debug
    paramString = '3.1213820	7.9821730	1.5075463	1.2085196	0.5774108	1.4533108	1.1247945	3.1213820	2.5609104	0.6966326	312.1341346	0.2523244	0.0070975	3.1083847	-0.165799	32.557	0.286198	0.66'
    paramValues = paramString.split('\t')
    paramDict = {}
    for i in range(len(paramValues)):
        paramValues[i] = float(paramValues[i]) * random.uniform(0.9, 1.1)
    return paramValues

def readData():
    # Read Data in Format:
    # energy (eV)
    # cellSize (angstrom)
    # x y z Fx Fy Fz
    # ...
    positions = []
    forces = []
    energy = 0
    cellSize = 0
    with open(INPUT_FILE, 'r') as f:
        energy = float(f.readline())
        cellSize = float(f.readline())
        for line in f:
            line = line.split()
            if len(line) == 6:
                for i in range(6):
                    line[i] = float(line[i])
                positions.append([line[0], line[1], line[2]])
                forces.append([line[3], line[4], line[5]])
    nAtoms = len(positions)
    atoms = Atoms(
        ELEM + str(nAtoms),
        positions = positions,
        cell = (cellSize, cellSize, cellSize),
        pbc = (1, 1, 1),
    )
    ghosts = np.zeros(nAtoms, dtype = 'int8')
    calc = KIMCalculator(MODEL)
    atoms.set_calculator(calc)
    calc.set_ghosts(ghosts)
    return atoms, forces, energy
    
def getResiduals(params, paramNames, atoms, forces, energy):
    nAtoms = atoms.get_number_of_atoms()
    pkim = atoms.calc.pkim
    for i in range(len(paramNames)):
        p = km.KIM_API_get_data_double(pkim, paramNames[i])
        p[0] = params[i]
        # print p
    km.KIM_API_model_reinit(pkim)
    tmpForces = atoms.get_forces()
    sumForce = sum(np.linalg.norm(tmpForces[i] - forces[i]) for i in range(nAtoms))
    print 'Force:', sumForce
    return sumForce * FORCE_WEIGHT + energy * ENERGY_WEIGHT

def main():
    paramNames = getParamNames()
    print paramNames
    initialParams = readParams()
    print initialParams
    atoms, forces, energy = readData()
    getResiduals(initialParams, paramNames, atoms, forces, energy)
    # res = optimize.fmin(getResiduals, initialParams, args = (paramNames, atoms, forces, energy), full_output = 1)
    # print res
    
main()