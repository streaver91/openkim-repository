"""
Fitting EPID Potential With Python Minimization

Date: 2015/09/12
Author: Junhao Li
"""
# ASE
import ase
from ase import Atoms
from ase.lattice import bulk

# KIM
from kimcalculator import KIMCalculator
import kimservice as km

# Others
import numpy as np
from scipy import optimize
import ctypes
import os
from copy import copy
import sys
import math
import random

# Configs
LATTICE = 'diamond'
ELEM = 'Si'
MODEL_DIR = '/home/openkim/openkim-repository/mo'
MODEL = 'EDIP_BOP_Bazant_Kaxiras_Si__MO_958932894036_001'
MODEL = 'Three_Body_Stillinger_Weber_Balamane_Si__MO_113686039439_001'
# MODEL = 'LennardJones612_UniversalShifted__MO_959249795837_001'
INPUT_FILE = 'forces_small.txt'
PARAMS_FILE = 'params.txt'
FORCE_WEIGHT = 0.9
ENERGY_WEIGHT = 0.1
MAX_ITER = 1e-5
FTOL = 1e-3

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
    paramString = '3.1213820	7.049556277	0.6022245584	4	0	1.80	21.0	1.2	2.0951	2.315	-0.3333333'
    paramValues = paramString.split('\t')
    paramDict = {}
    for i in range(len(paramValues)):
        # paramValues[i] = float(paramValues[i])
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
    return energy, cellSize, positions, forces
    
def getResiduals(params, paramNames, atoms, forces, energy):
    nAtoms = atoms.get_number_of_atoms()
    dummyforces = atoms.get_forces()
    for i in range(len(paramNames)):
        p = km.KIM_API_get_data_double(atoms.calc.pkim, paramNames[i])
        p[:] = params[i]
        # print paramNames[i], p
    km.KIM_API_model_reinit(atoms.calc.pkim)
    tmpForces = atoms.get_forces()
    # sumForce = sum(np.linalg.norm(tmpForces[i]) for i in range(nAtoms))
    sumForce = sum(np.linalg.norm(tmpForces[i] - forces[i]) for i in range(nAtoms))
    print 'Force Total:', sumForce
    return sumForce * FORCE_WEIGHT + energy * ENERGY_WEIGHT

def buildAtoms(cellSize, positions):
    nAtoms = len(positions)
    atoms = Atoms(
        ELEM + str(nAtoms),
        positions = positions,
        cell = (cellSize, cellSize, cellSize),
        pbc = (1, 1, 1),
    )
    calc = KIMCalculator(MODEL, check_before_update = True)
    ghosts = np.zeros(nAtoms, dtype = 'int8')
    atoms.set_calculator(calc)
    calc.set_ghosts(ghosts)
    return atoms
    
def main():
    paramNames = getParamNames()
    print paramNames
    initialParams = readParams()
    print initialParams
    energy, cellSize, positions, forces = readData()
    atoms = buildAtoms(cellSize, positions)
    getResiduals(initialParams, paramNames, atoms, forces, energy)
    res = optimize.fmin(
        getResiduals,
        initialParams, 
        args = (paramNames, atoms, forces, energy), 
        full_output = 1, 
        maxfun = MAX_ITER, 
        maxiter = MAX_ITER, 
        ftol = FTOL
    )
    print res
    
main()