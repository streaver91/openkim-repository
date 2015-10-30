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
CUR_METHOD = 'lm'
CUR_DATASET = 'large'
INPUT_FILE = 'forces_' + CUR_DATASET + '.txt'
OUTPUT_FILE = 'res_' + CUR_METHOD + '_' + CUR_DATASET + '.txt'
PARAMS_FILE = 'EDIPGenParams.txt'
METHODS = {
    'simplex': optimize.fmin,
    'powell': optimize.fmin_powell,
    'cg': optimize.fmin_cg,
    'bfgs': optimize.fmin_bfgs,
    'lm': optimize.leastsq,
}
FORCE_WEIGHT = 1.0
ENERGY_WEIGHT = 1.0
MAX_ITER = 1e4
INF = 1e20
TOL = 1e-10

output = []
funCalls = -3

def getParamNames():
    paramNames = []
    desFile = os.path.join(MODEL_DIR, MODEL, 'descriptor.kim')
    with open(desFile, 'r') as f:
        for line in f:
            if 'PARAM_FREE' in line:
                paramName = line.split()[0]
                paramNames.append(paramName)
    del paramNames[0]  # drop cutoff
    return paramNames

def readParamSets():
    paramSets = []
    f = open(PARAMS_FILE, 'r')
    for line in f:
        line = line.split(',')
        for i in range(len(line)):
            line[i] = float(line[i])
        paramSets.append(line)
    return paramSets

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

def getResiduals(params, paramNames, atoms, forces, energy, setId, output, returnVector = False):
    outputLength = len(output)
    if (not returnVector) and outputLength > MAX_ITER:
        return INF
    nAtoms = atoms.get_number_of_atoms()
    # npForces = forces.copy()
    
    # Change parameters and get forces
    for i in range(len(paramNames)):
        p = km.KIM_API_get_data_double(atoms.calc.pkim, paramNames[i])
        p[:] = params[i]
    km.KIM_API_model_reinit(atoms.calc.pkim)
    atoms.calc.update(atoms)
    tmpForces = atoms.get_forces()
    
    # Output residual
    # diff = np.apply_along_axis(np.linalg.norm, 1, tmpForces - npForces);
    diff = np.sqrt(((tmpForces - forces)**2).sum(axis = -1))
    diffEnergy = (energy - atoms.get_potential_energy()) / nAtoms
    diff = np.append(diff, diffEnergy)
    residual = (diff**2).sum()
    output.append(residual)
    
    if np.isnan(residual):
        return INF
    if outputLength % 100 == 0:
        print '#', outputLength, '\tSet:', setId, '\tResidual:', residual
    if returnVector:
        return diff
    else:
        return residual
    
def buildAtoms(cellSize, positions):
    nAtoms = len(positions)
    atoms = Atoms(
        ELEM + str(nAtoms),
        positions = positions,
        cell = (cellSize, cellSize, cellSize),
        pbc = (1, 1, 1),
    )
    calc = KIMCalculator(
        MODEL, 
        check_before_update = False, 
        updatenbl = False,
        manualupdate = True
    )
    atoms.set_calculator(calc)
    return atoms
    
def main():
# if True:
    paramNames = getParamNames()
    print paramNames
    paramSets = readParamSets()
    # paramSets = [paramSets[0]]
    print paramSets
    energy, cellSize, positions, forces = readData()
    atoms = buildAtoms(cellSize, positions)
    forces = np.array(forces)
    # 1/0
    
    # Warm up
    getResiduals(paramSets[0], paramNames, atoms, forces, energy, -1, [])
    getResiduals(paramSets[0], paramNames, atoms, forces, energy, -1, [])
    getResiduals(paramSets[0], paramNames, atoms, forces, energy, -1, [])
    
    print '--------------'
    output = []
    
    for i in range(len(paramSets)):
        paramSet = paramSets[i]
        output.append([])
        optimizer = METHODS[CUR_METHOD]
        if CUR_METHOD == 'simplex' or CUR_METHOD == 'powell':
            optimizer(
                getResiduals,
                paramSet,
                args = (paramNames, atoms, forces, energy, i, output[i]), 
                full_output = 1,
                maxfun = MAX_ITER,
                maxiter = MAX_ITER,
                ftol = TOL,
            )
        elif CUR_METHOD == 'lm':
            optimizer(
                getResiduals,
                paramSet,
                args = (paramNames, atoms, forces, energy, i, output[i], True), 
                full_output = 1,
                maxfev = int(MAX_ITER),
                # ftol = TOL,
            )
        else:
            optimizer(
                getResiduals,
                paramSet, 
                args = (paramNames, atoms, forces, energy, i, output[i]), 
                full_output = 1, 
                maxiter = MAX_ITER,
                gtol = TOL,
            )
        initResidual = getResiduals(paramSet, paramNames, atoms, forces, energy, i, [])
        curMin = initResidual
        # for j in range(len(output[i])):
        outputLength = len(output[i])
        for j in range(int(MAX_ITER)):
            if j >= outputLength:
                output[i].append(curMin)
                continue
            if output[i][j] < curMin:
                curMin = output[i][j]
            else:
                output[i][j] = curMin
        print '========='
    output = sorted(output, key = lambda set: set[0])
    outputT = []
    # print output
    for i in range(int(MAX_ITER)):
        outputT.append([])
        outputT[i].append(str(i))
        for j in range(len(paramSets)):
            # print i, j
            # print len(output[j])
            outputT[i].append(str(output[j][i]))
        outputT[i] = '\t'.join(outputT[i])
    outputT = '\r\n'.join(outputT)
    f = open(OUTPUT_FILE, 'w')
    f.write(outputT)
    print 'Data saved to', OUTPUT_FILE
    
main()