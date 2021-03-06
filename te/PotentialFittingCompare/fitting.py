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
# CUR_METHOD = 'lm'
CUR_METHOD = raw_input('Method: ')
CUR_PARAMSET = int(raw_input('Param Set (-1 for all): '))
# CUR_PARAMSET = 3
CUR_DATASET = 'large'
# INPUT_FILE = 'forces_' + CUR_DATASET + '.txt'
INPUT_FILE = 'forces.txt'
OUTPUT_FILE = 'res_' + CUR_METHOD + '_' + CUR_DATASET + '_GeStart.txt'
PARAMS_FILE = 'EDIPGenParams_Ge.txt'
MAX_DIFF = 0.2
METHODS = {
    'simplex': optimize.fmin,
    'powell': optimize.fmin_powell,
    'cg': optimize.fmin_cg,
    'bfgs': optimize.fmin_bfgs,
    'lm': optimize.leastsq,
    'lm2': optimize.leastsq,
    'ga': optimize.differential_evolution,
    'bh': optimize.basinhopping
}
FORCE_WEIGHT = 1.0
ENERGY_WEIGHT = 1.0
MAX_ITER = 2000
INF = 1e20
TOL = 1e-20
EPS = 1e-10
DX_PERCENT = 0.00001
EXCLUDED = [0, 7, 8, 14, 15, 16, 17]

output = []
funCalls = -3
curResidualTmp = 0

def getParamNames():
    paramNames = []
    desFile = os.path.join(MODEL_DIR, MODEL, 'descriptor.kim')
    with open(desFile, 'r') as f:
        for line in f:
            if 'PARAM_FREE' in line:
                paramName = line.split()[0]
                paramNames.append(paramName)
    
    lenExcluded = len(EXCLUDED)
    for i in range(lenExcluded):
        del paramNames[EXCLUDED[lenExcluded - i - 1]]
    return paramNames

def readParamSets():
    paramSets = []
    f = open(PARAMS_FILE, 'r')
    for line in f:
        line = line.split(',')
        for i in range(len(line)):
            line[i] = float(line[i])
        lenExcluded = len(EXCLUDED)
        for i in range(lenExcluded):
            del line[EXCLUDED[lenExcluded - i - 1]]
        print np.array(line)
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

def getResiduals(
        params,
        paramNames,
        atoms,
        forces,
        energy,
        setId,
        output,
        paramVals,
        initialSet = None,
        returnVector = False,
        driftPenalty = 0.0,
    ):
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
    # diff = np.sqrt(((tmpForces - forces)**2).sum(axis = -1))
    diff = np.reshape(tmpForces - forces, -1)
    # print diff
    # sys.exit(0)
    diffEnergy = (energy - atoms.get_potential_energy()) / nAtoms
    diff = np.append(diff, diffEnergy)

    residual = (diff**2).sum()
    output.append(residual)
    realResidual = residual
    curResidualTmp = residual
    if driftPenalty > EPS:

        # print '.....'
        # print np.array(params)
        # print np.array(initialSet)
        # print (np.array(params) - np.array(params))
        # print (np.array(initialSet) + np.array(params))
        # print np.array(initialSet).shape
        # print np.array(params).shape
        # print (np.array(params) - np.array(initialSet)) / (EPS + np.array(initialSet))
        diff = np.append(diff, (params - initialSet) / (EPS + initialSet) * driftPenalty)
        # diff = np.append(diff, (params - initialSet) / (EPS + initialSet) * np.sqrt(curResidualTmp))
        #print params - initialSet
        #print EPS + initialSet
        #print (params - initialSet) / (EPS + initialSet)
        #if len(output) > 4:
        #    sys.exit(0)
        residual = (diff**2).sum()

    # Store Intermediate parameter values
    paramVals.append(copy(params))
    # np.append(paramVals, params)

    if np.isnan(residual):
        return INF
    if outputLength % 100 == 0 and setId != -1:
        print '#', outputLength, '\tSet:', setId, '\tResidual:', realResidual
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

    # sys.exit(0)
    paramSets = readParamSets()

    if CUR_PARAMSET >= 0:
        paramSets = [paramSets[CUR_PARAMSET]]

    print paramSets
    energy, cellSize, positions, forces = readData()
    atoms = buildAtoms(cellSize, positions)
    forces = np.array(forces)

    paramValOrigin = []
    for i in range(len(paramNames)):
        p = km.KIM_API_get_data_double(atoms.calc.pkim, paramNames[i])
        paramValOrigin.append(p[0])
    print paramValOrigin
    # sys.exit(0)
    # 1/0

    # Warm up
    getResiduals(paramSets[0], paramNames, atoms, forces, energy, -1, paramValOrigin, [], [])
    getResiduals(paramSets[0], paramNames, atoms, forces, energy, -1, paramValOrigin, [], [])
    getResiduals(paramSets[0], paramNames, atoms, forces, energy, -1, paramValOrigin, [], [])

    print '--------------'
    output = []
    paramVals = []
    # paramVals = np.array([])
    for i in range(len(paramSets)):
        paramSet = np.array(paramSets[i])
        output.append([])
        paramVals.append([])
        # np.append(paramVals, np.array([]))
        optimizer = METHODS[CUR_METHOD]
        if CUR_METHOD == 'simplex' or CUR_METHOD == 'powell':
            res = optimizer(
                getResiduals,
                paramSet,
                args = (paramNames, atoms, forces, energy, i, output[i], paramVals[i]),
                full_output = 1,
                maxfun = MAX_ITER,
                maxiter = MAX_ITER,
                ftol = TOL,
            )
        elif CUR_METHOD == 'lm':
            res = optimizer(
                getResiduals,
                paramSet,
                args = (paramNames, atoms, forces, energy, i, output[i], paramVals[i], paramValOrigin, True),
                full_output = 1,
                maxfev = int(MAX_ITER),
                xtol = 0.000001
                # ftol = TOL,
            )
        elif CUR_METHOD == 'lm2':
            penaltyCoef = 1e6
            # dampingCoef = 0.5
            # terminalCoef = 1
            curResidual = getResiduals(paramSet, paramNames, atoms, forces, energy, -1, [], [])
            curResidualTmp = curResidual
            newParamSet = paramSet.copy()
            penaltyCoef = np.sqrt(curResidual)
            initialSet = paramSet.copy()
            while len(output[i]) < int(MAX_ITER):
                res = optimizer(
                    getResiduals,
                    newParamSet,
                    args = (
                        paramNames,
                        atoms,
                        forces,
                        energy,
                        i,
                        output[i],
                        paramVals[i],
                        initialSet,
                        True,
                        penaltyCoef,
                    ),
                    full_output = 1,
                    maxfev = 200,
                    xtol = 1e-10
                )
                # print len(output[i])
                newParamSet = paramVals[i][-1]
                print 'Current Drift Penalty: ', penaltyCoef
                curResidual = getResiduals(newParamSet, paramNames, atoms, forces, energy, -1, [], [])
                penaltyCoef = np.sqrt(curResidual)
            
            paramSet = initialSet
                # print 'Current Residual: ', getResiduals(paramSet, paramNames, atoms, forces, energy, -1, [], [])
                # paramSet = paramVals[i][-1]
        elif CUR_METHOD == 'lm3':
            # local opt on continuously depent variables
            # global opt on other variables
            print 'Hello lm3'
        elif CUR_METHOD == 'cg' or CUR_METHOD == 'bfgs':
            res = optimizer(
                getResiduals,
                paramSet,
                args = (paramNames, atoms, forces, energy, i, output[i], paramVals[i]),
                full_output = 1,
                maxiter = MAX_ITER,
                gtol = TOL,
            )
        elif CUR_METHOD == 'ga':
            bounds = []
            for j in range(len(paramNames)):
                bounds.append((paramSet[j] * (1 - MAX_DIFF), paramSet[j] * (1 + MAX_DIFF)))
            res = optimizer(
                getResiduals,
                bounds,
                args = (paramNames, atoms, forces, energy, i, output[i], paramVals[i]),
                tol = TOL,
                strategy = 'best1exp',
                disp = True,
            )
            print res
        elif CUR_METHOD == 'bh':
            res = optimizer(
                getResiduals,
                args = (paramNames, atoms, forces, energy, i, output[i], paramVals[i]),
                tol = TOL,
                strategy = 'best1exp',
                disp = True,
            )
        initResidual = getResiduals(paramSet, paramNames, atoms, forces, energy, -1, [], [])
        curMin = initResidual
        # for j in range(len(output[i])):
        outputLength = len(output[i])
        print "Length: ", str(len(output[i]))
        for j in range(int(MAX_ITER)):
            if j >= outputLength:
                output[i].append(curMin)
                continue
            if output[i][j] < curMin:
                curMin = output[i][j]
            else:
                output[i][j] = curMin
        # print output[i][0], output[i][1], output[i][2]
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
    
    if CUR_PARAMSET == -1:
        f = open(OUTPUT_FILE, 'w')
        f.write(outputT)
        print 'Residual Data Saved To', OUTPUT_FILE
    
    f = open('paramVals.txt', 'w')
    paramVal = np.array(paramVals[0])
    # print paramVal
    output2 = []
    # print paramVal.shape

    for i in range(paramVal.shape[0]):
        tmp = [str(i)]
        for j in range(paramVal.shape[1]):
            tmp.append(str(paramVal[i][j] / paramValOrigin[j]))
        output2.append('\t'.join(tmp))
    # print output2
    f.write('\r\n'.join(output2))

    # Obtain Jacobian
    # print np.array(paramVals).shape
    J = []
    HDiag = []
    for i in range(len(paramNames)):
        finalParams = copy(paramVals[0][-1])
        residualOrigin = getResiduals(finalParams, paramNames, atoms, forces, energy, -1, [], [])
        finalParams[i] *= 1 + DX_PERCENT
        residualPlus = getResiduals(finalParams, paramNames, atoms, forces, energy, -1, [], [])
        print residualPlus
        finalParams = copy(paramVals[0][-1])
        finalParams[i] *= 1 - DX_PERCENT
        residualMinus = getResiduals(finalParams, paramNames, atoms, forces, energy, -1, [], [])
        divI = (residualPlus - residualMinus) / (finalParams[i] * DX_PERCENT * 2)
        div2I = (residualPlus + residualMinus - 2 * residualOrigin) / (finalParams[i] * DX_PERCENT)**2
        print residualPlus, residualMinus, finalParams[i], finalParams[i] * DX_PERCENT * 2, divI

        J.append(divI)
        HDiag.append(div2I)

    print 'Jacobian: '
    for i in range(len(paramNames)):
        print paramNames[i], ':', finalParams[i], J[i], HDiag[i]



main()
