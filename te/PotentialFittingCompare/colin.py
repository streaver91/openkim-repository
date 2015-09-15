"""
    potfitToKim.py

    This script include files for building ASE simulation classes with openKIM
    models/tests and PotFit's input files.

    author: Colin Clement
    date: Jan 29 2015
"""

import numpy as np
import scipy as sp
import os
from copy import copy
import read_potfit as rpf
import ase
from kimcalculator import KIMCalculator
import kimservice as km

MODELDIR = "/home/openkim/openkim-repository/mo"

def makeASEAtoms(pfparamfile):
    potfit_params, configlist = rpf.load_potfit_paramfile(pfparamfile)
    AtomsClassList = []
    for config in configlist:
        cell = np.c_[config['X_unit'], config['Y_unit'], config['Z_unit']]
        atomslist = []
        atom_configs = config['atom_configs']
        for species in atom_configs:
            try:
                symbol = potfit_params['#C'][int(species)]
            except KeyError:
                symbol = 'Ar' #Make Argon default
                pass
            for atom_pos in atom_configs[species]:
                atomslist += [ase.Atom(symbol, atom_pos[:3])]
        AtomsClassList += [ase.Atoms(atomslist, cell=cell, pbc=True)]
    return AtomsClassList, potfit_params, configlist

def getModelParams(modelname):
    descriptfile = os.path.join(MODELDIR,modelname,'descriptor.kim')
    model_params = []
    with open(descriptfile, 'r') as infile:
        for line in infile:
            if "PARAM_FREE" in line:
                model_params += [line.split()[0]]
    return model_params

class AtomsLeastSquare(object):
    def __init__(self, potfit_paramfile, modelname, param_map = None, 
                use_energies = False):
        """
        potfit_paramfile : (str) location of potfit param file
        modelname : (str) KIM model name
        (optional)
        param_map : dictionary with keys of the potfit parameter names and
            values of the kim model parameter names from descriptor.kim
            If not supplied, will try to match them up
        use_energies : (bool) If true this object will use cohesive energies in
            its least square fitting
        """
        #------------------
        # Initialize objects
        #------------------
        self.Atomslist, self.potfit_params, self.configlist = \
            makeASEAtoms(potfit_paramfile)
        self.model_params = getModelParams(modelname)
        self.calclist = [KIMCalculator(modelname) for a in self.Atomslist]
        for atoms, calc in zip(self.Atomslist, self.calclist):
            atoms.set_calculator(calc)
        self.use_energies = use_energies
        #---------------------
        # Parameter organization
        #---------------------
        #Set and forget cutoff!
        if reduce(lambda x, y: x or y, map(lambda x: 'cutoff' in x,
            self.model_params)):
            self.set_param({'PARAM_FREE_cutoff': self.potfit_params['cutoff']})
            self.model_params.remove('PARAM_FREE_cutoff')
        if param_map is None:
            self.param_map = {}
            for kimp in self.model_params:
                for pfp in self.potfit_params:
                    if pfp in kimp:
                        self.param_map[kimp] = pfp
        else:
            self.param_map = param_map
        self.initial_params = []
        #Use potfit's initial parameter choice
        for name in self.model_params:
            potfitname = self.param_map[name]
            self.initial_params += [self.potfit_params[potfitname][0]]
        self.current_params = copy(self.initial_params)
        #-----------------------
        #Reference data formating
        #-----------------------
        self.ref_forces = []
        self.ref_energies = []
        for cfg in self.configlist:
            force = [] 
            for spec in cfg['atom_configs']:
                force += cfg['atom_configs'][spec]
            self.ref_energies += [cfg['coh_eng']]
            self.ref_forces += [np.array(force)[:,3:].copy()]

    def _update_calc(self):
        for atoms in self.Atomslist:
            atoms.set_calculator(self.calc)

    def set_param(self, paramdict):
        """
        paramdict : (dict) key - element from self.model_params
                           value - parameter value
        """
        for calc in self.calclist:
            for name in paramdict:
                p = km.KIM_API_get_data_double(calc.pkim,name)    
                p[0] = paramdict[name]
                try:
                    self.current_params[self.model_params.index(name)] = paramdict[name]
                except AttributeError:
                    pass #Not yet defined in __init__
            km.KIM_API_model_reinit(calc.pkim)

    def get_forces(self):
        return [a.get_forces() for a in self.Atomslist]

    def get_energies(self):
        return [a.get_total_energy() for a in self.Atomslist]

    def evaluate_residuals(self, paramlist):
        """
        paramlist : (list) of parameters to fit, ordered by the order they
        appear in self.param_map
        """
        paramdict = {name: val for val, name in zip(paramlist,
            self.model_params)} 
        self.set_param(paramdict)
        res = np.array([])
        forces = self.get_forces()
        energies = self.get_energies()
        for calc_f, ref_f, calc_en, ref_en in zip(forces, self.ref_forces,
                energies, self.ref_energies):
            if self.use_energies:
                res = np.r_[res, (calc_f - ref_f).ravel(),
                        np.array([calc_en-ref_en])]
            else:
                res = np.r_[res, (calc_f - ref_f).ravel()]
        return res

    def evaluate_cost(self, paramlist):
        res = self.evaluate_residuals(paramlist)
        return 0.5*np.sum(res**2)

    def fit(self):
        self.sol = sp.optimize.leastsq(self.evaluate_residuals,
                self.initial_params, full_output=1)
        self.best_fit = self.sol[0]
