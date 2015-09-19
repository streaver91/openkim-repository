from kimcalculator import KIMCalculator
import kimservice as ks
from ase.lattice import bulk
import numpy as np

model = 'Three_Body_Stillinger_Weber_Si__MO_405512056662_001'
param = 'PARAM_FREE_epsilon'

atoms = bulk('Si', 'diamond', cubic=True)
atoms.set_calculator(KIMCalculator(model, check_before_update=True))
atoms.set_positions(atoms.get_positions() + np.random.rand(*atoms.get_positions().shape))

forces0 = atoms.get_forces()

t = ks.KIM_API_get_data_double(atoms.calc.pkim, param)
# t[:] += 1e-3

# atoms.set_positions(atoms.get_positions() + 1e-14)

forces1 = atoms.get_forces()

print forces1-forces0
