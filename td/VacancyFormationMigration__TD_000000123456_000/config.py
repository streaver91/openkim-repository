# Configuration File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Date: 2015/10/30
# Author: Junhao Li <streaver91@gmail.com>

# Parameters for Production
FIRE_LOG = 'fire.log'
FIRE_MAX_STEPS = 500
FIRE_TOL = 1e-3 # absolute
FMIN_FTOL = 1e-10 # relative
CELL_SIZE_MIN = 3
CELL_SIZE_MAX = 5
DYNAMIC_CELL_SIZE = True # Determine supercell size according to lattice structure
MDMIN_TOL = 1e-3 # absolute
MDMIN_MAX_STEPS = 200
NEB_POINTS = 20
UNCERT_STEPS = 20
EPS = 1e-10

# Parameters for Debugging
# FIRE_MAX_STEPS = 200
# FIRE_TOL = 1e-2 # absolute

# Extrapolation Parameters
FITS_CNT = [2, 3] # Number of data points used for each fitting
FITS_ORDERS = [
    [0, 3],
    [0, 3],
] # Number of orders included in each fitting
# Fit Results Used (Corresponding to the above)
FITS_VFE_VALUE = 0 # Vacancy Formation Energy
FITS_VFE_UNCERT = [1]
FITS_VME_VALUE = 0 # Vacancy Migration Energy
FITS_VME_UNCERT = [1]

# Output Configuration
KEY_SOURCE_VALUE = 'source-value'
KEY_SOURCE_UNIT = 'source-unit'
KEY_SOURCE_UNCERT = 'source-std-uncert-value'
UNIT_ENERGY = 'eV'
UNIT_LENGTH = 'angstrom'
UNIT_ANGLE = 'degree'
UNIT_PRESSURE = 'GPa'
UNIT_VOLUME = UNIT_LENGTH + '^3'
SPACE_GROUPS = {
    'fcc': 'Fm-3m',
    'bcc': 'Im-3m', 
    'sc': 'Pm-3m', 
    'diamond': 'Fd-3m', 
    'hcp': 'P63/mmc',
}
WYCKOFF_CODES = {
    'fcc': ['4a'], 
    'bcc': ['2a'], 
    'sc': ['1a'], 
    'diamond': ['8a'], 
    'hcp': ['2d'],
}
WYCKOFF_SITES = {
    'fcc': [[0.0, 0.0, 0.0]],
    'bcc': [[0.0, 0.0, 0.0]],
    'sc': [[0.0, 0.0, 0.0]],
    'diamond': [[0.0, 0.0, 0.0]],
    'hcp': [[2.0 / 3.0, 1.0 / 3.0, 0.25]],
}
VFE_PROP_ID = 'tag:staff@noreply.openkim.org,2015-07-28:property/monovacancy-neutral-relaxed-formation-potential-energy-crystal-npt'
VME_PROP_ID = 'tag:staff@noreply.openkim.org,2015-09-16:property/monovacancy-neutral-migration-energy-crystal-npt'