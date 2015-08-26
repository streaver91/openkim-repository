import json
from ase.data import chemical_symbols
from ase.data import reference_states
from random import randint

with open('test_generator.json', 'w') as f:
    for element in chemical_symbols:
        kimnum = '{:012d}'.format(randint(0,10**12-1))
        f.write(json.dumps({
            'symbol': element,
            'lattice': 'hcp',
            'kimnum': kimnum,
            'version': '001',
        }) + '\n')
