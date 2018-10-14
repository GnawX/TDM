from ase import atoms, io
from ase.units import Ang, Bohr, _amu, _me
import numpy as np

es = io.read('end.xyz')
gs = io.read('start.xyz')

re = es.get_positions()  # in unit of Ang
rg = gs.get_positions()

dr = re - rg
m = es.get_masses()   # in unit of amu

dq = np.sqrt(np.dot(m, np.square(dr)).sum())  # in unit of sqrt(amu)*Ang

print dq
