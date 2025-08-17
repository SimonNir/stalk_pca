#!/usr/bin/env python3

from pyscf import dft
from numpy import savetxt


### generated system text ###
from pyscf import gto as gto_loc
mol = gto_loc.Mole()
mol.verbose  = 4
mol.atom     = '''
               C    0.37907653   0.00000076  -0.00000000
               H    0.00125712   1.02727264  -0.00000000
               H    0.00125664  -0.51363518   0.88964353
               H    0.00125664  -0.51363518  -0.88964353
               F   -2.19555408   0.00000373   0.00000000
               F    1.81270714   0.00000288  -0.00000000
               '''
mol.basis    = 'ccpvtz'
mol.unit     = 'A'
mol.ecp      = 'ccecp'
mol.charge   = -1
mol.spin     = 0
mol.symmetry = False
mol.cart     = True
mol.build()
### end generated system text ###



### generated calculation text ###
mf = dft.RKS(mol)
mf.xc = 'pbe'
e_scf = mf.kernel()
### end generated calculation text ###

savetxt('energy.dat', [[e_scf, 0.0]])