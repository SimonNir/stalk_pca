#!/usr/bin/env python3

from pyscf import dft
from numpy import savetxt


### generated system text ###
from pyscf import gto as gto_loc
mol = gto_loc.Mole()
mol.verbose  = 4
mol.atom     = '''
               C    0.37988314   0.00000229   0.00000000
               H    0.00126678   1.02598545  -0.00000000
               H    0.00126630  -0.51298894   0.88852731
               H    0.00126630  -0.51298894  -0.88852731
               F   -2.19606132   0.00000032  -0.00000000
               F    1.81237880  -0.00000054   0.00000000
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