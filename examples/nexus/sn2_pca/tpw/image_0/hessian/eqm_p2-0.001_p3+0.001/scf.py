#!/usr/bin/env python3

from pyscf import dft
from numpy import savetxt


### generated system text ###
from pyscf import gto as gto_loc
mol = gto_loc.Mole()
mol.verbose  = 4
mol.atom     = '''
               C    0.37864908   0.00000076  -0.00000000
               H    0.00166499   1.02704407  -0.00000000
               H    0.00166491  -0.51352088   0.88944558
               H    0.00166491  -0.51352088  -0.88944558
               F   -2.19660189   0.00000382   0.00000000
               F    1.81295800   0.00000274  -0.00000000
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