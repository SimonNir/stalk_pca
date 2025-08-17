#!/usr/bin/env python3

from pyscf import dft
from numpy import savetxt


### generated system text ###
from pyscf import gto as gto_loc
mol = gto_loc.Mole()
mol.verbose  = 4
mol.atom     = '''
               C    0.37950398   0.00000154   0.00000000
               H    0.00084925   1.02628620  -0.00000000
               H    0.00084838  -0.51314061   0.88878850
               H    0.00084838  -0.51314061  -0.88878850
               F   -2.19617432   0.00000188  -0.00000000
               F    1.81412434   0.00000125   0.00000000
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