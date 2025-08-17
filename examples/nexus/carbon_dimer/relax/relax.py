#!/usr/bin/env python3

from pyscf import dft


### generated system text ###
from pyscf import gto as gto_loc
mol = gto_loc.Mole()
mol.verbose  = 4
mol.atom     = '''
               C    0.00000000   0.00000000  -0.77000000
               C    0.00000000   0.00000000   0.77000000
               '''
mol.basis    = 'ccpvdz'
mol.unit     = 'A'
mol.ecp      = 'ccecp'
mol.charge   = 0
mol.spin     = 4
mol.symmetry = False
mol.build()
### end generated system text ###



### generated calculation text ###
mf = dft.RKS(mol)
mf.xc = 'pbe'
e_scf = mf.kernel()
### end generated calculation text ###

from pyscf.geomopt.geometric_solver import optimize
mol_eq = optimize(mf, maxsteps=100)
from pyscf.gto.mole import tofile
# Write to external file
tofile(mol_eq, 'relax.xyz', format='xyz')
# Write to output file 
print(mol_eq.atom_coords())
