#!/usr/bin/env python3

from ctypes import Structure
from numpy import sin, cos, ndarray, array, load

from pyscf import dft
from pyscf import gto
from pyscf import grad
from pyscf.geomopt.geometric_solver import optimize
from pyscf.gto.mole import tofile

from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write
from ase.data import chemical_symbols

from stalk.io.FilesLoader import FilesLoader
from stalk.io.XyzGeometry import XyzGeometry
from stalk import ParameterStructure
from stalk.nexus.NexusGeometry import NexusGeometry
from stalk.nexus.NexusPes import NexusPes
from stalk.nexus.QmcPes import QmcPes
from stalk.params.PesFunction import PesFunction
from stalk.util import EffectiveVariance

from nexus import generate_pyscf, generate_qmcpack, job, obj
from nexus import generate_physical_system, generate_convert4qmc
from structure import Structure

# This requires the following job arguments to be defined in local nxs.py
from nxs import pyscfjob, optjob, dmcjob, pwscfjob, p2qjob

# Pseudos (execute download_pseudos.sh in the working directory)
qmcpseudos = ['F.ccECP.xml', 'C.ccECP.xml', 'H.ccECP.xml']

# Forward mapping: produce parameter values from an array of atomic positions
_network = load('network_parameters.npz')
_components   = _network['components']      # (n_pcs, n_features)
_scaler_mean  = _network['scaler_mean']     # (n_features,)

try: 
    _scaler_scale = _network['scaler_scale']    # (n_features,) 
    if _scaler_scale is None: 
        _scaler_scale = 1
except: 
    print("unable to import scale; defaulting to 1 ")
    _scaler_scale = 1
    
_pca_mean     = _network.get('pca_mean', 0)  # (1, n_features) or 0

# Dimensions
_n_pcs, _n_features = _components.shape
_n_atoms = _n_features // 3


def forward(pos):
    """pos: (n_atoms,3) → params: (n_pcs,)"""
    flat = pos.reshape(1, -1)
    x = (flat - _scaler_mean) / _scaler_scale
    x = x - _pca_mean
    params = x @ _components.T
    return params.ravel()


def backward(params):
    """params: (n_pcs,) → pos: (n_atoms,3)"""
    p = params.reshape(1, -1)
    x = p @ _components
    x = x + _pca_mean
    flat = x * _scaler_scale + _scaler_mean
    return flat.reshape(_n_atoms, 3)


def kernel_pyscf(positions, elem):
    atom = []
    for el, pos in zip(elem, positions):
        atom.append([el, tuple(pos)])
    # end for
    mol = gto.Mole()
    mol.atom = atom
    mol.verbose = 3
    mol.basis = 'ccpvtz'
    mol.unit = 'A'
    mol.ecp = 'ccecp'
    mol.charge = -1
    mol.spin = 0
    mol.symmetry = False
    mol.cart = True
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = 'pbe'
    return mf
# end def


def relax_pyscf(structure: ParameterStructure, outfile='relax.xyz'):
    # Run mean-field calculation
    mf = kernel_pyscf(structure.pos, structure.elem)
    mf.kernel()
    energy = mf.e_tot

    # Geometry optimization
    mol_eq = optimize(mf, maxsteps=100)

    # Convert to ASE Atoms object
    symbols = [mol_eq.atom_symbol(i) for i in range(mol_eq.natm)]
    positions = mol_eq.atom_coords()
    atoms = Atoms(symbols=symbols, positions=positions)

    # Save XYZ with energy in comment
    write(outfile, atoms, comment=f"energy={energy}")
# end def

HA_TO_EV    =  27.2114   # ev/Ha
HA_PER_BOHR_TO_EV_PER_A    =  27.2114 / 0.529177 # (Ha/Bohr) * (eV/Ha) * 1/(A/Bohr) = eV/A
class PySCF(Calculator):    
    implemented_properties = ['energy', 'forces']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
    # end def

    def calculate(
        self,
        atoms=None,
        properties=['energy'],
        system_changes=all_changes
    ):
        Calculator.calculate(self, atoms, properties, system_changes)

        # Get positions and symbols from the Atoms object
        positions = self.atoms.get_positions()
        elem = self.atoms.get_chemical_symbols()

        # Calculate energy and forces here using your custom method
        energy = self.compute_energy(positions, elem)
        forces = self.compute_forces(positions, elem)

        # Set results to be used by ASE
        self.results = {'energy': energy, 'forces': forces}
    # end def

    def compute_energy(self, positions, elem):
        mf = kernel_pyscf(positions, elem)
        return mf.kernel()*HA_TO_EV
    # end def

    def compute_forces(self, positions, elem):
        mf = kernel_pyscf(positions, elem)
        mf.kernel()
        mf_grad = grad.RKS(mf)
        forces = -mf_grad.kernel()
        return forces*HA_PER_BOHR_TO_EV_PER_A
    # end def

# end class


def neb_image(structure: ParameterStructure):
    atoms = Atoms(structure.elem, positions=structure.pos)
    atoms.calc = PySCF()
    # Fix C atom
    atoms.set_constraint(FixAtoms(indices=[0]))
    return atoms
# end def


# Let us define an SCF PES job that is consistent with the earlier relaxation
def scf_pes_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        net_charge=-1,
        F=7,
        C=4,
        H=1,
    )
    scf = generate_pyscf(
        template='../pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path,
        mole=obj(
            verbose=4,
            ecp='ccecp',
            basis='ccpvtz',
            charge=-1,
            symmetry=False,
            cart=True
        ),
    )
    return [scf]
# end def


# Hessian based on the structural mappings
pes_pyscf = NexusPes(
    func=PesFunction(scf_pes_job),
    loader=FilesLoader({'suffix': 'energy.dat'})
)


# 4-5) Stochastic: Line-search
# return a simple 4-item DMC workflow for Nexus:
#   1: run SCF with norm-conserving ECPs to get orbitals
#   2: convert for QMCPACK
#   3: Optimize 2-body Jastrow coefficients
#   4: Run DMC with enough steps/block to meet the target errorbar sigma
#     it is important to first characterize the DMC run into var_eff
def dmc_pes_job(
    structure: Structure,
    path,
    sigma=None,
    samples=10,
    var_eff=None,
    **kwargs
):
    # Estimate the relative number of samples needed
    if isinstance(var_eff, EffectiveVariance):
        dmcsteps = var_eff.get_samples(sigma)
        # Optionally, clip steps to keep within sane range *may result in more ls iter required*
        # from numpy import clip
        # dmcsteps = clip(dmcsteps, a_min=10, a_max=1000)
    else:
        dmcsteps = samples
    # end if

    # Center the structure for QMCPACK
    system = generate_physical_system(
        structure=structure,
        net_charge=-1, 
        F=7,
        C=4,
        H=1,
    )
    scf = generate_pyscf(
        template='../pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path + 'scf',
        mole=obj(
            spin=4,
            verbose=4,
            ecp='ccecp',
            basis='ccpvtz',  # Use larger basis to promote QMC performance
            charge=-1,
            symmetry=False,
        ),
        save_qmc=True,
    )
    c4q = generate_convert4qmc(
        identifier='c4q',
        path=path + 'scf',
        job=job(cores=1),
        dependencies=(scf, 'orbitals'),
    )
    opt = generate_qmcpack(
        system=system,
        path=path + 'opt',
        job=job(**optjob),
        dependencies=[(c4q, 'orbitals')],
        cycles=8,
        identifier='opt',
        qmc='opt',
        input_type='basic',
        pseudos=qmcpseudos,
        J2=True,
        J1_size=6,
        J1_rcut=6.0,
        J2_size=8,
        J2_rcut=8.0,
        minmethod='oneshift',
        blocks=200,
        substeps=2,
        samples=100000,
        minwalkers=0.1,
    )
    dmc = generate_qmcpack(
        system=system,
        path=path + 'dmc',
        job=job(**dmcjob),
        dependencies=[(c4q, 'orbitals'), (opt, 'jastrow')],
        steps=dmcsteps,
        identifier='dmc',
        qmc='dmc',
        input_type='basic',
        pseudos=qmcpseudos,
        jastrows=[],
        walkers_per_rank=128,
        blocks=200,
        timestep=0.01,
        ntimesteps=1,
    )
    # Store the relative samples for printout
    dmc.samples = dmcsteps
    return [scf, c4q, opt, dmc]
# end def


# Configure a job generator and loader for the DMC PES
# -the suffix points to the correct nexus analyzer
# -the qmc_idx points to the correct QMC series (0: VMC; 1: DMC)
pes_dmc = NexusPes(
    dmc_pes_job,
    loader=QmcPes({'suffix': '/dmc/dmc.in.xml', 'qmc_idx': 1})
)
