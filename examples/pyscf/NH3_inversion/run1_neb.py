#!/usr/bin/env python3

from os import makedirs
import os

from ase.mep import NEB
from ase.optimize import BFGS
from ase import io

from params import neb_image, HA_PER_BOHR_TO_EV_PER_A, HA_TO_EV
from run0_relax_a import structure_relax as structure_a
from run0_relax_b import structure_relax as structure_b
from stalk.params.util import interpolate_params


# Generate base directory for NEB
basedir = 'neb'
makedirs(basedir, exist_ok=True)

# number of intermediate images
n_images = 3

traj_init = interpolate_params(structure_a, structure_b, n_images)
# Or, to be safer (against discontinuities, etc) but possibly less accurate, interpolate in cartesian

images = [neb_image(structure) for structure in traj_init]

# Try to load from disk
try:
    traj_neb = []
    for i in range(n_images + 2):
        xyz_file = f'{basedir}/image{i}.xyz'
        image = io.read(xyz_file)
        from params import HA_TO_EV
        traj_neb.append(structure_a.copy(pos=image.positions, surrogate_energy=image.get_total_energy()/HA_TO_EV))
    # end for
except FileNotFoundError:
    neb = NEB(images, climb=True, remove_rotation_and_translation=True)
    opt = BFGS(neb, trajectory=f'{basedir}/cineb.traj', logfile=f'{basedir}/cineb.log')
    opt.run(fmax=0.01)

    # write to xyz for future pca 
    io.write(f'{basedir}/cineb.xyz', io.read(f'{basedir}/cineb.traj@:'))
    if os.path.exists(f'{basedir}/cineb.traj'):
        os.remove(f'{basedir}/cineb.traj')

    positions = opt.atoms.get_positions().reshape(-1, *structure_a.pos.shape)
    traj_neb = []
    for i, image in enumerate(opt.atoms.images):
        traj_neb.append(structure_a.copy(pos=image.positions, surrogate_energy=image.get_total_energy()/HA_TO_EV))
        xyz_file = f'{basedir}/image{i}.xyz'
        io.write(xyz_file, image)
    # end for
# end try
