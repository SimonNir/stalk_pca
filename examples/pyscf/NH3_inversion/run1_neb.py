#!/usr/bin/env python3

from os import makedirs

from ase.mep import NEB
from ase.optimize import BFGS
from ase import io

from params import neb_image
from run0_relax_a import structure_relax as structure_a
from run0_relax_b import structure_relax as structure_b
from stalk.params.util import interpolate_params


# Generate base directory for NEB
basedir = 'neb'
makedirs(basedir, exist_ok=True)

# number of intermediate images
n_images = 3

traj_init = interpolate_params(structure_a, structure_b, n_images)
images = [neb_image(structure) for structure in traj_init]

# Try to load from disk
try:
    traj_neb = []
    for i in range(n_images + 2):
        xyz_file = f'{basedir}/image{i}.xyz'
        image = io.read(xyz_file)
        traj_neb.append(structure_a.copy(pos=image.positions))
    # end for
except FileNotFoundError:
    neb = NEB(images, climb=True)
    opt = BFGS(neb)
    opt.run(fmax=0.01)
    positions = opt.atoms.get_positions().reshape(-1, *structure_a.pos.shape)
    traj_neb = []
    for i, image in enumerate(opt.atoms.images):
        traj_neb.append(structure_a.copy(pos=image.positions))
        xyz_file = f'{basedir}/image{i}.xyz'
        io.write(xyz_file, image)
    # end for
# end try
