#!/usr/bin/env python

from stalk.lsi.TransitionPathway import TransitionPathway
from stalk.params.ParameterStructure import ParameterStructure
from ase.io import read 
from params import pes, forward, backward

traj_mep_aligned = read('neb_mep_aligned.xyz@:')

traj_neb = []
for i, image in enumerate(traj_mep_aligned):
    try:
        surrogate_energy = image.get_total_energy()
    except RuntimeError:
        print(f"Warning: no surrogate energy found for image {i}")
        surrogate_energy = None

    traj_neb.append(
        ParameterStructure(
            forward=forward,
            backward=backward,
            pos=image.positions.reshape(-1, 3),
            elem=image.get_chemical_symbols(),
            surrogate_energy=surrogate_energy,
            units='A',
            label=f"image{i}",
        )
    )

basedir = 'tpw'

tpw = TransitionPathway(
    path=basedir,
    all_images=traj_neb,
    active_indices=range(len(traj_neb))
)
tpw.calculate_hessians(pes=pes)
