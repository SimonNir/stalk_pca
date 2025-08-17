#!/usr/bin/env python

from stalk.lsi.TransitionPathway import TransitionPathway

from params import pes
from run1_neb import traj_neb

basedir = 'tpw'

tpw = TransitionPathway(
    path=basedir,
    all_images=traj_neb,
    active_indices=range(len(traj_neb))
)
tpw.calculate_hessians(pes=pes)
