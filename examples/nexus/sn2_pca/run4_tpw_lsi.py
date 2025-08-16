#!/usr/bin/env python

from matplotlib import pyplot as plt

from params import pes_dmc
from run3_tpw_surrogate import tpw

tpw.run_linesearches(
    num_iter=3,
    path='dmc_lsi',
    pes=pes_dmc,
    # add_sigma=True,
)

params_init, params_init_err = tpw.pathway_init
params_final, params_final_err = tpw.pathway_final
plt.errorbar(params_init[:, 1], params_init[:, 0], params_init_err[:, 0], xerr=params_init_err[:, 1], color='b', label='PBE pathway')
plt.errorbar(params_final[:, 1], params_final[:, 0], params_final_err[:, 0], xerr=params_final_err[:, 1], color='r', label='LDA+noise pathway')
plt.legend()
plt.show()
