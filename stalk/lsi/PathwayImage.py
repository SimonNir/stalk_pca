#!/usr/bin/env python3
'''LineSearchIteration class for treating iteration of subsequent parallel linesearches'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from copy import copy, deepcopy
from os import makedirs
from numpy import savetxt, loadtxt, array
from functools import partial

from numpy import ndarray, zeros
from stalk.lsi.LineSearchIteration import LineSearchIteration
from stalk.nexus import NexusPes, NexusStructure
from stalk.params import ParameterStructure
from stalk.params.ParameterHessian import ParameterHessian
from stalk.params.ParameterSet import ParameterSet
from stalk.params.PesFunction import PesFunction
from stalk.pls.TargetParallelLineSearch import TargetParallelLineSearch
from stalk.util.util import get_fraction_error, orthogonal_subspace_basis


class PathwayImage():
    _path = None  # Path is set upon Hessian calculation
    _lsi: LineSearchIteration = None  # Line-search iteration
    _structure: ParameterSet = None
    _hessian: ParameterHessian = None
    _tangent = None
    _subspace = None
    _surrogate = None
    _image_type: str | None = None  # 'eqm' | 'saddle' | 'intermediate'

    def __init__(
        self,
        structure: ParameterSet,
        image_type: str | None = None, 
        surrogate_energy: float | None = None,
    ):
        self._structure = structure
        self._image_type = image_type

        if surrogate_energy is not None:
            self.surrogate_energy = surrogate_energy
        elif structure.surrogate_energy is not None:
            self.surrogate_energy = structure.surrogate_energy
        else:
            self.surrogate_energy = None
    # end def

    @property
    def lsi(self):
        return self._lsi
    # end def

    @property
    def structure(self):
        return self._structure
    # end def

    @property
    def hessian(self):
        return self._hessian
    # end def

    @property
    def surrogate(self):
        return self._surrogate
    # end def

    @property # DEPRECATED; retained for compatibility
    def reaction_coordinate(self):
        return None 
    # end def

    @property
    def image_type(self):
        return self._image_type
    # end def

    @image_type.setter
    def image_type(self, val: str):
        if val not in ('eqm', 'saddle', 'intermediate'):
            raise ValueError("image_type must be 'eqm', 'saddle', or 'intermediate'")
        self._image_type = val
        # end if 
    # end def

    @property
    def structure_init(self):
        if self._tangent is None:
            return self.lsi.structure_init
        else:
            return extend_structure(self.structure, self.lsi.structure_init, self._subspace)
        # end if
    # end def

    @property
    def structure_final(self):
        if self._tangent is None:
            return self.lsi.structure_final
        else:
            return extend_structure_errors(
                self.structure,
                self.lsi.structure_final,
                self._subspace
            )
        # end if
    # end def

    def calculate_hessian(
        self,
        tangent,
        pes: PesFunction,
        path='',
        **hessian_args  # dp=0.001, dpos_mode=False, structure=None
    ):
        hessian_file = f'{path}/hessian.dat'
        makedirs(path, exist_ok=True)
        if self.image_type in ('eqm', 'saddle'):
            # eqm and saddle are calculated in full parametric space
            hessian = ParameterHessian(structure=self.structure)
            subspace = None
            pes_comp = pes
        elif tangent is not None:
            # Calculate orthogonal subspace exluding the tangent direction
            subspace = orthogonal_subspace_basis(tangent)
            # Make a copy of the PesFunction wrapper, then replace the func
            pes_comp = copy(pes)
            pes_comp.func = partial(extended_pes, self.structure, subspace, pes)

             # Create the appropriate subspace structure type based on the original structure type
            if isinstance(self.structure, NexusStructure):
                structure_sub = NexusStructure(params=zeros(len(subspace)))
            else:
                structure_sub = ParameterSet(params=zeros(len(subspace)))

            hessian = ParameterHessian(structure=structure_sub, require_consistent=False)
            # In the future, we can add subspace forward and backward functions to allow 
            # subspace hessian consistency checks
            # After about 4 hours spent attempting to do this, SDN decided to just skip it for now. 
            # If the full forward and backward are consistent, then the subspace equivalents will be, so 
            # checking here is not really necessary given that the eqm and saddle checks will already verify the 
            # overall parameterization validity
        else:
            raise ValueError("tangent cannot be None for intermediate images")
        # end if
        try:
            hessian_array = loadtxt(hessian_file, ndmin=2)
            hessian.init_hessian_array(hessian_array)
        except FileNotFoundError:
            hessian.compute_fdiff(
                pes=pes_comp,
                path=path,
                **hessian_args
            )
            savetxt(hessian_file, hessian.hessian)
        # end try
        self._path = path
        self._subspace = subspace
        self._tangent = tangent
        self._hessian = hessian
    # end def

    def generate_surrogate(
        self,
        pes: PesFunction = None,
        overwrite=False,
        **surrogate_args
    ):
        if self.image_type in ('eqm', 'saddle'):
            pes_sub = pes
        elif self._tangent is not None:
            # Make a copy of the PesFunction wrapper, then replace the func
            pes_sub = copy(pes)
            pes_sub.func = partial(extended_pes, self.structure, self._subspace, pes)
        else:
            raise ValueError("tangent cannot be None for intermediate images")
        # end if

        hessian = self.hessian
        eigvals = hessian.lambdas
        sgn_list = []
        neg_eig_indices = [i for i, val in enumerate(eigvals) if val < 0]
        if len(neg_eig_indices) > 0:
            if self.image_type != 'saddle':
                raise ValueError("Negative Hessian eigenvalue found for non-saddle image. This is not allowed.")
            if len(neg_eig_indices) > 1:
                raise ValueError("Multiple negative Hessian eigenvalues found for saddle image. This should not happen for a true first-order saddle at surrogate level.")
            # For saddle: maximize along negative, minimize along positive
            for i, val in enumerate(eigvals):
                sgn_list.append(-1 if i in neg_eig_indices else 1)
        else:
            # eqm or intermediate: minimize along all
            sgn_list = [1] * len(eigvals)

        path = '{}/surrogate'.format(self._path)
        surrogate = TargetParallelLineSearch(
            load='data.p',
            path=path,
            hessian=self.hessian,
            pes=pes_sub,
            sgn_list=sgn_list,
            **surrogate_args
        )
        surrogate.write_to_disk(overwrite=overwrite)
        self._surrogate = surrogate
    # end def

    def optimize_surrogate(
        self,
        overwrite=True,
        **optimize_args
    ):
        self.surrogate.bracket_target_biases()
        self.surrogate.optimize(**optimize_args)
        self.surrogate.write_to_disk(overwrite=overwrite)
    # end def

    def run_linesearch(
        self,
        num_iter=3,
        path='lsi',
        pes: PesFunction = None,
        add_sigma=False,
        **lsi_args
    ):
        if self.image_type in ('eqm', 'saddle'):
            pes_comp = pes
        elif self._tangent is not None:
            # Make a copy of the PesFunction wrapper, then replace the func
            pes_comp = copy(pes)
            pes_comp.func = partial(extended_pes, self.structure, self._subspace, pes)
        else:
            raise ValueError("tangent cannot be None for intermediate images")
        # end if
        lsi = LineSearchIteration(
            path=self._path + path,
            surrogate=self.surrogate,
            pes=pes_comp,
            **lsi_args
        )
        for i in range(num_iter):
            lsi.propagate(i, add_sigma=add_sigma)
        # end for
        self._lsi = lsi
    # end def

    # DEPRECATED
    def __lt__(self, other):
        return (hasattr(other, 'reaction_coordinate') and
                self.reaction_coordinate > other.reaction_coordinate)
    # end def

# end class


def extend_structure(structure0: ParameterSet, structure_sub: ParameterSet, subspace):
    # If structure0 is a NexusStructure, preserve that type
    structure = structure0.copy(label=structure_sub.label)
    structure.shift_params(structure_sub.params @ subspace)
    if isinstance(structure0, NexusStructure):
        assert isinstance(structure, NexusStructure)
    return structure
# end def


def extend_structure_errors(
    structure0: ParameterSet,
    structure_sub: ParameterSet,
    subspace,
    N=200,
    fraction=0.025
):
    ps_sub = structure_sub.get_params_distribution(N=N)
    ps = []
    for p_sub in ps_sub:
        ps.append(structure0.params + p_sub @ subspace)
    # end for
    ps = array(ps).T
    params_err = [get_fraction_error(p, fraction=fraction)[1] for p in ps]
    structure = structure0.copy(label=structure_sub.label)
    structure.set_params(
        structure0.params + structure_sub.params @ subspace,
        params_err=params_err
    )
    structure.value = structure_sub.value
    structure.error = structure_sub.error
    return structure
# end def


def extended_pes(
    structure: ParameterSet,
    subspace: ndarray,
    pes: PesFunction,
    structure_sub: ParameterSet,
    **kwargs
):
    if isinstance(pes, NexusPes):
        assert isinstance(structure, NexusStructure)
    new_structure = extend_structure(structure, structure_sub, subspace)
    return pes.func(new_structure, **kwargs)
# end def
