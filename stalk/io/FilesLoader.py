#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import loadtxt

from stalk.io.PesLoader import PesLoader
from stalk.params.PesResult import PesResult


class FilesLoader(PesLoader):

    def __init__(self, func, args={}):
        self._func = None
        self.args = args
    # end def

    def _load(self, structure, suffix='energy.dat', **kwargs):
        # Use the job_path attribute of the structure instead of its string representation
        if hasattr(structure, 'job_path'):
            path = structure.job_path
        else:
            path = ''
        value, error = loadtxt(f'{path}{suffix}')
        return PesResult(value, error)
    # end def

# end class
