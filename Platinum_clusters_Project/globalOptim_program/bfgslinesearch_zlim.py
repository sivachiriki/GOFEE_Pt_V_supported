from __future__ import print_function

# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

import time
import numpy as np
from numpy import eye, absolute, sqrt, isinf
from ase.utils.linesearch import LineSearch
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.constraints import FixAtoms
from ase.utils import basestring


# These have been copied from Numeric's MLab.py
# I don't think they made the transition to scipy_core

# Modified from scipy_optimize
abs = absolute
pymin = min
pymax = max
__version__ = '0.1'


class BFGSLineSearch_zlim(BFGSLineSearch):
    def __init__(self, atoms, restart=None, logfile='-', maxstep=.2,
                 trajectory=None, c1=0.23, c2=0.46, alpha=10.0, stpmax=50.0,
                 master=None, force_consistent=None, zlim=None):
        """
        add zlim to BFGSLineSearch:

        zlim = [z_min, z_max]
        """
        
        self.zlim = zlim

        BFGSLineSearch.__init__(self, atoms, restart=restart, logfile=logfile, maxstep=maxstep,
                                trajectory=trajectory, c1=c1, c2=c2, alpha=alpha, stpmax=stpmax,
                                master=master, force_consistent=force_consistent)
        
    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()
        if hasattr(self.atoms, 'get_curvature'):
            return ((forces**2).sum(axis=1).max() < self.fmax**2 and
                    self.atoms.get_curvature() < 0.0)
        if self.zlim is not None:
            # get indices of non-fixed atoms
            indices = np.arange(self.atoms.get_number_of_atoms())
            for constraint in self.atoms.constraints:
                if isinstance(constraint, FixAtoms):
                    indices_fixed = constraint.get_indices()
                    indices = np.delete(np.arange(self.atoms.get_number_of_atoms()), indices_fixed)
            
            pos_z = self.atoms.positions[indices,2]
            return np.any(pos_z < self.zlim[0]) or np.any(pos_z > self.zlim[1])
        return (forces**2).sum(axis=1).max() < self.fmax**2

