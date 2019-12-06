from ase.calculators.dftb import Dftb
import numpy as np
from scipy.spatial.distance import cdist
from time import time
from custom_calculators import krr_calculator
from ase.io import read, write, Trajectory
from ase.io.trajectory import TrajectoryWriter
from ase.optimize import BFGS, BFGSLineSearch
from bfgslinesearch_zlim import BFGSLineSearch_zlim
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.calculators.singlepoint import SinglePointCalculator
from ase.ga.relax_attaches import VariansBreak
from ase.constraints import FixAtoms

from population import population
from parallel_utilities import sync_atoms
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer,PW
from gpaw import extra_parameters
extra_parameters['blacs'] = True
from gpaw.utilities import h2gpts
import traceback
import sys
import os
import ase.parallel as mpi
world = mpi.world

def relaxGPAW(structure, label, forcemax=0.1, niter_max=1, steps=10):
    # Set calculator 
    # loop a number of times to capture if minimization stops with high force
    # due to the VariansBreak calls
    niter = 0

    # If the structure is already fully relaxed just return it
    traj = Trajectory(label+'_lcao.traj','w', structure)
    while (structure.get_forces()**2).sum(axis = 1).max()**0.5 > forcemax and niter < niter_max:
        dyn = BFGS(structure,
                   logfile=label+'.log')
        vb = VariansBreak(structure, dyn, min_stdev = 0.01, N = 15)
        dyn.attach(traj)
        dyn.attach(vb)
        dyn.run(fmax = forcemax, steps = steps)
        niter += 1
    #print('relaxgpaw over',flush=True)
    return structure

calc = GPAW(mode=PW(500),xc='PBE',
            basis='dzp',
            kpts = (3,3,1))

traj = Trajectory('V2O5_H2OTiO2_101_DFTrelaxed.traj','w')
name ='V2O5_H2OTiO2_101surface_gm'
data = read('V2O5_H2O_TiO2_101_BE2s.traj@:')
for i in range(0,len(data)):
    name ='V2O5H2OTiO2_101surface_isomer_{}'.format(i)
    a=data[i]
    a.set_calculator(calc)
    a_relaxed = relaxGPAW(a, name, forcemax=0.01, niter_max=2,steps=100)
    traj.write(a_relaxed)
