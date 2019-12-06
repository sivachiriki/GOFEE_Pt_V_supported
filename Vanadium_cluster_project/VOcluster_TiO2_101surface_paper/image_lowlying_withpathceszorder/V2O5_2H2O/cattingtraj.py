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

data1=read('V2O5_2H2O_TiO2_101_DFTlowlyingsorted.traj@:')
for j in range(0,len(data1)):
    print(data1[j].get_potential_energy())

numbers = [0,1,2,6,7,8,9,13,19,20,22,25,37,39,44,46,57,65,67,70,73,83]
traj = Trajectory('V2O5_2H2O_TiO2_101_DFTefxed_rhighlevel.traj','w')
for j,a in enumerate(numbers):
    atoms = data1[a]
    traj.write(atoms)
