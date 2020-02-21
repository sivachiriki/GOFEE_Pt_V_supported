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

data1=read('Pt7O2_CO_Al2O3_KRRfund9l_DFTrelaxedsorted.traj@:')
traj = Trajectory('Pt7O2CO_Al2O3_KRRfund9l_DFTrelaxedsorted.traj','w')
#numbers =[0,5,7,8,9,11,12,13,16,18,21,24,35,37,43,44,46] # Oxide
numbers =[0,2,4,6,7,9,11,14,15,16,17,18,19,20,21,22,27,28,29,30,31,37,38,39,40,42] # CO
for a,j in enumerate(numbers):
#for j in range(len(data1)):
 #   if (j != 14):
    traj.write(data1[j])
       #print(data1[j].get_potential_energy())
