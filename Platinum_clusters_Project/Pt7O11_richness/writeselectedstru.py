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

data1=read('Pt7O11_Al2O3_KRRfund9l_DFTrelaxedsorted.traj@:')
traj = Trajectory('Pt7O11_Al2O3_KRRfund9l_DFTrelaxedsorted1.traj','w')
numbers =[0,2,5,6,10,14]
#for a,j in enumerate(numbers):
for j in range(len(data1)):
    if (j != 3):
       traj.write(data1[j])
       #print(data1[j].get_potential_energy())
