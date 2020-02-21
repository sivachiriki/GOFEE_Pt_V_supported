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

data1=read('Pt7O4CO_Al2O3_KRRfund9l_DFTrelaxedsorted.traj@:')
traj = Trajectory('Pt7O4CO_Al2O3_KRRfund9l_DFTrelaxedsorted.traj','w')
#numbers =[0,1,9,11,13,15,18,19,20,21,24,29,30,32,34,35,36,40,44,45,51,53]
numbers =[0,2,4,6,8,9,10,11,15,17,18,19,21,33,34]
print(len(numbers))
for a,j in enumerate(numbers):
#for j in range(len(data1)):
#    if (j != 9 and j !=12 and j !=16 ):
    traj.write(data1[j])
       #print(data1[j].get_potential_energy())
