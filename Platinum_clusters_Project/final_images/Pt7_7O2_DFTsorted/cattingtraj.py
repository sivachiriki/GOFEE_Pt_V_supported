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

traj = Trajectory('Pt7_7O2_Al2O3_KRRfund9l_DFTrelaxed.traj','w')
for i in range(1,12):
    data =[]
    stru = '../Pt7_7O2_{}/Pt7_7O2_Al2O3_KRRfund9l_DFTrelaxed.traj'.format(i)
    data = read(stru+'@:')
    data.sort(key=lambda x: x.get_potential_energy())
    for j in range(0,len(data)):
        traj.write(data[j])

data1=read('Pt7_7O2_Al2O3_KRRfund9l_DFTrelaxed.traj@:')
data1.sort(key=lambda x: x.get_potential_energy())
traj = Trajectory('Pt7_7O2_Al2O3_KRRfund9l_DFTrelaxedsorted.traj','w')
for j in range(2,len(data1)):
    traj.write(data1[j])
    print(data1[j].get_potential_energy())
