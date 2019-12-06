import numpy as np

from ase.optimize import BFGS, BFGSLineSearch
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer, PW
from gpaw import extra_parameters
extra_parameters['blacs'] = True
from gpaw.utilities import h2gpts
from ase.ga.relax_attaches import VariansBreak
from ase.calculators.dftb import Dftb

#from kernels import RBF, ConstantKernel
from gaussComparator import gaussComparator
from featureCalculators_multi.angular_fingerprintFeature_cy import Angular_Fingerprint
from delta_functions_multi.delta import delta as deltaFunc
from GPR import GPR

from ase.io import read, write
from ase.ga.offspring_creator import OperationSelector
from standardmutations_mkb import RattleMutation, PermutationMutation


from prepare_startGenerator import prepare_startGenerator

import sys

from globalOptim import globalOptim

### Load startGenerator + set up mutations ###

sg = prepare_startGenerator()
atom_numbers_to_optimize = sg.atom_numbers
n_to_optimize = len(atom_numbers_to_optimize)
blmin = sg.blmin
slab = sg.slab

mutationSelector = OperationSelector([0.2, 0.6, 0.2],
                                     [sg,
                                      RattleMutation(blmin, n_to_optimize,
                                                     rattle_strength=2, rattle_prop=0.20, min_z=9.0,
                                                     descriptor='RattleMutation_large'),
                                      RattleMutation(blmin, n_to_optimize,
                                                     rattle_strength=1.2, rattle_prop=0.5,
                                                     descriptor='RattleMutation_small')])

### Set up feature ###

# Template structure
a = sg.get_new_candidate()

# Radial part
Rc1 = 6
binwidth1 = 0.2
sigma1 = 0.2

# Angular part
Rc2 = 4
Nbins2 = 30
sigma2 = 0.2
gamma = 2

# Radial/angular weighting
eta = 20
use_angular = True

# Initialize feature
featureCalculator = Angular_Fingerprint(a, Rc1=Rc1, Rc2=Rc2, binwidth1=binwidth1, Nbins2=Nbins2, sigma1=sigma1, sigma2=sigma2, gamma=gamma, eta=eta, use_angular=use_angular)

### Set up calculator ###
calc = GPAW(mode=PW(400),xc='PBE',
            basis='dzp',
            kpts = (1,1,1))

### Set up KRR-model ###
reg = [1e-5]
sigma_a = [3]
sigma_b = [0.3]
kwargs1 = {'sigma': sigma_a[0]}
kwargs2 = {'sigma': sigma_b[0]}
k1 = gaussComparator(featureCalculator=featureCalculator, **kwargs1)
k2 = gaussComparator(featureCalculator=featureCalculator, **kwargs2)

comparator = 0.99*k1 + 0.01*k2

def bias_func(x): return np.mean(x)
delta = deltaFunc(atoms=a, rcut=6)
gpr = GPR(comparator=comparator,
          featureCalculator=featureCalculator,
          delta_function=delta,
          bias_func=bias_func,
          regGS=reg)

### Savefile setup ###

savefiles_path = sys.argv[1]
try:
    run_num = sys.argv[2]
except IndexError:
    run_num = ''
savefiles_namebase = savefiles_path + 'global' + run_num + '_' 


### Load start population if supplied ###

try:
    start_pop_file = sys.argv[3]
    start_pop = read(start_pop_file, index=':')
except IndexError:
    start_pop = None


### Initialize and run search ###
    
optimizer = globalOptim(traj_namebase=savefiles_namebase,
                        MLmodel=gpr,
                        startGenerator=sg,
                        mutationSelector=mutationSelector,
                        calc=calc,
                        startStructures=start_pop,
                        slab=slab,
                        population_size=10,
                        kappa=2,
                        Niter=1000,
                        Ninit=2,
                        zlim=[9.0, 15],
                        dualPoint=True,
                        errorRelax=True,
                        relaxFinalPop=False,
                        max_similarity_dist=0.01,
                        Nuse_pop_as_candidates=3)

optimizer.runOptimizer()
