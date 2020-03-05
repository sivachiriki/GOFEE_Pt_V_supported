import numpy as np
from ase.neb import NEB
from ase.optimize import MDMin, BFGS
from ase.io import read, write

from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer,PW
from gpaw import extra_parameters
extra_parameters['blacs'] = True
from gpaw.utilities import h2gpts

from ase.parallel import rank, size

initial = read('2VO3Htrans_initi_TiO21013by1cell2layer_DFTrelxd.traj', index='0')
final = read('2VO3Htrans_final_TiO21014by1cell2layer_DFTrelxd.traj', index='0')

Nimages = 5
n = size // Nimages      # number of cpu's per image
j = 1 + rank // n  # my image number
assert Nimages * n == size

data =read('2VO3Htransmodifiednebimages2.traj@:')
images =[initial]
# Make moving images + set calculator
for i in range(Nimages):
    ranks = range(i * n, (i + 1) * n)
    image =data[i].copy()
    #image = initial.copy()

    if rank in ranks:
        calc = GPAW(mode=PW(500),xc='PBE',
            basis='dzp',
            maxiter=99,
            kpts = (2,2,1),
            txt='neb{}.txt'.format(j),
            communicator=ranks)
        image.set_calculator(calc)
    images.append(image)

images.append(final)

# Initialize neb object
neb = NEB(images, k=0.5, parallel=True, method='eb', climb=False)
#neb.interpolate()

# Converge path roughly before invoking the CI method 
qn = BFGS(neb, trajectory='neb_rough.traj')
qn.run(fmax=0.05)

# Turn on CI
neb.climb = True

# Fully converge the path including an image climbing to the saddle point
qn = BFGS(neb, trajectory='neb_CI.traj')
qn.run(fmax=0.05)

write('neb_done.traj', images)
