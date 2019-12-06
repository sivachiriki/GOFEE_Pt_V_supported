from ase.io import read
from ase.calculators.vasp import Vasp
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer,PW
from ase.vibrations import Infrared
water = read('O2_molecule_DFTrelaxed_kpts331.traj')  # read pre-relaxed structure of water

calc = GPAW(mode=PW(500),xc='PBE',
            kpts=(3,3,1),
            symmetry={'point_group': False},
            basis='dzp')

water.set_calculator(calc)
ir = Infrared(water)
ir.run()
ir.summary()
