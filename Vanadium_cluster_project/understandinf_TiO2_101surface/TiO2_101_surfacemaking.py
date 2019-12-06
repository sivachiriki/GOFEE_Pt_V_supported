from sys import argv
from ase.spacegroup import crystal
from ase.build import surface, add_vacuum
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.visualize import view
from ase.io import write,read

slab = read('TiO2.cif')
slab = surface(slab, (1,0,1), 3, vacuum=6)
slab = slab*(1, 2, 1)
#extra_atoms = []
#for atom in slab:
#    if atom.z >= 12.50:
#        extra_atoms.append(atom.index)
#del slab[extra_atoms]
#extra_atoms = []
#for atom in slab:
#    if atom.z <= 12.05:
#        extra_atoms.append(atom.index)
#del slab[extra_atoms]
#slab.positions[:, 2] -= 6.0
write('TiO2_101surface.traj',slab)
view(slab)

