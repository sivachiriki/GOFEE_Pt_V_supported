from sys import argv
from ase.spacegroup import crystal
from ase.build import surface, add_vacuum
from ase.utils.geometry import sort
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.visualize import view
from ase.io import write

calculate = False or len(argv) > 1

a = 5.116
b = 5.238
c = 5.273
coordinate = [
    (0.27708,   0.04171,   0.20995),
    (0.07604,   0.34747,   0.33143),
    (0.44632,   0.75893,   0.48151),
]
slab = crystal(['Zr', 'O', 'O'], basis=coordinate, spacegroup=14, cellpar=[a, b, c, 90, 99.13, 90])

slab111 = surface(slab, (-1, 1, 1), 5, vacuum=6)
slab111 = slab111*(2, 2, 1)
layer = 3
extra_atoms = []
if layer == 3: 
    for atom in slab111:
        if atom.z > 20 or atom.z < 10.5:
            extra_atoms.append(atom.index)
    del slab111[extra_atoms]
elif layer == 4: 
    for atom in slab111:
        if atom.z > 20 or atom.z < 7.5:
            extra_atoms.append(atom.index)
    del slab111[extra_atoms]

slab111.center()
slab111.center(axis=2, vacuum=6)

fixed_atoms = FixAtoms(indices=[atom.index for atom in slab111 if atom.z < 12])
slab111.set_constraint(fixed_atoms) # Fix bottom two

add_vacuum(slab111, 5)

slab111 = sort(slab111)

if not calculate:
    view(slab111)
    quit()

calculator = Vasp(xc='PBE', kpts=(2, 2, 1), gamma=True, npar=1, setups={'Zr': '_sv'},
                encut=400, algo='Fast',
                ismear=0, sigma=0.05,
                ibrion=2, nsw=400, ediffg=-0.02,
                lreal=False, lcharg=False, lwave=False
                )

slab111.set_calculator(calculator)
final_energy = slab111.get_potential_energy()
#print 'ZrO2-111-2x2-1-2 : {}'.format(final_energy)
