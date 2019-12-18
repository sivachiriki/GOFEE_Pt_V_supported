from ase.io import read, write
from ase.visualize import view

data=read('V2O5_TiO2_101_1by4cell_GMDFTrelaxed.traj@:')
image =data[0]
write('V2O5_TiO2_101_1by4cell_GMDFTrelaxed.xyz',image,vec_cell=True)
